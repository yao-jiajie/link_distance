# 前端 patch 合并与容器结构优化

这是原优化顺序的第七项。边界 patch 合并和支撑平面注册改用连续数组；在保持原有
矩形划分、编号和来源信息的前提下，减少有序树节点分配和间接访存。

## 实现与不变量

`src/octree_boundary.cpp` 保留原有的外露面提取和邻居二分查询。合并阶段先统计六个
朝向的面数，再为每个朝向一次性预留紧凑记录数组。记录保存整数 `plane/v/u` 和原始
face ID，按 `(plane,v,u)` 排序；朝向按 `(axis,sign)` 遍历。因此扫描顺序与旧版
`map<PlaneKey,map<RowPosition,face_id>>` 完全一致，公开的 `faces` 数组保持原顺序。

每行的连续区间写入可复用的 `previous/current` 数组。两个有序区间序列通过单调游标
匹配；只有端点完全一致、且前一个矩形的 `v1` 等于当前行号时才延长矩形。断行、空洞、
区间变宽或变窄以及相反法向都不会被跨越。来源 face ID 按原顺序直接追加到 patch，
不再为每行创建临时来源列表，也不再为每个行区间分配 map 节点。

合并仍为 `O(F log F)` 时间、`O(F)` 辅助存储；排序后的扫描为线性时间。它不扫描
包围盒中的空体素，不尝试最少矩形划分，也不改变 patch 数量。64 位平台上的排序记录
包含三个 64 位整数和一个 `size_t`；实际峰值 RSS 和分配次数未单独测量。

`src/octree_boundary_tree.cpp` 的支撑平面注册先按轴收集、排序并去重整数坐标，再分配
按 `(axis,coordinate)` 排列的支撑数组。随后按输入 patch ID 顺序二分查找唯一坐标，
追加来源并合并法向位掩码。因此即使调用者打乱 patch 顺序，也与旧版注册结果一致。
持久化坐标数组重新按唯一坐标数量分配，避免缓存快照保留按 patch 数量扩张的容量。
这一步的时间上界为 `O(P log P + P log S)`，其中 `P` 是 patch 数，`S` 是唯一支撑数；
收益来自整数排序和连续存储，测试中也包含大量 patch 共用少量支撑面的场景。

公共 API、CLI 和帧缓存配置均保持兼容。direct、flat/local、Octree 边界路径以及帧缓存
未命中后的重建会自动使用新实现；已经命中的帧缓存无需重新构建这部分几何。
边界独立证书、全部体/面/边/顶点 strata、工作量上限、树构造和 LSE 求值逻辑保留。

## 正式同进程基准

2026-09-09，本机 GCC Release，`h=0.04 m`。每个场景在同一进程内执行两轮预热，
随后交替顺序运行旧 map 实现和生产实现 31 对，取中位数。计时期间没有并行运行本次
构建或测试任务。每对完成后逐项比较结果，并分别验证两份边界的所有外露面和 patch 覆盖。
旧版算法冻结在 `tests/boundary_frontend_reference.hpp`，只供测试及基准使用。

| 场景 | 体素 / 面 / patch | 旧 / 新合并 ms | 合并下降 | 旧 / 新前端 ms | 前端下降 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 相机 | 341 / 782 / 318 | 0.150478 / 0.071856 | 52.2% | 0.330670 / 0.215106 | 34.9% |
| 大平板 | 9,216 / 18,816 / 6 | 1.891379 / 1.330244 | 29.7% | 5.542181 / 4.744416 | 14.4% |
| 实心立方体 | 13,824 / 3,456 / 6 | 0.305016 / 0.150703 | 50.6% | 5.960507 / 5.646625 | 5.3% |
| 棋盘格 | 4,000 / 24,000 / 24,000 | 4.821680 / 2.923857 | 39.4% | 7.565468 / 5.077064 | 32.9% |
| 稀疏分离体 | 30 / 180 / 180 | 0.037384 / 0.010413 | 72.1% | 0.090942 / 0.039920 | 56.1% |

相机的支撑注册从 `0.027142 ms` 降至 `0.019044 ms`，下降约 `29.8%`。
这里的“前端”包含外露面提取、patch 合并、支撑注册以及这些步骤临时容器的释放；
不包含体素化、树构造、证书、查询、最终结果析构或文件 I/O。内部 `merge_ms` 与历史
实现一样，在临时容器析构前取值，故不能直接用内部子项之和代替前端墙钟时间。

14 个场景的前端中位数均下降。实心体素块仍以未修改的邻居查询为主，完整前端收益较小；
微小场景的亚微秒数值对计时开销和系统噪声敏感。本项的收益不能套用于完整流水线或查询吞吐。
完整表格、原始逐轮时间及源码/输入/二进制哈希见
[正式基准报告](../output/boundary_frontend_benchmark_final/report.md)。

## 端到端对照与验证

另保存本次修改前的 `build_halfspace_tree`，对 L 形、闭合空腔和相机分别预热一对，
再交替启动新旧程序 7 对，使用 direct DAG、1 mm LSE 标量误差预算、1000 个随机验证样本
以及原有全 strata 和距离审计。每对运行后核对以下结果：

- 原始/简化树、网格、全部 patch 和来源、支撑平面、分割记录、占据键及发布入口逐字节一致。
- 符号诊断、验证报告、LSE 配置、距离采样 CSV 中的值和梯度逐字节一致。
- 距离审计结果（包括零面射线）一致，排除运行时间及该对照工具已有的缓存计数剔除项。
- 另外逐项比较 `benchmark.json` 的全部非耗时字段，因此缓存计数、工作量、strata、探针、
  梯度检查和误差统计也必须完全一致。

相机 `geometry_core_ms` 从 `0.898292` 降至 `0.777104`。流水线 `total_ms` 中位数为
`193.378904 / 190.496590`，主要由未修改的验证耗时构成；它包含构建及验证，排除输入/输出
I/O 和可选查询基准。该总耗时差含系统噪声，不能全部归因于前端优化。

Release 全部 **37/37 CTest** 通过。新增差分测试覆盖 255 种 `2×2×2` 非空占据组合、
100 组随机体素、三个轴上的区间分裂/合并和断行、薄板、空腔、边/顶点接触、稀疏大跨度、
合法 `int64_t` 极端邻居坐标、owner 元数据、打乱的 patch 输入，以及非法键和轴号。
原有边界变异拒绝测试和所有流水线回归继续通过。

ASan/UBSan Debug 下，`boundary_frontend`、`octree_boundary`、`octree_boundary_tree`、
`boundary_direct_tree`、`boundary_direct_frame_cache`、`boundary_smooth_distance` 六项核心
测试通过；相机同进程前端基准另通过 sanitizer 检查。使用 `ASAN_OPTIONS=detect_leaks=0`
和 `UBSAN_OPTIONS=halt_on_error=1`。

## 复现

```bash
cmake -S . -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
ctest --test-dir build --output-on-failure
python3 scripts/benchmark_boundary_frontend.py \
  output/boundary_frontend_comparison_new --repeats 31
```

若另有修改前的可执行文件，可追加端到端对照：

```bash
python3 scripts/benchmark_boundary_frontend.py \
  output/boundary_frontend_with_pipeline_new --repeats 31 \
  --baseline-binary /path/to/pre-change/build_halfspace_tree --pipeline-repeats 7
```

输出目录必须不存在。默认基准直接包含旧算法，无需旧版可执行文件；可选端到端对照仅
对所选场景中的 `l_voxels`、`closed_cavity` 和 `camera` 执行。已有历史基准产物保留，
今后运行历史脚本时，其共享边界前端会采用当前实现。
