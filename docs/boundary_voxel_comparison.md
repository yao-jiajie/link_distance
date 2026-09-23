# Flat 删除 Octree：实现、对照协议与结果

本次只删除 **flat 正常路径的 Octree 层级构建与依赖**，保留占据体素集合。
旧路径保留作为基准，默认仍为 `octree`；不删除全工程 Octree，也不改变 MIN/MAX 后端。

## 1. 前后流程

```text
旧 octree：
点云 → 体素索引/排序 → Octree 节点 → 精确 rect 分组
     → 暴露面/patches → 坐标压缩/sweep → 自由盒 → flat 树

新 voxels：
点云 → 体素索引/排序去重
     → 暴露面/patches → 坐标压缩/sweep → 自由盒 → 同一 flat 树
```

新路径没有 corner 点集，没有全局 merge，没有新的几何近似。
仅保存排序去重的整数体素索引、分辨率和范围。不分配节点、子指针或点 ID 列表。
为了保持随机查询完全相同，计算一个与旧根盒相同的采样范围；这是一个盒子，不是八叉树。

`OctreeBox` 等现有几何类型名称以及 `--pipeline octree` 分发名继续复用，
不代表运行时创建了 Octree。`octree_built` / `octree_nodes` 明确输出真实情况。

正常 `voxels` 路径的构建、验证、查询和导出均不创建 Octree。
若 sweep 超限，才从原输入懒构建原 Octree/rect，输出相同的旧障碍树；这部分耗时计入核心总时间。
因此回退不是无 Octree 路径，也不保证更快；此时两份占据存储都保留。

## 2. 查询、验证和兼容性

新 `SparseVoxelOccupancy::contains()` 判断的是**闭占据体素并集**，不是单次 floor 的半开归属。
先检查 bbox，再按实际导出的 `double(index)*voxel_size` 平面二分定位。
每轴最多两个候选，面/边/顶点处最多检查 8 个体素索引；不使用 epsilon。
避免在最大支持索引的上边界调用超范围的点索引函数。

flat 的标量及查询规则不变：负值非占据，正值占据内部；零值还须查询原始占据。
本改动**不消除人工零值，也不改变树节点数量**。它不构成 CBF、距离或梯度认证。

新旧路径共用：

- 暴露面和矩形 patches 的提取/整数覆盖验证；
- 支撑平面注册、X/Y/Z/best-axis 扫描、有限自由盒精确合并；
- 原始/简化树构建、半空间来源验证；
- 全维符号诊断，包括无界区间、人工零面/边/顶点。

新路径用 raw point→voxel key 和闭盒包含验证取代 Octree 节点证书。
`--boundary-validation-reference voxels` 可让**两种模式**都用同一个细体素 AABB 树作独立验证参考，
避免把“旧 rect 参考树”和“新细体素参考树”的不同工作量混进对比。
该参考只在验证阶段构建，不进入核心时间；不创建 Octree。
默认 `native` 时旧模式仍用原 rect 参考，新模式用细体素参考，因此完整验证耗时不宜直接混比。

导出变化：

- MIN/MAX raw/simplified 树、自由盒、支撑平面和 `boundary.off` 与旧路径逐字节一致。
- `strict_free_space.json` 增加 `occupancy_backend` 和 `occupancy_schema`。
- 新 `occupancy_reference.json`：schema 2、`representation=sorted_voxel_indices`、
  `root_id=null`、`nodes=[]`，保留 `voxel_size` 和 `occupied_voxel_indices`。
- 新边界 JSON 的 `face_owner_kind=voxel`；遗留字段 `convex_cell_id(s)` 指体素 ID，
  不再指旧 rect 分组。面几何、方向和 patch 来源未变，不能拿该字段反推旧凸块编号。
- 旧 Octree 占据格式保持；消费者不能对 schema 2 继续从节点 0 开始遍历。

`voxels` 当前只支持 flat；与 local 混用直接报错。
`--max-depth` 不约束没有层级的新正常路径，只有触发旧 Octree 回退时才适用；
回退深度不足会报错，不暗中改变体素分辨率。

## 3. 使用

```bash
# 删除前：Octree + rect
build/build_halfspace_tree data/realistic_sparse_camera.xyz output/flat_before_new \
  --pipeline octree --tree-source boundary --voxel-size 0.04 \
  --boundary-sweep y --boundary-expression flat --boundary-occupancy octree \
  --boundary-validation-reference voxels --boundary-diagnostics true

# 删除后：直接占据体素集合
build/build_halfspace_tree data/realistic_sparse_camera.xyz output/flat_after_new \
  --pipeline octree --tree-source boundary --voxel-size 0.04 \
  --boundary-sweep y --boundary-expression flat --boundary-occupancy voxels \
  --boundary-validation-reference voxels --boundary-diagnostics true

# 完整固定轴配对基准，输出目录必须是新的
python3 scripts/benchmark_boundary_voxels.py output/boundary_voxel_comparison_new

# 相机较多轮确认
python3 scripts/benchmark_boundary_voxels.py output/boundary_voxel_camera_confirmation_new \
  --cases camera --repeats 31 --warmups 3 --query-repeats 7
```

核心实现：`include/rokae_demo/sparse_voxel_occupancy.hpp`、`src/sparse_voxel_occupancy.cpp`。
现有 boundary 提取/构建函数拆出共用数据接口，旧包装函数保持；CLI 增加占据来源选择和懒回退。
没有重写 sweep、树结构、求值器或 simplifier，旧 volume/rect、Alpha-tet、boundary local 均保留。

## 4. 对比协议

同一 Release 二进制、固定 Y、h=0.04 m；16 场景 × 两种占据后端，共 32 配置。
每场景 1 轮配对预热、7 轮正式配对，串行交替先后顺序；不在计时期间并行编译或测试。
查询另运行 3 轮配对进程，每项每轮取 7 批 × 5 次循环的中位数。
构建样本每次启动新 CLI 进程；预热不能摊销进程内静态表的首次初始化。
这是本机批处理构建对比，不是长期运行 ROS 节点的稳态延迟或硬实时上界。

两模式使用相同细体素验证树、原始点/细盒/patch/压缩单元探针、1000 个随机验证点和完整符号诊断。
查询包括 1024 个固定种子随机点、真实边界点，以及从自由盒面中心筛选的人工零接缝点。
三组查询点写入 `query_samples.json`，逐字节比较；没有人工接缝时 seam 时间为 0，表示未测而不是零成本。

每轮比较 raw/simplified 树、自由盒、支撑注册表、边界网格逐字节一致，检查几何计数及符号诊断一致。
正常新路径必须 `octree_nodes=0`；回退样例必须明确标记已构建 Octree，不能算成删除成功。

`occupancy_storage_bytes` 是占据对象及其 vector capacity 的持有存储估算：
包括旧节点/叶子索引/点 ID，或新体素索引数组；回退时包括两者。
**不包含 allocator 开销、rect 分组、临时数组、整棵表达树，不是进程 RSS 或峰值内存。**

## 5. 实测记录

首轮完整矩阵已保存于 [32 配置报告](../output/boundary_voxel_comparison/report.md)、
[summary](../output/boundary_voxel_comparison/summary.json)、
[manifest](../output/boundary_voxel_comparison/manifest.json)。不覆盖历史 local/XYZ 基准。

相机输入 1655 点、341 占据体素。首轮核心时间中位数为 1.555 → 1.457 ms，
持有占据存储估算为 107.02 → 38.97 KiB，Octree 节点 456 → 0；最终树始终为 721 节点。
预处理减少不等于总时间按相同比例减少：后续边界/扫描/树构建仍占主要部分。
完整验证在本轮启用了细体素参考和符号诊断，不能直接与历史不同配置的总时间比较。

首轮 `sparse_disconnected` 触发原回退：核心中位数 0.529 → 0.573 ms，
新模式仍建 204 个 Octree 节点，并额外保留稀疏体素集合。这是预期的负面结果，不删除或隐藏。

相机单次波动较大，因此追加 3 轮配对预热、31 轮正式配对、7 轮查询进程确认。
完整结果保留于 [相机确认报告](../output/boundary_voxel_camera_confirmation/report.md)、
[相机统计](../output/boundary_voxel_camera_confirmation/summary.json)，没有替换首轮矩阵。

| 相机确认指标 | 原 Octree | 无 Octree |
| --- | ---: | ---: |
| 占据索引/结构构建中位数 ms | 0.303612 | 0.172824 |
| 正常路径 rect 分组中位数 ms | 0.074664 | 0 |
| 核心构建中位数 ms | 1.568554 | 1.379067 |
| 核心单次 min–max ms | 1.461277–1.866237 | 1.259369–1.641217 |
| 构建＋完整验证中位数 ms | 259.878 | 260.286 |
| 占据持有存储估算 KiB | 107.02 | 38.97 |
| Octree 节点 | 456 | 0 |
| 简化树节点 / 平面引用 | 721 / 618 | 721 / 618 |
| 随机点严格查询 μs | 3.761 | 3.652 |
| 人工零接缝严格查询 μs | 3.834 | 3.998 |
| 真实边界占据查询本身 μs | 0.079 | 0.174 |
| 人工零面 / 零边 / 零顶点 | 1181 / 1725 / 645 | 1181 / 1725 / 645 |

核心中位数减少约 **0.189 ms / 12.1%**；占据持有存储估算减少约 **63.6%**。
这些分项各自取中位数，不应要求中位数相加严格等于总中位数。
随机查询的零值守卫调用比例为 0，因此该项的小幅变化不能归因于新占据查询更快。
369 个人工零接缝查询点完全相同；这组结果显示新闭体素查询可能稍慢。
真实边界占据查询本身也比层级查询慢，不能宣称查询全面提速。
完整验证仍占绝大多数总耗时；删除 Octree 没有解决这部分成本。

结论：此修改主要节省正常路径的结构构建和占据存储，几何及标量输出保持不变。
回退场景会补建旧结构，可能更慢；保留 opt-in，未切换默认后端。

## 6. 回归记录

Release 全工程编译成功，20/20 CTest 通过，包括：

- 新闭体素查询对照原 Octree 与独立细盒穷举；覆盖面、边、顶点、两侧 `nextafter`、
  负/小数坐标、接近 ±2^50 的支持索引、重复点及随机体素。
- X/Y/Z/best-axis 的输出树逐字节相同、导出后独立解释 schema 2。
- 三种 sweep 工作上限及懒回退、正常路径 `max-depth=0`、非法参数和 local 混用拒绝。
- 配对基准回归，保留人工接缝与回退负面结果；旧 local/XYZ/volume/Alpha-tet 测试继续通过。

相机旧路径 raw/simplified 树还与历史 local/XYZ 基准逐字节一致。
两轮完整基准及当前 Release 可执行文件的 SHA-256 均为
`4c766d30076710922f793db272a401725f65562175ce7e3180527f2ab62ecfd1`。

另以 Debug、`-fsanitize=address,undefined -fno-omit-frame-pointer` 编译并运行
`sparse_voxel_occupancy`、`boundary_local_tree`、`octree_boundary_tree`，3/3 通过。
运行设置 `ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1`；未进行内存泄漏检测。
