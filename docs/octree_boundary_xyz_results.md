# Boundary XYZ / best-axis：完整五组对照

已补齐 X/Y/Z 扫描和 best-axis，默认 volume/rect、Alpha-tet 与共享树结构保持不变。
**三轴都能完整表示同一个非占据域；换轴改变的是分块及树复杂度，不是几何完整性。**
这里的完整是相对于当前闭占据体素并集，不代表稀疏观测已经完整重建真实障碍物。

## 实现及选择规则

- `--tree-source boundary --boundary-sweep x|y|z|best-axis`，默认仍为 x。
- 单轴扫描采用循环横轴 X→YZ、Y→ZX、Z→XY；后续盒合并固定按世界 XYZ 顺序。
- best-axis 实际构建三个候选，按 `(fallback, tree_nodes, plane_references, finite_prisms, axis)` 选择。
  成功候选优先于超限回退，完全相同则按 X/Y/Z 排序。全部超限时返回原 rect 障碍树。
- 共同预处理只做一次；三个候选的完整构建和选择时间**全部计入**核心耗时。
  每个候选都做整数覆盖证书、原始/简化树求值和严格分类验证。
- `axis_candidates.json` 保存各轴规模、构建分项、回退原因、验证状态与最终选择。
  `benchmark.json` 的阶段耗时为所有候选之和，几何计数描述最终选中结果。
- best-axis 不自动与 volume 择优，也不保证全局最小树或最快构建。

## 实验设置

2026-09-08，本机 Release，GCC 11.4，CGAL 5.6.3；单进程串行，无并行三轴构建。
体素边长 h=0.04 m。15 个场景 × 5 种模式，共 **75 个配置**；每配置预热 1 次、正式 5 次，
五种模式轮换顺序，报告中位数。计时期间没有并行编译或测试任务。
模拟/合成场景包括实心、相邻块、L/U、阶梯、闭腔、窄通道、薄壁、不连通小块、随机簇、
U/椭球表面和模拟相机点云；没有真实传感器扫描、CBF 梯度或 QP 基准。

查询另起进程，在相同的 1024 个固定种子查询点上比较严格非占据判定。
每个方向与 volume 在同一进程配对；7 批 × 每批 5 轮，取每查询中位数，包含零值守卫。
表中 volume 查询时间为四个配对基线估计的中位数；原始配对结果完整保留。
查询计时有系统噪声，相同树的细小时间差不应当作算法收益。

`core` = 占据构建 + 原 rect 分区 + 边界提取 + 全部候选构建 + 选择。
`构建+验证` 额外包含完整证书、参考树和全部探针，不含文件 I/O、可选查询基准。
进程墙钟时间另存。**核心构建毫秒数不是一次 CLI 完整执行时间。**

## 模拟相机：Y 树最小，但 best-axis 构建明显更贵

1655 个原始点，341 个占据体素，782 个暴露单位面，318 个合并 patch。
34 个无向支撑平面，压缩间隔 X/Y/Z = 10/8/13，总共 1040 个开单元。

| 模式 | 树节点 | 平面引用 | 核心构建 ms | 构建+验证 ms | 严格查询 μs |
| --- | ---: | ---: | ---: | ---: | ---: |
| volume rect | 1051 | 900 | 1.681 | 84.166 | 5.137 |
| boundary X | 763 | 654 | 1.596 | 146.132 | 3.737 |
| boundary Y | 721 | 618 | 1.543 | 144.196 | 3.312 |
| boundary Z | 861 | 738 | 1.651 | 150.705 | 4.178 |
| best-axis → Y | 721 | 618 | 3.723 | 363.564 | 3.277 |

Y 相对原 volume 树减少 **31.4%** 节点，相对 X 再减少 **5.5%**。
有限自由盒 X/Y/Z 分别为 108/102/122；原 volume 是 150 个占据盒，两者计数语义不同。
自由盒计数不包含 bbox 外部六个无界半空间分支。
唯一平面定义数 volume/X/Y/Z 为 62/66/67/68：树引用减少不代表平面字典也减少。

best-axis 三候选构建中位数约 2.744 ms，选择约 0.318 ms，选中候选自身构建约 0.809 ms。
这些分项不包含共同预处理，且分别取中位数，不应强求逐项中位数精确相加。
不能把 0.809 ms 当作 best-axis 核心时间；实际是 3.723 ms，约为直接 Y 的 2.41 倍。
当前选择器序列化三个候选统计节点数，选择费用也没有隐藏。

best-axis 的完整整数证书中位数约 5.044 ms；相机执行 7718 次探针（位置可能重复）、三个候选共 23154 次分类检查。
最终 Y 的验证零值守卫调用 3084 次，三个候选合计 9245 次。
均匀随机查询零值守卫比例为 0，不代表机器人轨迹不会遇到边界；另有边界点查询压力测试。
默认完整验证保留，是整次执行的主要时间开销。

## 跨场景结果：没有通用最优方向，也不是处处优于 volume

| 场景 | volume 节点 | X / Y / Z 节点 | 选择 | volume 核心 ms | best-axis 核心 ms |
| --- | ---: | --- | --- | ---: | ---: |
| 单体素 | 7 | 7 / 7 / 7 | X | 0.035 | 0.078 |
| 实心 3³ | 57 | 7 / 7 / 7 | X | 0.150 | 0.134 |
| 闭腔 | 71 | 14 / 14 / 14 | X | 0.231 | 0.241 |
| 两个不连通体素 | 15 | 42 / 42 / 35 | Z | 0.078 | 0.257 |
| 随机体素簇 | 85 | 126 / 126 / 112 | Z | 0.200 | 0.658 |
| U 表面 256 点 | 561 | 329 / 392 / 343 | X | 0.804 | 1.721 |
| 椭球表面 256 点 | 827 | 938 / 924 / 1064 | Y | 1.164 | 4.012 |
| 分散 30 体素 | 211 | 211 / 211 / 211（回退） | X | 0.373 | 1.119 |

椭球即使挑最小方向，仍然比 volume 多 97 个节点；两个不连通体素、随机小簇也明显不适合这一表达。
分散样例各轴都触发 `compressed_cell_limit`，完整回退到原 rect，无漏包，但前置工作和三次回退树构建仍有成本。
成功的 boundary 树不会仅因比 volume 大而被悄悄替换，负结果完整保留。

## 验证与兼容性

- 全部 75 个配置通过。各轴及 best-axis 全部候选的原始点误判自由、探针误判自由/漏判自由、
  原始/简化树标量不一致计数均为 0。
- 每个成功候选通过真实边界来源、压缩单元均匀性、完整自由分区覆盖及树系数的整数证书。
  回退候选按原 rect 树验证，不伪装成 boundary-only 来源证书。
- 小例穷举面、边、顶点和人工接缝；独立 Python 消费导出资源重建严格分类。
- 覆盖三轴平局、部分方向超限、全轴超限、负/小数坐标、输入逆序和非法轴参数。
- best-axis 的最终树与独立构建的选中方向逐字节一致；75 配置基准对此自动检查。
- 全量 Release 构建成功，**14/14 CTest 通过**，包括新五组基准回归测试。
- 新三轴核心测试通过 AddressSanitizer + UndefinedBehaviorSanitizer；容器中关闭 LeakSanitizer，未做泄漏检查。
- 新 X 相机树与历史 X 原型逐字节一致。默认 volume 相机树与历史 X 基准中的 volume、
  以及 Phase 4 h=0.04 m 留存树均逐字节一致。

## 使用边界与结论

当前几何与五组基准已补齐。针对这个相机样例，直接 Y 比每帧 best-axis 更经济；
U 样例却选 X，因此没有把全局默认方向改成 Y，也没有把 boundary 替换成默认后端。
best-axis 适合评估树规模收益，不应仅凭其名字理解为最快实时模式。

新树仍只表示 `closure(R³\V)`。严格非占据分类必须遵守 `strict_free_space.json`：
负值自由，正值占据，零值查询**同一帧**闭占据 Octree；回退时使用原障碍树 `root>0`。
分类一致不等于跨轴 MIN/MAX 标量或梯度一致，人工接缝仍可有零值。
**不能直接把这套分类适配器当作已验证的 CBF 标量/梯度，也不能把未占据等同于已观测安全自由。**

## 复现和文件

```bash
build/build_halfspace_tree data/realistic_sparse_camera.xyz output/camera_boundary_y_new \
  --pipeline octree --tree-source boundary --voxel-size 0.04 --boundary-sweep y

python3 scripts/benchmark_boundary_tree.py output/boundary_xyz_benchmark_new
ctest --test-dir build --output-on-failure
```

- [完整 75 配置表](../output/boundary_xyz_benchmark/report.md)
- [中位数统计](../output/boundary_xyz_benchmark/summary.json)、[输入/二进制哈希](../output/boundary_xyz_benchmark/manifest.json)
- [相机三候选](../output/boundary_xyz_benchmark/camera/best-axis/axis_candidates.json)
- [相机各轮完整记录](../output/boundary_xyz_benchmark/camera/best-axis/trials.json)
- [相机严格查询入口](../output/boundary_xyz_benchmark/camera/best-axis/strict_free_space.json)
- [实现方案与导出约定](octree_boundary_tree_plan.md)、[历史 X-only 结果](octree_boundary_tree_results.md)

单次目录下的 `benchmark.json`、`axis_candidates.json` 对应最后一轮，表格和 `summary.json` 则来自五轮中位数。
历史输出未删除，基准脚本要求新的输出目录，避免覆盖既有实验记录。
