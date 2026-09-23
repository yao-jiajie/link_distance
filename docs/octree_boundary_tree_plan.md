# 修正后的 boundary XYZ / best-axis 方案

后续可选 `--boundary-expression local` 的构建规则、全维符号诊断和失败反例见
[局部构建结果](boundary_local_tree_results.md)。该实验未完成通用严格零水平集目标，不替换本文默认 flat 路径。

flat 新增可选 `--boundary-occupancy voxels`，正常路径不建 Octree，超限才懒构建原回退。
新旧输出兼容性、占据格式与实测对照见 [删除 Octree 对比](boundary_voxel_comparison.md)。

范围：保留默认 `octree/rect`，新增有上限的 X/Y/Z 扫描、自由矩形分组、XYZ 相邻盒合并、
共享树输出与验证。支持单轴和 best-axis；不实现全局最优分解或新的几何内核。

## 严格补集接口

V 是闭占据体素并集。有限闭半空间的 MIN/MAX 树不能以 `<=0` 表示开集 R³\V。
因此不修改现有树语义，也不声称新树独立表示严格补集：

- 新树 `q` 表示 `closure(R³\V)`，仍以 `q<=0` 表示其自身集合。
- 新接口 `boundaryStrictFree()`：q<0 返回 true；q>0 返回 false；q==0 时返回 `!octree.contains(p)`。
- 零值处的 Octree 查询同时排除障碍物边界、保留自由分块的人工接缝。
- 原始点和占据体素的面、边、顶点属于 V，严格自由查询必须返回 false。
- 不使用 epsilon。非占据域不是传感器确认的自由空间，也不是已完成的 CBF 安全认证。
- 工作上限触发时输出原 rect 障碍树，严格自由查询改为 `g>0`，不改变 fallback 的集合语义。

输出必须区分闭包树和严格查询资源。消费者必须加载严格查询 manifest 及占据参考，
不能只对闭包树调用 `evaluateInside()` 来判断严格自由。F 值不是欧氏距离；本原型不实现 CBF 梯度。

## 几何与来源

直接复用原 `extractOctreeBoundary()`。平面注册表是 `(axis, integer coordinate)`，
保留来源 patch 和原始外法向；半空间叶子引用 `(plane_id, coefficient sign)`，允许翻转保留侧。
每轴压缩坐标仅来自该轴法向的真实暴露面，不使用任意网格细分位置。
未引入新平面位置不代表不存在内部切割，单独记录叶子来源，不能以 patch 数代替树大小。

## 三轴扫描与选择

扫描轴 `s` 的局部横轴为 `u=(s+1)%3, v=(s+2)%3`：X→YZ，Y→ZX，Z→XY。
每个事件平面按有向 patch 更新 UV 占据 mask。负向面进入，正向面退出；必须验证原状态，
不使用 XOR 或浮点采样。每个开 slab 对非占据 mask 做两种转置行区间合并，
取矩形较少者，平局固定选局部 v 方向 runs。有限矩形挤出成闭棱柱，再映射回世界 XYZ。
重复**世界坐标 X/Y/Z 固定顺序**的相邻同截面合并，最多指定轮数；达到轮数上限保留当前精确分区，不放松几何。
bbox 外部用六个单半空间分支，不伪造无限坐标 AABB。

`--boundary-sweep x|y|z|best-axis`，默认仍为 x。best-axis 串行构建三个方向，
按 `(是否回退, 简化树节点数, 平面引用数, 有限自由棱柱数, 轴序 X/Y/Z)` 选最小值。
成功的 boundary 候选优先于超限回退；全部超限则使用 X 候选的原 rect 树。
它不与 volume 的树大小自动择优，也不按耗时挑轴，更不保证全局最小树。
所有候选都经过完整整数证书与严格分类探针检查，不能仅验证最终赢家。

默认上限：16384 个压缩三维开单元、1048576 次事件格更新、4096 个合并前自由棱柱、4 轮合并。
前三项超限回退原 rect；合并未收敛只停止进一步压缩，输出是否收敛。
这些是工作量上限，不是毫秒截止时间；提取边界和构建 fallback 本身也要花时间。

## 证明与检查

1. 原始暴露单位面与 patch 精确覆盖及朝向，复用既有整数验证器。
2. 按每个压缩单元内的源体素数量与整数体积验证均匀占据，不靠采样猜测。
3. 验证自由棱柱在所有压缩三维开单元上无遗漏、不覆盖占据、不重叠。
4. 输出集合是有限闭棱柱/半空间并集，以上开单元覆盖及坐标来源确定其闭包。
5. 精确比对重建的树系数和结构，原始/简化树求值一致。
6. 严格查询用零值守卫修复面、边、顶点归属。单元测试额外穷举小例的面/边/顶点及内部接缝。

分别记录核心构建、证书、补充探针、原始/简化闭包树查询、
严格自由查询、零值守卫调用比例。默认 AABB 和非 AABB hull 的行为不变。

## 使用与产物

```bash
build/build_halfspace_tree data/realistic_sparse_camera.xyz output/camera_boundary_x \
  --pipeline octree --tree-source boundary --voxel-size 0.04 \
  --boundary-sweep x --boundary-max-cells 16384 \
  --boundary-max-events 1048576 --boundary-max-prisms 4096 \
  --boundary-merge-passes 4 --validation-samples 1000

build/build_halfspace_tree data/realistic_sparse_camera.xyz output/camera_boundary_best \
  --pipeline octree --tree-source boundary --voxel-size 0.04 \
  --boundary-sweep best-axis --validation-samples 1000

python3 scripts/benchmark_boundary_tree.py output/boundary_xyz_benchmark_new
```

不传 `--tree-source` 或指定 `volume` 保持原行为。`boundary` 不接受 hull、其他 pruning 模式、
关闭精确边界或无效 sweep 值。不默认做查询基准；需要时增加 `--query-benchmark true`。

基准默认包含 15 个场景、volume/X/Y/Z/best-axis 五组；轮换执行顺序，预热 1 次、测量 5 次。
`--cases camera single` 可缩小矩阵。`best-axis` 的核心时间包含共同预处理一次、
三个候选完整构建及选择；`all_candidates_ms` 与 `selection_ms` 单列。
构建分项时间是**所有尝试候选之和**，并非只记录赢家。
`selected_candidate_build_ms` 仅供诊断，不能当作 best-axis 总耗时。
计数项默认描述最终选中结果；`candidate_count`、`candidate_point_checks` 描述全部候选工作。
完整验证、I/O、查询基准另列，`total_ms` 包含构建和完整验证但不含 I/O/查询基准。
三轴集合相同但 MIN/MAX 标量不必相同，验证跨轴**分类**一致，不声称跨轴距离值或梯度一致。

产物：

- `strict_free_space.json`：入口 manifest，定义严格查询的分支规则与资源路径。
- `axis_candidates.json`：请求模式、选中轴、确定性选择规则、逐候选规模/耗时/回退原因/验证状态。
- `free_closure_tree_raw.json` / `free_closure_tree_simplified.json`：既有树 schema，集合是非占据域闭包。
- `occupancy_reference.json`：同一帧的原始闭占据 Octree，零值时需访问所有接触该点的子节点，不能只 floor 到一个体素。
- `boundary_supports.json`：真实轴向支撑平面、来源 patch、原始外法向标志、最终叶子保留侧。
- `free_prisms.json`：有限自由棱柱的网格坐标，额外六个无界半空间不作为有限盒子计数。
- `boundary_patches.json` / `boundary.off`：原始**障碍物**的边界；不是自由空间体网格。
- `validation_report.json` / `benchmark.json`：证书、补充探针、工作上限、分项时间和 optional 查询计时。

回退时 manifest 指向 `obstacle_tree_raw.json` / `obstacle_tree_simplified.json`，严格查询为 `root>0`。
原障碍树的内部平面不满足 boundary-only 来源约束，因此回退时 `boundary_leaf_provenance_applicable=0`，
不把它伪装成通过来源约束的 boundary tree。输出目录可能保留较早运行的其他文件，消费者必须以当前 manifest 为入口。
树与占据参考必须来自同一次构建；在线替换应整体更新，不得混用两帧资源。
零值守卫只修正集合分类，不会消除树标量在人工接缝上的零水平面。
回退也有不同的标量符号约定，必须通过 manifest/查询接口解释。
因此本原型不能不经分析就直接替换 CBF 的标量函数或梯度输入。

## 文件与阶段

新增 `include/rokae_demo/octree_boundary_tree.hpp`、`src/octree_boundary_tree.cpp`、
`src/octree_boundary_tree_pipeline.cpp`。只在原 octree dispatcher 增加路由与默认 volume 参数解析，
没有修改原树数据结构、求值器、rect/hull 算法或 Alpha-tet。

其他新增文件：

- `tests/octree_boundary_tree_test.cpp`
- `tests/octree_boundary_tree_pipeline_test.py`
- `scripts/benchmark_boundary_tree.py`
- `docs/octree_boundary_tree_plan.md`、`docs/octree_boundary_tree_results.md`

修改已有文件：`CMakeLists.txt`（构建/测试目标）、`src/octree_pipeline.cpp`（分发/参数）、`README.md`（入口说明）。

已完成注册表、XYZ 扫描、best-axis、有限轮合并、树/严格查询适配、导出和测试。
历史 X 原型结论保留于 [X 结果报告](octree_boundary_tree_results.md)，
完整五组对照见 [XYZ/best-axis 结果](octree_boundary_xyz_results.md)。默认 AABB、Alpha-tet 均未替换。
