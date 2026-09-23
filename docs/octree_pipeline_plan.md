# Octree 前端方案与分阶段验收

## 表达对象与保证

输入为米制 XYZ 点云。世界坐标原点固定为体素网格原点，分辨率为
`voxel_size`。设原始点集为 P、occupied voxel 的闭 AABB 并集为 V、
输出凸块并集为 C。目标是不变量 `P ⊆ V ⊆ C`，边界视为 inside。
这是对观测点和占据体素的保证，不是对未观测的真实障碍物内部的保证。
稀疏表面点云中的空体素属于 unknown，不能无依据地解释为 free space。

新前端与 `alpha-tet`、`exact-nef` 共存，不调用 mesh convex decomposition。
共享 ConvexCluster、半空间、树构造/化简/序列化和凸块验证；不修改
MIN/MAX 树数据结构，也不改变 Alpha Wrap 的 tetra/aggregation 算法。
`F(x)=min_j max_i(n_ji·x-d_ji)`，`F<=0` 表示占据；这不是光滑距离场，
也不等于已经实现机器人 CBF 控制器。

## 对原提案的修正

- AABB 始终凸，AABB 阶段的 occupancy error 与所提 concavity error
  都等于 `1-volume(V∩node)/volume(node)`，不能当作两个独立指标。
  将其称为 fill error，后续接受条件用 `<=`，零阈值允许无膨胀合并。
- 连通的 unknown 区域不能自动作为 free-space intrusion。后续若要保留
  U/L 凹槽，需区分已知自由空间、unknown 和 occupied，并补充几何长度
  误差上限；不能仅依赖体积比例。没有传感器射线/位姿时不宣称此保证。
- Phase 2 优先自底向上恢复完全占满的 parent；近似恢复另设体积膨胀
  和最大填充距离预算，并以凸块/半空间引用数量下降为收益。
  不做全局 O(N²) 合并，也不重复执行同一个 top-down/bottom-up 判据。
- 删除内部面、合并共面多边形只改变边界表示。非凸边界的所有平面不能
  直接做 MAX，否则会改变几何甚至漏包。Phase 3 必须保留凸块划分，
  或对候选凸包另做保守覆盖与误差检查。
- Octree 和 Alpha Wrap 几何本来可能不同，不能要求两者逐点判定一致。
  一致性要求是各自源几何与各自树、以及树化简前后的一致性。

## Phase 1：精确占据基线（已完成，以 --octree-pruning none 保留）

1. 使用现有 `voxelIndex`：负坐标 floor 语义、有限数检查、网格面边界
   一致性。occupied index 排序去重；不生成给 Alpha Wrap 的 corner 点云。
2. 以最小 occupied index 为 root 下角，边长为覆盖全部 occupied index
   所需的最小 2 的幂。root 的所有子节点沿整数网格二分；只存 occupied
   分支。跨帧 root 可变化，Phase 1 不承诺增量更新或跨帧稳定节点 ID。
3. 叶节点尺寸严格为一个体素。`max_depth` 不足时报错，不偷偷增大体素。
   点索引只保存在叶节点，内部节点保留计数，避免每层复制所有点索引。
4. 一个 occupied leaf 输出一个 AABB ConvexCluster，解析生成 6 个支撑
   半空间。AABB 已经是最小面表示，局部 plane reduction 为恒等操作；
   全局相同有向平面复用 leaf 定义，但不合并互相相邻的凸块。
5. 复用已有 MIN/MAX 树和化简器。不同凸块对同一平面的引用仍单独计数，
   不能把 unique plane 数量误报为树规模。
6. Phase 1 不执行 pruning/refinement，故 `C=V`，以 V 为参考的体积膨胀
   为零；不承诺此时的树比 Alpha Wrap 更小。为后续优化提供正确性基线。

Phase 1 接口：`--pipeline octree --voxel-size <m> --octree-pruning none [--max-depth 20]
[--validation-samples 1000]`。Phase 1 不接受未实现的误差阈值或 Alpha Wrap
参数，以免参数被静默忽略。单点/共面/共线输入有效，空输入或非有限坐标失败。

## 验证与计时

- 完整检查每个 raw point 的 voxel 覆盖和最终 raw/simplified tree 覆盖。
- 检查每个 AABB 的全部 8 个角满足支撑平面，且质心严格在内部。
- 验证 occupied leaf 与凸块一一对应，覆盖所有 occupied voxel，内部不交叠；
  加上解析六平面的构造不变量，提供整体体素覆盖保证，不仅靠随机采样。
- 随机采样及全部体素角点/质心：occupied voxel 几何、凸块几何、raw tree、
  simplified tree 判定一致。失败返回非零，写出 validation_report.json。
- 不同后端的 benchmark 使用同一原始输入和米制分辨率，分别报告几何规模、
  构建时间、验证时间、文件 I/O；随机样本数一致。先比较 Phase 1 基线，
  完整质量/速度权衡属于 Phase 4，不宣称硬实时延迟上界。

统计至少包括 raw_points、occupied_voxels、octree_nodes、leaf_nodes、
convex_cells、planes_before_reduce、planes_after_reduce、unique_planes、
tree_nodes、voxel_resolution、split_time、pruning_time、plane_time、
tree_time、total_time、validation_time、raw_points_outside_tree、volume_inflation。
时间单位 ms，total_time 包括占据索引、树构建和验证，不含输入/输出文件 I/O；
另列 core_runtime_ms 与 input/output I/O。volume 单位 m³。
当前 pruning_time 包含节点选择集合的构造，none 模式也有小的非零开销；
剪枝实际合并次数另列 pruning_merges，none 模式为 0。

## 后续阶段

- Phase 2：全满 parent 精确 pruning，之后评估有显式误差预算的近似 pruning。
- Phase 3：精确外边界提取与共面矩形面片合并（见下文），保留完整凸块支撑面，
  不把非凸边界直接变成 AND。非轴对齐近似凸块不在本阶段实现范围内。
- Phase 4：真实与合成稀疏数据、U/L/薄壁/噪声场景下，对比 Alpha Wrap 的
  时间、点包含率、体积膨胀、树大小和树求值开销。

### Phase 4 实施状态

已实现多场景评估脚本、确定性合成表面点生成器和现有 C++ 树求值器的只读计时工具。
默认矩阵为 5 种形状 × 2 档点数 × 2 档噪声 × 2 档体素尺寸，以及模拟相机点云
的两种尺寸，共 42 组；主对比为 Alpha-tet / Octree rect，其余剪枝模式可选。
参见 [评估协议](octree_phase4_protocol.md)。

分别报告原始点覆盖、体素探针、公共随机查询、理想表面/实体内部和 L/U 凹槽采样，
不把原始点覆盖等同于未观测实体的保守覆盖。原有几何/树后端未改变。
仓库相机点云是模拟数据；支持追加 `--measured-input`，实测扫描数据的验收仍待提供输入，
不能把本轮合成/模拟测试宣称为真实机器人部署验证。误差驱动近似凸块依旧未实现。
本轮 42 组本机统计及发现的问题见 [Phase 4 结果](octree_phase4_results.md)。

### Phase 2 实施状态

已实现精确全满 parent pruning（`--octree-pruning exact`），以及在其基础上
新增的局部满长方体分组（当前默认 `--octree-pruning rect`）。支持 `none`
回退到 Phase 1。未实现近似填充；可选边界 refinement 见 Phase 3。
使用反向节点顺序进行 bottom-up 全满检查，再线性生成不重叠的选中节点集合与
细体素归属；每次恢复由 8 块/48 个半空间引用变成 1 块/6 个引用。
输出层级删除被替代的 descendants，内存保留原始 fine occupancy 作为验证基准。
节点容量检查避免直接计算超大 root 的 w³，从而避免整数溢出。

增强 Phase 2：对没有被更大满立方体覆盖的最细 2×2×2 sibling group，
枚举 27 种满长方体，使用 256 状态小表求最少不重叠分组。支持 2×1×1、
2×2×1 等非立方体 AABB；保留原先可恢复的更大满立方体。不跨 parent
合并，也不声称全局最优。按每个源体素唯一归属、整数容量和局部 mask
校验保证占据并集不变。长方体凸块与立方体空间索引分开输出，不伪造 Octree 节点。

注意：剪枝保持零阈值 inside/outside 判定，不保证 MIN/MAX 函数数值完全不变。
跨 pruning 验证判定一致，并报告随机点和恢复父块中心上的数值变化；同一剪枝结果的
raw/simplified tree 仍需数值一致。该数值变化不作为失败，也不是全局误差上界。

### Phase 3 实施范围：精确边界表示

接口 `--octree-boundary exact|none`，默认 `none` 保持实时构建路径不增加边界开销。
可以与三种 pruning 模式任意组合。启用时输出 phase=3，原 pruning_mode 另存。

1. 对每个 fine occupied voxel 的 6 个有向面查询相邻体素。相邻体素存在则删除
   两侧面，不论两体素是否属于同一个凸块；保留内腔壁。用稀疏索引，不扫描整个 bbox。
2. 按轴向、法向符号、整数平面位置分组，同一平面逐行形成连续区间，只把相邻行
   完全相同的区间扩展为矩形。得到精确、互不重叠的矩形多边形铺砌，不跨空洞或缺口。
   这是确定性共面合并，不承诺最少面片；带洞/L 型平面保留为多个矩形，不伪造填洞多边形。
3. `boundary_patches.json` 保存有向矩形的整数/米制角点、单位法向、平面 offset、
   source_face_ids、convex_cell_ids。边界平面仅是局部有向平面，不是整体支撑集。
4. `boundary.off` 为外露细网格面的三角化，用于可视化；不使用合并矩形直接拼网格，
   避免 T 形接缝。仅边/点接触仍可能非流形；不保证适合作为实体 mesh 布尔输入。
5. 不改 ConvexCluster、六面提取、plane reduction、树结构/求值器。当前 AABB
   的 6 面已经不可再减，边界面片减少不等于半空间引用/树节点减少。不能从凸块内删除
   一个面，只因它在整体内部；否则该块可能向外无限延伸，改变整棵树的几何。
6. 完整验证所有预期外露面的唯一性、来源、方向，并用整数矩形容量与唯一来源集合
   证明每个 patch 全满、所有源面恰好覆盖一次。该证明不依赖随机点数。原始点、细体素
   与树的一致性检查全部保留；集成测试另外检查 Phase 3 开关前后树/凸块逐字节相同。

边界阶段复杂度 O(M log M + F log F)、空间 O(M+F)，M 为 occupied voxel 数，
F 为外露面数（F≤6M）。不做全局 O(M²) 合并，也不调用 mesh convex decomposition。
`boundary_time` 计入 core；其 extraction/merge 子项分别输出。完整边界证书检查计入
validation，boundary JSON/OFF 准备和写出计入 output I/O。没有虚报 tree 收益。

### 近似恢复的约束定义（未实现）

对边长为 w 个细体素的 parent，`E_fill = 1-N_occ/w³`；
最终互不重叠的输出节点集合 S 的总体积膨胀为
`sum(w_s³)/N_occ_total - 1`。近似恢复必须同时满足局部 fill error 与
全局 volume inflation 预算，不能把局部阈值误当作全局误差上限。
精确模式要求 `N_occ=w³`，因此 C 不变。

几何误差采用长度量 `delta(C)=sup_{x∈C\V} dist(x,V)`，单位 m；
它与体积比例不同，能约束窄而深的填充。不能把随机采样距离当作严格
上界：一种可验证的保守上界是对每个拟填充细体素取中心到 V 的距离，
再加半个体素对角线；该做法可能昂贵，只作为后续正确性基线。
其更快的保守层级上界需单独验证，不包含在 Phase 1 中。

若有独立的已知自由空间集合 Q，则另外要求 `C∩Q=∅`（或明确声明非零
侵入预算）；只有 XYZ 时 Q 不可推断。长度/体积限制也不能独自保证
所有凹槽拓扑不变，严格拓扑保留需禁填对应区域，或保持精确模式。
最终接受 parent 还要求替换后的凸块数/半空间引用数减少。
