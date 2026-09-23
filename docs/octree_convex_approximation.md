# 误差驱动的非 AABB 层级凸近似

新增可选 `--octree-approx hull`；默认仍为已有的精确 `rect` 路径。
Alpha-tet、Nef benchmark、共享 `ConvexCluster` 和 MIN–MAX 树结构保留。
本模式已经生成真实带斜面的凸多面体，不是把一组体素简单替换成外接 AABB。

## 流程及保证范围

```text
原始 XYZ（m） → 占据体素 V → 稀疏 Octree → 现有精确 rect 初始分块
  → 自底向上尝试节点内体素角点凸包
  → 检查实际导出半空间几何的体积、填充距离、平面数量
  → 接受候选 / 保留子块
  → ConvexCluster → 既有 MIN(MAX(halfspaces)) → 完整验证
```

参照集合 V 是原始占据体素的闭集并集，不是未知障碍物完整实体。
对最终半空间集合 C，要求 `raw ⊆ V ⊆ C`，允许在明确预算内填充 `C \ V`。
不把未占据位置解释成已经由传感器确认的自由空间；没有射线、法向或自由空间输入。
不能据此声称完成了稀疏观测的实体重建、机器人安全控制器或硬实时实现。

## 两项独立误差预算

1. 体积膨胀率：`(vol(C) - vol(V)) / vol(V) ≤ max_volume_inflation`。
   同时约束每个被接受节点的局部体积膨胀，防止全局预算集中填充一个局部凹槽。
2. 米制填充距离：`sup(x ∈ C) dist(x, V) ≤ max_fill_distance`。
   它衡量新增占据区域离源体素有多远，不是只检查原始点、凸包顶点或随机样本。

距离上界通过有理数覆盖盒细分得到：对盒 B，使用
`min_A max_(x ∈ B) dist(x,A)`，其中 A 遍历原始精确分块。
这一量是到原始体素并集距离的安全上界。证明盒子在候选外侧时可跳过；
若找到候选内部超出预算的中心点则拒绝；否则沿最长轴二分。
达到盒数限制仍不能证明通过时，拒绝本次候选，不放宽阈值。
参考 AABB 可以来自不同原始分块，界可能偏松，因此可能拒绝实际可行的候选。
小的 L 型填充可以被允许；足够宽的 U 型槽会被距离预算阻止，
但没有“任意 U/L 形槽永远保持”的无条件保证。

## 半空间数值处理

候选凸包输入为节点内全部原始占据体素角点，而非体素中心或点云样本本身。
每个理想凸包法向转换为 double 后，重新计算该法向对所有源角点的支撑值，
并向外舍入 offset，留出浮点点积余量。再添加未放大的轴向边界约束，
将几何限制在节点及源角点紧包围盒中。

体积和距离预算作用于这些 **实际 double 系数定义的半空间交集**，
不是舍入前的理想凸包。使用 CGAL 精确构造内核重建交集和计算体积。
理想凸包体积仅用于提前拒绝必然超预算的候选，不能用作最终验收体积。
报告值从精确数值统一转换，避免 Lazy exact 的不同表达式产生不同浮点统计值。

约简只删除精确冗余平面/相同朝向的精确重复支撑面，不做近法向拟合。
新路径以完全相同的法向分量和 offset 复用全局叶子，
不使用旧后端允许的近法向叶子匹配；旧后端的默认行为未改变。
验证重新构造约简后交集，并检查其所有顶点满足被删除的约束。
源体素所有角点满足每个支撑半空间，因此整个源体素被覆盖。

CGAL 接口参考：[Convex Hull 3 / halfspace intersection](https://doc.cgal.org/5.6.2/Convex_hull_3/group__PkgConvexHull3Functions.html)。
实现实际使用工程锁定的 5.6.3 本地头文件，不依赖升级后的 API。

## 层级选择与限制

- 仅在已有节点内 bottom-up 恢复；没有全局两两合并或传统 mesh convex decomposition。
- 至少替代两个当前块，且约简后平面引用数严格下降，才接受新块。
- 保留未通过的子块；每个候选的角点和误差参照始终来自原始体素。
  因此父候选可能收回子近似新增的占据区域，但不会收回原始体素。
- 节点来源与轴向约束保证最终各区域内部不交叠，可将体积相加。
- 固定最细体素分辨率，按误差选择最终多尺度凸块；不是空间自适应 Alpha Wrap。
- 这不是全局最优凸分解，也不是按所有可能方向寻找最优切分。
- 工作量限制是每候选的源体素数和覆盖盒数限制，不是毫秒截止时间或整帧复杂度上界。

## 使用

L 形三体素可复现实例，以下参数只是测试预算，不是机器人安全参数推荐：

```bash
build/build_halfspace_tree tests/data/convex_l.xyz output/l_hull \
  --pipeline octree --voxel-size 0.01 --octree-approx hull \
  --max-volume-inflation 0.2 --max-fill-distance 0.008 \
  --convex-max-voxels 64 --convex-test-max-boxes 512
```

两个误差参数必须显式给出，均可为零；零预算只允许几何不变的恢复。
默认工作量限制为 128 个源体素、4096 个覆盖盒。
`--octree-pruning rect` 可省略，不接受与其他 pruning 模式混用。
`--octree-boundary exact` 在 hull 模式下报错，避免将源占据体素边界误当成近似几何边界。
普通 `--octree-approx none` 与不传此参数等价。

## 输出与验证

- `convex_clusters.json`：实际支撑平面、源体素索引、节点来源、非 AABB 类型、体积、填充距离上界。
- `convex_cells.off`：每个最终凸块的三角化显示网格，保留块间接缝。
  它不是布尔并集的全局外边界；double 网格顶点仅用于显示，半空间是权威几何。
- `octree.json`：完整细占据层级及源体素到最终块映射，标记 `fine_occupancy_reference`，
  不伪装成仅由 AABB 构成的最终树切割。
- `logic_tree_raw.json` / `logic_tree_simplified.json`：既有树 schema，pipeline 标记为 `octree-hull`。
- `validation_report.json`：重新检查来源完整覆盖、实际几何和误差证书；
  再测试全部原始点、源角点、细体素中心及随机点。
- `benchmark.json` / CSV：原始/最终规模、候选接受/拒绝原因、体积/距离预算与实际值、
  核心构建、完整验证和 I/O 的分开计时。

原始树和简化树要求标量与 inside 判断一致；原始精确占据与近似几何只要求不漏包，
新增占据样本单独统计。不同分块的 F 值不一定相等，F 也不等于欧氏距离。
完整几何证书不是靠随机采样通过推断出来的；随机采样是补充实现检查。

## 测试和对比

首轮实测分析见 [octree_convex_results.md](octree_convex_results.md)：树有压缩，但核心构建显著慢于 rect。

```bash
cmake --build build -j1
ctest --test-dir build --output-on-failure
python3 scripts/benchmark_octree_convex.py output/octree_convex_benchmark
```

脚本对三体素 L、256 点 U/椭球表面、1655 点模拟相机数据，
比较 rect、两档距离预算的 hull 和 Alpha-tet。默认预热 1 次、正式 3 次，
串行轮换顺序，记录输入/二进制哈希、每次命令、时间和公共查询结果。
两档 hull 均使用局部/全局体积预算 20%，距离预算分别为 `0.5h` / `0.8h`，
每候选上限 64 个体素、512 个覆盖盒。Alpha-tet 为 alpha=3h、offset=h、exact reduction。
Alpha-tet 没有本模式的同一套体积/距离保证，不应只比树大小就判定优劣。

新增测试覆盖：非 AABB L 棱柱、体积/距离拒绝、零误差恢复、计算上限安全回退、
负坐标和非二进制精确的网格尺寸、U 槽、随机占据掩码、输入次序、
损坏的元数据/平面验证失败、近法向叶子不被错误合并、导出树独立求值及旧默认路径兼容。
