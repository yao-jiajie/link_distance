# Octree Phase 2：精确层级剪枝

本文为仅全满父立方体剪枝的历史结果。当前增强版默认是
`--octree-pruning rect`；显式 `--octree-pruning exact` 可复现本文的几何。

本次历史测量时已实现并默认启用 `--octree-pruning exact`。只恢复完全占满的父节点，
不填充任何未占据体素；近似误差预算、近似合并与 Phase 3 边界 refinement
尚未实现。`--octree-pruning none` 保留 Phase 1。

## 实现

在 `src/octree_occupancy.cpp` 中新增 `pruneOctree` 和
`validateOctreePruning`，复用原有 AABB 半空间提取与共享树后端。

1. 沿反向创建顺序 bottom-up 检查 full 状态：叶体素为 full，父节点仅在
   八个子节点都存在且 full 时才能恢复。每次恢复从 8 块 / 48 个半空间引用
   降为 1 块 / 6 个引用；不做全局邻居两两合并。
2. 线性生成层级上的不相交选中节点集合和 fine voxel → convex cell 归属。
   保留按首个细体素的字典序编号，原始点排列变化不改变输出几何/树编号。
3. `octree.json` 只输出有效节点，被恢复父节点没有 children，其 descendants
   不再输出。节点 ID 保留原值，因此可能不连续，必须按 id 而非数组下标查找。
4. 内存中仍保留未剪枝占据树作独立验证基准，不宣称已释放这部分内存。
   输出 `occupied_voxel_indices`、`source_voxel_ids` 及原始点归属可追溯覆盖来源。
5. 剪枝复杂度 O(octree_nodes + occupied_voxels)。大跨度稀疏根节点不会触发
   dense grid 分配；容量验证也不会直接计算可能溢出的 width³。

## 保守性与 CBF 数值语义

结构证书验证每个源占据体素恰好归属一个 full 输出节点，没有父子重叠、遗漏
或错误填充。再检查每个凸块的六个解析面与共享树 leaf 一致，以及全部角点满足
支撑平面。因此 `P_raw ⊆ V = C`，体积膨胀相对 V 为 0；未观测真实物体内部
仍不在该保证范围内。

完整 raw point containment 和所有 fine voxel 角点/中心检查不会因剪枝减少。
随机验证开启时还在验证阶段构造未剪枝树，在随机点和恢复父节点的中心比较
剪枝前后 inside/outside 判定。未剪枝参考树的构建不计入核心构建时间。

**剪枝的判定等价不代表函数值等价。** 被删除的内部网格面原先可能令 F=0，
合并之后同一位置可有 F<0。测试中 [-2,2]³ 的 64 个细体素合成一个 AABB，
中心 F 从 0 变成 -2。当前相机数据恢复父块中心处观察到的变化约 0.03。
这不是漏包，也不是全局数值误差上界；若直接将 F 用于 CBF，不能假定梯度、
数值余量和控制行为保持不变。仍未修改已有 MIN/MAX 数据结构或评价器。

## 验证结果

- 5 项 CTest 全部通过，包括新增 `octree_pruning`。
- 满 4×4×4 块：64 → 1 凸块，449 → 7 树节点，递归恢复 9 次。
- 同一块缺一个体素：63 → 14 凸块，仅恢复其余 7 个满父块。
- U/L 形、内部空腔、单点、重复点、乱序点、米制小数网格、随机占据子集、
  超大空跨度测试通过；故意破坏 cut、归属或合并计数会被拒绝。
- Phase 1 回退的 raw/simplified 树与修改前逐字节一致；Alpha Wrap 的 mesh
  与 simplified tree 也与历史基准逐字节一致。
- 新剪枝模块的 AddressSanitizer / UndefinedBehaviorSanitizer 测试通过。
  因执行环境的 ptrace 限制，LeakSanitizer 以 `detect_leaks=0` 禁用，不宣称
  已完成内存泄漏检测。

## 本机实测

Release，输入 `data/realistic_sparse_camera.xyz`，1655 个点，分辨率 0.03 m。
每条路线预热一次、重复七次、轮换运行顺序。随机样本数 1000。
Alpha-tet 使用 alpha=0.09 m、offset=0.03 m、exact plane reduction。

| 指标 | Alpha-tet | Octree none | Octree exact |
| --- | ---: | ---: | ---: |
| occupied voxels | 543 | 543 | 543 |
| 凸块 | 101 | 543 | 529 |
| 半空间引用 | 524 | 3258 | 3174 |
| unique planes | 524 | 84 | 83 |
| 树节点 | 626 | 3802 | 3704 |
| 核心构建中位数 | 11.57 ms | 5.31 ms | 5.37 ms |
| 验证中位数 | 30.89 ms | 377.73 ms | 381.01 ms |
| 构建＋验证中位数（不含 I/O） | 42.46 ms | 383.47 ms | 386.38 ms |
| 进程 wall 中位数（含 I/O） | 61.32 ms | 399.99 ms | 401.96 ms |
| raw points outside tree | 0 | 0 | 0 |

统计列分别取中位数，不保证分项中位数之和等于总量中位数。
Octree exact 只找到 2 次八子块恢复，有效 octree 节点由 767 降至 751，
输出树缩小约 2.58%。剪枝本身约 0.02 ms；这组测量没有显示明显的核心提速。

Octree 的 fine-grid 验证覆盖保持不变：none 为 7542 个测试点，exact 为
7544 个（额外两个恢复父块中心）；exact 还对 1002 个点评价未剪枝参考树。
因此验证成本不能仅按缩树比例推断，也不能与 Alpha-tet 验证直接等工作量比较。
相对 occupied union 的 volume inflation 为 0；两种分块的双精度体积求和
差约 6.94e-18 m³，单独记录为 `volume_roundoff_m3`，不视为几何膨胀。

结论：满体素区域可以高效压缩，但稀疏表面数据只有很少的完整八子块，精确
剪枝的收益有限。目前仍不能宣称达到实时总延迟目标。若继续降低树规模，需要
另行评估带明确误差预算的近似策略，不能直接把未占据体素填满。

## 复现

```bash
build/build_halfspace_tree data/realistic_sparse_camera.xyz \
  output/camera_octree_phase2 --pipeline octree --voxel-size 0.03 \
  --octree-pruning exact --validation-samples 1000

python3 scripts/benchmark_octree_phase1.py data/realistic_sparse_camera.xyz \
  output/octree_phase2_benchmark --voxel-size 0.03 --repeats 7 \
  --validation-samples 1000 --compare-pruning

ctest --test-dir build --output-on-failure
```

完整数据位于 `output/octree_phase2_benchmark/summary.json`、`summary.csv`，
三条路线最后一次输出分别在 `alpha-tet/`、`octree-none/`、`octree-exact/`。
