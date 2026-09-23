# Phase 2 增强：局部精确长方体合并

已实现，默认模式为 `--octree-pruning rect`。保留 `exact`（仅满立方体恢复）
和 `none`（逐体素），没有执行 Phase 3，也没有启用近似填充。

## 实现范围

新增 `include/rokae_demo/octree_rectangular_pruning.hpp` 和
`src/octree_rectangular_pruning.cpp`，复用原来的精确立方体剪枝、AABB 六面
提取、ConvexCluster、MIN/MAX 树与验证框架，不改树数据结构/求值器。

先执行原来的 full-parent recovery，保留更大的满立方体。对剩余的最细
2×2×2 sibling group，枚举 27 种矩形体候选（含单体素、2×1×1、2×2×1、
2×2×2）。以占据 mask 为状态，在 256 状态小表中求不重叠分组的最少块数：
每次选一个包含首个 occupied bit 且完全在 occupied mask 内的候选，再处理余集。
每组是局部最优，不跨 parent 合并，也不声称全局最优。

候选必须全满，不扩大 occupied union，不依赖采样估计决定是否合并。
分组遍历复杂度 O(nodes + occupied voxels)，小表大小固定，不做全局 O(N²) 合并。
最终按首个源体素字典序编号；原始点乱序不改变几何或树编号。

## 正确性与输出语义

验证每个源体素恰好归属一个输出块、整数索引位于块内，每块收到的唯一体素数量
恰好等于其整数容量。配合六个解析支撑面的检查，证明无遗漏、无额外填充且内部
不重叠，即 `P_raw ⊆ occupied voxel union = output geometry`。
保证不涉及未观测的真实障碍物内部。

`octree.json` 仍是精确立方体剪枝后的空间索引。非立方体不能假装成为八叉树
节点，故真实凸块单独保存在 `convex_clusters.json`：

- `extent_voxels=[nx,ny,nz]` 是三个方向的细体素跨度。
- 非立方体的 `width_voxels=null`、`octree_node=-1`，另带 `source_parent_node`
  和 `sibling_mask`；满立方体保留原节点 ID。
- `source_voxel_ids` 对应 `octree.json` 的 `occupied_voxel_indices`。
- 多个空间索引叶节点可指向同一 `convex_cell_id`；节点点索引列表仍恰好覆盖原始点
  一次，不因矩形合并重复列出点。不能用 octree 节点数量代替凸块数量。

内存保留 fine occupancy 作为验证基准，不宣称释放了该部分空间。
完整原始点与所有 fine voxel 角点/中心检查继续保留；随机验证另比较未剪枝树
和所有合并块中心。内外判定必须一致，但 MIN/MAX 函数值可能改变，当前观察到
最大变化约 0.03；该值不是全局上界，不保证 CBF 的梯度/控制行为不变。

## 本机实测

Release，`data/realistic_sparse_camera.xyz`，1655 点，voxel=0.03 m。
每条路线预热 1 次、测量 7 次、轮换顺序；随机验证样本数均为 1000。
Alpha-tet 使用 alpha=0.09 m、offset=0.03 m、exact plane reduction。

| 指标 | Alpha-tet | Octree none | Octree exact | Octree rect |
| --- | ---: | ---: | ---: | ---: |
| 凸块 | 101 | 543 | 529 | 270 |
| 半空间引用 | 524 | 3258 | 3174 | 1620 |
| unique planes | 524 | 84 | 83 | 80 |
| 树节点 | 626 | 3802 | 3704 | 1891 |
| 核心构建中位数 | 11.07 ms | 5.33 ms | 5.27 ms | 2.92 ms |
| 验证中位数 | 29.21 ms | 358.42 ms | 412.16 ms | 202.30 ms |
| 构建＋验证中位数，不含 I/O | 40.24 ms | 364.02 ms | 417.75 ms | 205.36 ms |
| 进程 wall 中位数，含 I/O | 58.70 ms | 379.74 ms | 434.17 ms | 217.01 ms |
| raw points outside tree | 0 | 0 | 0 | 0 |

各列中位数单独计算，不要求分项之和等于总量中位数。本次 rect 相比 exact：
树节点下降约 48.95%，构建时间下降约 44.59%。这是本机本数据的实测，不是硬实时
上界。总延迟仍约 205 ms，不能把 2.92 ms 的构建时间当成含验证的实时延迟。

当前共 543 个占据体素，原立方体剪枝恢复 2 个父块。新增分组影响 130 个 sibling
group，产生 161 个非立方体块（112 个双体素条块、49 个四体素薄块），额外减少
259 个凸块，准确复现先前只读试算的 270 块结果。立方体空间索引仍为 751 个节点、
529 个索引叶节点；270 是最终凸块数，不是空间索引叶节点数。

三个 Octree 模式的 occupied union 完全相同，体积约 0.014661 m³，体积膨胀为 0。
rect 的双精度体积求和差约 1.60e-16 m³，另记 `volume_roundoff_m3`。
rect 验证点数 7705：1655 raw + 543×9 fine-grid probes + 1000 random + 163 merged centers。
其中 1163 个点额外评价未剪枝参考树；不与其他后端的 validation_ms 作等工作量比较。

## 测试

- 6 项 CTest 通过；集成测试固化了相机数据的 270 块 / 1891 节点回归。
- 穷举全部 256 个 sibling mask，使用独立的矩形子集枚举核对最少块数，并检查
  完整覆盖、无重叠。非空 mask 均通过真实 AABB/tree 构造与边界/空体素探测。
- 保留 64 体素满块→1 块；支持条块、薄片、U/L、空腔、多个 parent、重复/乱序点、
  负坐标、小数分辨率、随机子集和超大稀疏跨度。
- 故意破坏体素归属、mask、box、容量或合并计数会被拒绝。
- `none` 和 `exact` 的 raw/simplified 树与修改前逐字节一致；Alpha-tet mesh 和
  simplified tree 与旧基准逐字节一致。
- 新模块通过 AddressSanitizer / UndefinedBehaviorSanitizer；执行环境限制下
  使用 `detect_leaks=0`，未宣称完成 LeakSanitizer 检查。

## 复现

```bash
build/build_halfspace_tree data/realistic_sparse_camera.xyz \
  output/camera_octree_rect --pipeline octree --voxel-size 0.03 \
  --octree-pruning rect --validation-samples 1000

python3 scripts/benchmark_octree_phase1.py data/realistic_sparse_camera.xyz \
  output/octree_rect_benchmark --voxel-size 0.03 --repeats 7 \
  --validation-samples 1000 --compare-pruning

ctest --test-dir build --output-on-failure
```

完整测量保存在 `output/octree_rect_benchmark/summary.json` / `summary.csv`。
最后一次矩形模式的实际几何、表达树和验证报告位于 `output/octree_rect_benchmark/octree-rect/`。
