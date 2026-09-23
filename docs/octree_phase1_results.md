# Octree Phase 1 实施与基线结果

本文为 Phase 1 当时的测量记录。当前默认已启用 Phase 2 精确剪枝；
复现本阶段几何请显式使用 `--octree-pruning none`，原 benchmark 脚本已固定该模式。

已完成：方案细化、稀疏 octree 占据模块、AABB 半空间输出、共享树后端接入、
完整原始点包含检查和基线 benchmark。没有实现 Phase 2 pruning、近似填充或
Phase 3 boundary refinement。

## 代码边界

- `include/rokae_demo/octree_occupancy.hpp` / `src/octree_occupancy.cpp`：
  复用原有 voxelIndex，建立排序去重的 occupied keys、整数网格 dyadic root、
  occupied-only 节点和叶点索引；解析提取每个叶体素的 6 个 AABB 支撑面。
- `include/rokae_demo/convex_cluster.hpp` / `src/convex_cluster.cpp`：
  从 alpha-tet 提取原有 ConvexCluster、Plane、树构造/导出、凸块验证。
  没有重新设计 tetra aggregation 或 MIN/MAX 树；树重复构建时清空旧 leaf IDs。
- `src/octree_pipeline.cpp`：严格 XYZ/参数检查、原始点与体素覆盖证书、
  raw/simplified tree 验证、输出文件及分阶段计时。
- `tests/octree_occupancy_test.cpp` / `tests/octree_pipeline_test.py`：
  数值边界、负坐标、重复/乱序点、单点/共线/共面、U/L 凹槽、深度不足、
  大空跨度稀疏存储、无效输入、导出 JSON 的独立求值与故意破坏结构的拒绝测试。

## 实测

本机 Release (`-O3 -DNDEBUG`)，输入 `data/realistic_sparse_camera.xyz`：
1655 点，voxel size = 0.03 m。每个 backend 预热 1 次、测量 7 次、交替运行，
随机验证样本数均为 1000。Alpha Wrap 使用 alpha = 0.09 m、offset = 0.03 m、
exact plane reduction。时间是本次测量，不是硬实时上界。

```bash
python3 scripts/benchmark_octree_phase1.py data/realistic_sparse_camera.xyz \
  output/octree_phase1_benchmark --voxel-size 0.03 --repeats 7 \
  --validation-samples 1000
```

| 指标 | Alpha-tet | Octree Phase 1 |
| --- | ---: | ---: |
| occupied voxels | 543 | 543 |
| convex cells | 101 | 543 |
| 半空间引用 | 524 | 3258 |
| unique directed planes | 524 | 84 |
| tree nodes | 626 | 3802 |
| 核心构建中位数 | 13.11 ms | 6.04 ms |
| 核心构建 min–max | 10.95–14.66 ms | 5.25–6.77 ms |
| 验证中位数 | 33.82 ms | 414.94 ms |
| 构建+验证中位数（不含 I/O） | 46.92 ms | 421.23 ms |
| 整个子进程 wall 中位数（含 I/O） | 67.29 ms | 437.46 ms |
| raw points outside tree | 0 | 0 |

各列中位数独立计算，因此分项中位数之和不一定等于总时间中位数。
Octree 的 767 个节点包含 543 个叶节点，required depth = 5。
Octree 分项中位数：occupancy/indexing + split 0.45 ms，AABB plane extraction
0.83 ms，共享 tree build + simplification + node serialization 4.76 ms。

这两种几何不是等价近似：Alpha Wrap 输出另一种连续包络，Octree 精确输出
占据体素并集。Octree 的体积为 0.014661 m³，相对于该占据体素并集的膨胀为 0；
不能把它解释成相对于真实障碍物的 0 误差。

## 结论与限制

当前 octree 前端本身较轻，核心构建约为 Alpha-tet 的 46%，但输出树约大 6.1 倍。
84 个 unique planes 不等于 84 个树节点，3258 次叶引用仍然存在。
Phase 1 不能称为已达到实时目标：默认完整验证需要遍历较大的树，此外还检查
每个体素的 8 个角和质心；总验证点为 7542，比 Alpha-tet 的验证工作量更大。
不能直接将两者 validation_ms 当成相同工作负载的速度对比。

Phase 2 的优化目标应是减少凸块/叶引用；精确全满 parent pruning 对稀疏表面
占据可能收益有限，需要实测后再决定是否引入明确误差预算的近似合并。
验证耗时也需另行优化，但不能通过取消完整 raw point containment 来掩盖问题。

共享后端提取后，基准输入的 Alpha Wrap mesh、tetra raw tree、cluster raw tree、
simplified tree 均与提取前逐字节相同；现有 voxel/alpha-tet 回归与新 Octree
CTest 共 4 项通过。

独立的 AddressSanitizer / UndefinedBehaviorSanitizer 模块测试通过。
当前执行环境不支持 LeakSanitizer 的 ptrace 检查，故以 `detect_leaks=0`
运行；不宣称已完成泄漏检测。旧 exact-nef 的 L-prism 回归也通过，输出
2 个凸块、9 个 unique leaves、15 个树节点，精确对称差为空。

完整本次数据和输出几何位于 `output/octree_phase1_benchmark/summary.json`、
`summary.csv`、`alpha-tet/`、`octree/`。
