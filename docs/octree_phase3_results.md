# Phase 3：精确外边界与共面面片

已实现可选 `--octree-boundary exact`。默认 `none`，沿用 Phase 2 rect 的实时候选
路径。此阶段优化的是边界表示，不是把非凸整体的边界平面直接变成 MAX。

新增模块 `include/rokae_demo/octree_boundary.hpp` / `src/octree_boundary.cpp`。
从 fine occupied union 删除内部共享面（包括跨凸块接口），按有向平面用稀疏行区间
合并连续矩形面片。保持内腔、L/U 缺口、薄壁和离散分量；没有增加占据体积。
不做两两全局 merge、不做 mesh convex decomposition，也不引入近似填充。
算法与验收约束见 [方案文档](octree_pipeline_plan.md)。

## 收益与限制

目前 ConvexCluster 都是 AABB，6 个支撑面不可再减。外边界去内部面只用于边界
输出，不删除凸块自身的任何约束，所以分组、半空间、MIN/MAX 树保持不变。
Phase 3 开关前后的凸块和 raw/simplified 树文件逐字节相同，标量场也没有变化。
没有实现非轴向凸近似，也没有宣称只做面片合并就能缩小树。

`boundary_patches.json` 输出精确的矩形多边形铺砌，包括整数/米制角点、有向法向、
offset 和面/凸块来源。带洞或 L 型区域分成多个矩形，不会被一个大多边形填掉。
该铺砌不承诺面片最少；面片平面也不是非凸整体的支撑半空间集合。

`boundary.off` 输出细网格外露面的有向三角形，用于可视化。保留细分以避免合并矩形
拼接造成 T 形接缝；边/点接触仍可能非流形，不能宣称是任意 mesh 布尔运算的有效输入。
OFF 面数不会与合并后面片数相同。单位统一 m，保守性范围是原始点与占据体素，
不包含没有观测到的真实障碍物内部。

## 本机实测

Release (`-O3 -DNDEBUG`)，1655 点相机样例，voxel=0.03 m，543 occupied voxels。
预热 1 次，测量 7 次，轮换后端顺序，随机验证各 1000 点。最终测量不与编译或测试
并行运行。Alpha-tet 使用 alpha=0.09 m、offset=0.03 m、exact plane reduction。

| 指标 | Alpha-tet | Octree rect | Octree rect + Phase 3 |
| --- | ---: | ---: | ---: |
| 凸块 | 101 | 270 | 270 |
| 半空间引用 | 524 | 1620 | 1620 |
| unique planes | 524 | 80 | 80 |
| 树节点 | 626 | 1891 | 1891 |
| 核心构建中位数 | 12.09 ms | 3.06 ms | 3.64 ms |
| 完整验证中位数 | 32.29 ms | 222.17 ms | 222.06 ms |
| 构建＋验证中位数，不含 I/O | 44.30 ms | 225.25 ms | 225.86 ms |
| 进程 wall 中位数，含 I/O | 63.63 ms | 237.67 ms | 246.48 ms |
| 原始点漏包 | 0 | 0 | 0 |

边界阶段的规模变化：3258 个体素有向面，删除 1862 个内部面（931 对），剩 1396 个
外露面，再合并为 639 个矩形面片。共面合并减少约 54.23% 的外露面片数量。
可视化 OFF 保留 2792 个三角形，不把此数误报为 639。

Phase 3 构建子阶段中位数：提取 0.2345 ms，面片合并 0.2764 ms，完整边界阶段
0.5467 ms（含中间容器清理等）。边界证书检查 0.6395 ms，边界文件准备与写出
7.9992 ms。前者分别计入 core / validation，后者计入 output I/O。
各列/子项中位数独立计算，不要求相加相等；完整验证的微小差异不代表优化收益。
Octree 还验证全部体素角点、质心和未剪枝参考树，不与 Alpha-tet 作等工作量验证比较。

两个 Octree 路线都无 false negative / false positive，体积膨胀为 0。
此阶段增加约 0.55 ms 边界构建工作，并未加速生成树；故默认关闭。
含验证总延迟约 226 ms，不能把 3.64 ms 核心构建耗时当作实时端到端延迟。
下一步若目标仍是缩树，应另评估凸块数量的减少或有显式误差预算的保守凸近似，
而不是删除 AABB 必需的支撑面。本阶段不提前实施这些更改。

## 验收

- 8 项 CTest 通过，保留所有原有 voxel / Octree / Alpha-tet 回归测试。
- 穷举非空 2×2×2 占据 mask（255 个），验证面片来源、完整覆盖与方向。
- 单体素、完整块、带内腔的壳、平面孔洞、L/U、薄片、边/点接触、负坐标、
  小数网格、重复/乱序、超大稀疏跨度和随机占据子集。
- 故意破坏面方向/来源、重复面、重复 patch 来源、缺失面片、扩张矩形、溢出跨度
  的验证失败；每个外露面恰好归属一个全满矩形面片的证书不依赖随机采样。
- Python 独立重建面集合与矩形覆盖，并检查 OFF 顶点、朝向、三角形数和有向体积。
- 三种 pruning 模式与 Phase 3 开关组合均验证；零随机点仍做完整证书和 raw point 检查。
- 相机样例的 ConvexCluster 和树前后逐字节相同；Alpha-tet wrap mesh 与树也与
  Phase 2 留存基准逐字节相同，后端算法未修改。
- 新模块的 AddressSanitizer / UndefinedBehaviorSanitizer 检查通过。
  沿用执行环境要求设置 `detect_leaks=0`，不宣称已完成泄漏检查。

## 复现与文件

```bash
build/build_halfspace_tree data/realistic_sparse_camera.xyz \
  output/camera_octree_boundary --pipeline octree --voxel-size 0.03 \
  --octree-pruning rect --octree-boundary exact --validation-samples 1000

python3 scripts/benchmark_octree_phase1.py data/realistic_sparse_camera.xyz \
  output/octree_phase3_benchmark --voxel-size 0.03 --repeats 7 \
  --validation-samples 1000 --compare-boundary

ctest --test-dir build --output-on-failure
```

实际结果：`output/octree_phase3_benchmark/octree-boundary/boundary.off`、
`boundary_patches.json`、`validation_report.json`，以及上级目录 `summary.json/csv`。
每个模式建议使用新目录；关闭边界阶段不删除已有文件，会提示旧边界文件没有更新，
当前运行应以 `boundary_mode` / `boundary_enabled` 为准。
