# 二元 LSE 专用内核

## 本轮范围

这是执行层结构共享之后的第二项优化，只针对二元 MIN/MAX 的查询内核。
不改变点云/体素/边界/支撑平面/canonical cells/原树/执行 DAG、LSE 定义或参数。
完整验证不删减，不引入指数截断、fast-math、倒数乘法或新的距离近似。

相机数据沿用结构共享后的 429 个执行节点：61 个叶节点、297 个二元节点、
71 个三元节点。二元路径减少循环、遍历分支和通用索引开销；三元及其他元数
仍走原通用内核。每次函数求值仍有 439 次 `exp`、368 次 `log1p`，
不是通过减少小权重项来加速。

## 实际实现

`SmoothTreeOptions` 增加第四个参数 `binary_kernel`，在一次查询开始时选择
编译出的通用或二元特化执行循环，避免在每个节点重复读取配置开关。
特化循环仅对 `count==2` 的内部节点使用展开计算，其余节点使用原实现。
没有修改最终 TreeRepresentation 的节点类型，也未生成新的几何切分。

二元节点：

1. 直接读取两个 child 值。
2. 使用与原有序扫描相同的严格大小比较选 anchor；并列时保留第一个出现项。
3. 按原有减法方向计算 gap，再计算 `e=exp(-beta*gap)`。
4. 保留原 `sum=0; sum+=e`、修正量和逐分量梯度加权/除法算式。
5. 继续做指数和、非负修正、上界、有限值及有限梯度检查。

```text
MAX2 correction = log1p(sum) / beta
MIN2 correction = log1p((1-sum)/(1+sum)) / beta
result = anchor + correction
gradient = (anchor_gradient + e*other_gradient) / (1+sum)
```

不使用 `abs(a-b)`、交换求和次序或重新组合归一化公式，不把重复 child reference
删掉。公式的实数保守性、自动 beta、路径误差界以及非精确 SDF 限制不变。
成功查询的值和梯度通过逐位对照及导出 CSV 比较；这不是对任意平台、数学库、
浮点环境的形式化认证，也不承诺数学库 errno/异常标志的完整行为一致性。

工作区大小不变，查询不分配内存。保守工作量预算沿用原来的节点/边估计，
不会因为内核改动而放宽验收限额。没有复制一份特化指令数组或增加每帧索引预处理。

## 使用方法

```bash
build/build_halfspace_tree data/realistic_sparse_camera.xyz output/camera_binary_new \
  --pipeline octree --tree-source boundary --boundary-expression direct \
  --voxel-size .04 --distance-field lse --lse-error .001 \
  --execution-sharing structural --lse-kernel binary --query-benchmark true
```

`--lse-kernel generic` 保留上一阶段的通用内核，用作对照。
当前 direct CLI 启用 LSE 后默认 `binary`；下面的实测并不保证所有场景/查询类型都更快。
`--lse-kernel` 只允许与 `--distance-field lse` 一起使用；非法值和混用报错。
不启用 LSE 时，硬树路径不变。

C++ API 保留旧默认，显式开启两个优化：

```cpp
HS::SmoothTree field(tree.raw, tree.leaves, {0.0, 0.001, true, true});
auto workspace = field.makeGradientWorkspace();
const auto sample = field.evaluateWithGradient(point, workspace);
```

新增指标：

- `smooth_binary_kernel`：是否启用二元内核。
- `smooth_binary_nodes`：实际进入二元路径的执行节点数，generic 模式为零。
- `smooth_generic_nodes`：仍进入通用路径的内部节点数，不含叶节点。

## 验证和 benchmark

单元对照覆盖 MIN/MAX、1/2/3/7 元、重复引用、结构共享开关、并列、signed zero、
最小非正规数、极大有限坐标、极端 beta、非法查询/工作区和叶值溢出。
混合元数的嵌套 DAG 与非轴向法向量覆盖梯度组合/抵消，函数值和梯度分量逐位比较。
原 255 个非空 2×2×2 占据模式增加二元内核与旧实现的符号及值比较。

复用上一阶段成对 benchmark，新增 `--comparison kernel`：

```bash
python3 scripts/benchmark_execution_sharing.py output/lse_kernel_comparison \
  --comparison kernel --errors .001 .005 .01
```

两种内核固定使用相同的 structural 执行共享，轮换 generic/binary 顺序。
核心计时含编译；完整验证、输入点和查询点相同，I/O 和额外查询 benchmark 独立。
每轮核对原树、几何、参数、逐点值/梯度 CSV 的字节一致性；完整距离诊断除时间外
完全相同，连 `work_used` 也不允许变化。

历史 `benchmark_smooth_distance.py` 以及默认 `--comparison sharing` 均显式锁定
`--lse-kernel generic`，避免把两个优化的收益混在一起。

## 修改文件

- `include/rokae_demo/halfspace_smooth_tree.hpp`：二元展开内核、选项和统计。
- `src/boundary_direct_pipeline.cpp`：CLI 和 benchmark 指标。
- 原 C++/CLI 回归测试、`tests/lse_kernel_benchmark_test.py`、`CMakeLists.txt`。
- `scripts/benchmark_execution_sharing.py`：复用成对工具对比 kernel。
- `scripts/benchmark_smooth_distance.py`：历史实验固定 generic，以及相关文档。

## 实测结果

2026-09-09，本机 Release，9 场景 × 3 预算 × 2 内核，共 54 组。
每组 1 次预热 + 3 次构建；查询为 3 次独立进程的中位数，每个进程内使用
7 组 × 5 轮中位数。两种模式顺序逐轮轮换。

[完整测量表](../output/lse_kernel_comparison/report.md)、
[统计与查询记录](../output/lse_kernel_comparison/summary.json)、
[可执行文件和输入指纹](../output/lse_kernel_comparison/manifest.json)。

相机场景固定相同的 429 节点结构共享程序：

| 平滑预算 | 内核 | 核心构建 ms | 构建+完整验证 ms | LSE 值 μs/点 | 值+梯度 μs/点 |
| --- | --- | ---: | ---: | ---: | ---: |
| 1 mm | generic | 1.297 | 281.556 | 12.684 | 13.632 |
| 1 mm | binary | 1.250 | 263.624 | 11.049 | 13.171 |
| 5 mm | generic | 1.212 | 320.318 | 15.387 | 17.049 |
| 5 mm | binary | 1.342 | 306.657 | 13.903 | 16.806 |
| 10 mm | generic | 1.211 | 317.692 | 15.593 | 17.384 |
| 10 mm | binary | 1.273 | 299.825 | 14.219 | 17.133 |

1 mm 档：值查询耗时减少约 **12.9%**，值+梯度仅减少约 **3.4%**。
5/10 mm 档值查询约减少 9.6%/8.8%；梯度仅约减少 1.4%，与计时波动较接近，
不能宣称梯度查询普遍明显提速。

此次没有构建算法或指令预处理优化，核心构建的正负差异不应解释为稳定的构建收益。
完整验证更快主要是重复执行平滑查询的成本下降，不是验证点或证书被删减。
原本的函数值+梯度接口仍只需调用一次，不能把表中两种查询时间相加当成必要成本。

小场景也不保证收益：1 mm 档 closed cavity 的值+梯度由约 0.501 变为 0.523 μs，
薄壁由约 0.228 变为 0.234 μs；存在数个百分点的退步。通用内核因此保留可选，
用户应按实际批量规模、查询类型和机器测量。默认 binary 并不意味着自动择优。

验证结果：

- 31 项 CTest 全部通过，包含独立二元内核端到端 benchmark。
- 7 项核心 ASan/UBSan 全部通过；ASan leak detection 关闭。
- 54 组完整验证全部通过，27 对原树/几何/参数/值与梯度 CSV 均逐字节相同。
- 所有距离误差、有限差分、射线诊断及工作量计数一致，仅计时变化。
- 相机 1 mm 档的优化输出与上一阶段 `execution_sharing_comparison` 中的
  原树、严格自由空间 manifest、函数定义和逐点 CSV 也逐字节一致。

结论：二元内核已完成，值查询有收益，梯度收益有限；欧氏距离精度未改变。
没有进一步采用近似数学函数或低权重截断。direct LSE 流水线现已统一复用一份
不可变 `CompiledTree`，后续优化重点转向验证缓存和多帧复用，而不是继续把循环
展开的收益外推。
