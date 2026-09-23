# Boundary local：受限构建规则、严格符号诊断与结果

本阶段已实现并验证人工零值诊断、面邻接图、局部重写、CLI 导出与固定单轴基准。
**没有完成任意输入的 `q=0 ⇔ x∈∂V` 构建算法。** L 形自由空腔通过，但两个隔开的体素
就能使当前局部模板停止并残留人工零值。这里是特定规则的失败反例，不是所有局部方法不可能成功的证明。

旧 volume/rect、Alpha-tet、boundary flat 均保留；默认行为不变。
没有引入 BSP、R-function、smoothing、Nef 或全局 Boolean 优化。

## 1. 实现位置及调用方式

新增：

- `include/rokae_demo/boundary_local_tree.hpp`
- `src/boundary_local_tree.cpp`
- `tests/boundary_local_tree_test.cpp`
- `tests/boundary_local_pipeline_test.py`
- `tests/boundary_local_benchmark_test.py`
- `scripts/benchmark_boundary_local.py`

修改 `src/octree_boundary_tree_pipeline.cpp` 接入可选模式；CMake 注册库源文件和测试。
共享 `ConvexCluster`、半空间叶子、MIN/MAX 节点结构、求值器和 simplifier 未改。

```bash
# 仅诊断旧 flat 标量，不重写
build/build_halfspace_tree data/realistic_sparse_camera.xyz output/camera_flat_sign_new \
  --pipeline octree --tree-source boundary --voxel-size 0.04 \
  --boundary-sweep y --boundary-expression flat --boundary-diagnostics true

# 实验性局部构建，自动启用完整符号诊断
build/build_halfspace_tree data/realistic_sparse_camera.xyz output/camera_local_new \
  --pipeline octree --tree-source boundary --voxel-size 0.04 \
  --boundary-sweep y --boundary-expression local

python3 scripts/benchmark_boundary_local.py output/boundary_local_benchmark_new
```

`--boundary-expression flat|local` 默认 flat。local 暂不接受 best-axis：本阶段固定一个方向，
避免用旧 flat 树规模为新的 local 表达挑选方向。
保留平面字典和每个叶子的边界支撑平面来源；没有全局删除平面定义。
原始 flat 参考树仍会构建并导出，其成本计入 local 的总核心时间。

## 2. 独立全维符号诊断

定义 V 为原始闭占据体素并集。目标是：

```text
R³ \ V        → q < 0
∂V            → q = 0
interior(V)   → q > 0
```

诊断使用已有整数平面坐标，把每轴分为相对开区间和单个平面位置，包括两端无界区间。
三个轴的笛卡尔积分为体、面、边、顶点 strata；不只抽样 patch 中心。
每个叶子的符号在单个 stratum 内固定，直接用整数顺序关系计算 -1/0/+1，再传播 MIN/MAX。

期望符号由邻接的三维开单元占据决定：全部占据为正，全部非占据为负，混合为真实边界零。
压缩占据单元均匀性通过源细体素数量与整数体积验证，非依赖中心采样猜测。
参考 occupancy 与 boundary tree 必须预先验证且来自同一帧。

原始构建树和简化树都接受独立符号检查。正确性**不依赖 simplifier 消缝**。
符号证书证明注册轴对齐平面在实数语义下的分类；世界坐标双精度求值另外做原始点、
边界和小例全 strata 验证。它不是任意数值误差、传感器噪声或机器人控制安全的认证。

统计定义：

- `artificial_zero_faces/edges/vertices`：严格非占据域内错误为零的对应维度 stratum 数；不是原始面片数。
- `unsafe_free_count`：非自由点被判为负值；`true_boundary_sign_errors`：真实边界不是零。
- `occupied_interior_zero_count`：占据内部错误为零，包含内部体素面、边、顶点。
- `strict_sign_certified`：完整枚举完成，且全部严格符号要求通过。
- `sign_certificate_applied`：是否确实执行完证书；跳过时零计数不能解释为“无错误”。

`boundary_sign_diagnostics.json` 同时保存 flat/local 结果、坐标与最多 24 个失败 witness。
每轴 stratum 编号 `2*k+1` 是第 k 个平面；偶数是开区间，0 和最后一个偶数是无界区间。
这样 witness 不需要虚构无穷远浮点坐标。

## 3. 面邻接与局部规则

邻接图包含：

- 有限自由棱柱之间的正面积共享面，支持部分矩形重叠；
- 有限自由棱柱与 bbox 外部六个无界半空间之间的正面积接口。

按支撑平面分桶检查两侧矩形相交，边/顶点接触不参与合并，但进入符号验证。
外部半空间只在**邻接搜索**时使用有限横向范围，几何始终保持无界。
候选按共享面积优先、固定坐标/编号打破平局。分桶内仍可能有二次比较，设置显式比较上限；
不声称这是无条件线性时间的邻接算法。

每个活动区域保存表达式、保守外包界和原始区域来源。合并后更新区域所有权，沿已有邻接重新尝试。
只接受两侧根 MAX 因子仍存在精确相反切割叶子 `s`、`-s` 的模板。
提取共同/支配外包约束 C 后，残余表达为 a、b。候选构造为：

```text
MAX(C,
    MIN(a,b),
    MIN(a,-s),
    MIN(s,b))
```

证明依据：

```text
MIN(MAX(a,s), MAX(b,-s))
= MAX(MIN(a,b), MIN(a,-s), MIN(s,b), MIN(s,-s))
```

最后一项 `MIN(s,-s)=-|s|` 恒不大于零，所以在这个 MAX 上下文删除它保持 `<=0` 集合。
但标量可以改变：在 `s=0,a<0,b<0` 时新残余值为负，消除人工接口零值。
若公共外约束 C 在该点也为零，零值仍可能保留，不能只依据合并次数宣布成功。

这不是简单删除两个分支中的正/反叶子。真实边界需要的同平面引用可保留在其他上下文。
构建中只做局部常量、相同叶子、同操作节点整理；不进行全局 Boolean 最小化。
随后调用原有 scalar-preserving simplifier；它不负责使零值消失。

## 4. 成功的小例：L 形自由空腔

在一个占据外壳中挖出：

```text
A = [-1,0] × [0,1] × [0,1]
B = [0,1]  × [0,2] × [0,1]
```

单元测试使用上述米坐标；基准按 0.04 缩放。这个例子是**自由区域呈 L 形**，
与基准中 `l_voxels` 的 L 形占据障碍不同，后者仍有外部接缝未消除。

三轴单元测试均验证：

- `(0,0.5,0.5)` 是自由内部：flat 为零，local 严格为负；
- `(-0.5,1.5,0.5)` 是凹口中的占据内部：local 严格为正；
- `(0,1.5,0.5)` 是真实障碍边界：local 仍为零；
- 全部有界和无界 strata 通过，而不只是这三个点。

固定 Y 的 `l_free_cavity` 基准：一个人工零面消失，零边/零顶点均为零；
树从 21 节点降至 17 节点；一次逻辑合并，`strict_sign_certified=true`。

## 5. 最小失败反例：两个分开的体素

米坐标定义：

```text
V = [0,1]³ ∪ ([2,3] × [0,1] × [0,1])
```

两个体素间仅有一个有限自由棱柱 `[1,2]×[0,1]×[0,1]`，它与四个外部半空间有面邻接。
图正确记录 4 个接口；当前规则接受 2 次合并，另 2 次因根切割因子不可匹配而失败。

合并后的表达在整理同操作节点后相当于：

```text
q = MIN(x, 3-x, y, 1-y, z, 1-z,
        MAX(1-x, x-2, -z, z-1))
```

它保留正确的闭包集合，但在 `p=(1.5,0.5,0)` 有 `q(p)=0`。
该点在两个占据体素之间，存在完全非占据的小邻域，因此应该严格为负。
沿 z 的接口项嵌入 MIN 分支，已经不是活动区域根 MAX 的直接因子；当前模板不继续深入重写。

最终还剩 2 个人工零面、4 个人工零边，`strict_sign_certified=false`。
本反例在 X/Y/Z 单元测试均保留，失败 witness 导出到 JSON。
单体素没有这个问题，因此两个体素已经足够暴露本实现的限制。
不能由此推断“纯 MIN/MAX 不能表达此场景”或“所有局部 factoring 都失败”。

## 6. 完整固定轴基准

本机 Release，GCC 11.4，CGAL 5.6.3。固定 Y，h=0.04 m；11 场景 × volume/flat/local 三模式，
共 33 配置，每配置预热 1 次、正式 5 次，串行轮换执行。
查询另进程使用相同 1024 个种子点、7 批 × 5 轮，保留配对 volume 基线。
本阶段收尾于 2026-09-09；原始各轮日志、输入/二进制哈希均保留，不覆盖历史五组基准。

模拟相机：1655 点，341 占据体素：

| 指标 | volume | flat Y | local Y |
| --- | ---: | ---: | ---: |
| 构建原始树唯一节点 | — | 721 | 727 |
| 简化树唯一节点 | 1051 | 721 | 1674 |
| 简化树递归展开节点引用 | — | 721 | 3244 |
| 原始树展开平面引用 | 900 | 618 | 2020 |
| 核心构建中位数 ms | 1.764 | 1.594 | 8.229 |
| 构建＋完整验证中位数 ms | 89.291 | 182.394 | 648.749 |
| 严格非占据查询 μs | 5.630 | 3.425 | 21.398 |
| 人工零面 | 未评估 | 1181 | 763 |
| 人工零边 | 未评估 | 1725 | 1282 |
| 人工零顶点 | 未评估 | 645 | 538 |
| 严格符号证书 | 未评估 | false | false |

local：356 个初始面邻接、82 次成功逻辑合并、287 次失败尝试。
失败尝试按活动区域版本和支撑平面去重；它不是初始邻接数减合并数。
`cancelled_interface_pairs` 统计成功删除非严格恒真切割子项的模板次数，不代表整个邻接图已经消缝。

相机 local **没有达到严格零水平集目标，并且构建/查询明显变慢**。
原始树 727 是按共享指针去重后的节点数；原求值器按引用递归访问，不缓存共享子树结果。
旧 simplifier 会重新构造分支，共享关系可能减少，因此“简化”后唯一节点反而达到 1674，
该分项时间中位数约 5.476 ms。不能只报原始唯一节点数或只报一次成功合并后的小子树。

核心时间包含 flat 参考构建、图构建、全部局部尝试和原结构简化。
完整符号证书、参考树和探针计入验证；I/O 与单独查询基准另计。
flat 本轮启用了符号诊断，因此完整验证时间不能直接与历史未启用诊断的 Y 行混用。
零值守卫仍包含在未认证候选的查询耗时内；不是机器人真实轨迹或硬实时延迟保证。

## 7. 输出、工作量限制和成功标准

local 输出 `free_local_tree_raw.json` / `free_local_tree_simplified.json`，并保留原
`free_closure_tree_*` 参考树；`strict_free_space.json` 是当前帧入口。

```text
local 严格符号认证通过：classifier=root_lt_zero，requires_zero_guard=false
local 集合安全但残余零值：保留原 zero guard，strict_sign_certified=false
原 flat 超限回退：obstacle_tree，classifier=root_gt_zero
local 证书未执行完/几何错误：CLI 报错，不发布未认证的 local 结果
```

`validation_report.json` 的 `passed=true` 表示集合/带适配器分类验证通过；
只有 `local_target_complete=true` 才表示此次 local 原始标量已满足全部严格符号要求。
在 flat 仅诊断模式中，即使诊断发现标量已经严格，仍不更改原 manifest 的 zero-guard 接口。
目录可能保留旧运行的其他产物；消费者必须读取当前成功运行的 manifest，不能猜文件名。

新增默认工作上限：

| 参数 | 默认值 | 含义 |
| --- | ---: | --- |
| `--boundary-max-local-nodes` | 16384 | 构建树递归展开节点上限 |
| `--boundary-max-local-attempts` | 4096 | 局部候选尝试上限 |
| `--boundary-max-adjacency-checks` | 1048576 | 分桶内面配对检查上限 |
| `--boundary-max-strata` | 200000 | 符号分层单元上限 |
| `--boundary-max-sign-ops` | 200000000 | 每次符号证书的估计树/叶子处理量上限 |

表达深度另限 128。局部构建超限只停止进一步重写/保留安全参考，不改几何；
非回退 local 若符号证书超限则报错，而不是把“没有发现错误”当成通过。
这些是有限工作量约束，不是端到端内存、CPU 或毫秒硬截止承诺。

## 8. 验证记录和本阶段结论

核心测试覆盖三轴、L/U/阶梯空腔、闭腔、负/小数坐标、12 组随机体素、逆序输入、
原始/简化标量一致、部分工作上限。另有独立构造的“只有人工零边/零顶点而没有人工零面”诊断用例。
Python 独立读取导出树、世界坐标细盒和边界 patches 检查符号与安全回退。
基准回归保留成功的 L 自由空腔和失败的对齐双体素，不以退出码或 `passed` 冒充严格符号成功。

2026-09-09 收尾复核：Release 全工程编译成功，`ctest --test-dir build --output-on-failure -j2`
共 **17/17 通过**。当前可执行文件 SHA-256 与完整基准 manifest 一致。
相机 volume / flat Y 简化树与历史 XYZ 基准逐字节一致，local 导出的原始 flat 参考树也与独立 flat 运行一致。
另以 Debug、`-fsanitize=address,undefined -fno-omit-frame-pointer` 独立编译核心测试并通过，
运行设置 `ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1`；本次未检查内存泄漏。

此阶段结论是：**受限局部规则可以改善部分形状，但当前不足以作为通用纯边界标量构建器或实时替代品。**
下一步若要求消除嵌套 MIN 上下文中的接缝，需要另行明确允许的规则与复杂度预算并重新论证，
不在本阶段自动扩展。即使某个样例严格符号通过，也不意味着欧氏距离、光滑梯度或 CBF/QP 已验证。

原始记录：

- [完整 33 配置表](../output/boundary_local_benchmark/report.md)
- [统计](../output/boundary_local_benchmark/summary.json)、[输入与二进制哈希](../output/boundary_local_benchmark/manifest.json)
- [相机诊断](../output/boundary_local_benchmark/camera/local/boundary_sign_diagnostics.json)
- [L 自由空腔严格分类入口](../output/boundary_local_benchmark/l_free_cavity/local/strict_free_space.json)
- [双体素最小失败诊断](../output/boundary_local_benchmark/aligned_gap/local/boundary_sign_diagnostics.json)

`benchmark.json` 对应最后一轮；表格是 5 轮中位数。禁止用单轮分项冒充表格中的总时间。
