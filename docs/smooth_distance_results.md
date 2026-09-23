# 表达树 → 保守 LSE 距离代理函数

## 范围与结论

新增的是现有 boundary-direct DAG 的执行适配层，提供函数值和解析梯度。
不重建体素、边界或表达树，不修改 MIN/MAX/LEAF 数据结构，不调用简化器。
后续添加了 [执行层相同子表达式共享](execution_sharing_results.md)，保留相同函数。
本页历史性能测量及 `benchmark_smooth_distance.py` 显式使用 `--execution-sharing none`
和 `--lse-kernel generic`。当前 direct CLI 使用 `structural` 和
[二元 LSE 内核](lse_binary_kernel_results.md)，不要把各阶段执行耗时混为同一版本。
原有 exact tree、严格符号证书和 `strict_free_space.json` 保持原样。
旧 flat、local、volume/rect 和 Alpha-tet 路径不变。

输出是以米为单位的**平滑有符号距离代理**，而不是经过距离重初始化的真实欧氏 SDF。
正值表示自由；负值覆盖原始占据区域，并可能覆盖原边界外的一层自由空间。
平滑后的零等值面通常向外移动，不再声明 `d=0 iff x∈∂V`。

距离基准是当前**闭合体素并集**的有限表面，不是点云点到点距离，也不是未观测实体的真实表面。
点云稀疏采样、体素化与真实物体之间的误差不在本次 SDF 误差中。

## 函数定义与保守性

原 boundary-direct 树为 `q`：自由区域负，占据内部正；硬距离代理为 `d0=-q`。
叶值保持 `nᵀx-b`，本路径法向量为单位轴向量。
对一个有 `k` 个子节点的内部节点，使用：

```text
MAX_upper(a) =  log(sum exp( beta * a_i)) / beta
MIN_upper(a) = -log(mean exp(-beta * a_i)) / beta

q_upper = 将 q 的内部节点替换为上述算子
d_lse   = -q_upper
```

两种算子都不小于相应硬最值，且对子输入单调。
从叶到根归纳，实数算术中 `q_upper >= q`，所以 `d_lse <= d0`。
因此 `d_lse > 0` 不会把占据点声明为自由；这是单侧保守的自由空间判据。
普通未经归一化的 soft-min 不具备上述上界，不用于这里的 MIN。

逐节点累计误差权重：

```text
W_leaf = 0
W_node = max(W_child) + log(child_count)
0 <= d0 - d_lse <= W_root / beta
```

`--lse-error eps` 自动由 DAG 路径权重选择 `beta`，并向上取相邻浮点值。
`eps` 控制的是平滑引入的**标量误差**，不是几何 SDF 误差或零面位移。
`beta` 的单位为 `1/m`，`1/beta` 是以米为单位的平滑尺度。

这里的保守性是实数公式的证明。实现使用 double 与标准数学库，检查有限值、
移位指数和、节点上界、采样保守关系及误差界；没有区间算术数学库认证。
采样中 `d_lse > d0` 或占据点 `d_lse > 0` 会直接报错，不用分类容差放行。
累计误差界和有限差分检查另有显式浮点舍入容差，不等于安全认证。

## 稳定求值与梯度

选取节点最值 `anchor`，只计算非正指数，避免正指数溢出。
除一个选中的 anchor 出现项以外，令 `t = Σ exp(-beta*gap)`。
对应修正量为：

```text
MAX: log1p(t) / beta
MIN: log1p((k-1-t)/(1+t)) / beta
```

二者均非负。MIN 在全部并列时修正量为零，避免相近对数相减带来的负修正。
不通过 `max(result, hard_value)` 之类硬截断伪造保守的“平滑”函数。
指数负溢出允许贡献成为零；非有限叶值、结果或梯度拒绝。

节点梯度为子梯度的指数权重凸组合；MIN 使用 `exp(-beta*a)`，MAX 使用
`exp(beta*a)` 的归一化权重，根节点再取负号。单位叶法向量使实数公式满足
`||∇d_lse|| <= 1`，但梯度可能很小甚至为零，**不保证 `||∇d||=1`**。
不把本实现称为已验证的 CBF 控制器，也不自动归一化梯度或除以其范数。

每个存储的 DAG 节点只算一次；每次出现的 child reference 都参与求和。
共享节点缓存不展开，也不消除重复引用，因为硬最值的幂等律、吸收律与任意
MIN 重组在 LSE 下不保持标量值。平滑结果绑定导出的树结构与 beta。
查询复用调用者工作区，无单次查询的堆分配；程序可共享，但线程各自使用工作区。

## C++ 接口与命令

```cpp
#include <rokae_demo/halfspace_smooth_tree.hpp>

namespace HS = rokae_demo::halfspace;
// existing_tree 为已经严格验证的 boundary-direct TreeRepresentation。
HS::SmoothTree field(existing_tree.raw, existing_tree.leaves, {0.0, 0.001});
auto scalar_work = field.makeValueWorkspace();
auto gradient_work = field.makeGradientWorkspace();
HS::Vec3 p{0.1, 0.2, 0.3};
double distance_proxy = field.evaluate(p, scalar_work);
auto sample = field.evaluateWithGradient(p, gradient_work);
// sample.value: m；sample.gradient: 对世界坐标 x,y,z 的导数。
```

```bash
build/build_halfspace_tree data/realistic_sparse_camera.xyz output/camera_lse_new \
  --pipeline octree --tree-source boundary --boundary-expression direct \
  --voxel-size 0.04 --distance-field lse --lse-error 0.001 \
  --query-benchmark true
```

也可将 `--lse-error 0.001` 替换为 `--lse-beta 10000`。二者必须且只能指定一个，
非正值、NaN、无穷值及混用报错。默认不启用 LSE；`--distance-field exact` 仅表示
保留原硬树，不表示精确欧氏距离。该原型只挂接在 boundary-direct 模式。

附加验证参数：

- `--sdf-max-work`：估算工作量上限，默认 100000000；超出拒绝，不跳过验证。
- `--sdf-gradient-samples`：有限差分点数上限，默认 128；每点最多三个方向。
- `--sdf-zero-rays`：边界面片中心外向射线数上限，默认 64。

后两者可以为零以显式不做该附加诊断；原硬树符号认证、原始点和距离探针检查仍执行。
过大的 beta 若使所有所请求有限差分步长不可表示，则该验证失败。
新数据请使用新输出目录：失败运行不会清除旧目录里的历史文件。

## 距离基准、验证和输出

`BoundaryDistanceReference` 对每个真实暴露矩形面片做坐标截断，计算查询点到
**有限矩形**的欧氏距离，再取最小值。用同一份原始闭体素占据判定赋符号：
外正、内负、真实边界零。它只作为独立基准，不进入 LSE 函数或梯度计算。
使用前由现有 exact 证书验证面片与占据快照一致。

探针包含全部原始点、所有占据体素中心和角点、每个面片中心及两侧偏移点、
包围盒外侧的确定性棱角探针、固定种子随机点。分别报告：

- `smoothing_error_*`：`|d_lse-d0|`。
- `hard_sdf_error_*`：`|d0-d_ref|`，即平滑前已有的距离畸变。
- `sdf_error_*`：`|d_lse-d_ref|`。
- `near_sdf_error_*`：`|d_ref| <= voxel_size` 子集的误差。
- 内外误差、真实边界值偏差、梯度范数、有限差分误差、错误自由分类数。

误差统计包括样本数、均值、RMSE、P95、最大值；均为采样而非全空间上界。
外向零面射线在两个体素边长内分 64 步寻找符号包围，然后二分 40 次。
遇到其他占据体素则标记 blocked；无符号包围则标记 unbracketed。
这只是第一个**采样到的**外向过零包围，不声称寻找全部零点，也不是全局 Hausdorff 距离。

输出文件：

| 文件 | 含义 |
| --- | --- |
| `smooth_distance.json` | 原树文件引用、beta、误差界、符号约定及非精确 SDF 声明 |
| `distance_validation.json` | 数值/梯度检查及三类距离误差、射线诊断 |
| `distance_samples.csv` | 每个采样点的硬代理、平滑代理、真实距离基准和解析梯度 |
| `benchmark.json` | 原指标加编译、平滑验证、值/梯度/参考距离查询耗时 |

`distance_validation.json` 的 `passed=true` 仅指指定数值检查通过，**不是距离精度达标**。
原树 JSON、边界网格、严格符号诊断及 `strict_free_space.json` 不变，不能把平滑函数
替换进旧 manifest 却继续声称严格原始零边界。

## 文件改动

- `include/rokae_demo/halfspace_smooth_tree.hpp`：LSE DAG 执行、误差界、解析梯度。
- `include/rokae_demo/halfspace_compiled_tree.hpp`：仅增加执行适配类的 friend 访问。
- `include/rokae_demo/boundary_distance_reference.hpp`：有限边界距离基准。
- `include/rokae_demo/boundary_smooth_distance.hpp` 与 `src/boundary_smooth_distance.cpp`：距离审计与导出。
- `src/boundary_direct_pipeline.cpp`：可选 CLI、编译、验证及独立查询计时。
- `tests/*smooth*`、`CMakeLists.txt`：单元、CLI、benchmark 回归测试。
- `scripts/benchmark_smooth_distance.py`：同树、同查询点、不同平滑预算对比。

## 已知距离限制

硬树零边界正确，并不意味着硬树数值是欧氏距离。
例如 `[0,0.12]^3 m` 立方体外的点 `(-0.01,-0.01,0.06) m`：
硬树代理距离为 `0.01 m`，真实距离为 `sqrt(2)*0.01 m`。
令 beta 趋于无穷只能逼近硬代理，不能消除这类原始距离畸变。
本阶段不另加距离重初始化、BSP、R-function、Nef 或全局 Boolean 优化来掩盖该限制。

## 实测结果

完整可复现结果由下列命令生成，构建与完整验证分列，IO 和额外查询 benchmark 不计入核心构建：

```bash
python3 scripts/benchmark_smooth_distance.py output/smooth_distance_benchmark
```

本机 Release，2026-09-09，h=0.04 m，9 个场景 × 3 档预算：
solid cuboid、L、U、stair、closed cavity、narrow corridor、disconnected、thin wall、camera。
每档 1 次预热 + 3 次构建取中位数；另一个进程测查询，7 组各 5 轮取中位数。
查询点完全相同；计时为本机样本，不是实时性保证，各档没有承诺按 beta 单调变化。

完整 [27 组测量表](../output/smooth_distance_benchmark/report.md)、
[原始统计](../output/smooth_distance_benchmark/summary.json) 和
[可执行文件及输入指纹](../output/smooth_distance_benchmark/manifest.json) 已保存。

相机：1655 点、341 occupied voxels、318 finite patches、524 DAG nodes、61 个方向半空间。
所有预算都沿用同一份树。查询访问 524 个存储节点而非展开后的 51402 个节点引用。

| 指标 | 1 mm 预算 | 5 mm 预算 | 10 mm 预算 |
| --- | ---: | ---: | ---: |
| beta（1/m） | 19879.3 | 3975.9 | 1987.9 |
| 实测最大平滑误差（mm） | 0.790 | 3.948 | 7.895 |
| 近边界真实距离 RMSE（mm） | 0.863 | 2.407 | 4.740 |
| 近边界真实距离最大误差（mm） | 9.661 | 12.064 | 15.125 |
| 全体采样真实距离 RMSE（mm） | 18.041 | 18.552 | 19.430 |
| 全体采样真实距离最大误差（mm） | 171.984 | 174.493 | 177.433 |
| 64 条射线最大外移（mm） | 0.743 | 3.717 | 7.434 |
| LSE 值查询（μs/点） | 14.850 | 20.588 | 18.764 |
| LSE 值+梯度（μs/点） | 16.711 | 24.165 | 23.513 |
| 核心构建，含执行编译（ms） | 1.106 | 1.209 | 1.291 |
| 构建+完整验证（ms） | 334.748 | 380.320 | 403.624 |

每档 6705 个距离探针，其中 5881 个满足 `|d_ref|<=0.04 m`。
原硬代理在同一批探针上已经有 **171.356 mm 最大误差、17.942 mm RMSE**。
以 1 mm 档最坏点 `(0.448223,0.410329,1.065369) m` 为例：
真实距离 `0.461685 m`、硬代理 `0.290329 m`、平滑代理 `0.289701 m`。
可见进一步增大 beta 也无法弥补原硬树的距离畸变。

查询对照：相机原硬 DAG 为 **1.595 μs/点**，1 mm 档有限面片参考距离为
**3.184 μs/点**；本场景 LSE 求值并不比直接参考距离更快。
参考距离只给值、不是到处可微的函数；LSE 的功能是提供平滑代理及解析梯度，
不能据此推断它在距离精度或计算速度上优于有限表面距离。

1 mm 档额外 LSE 编译/工作区准备约 `0.138 ms`；完整距离验证约 `306.387 ms`。
核心时间**不含**完整验证；当前 CLI 确实执行完整验证，不能将其每次运行称为 1.1 ms。
执行程序估算 18376 B；标量工作区 4192 B，梯度工作区 16768 B。
1 mm 档最小梯度范数约 0.0153，近边界有 35 个采样点范数小于 0.1，
不能假定该代理满足单位梯度或据此直接保证控制器可行性。

## 验证结果

- 29 项 CTest 全部通过；7 项核心 ASan/UBSan 测试全部通过（ASan leak detection 关闭）。
- 穷举 2×2×2 的全部 255 个非空占据模式，覆盖面/边/顶点接触等情况。
- 27 组 benchmark 的错误自由分类、原始点正距离、上界次序违反、累计误差界违反均为零。
- 所有组解析梯度有限差分检查通过，最大绝对差约 `1.283e-7`。
- 所有组所请求射线均找到过零包围，blocked/unbracketed 均为零；不推广为未采样位置的保证。
- 同一数据的 exact/LSE 树、网格、符号诊断、严格自由空间 manifest 逐字节一致。
  相机输出也与上一轮 DAG benchmark 对应文件逐字节一致。
- 硬树 artificial zero faces/edges/vertices 仍为零；这是原硬树证书，**不是平滑零面的证书**。
- CLI 回归从导出的 JSON 独立重算有限矩形距离、LSE 值、逐节点权重梯度及 DAG 路径误差界，
  逐点对照 CSV，并用外部求值器再做有限差分；不只依赖 C++ 内部验证。
  覆盖共享面、edge-touch、vertex-touch、L/U、腔体、负坐标和薄壁，
  检查绝对 beta/误差预算两种参数及工作量限额拒绝。
  显式关闭附加梯度/射线诊断后，原硬树证书和距离探针检查仍必须存在。
- 对相机 1 mm 档已有导出包的 6705 个探针另做 Python 独立复算：
  值与 C++ CSV 的最大差 `8.674e-19 m`，梯度分量最大差 `3.442e-15`；
  DAG 路径误差界也一致。这是实现一致性误差，不是相对真实边界的距离误差。

结论：函数值、解析梯度、单侧保守构造、误差分解与 benchmark 已完成。
**单纯 LSE 替换不能把现有 MIN/MAX 树变成全空间高精度欧氏 SDF**。
若后续目标是例如“距离真实边界 1 mm 精度”，需要另行约定有效查询带和精度指标，
再评估距离建模方案；本次未自行扩大算法范围。
