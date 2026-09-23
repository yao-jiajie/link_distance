# 执行层相同子表达式共享

## 本轮范围

这是“点云到函数”优化清单的第一项，不是重构几何前端或进一步平滑。
第二项 [二元 LSE 内核](lse_binary_kernel_results.md) 已另行完成；本页记录仍为
通用内核的共享前后比较。对应 benchmark 默认固定 `--lse-kernel generic`。
保持 raw point cloud、occupied voxels、暴露面、patches、支撑平面、canonical cells、
Phi、原 TreeRepresentation、MIN/MAX/LEAF 和 LSE 公式全部不变。
不删验证，不重排子引用，不调用简化器，不开启 fast-math。

原 DAG 仅通过相同指针共享计算，不同指针的同构表达式仍重复执行。
新增执行层结构匹配：叶节点按同一 leaf ID；内部节点按
`operator + 有序的已映射 child IDs`。只复用完全相同表达式的执行结果。

相机树源节点仍为 524；结构分析与实际编译均得到 429 个执行节点。
每次 LSE 查询的指数调用由 543 变成 439，对数调用由 463 变成 368。
改变的是执行工作量，不是输出函数或真实距离误差。

## 实现细节与边界

`CompiledTree` 先完整检查原始输入的节点编号、算子、叶索引、有限平面系数、
拓扑顺序、根、环/前向引用和可达性，再开始共享。不能靠合并隐藏非法节点。

编译遍历保持 postorder，记录原节点到执行节点的映射。
哈希仅用于选择候选桶，之后完整比较算子、leaf ID、子引用数量和顺序；
哈希碰撞不会导致错误合并。已有指令和边缓冲区原地压紧，保留容量，避免额外缩容复制。
源树及其序列化结果不修改；编译临时哈希表在查询前释放。

以下行为明确禁止：

- 把 `MAX(A,A)` 变成 `A`：LSE 有 `log(2)/beta` 偏移，重复出现必须保留。
- 把不同顺序的 children 合并或排序：不改变浮点计算顺序、并列选择和 signed zero。
- 根据近似相等的平面系数合并不同 leaf ID。
- 应用硬树吸收律、幂等律、任意展平或重新结合 LSE。

相同算子按同一顺序接收逐位相同的子值及子梯度，再执行同一算式，因此结果可复用。
误差权重 `W`、自动 beta、理论误差界和 Lipschitz 界保持不变。
审计的舍入容差继续按原始节点数计算，不因执行节点减少而改变验收标准。
操作预算使用实际执行节点/边计数；原构树分配上限、展开上限和全部符号证书上限仍保留。

查询循环未加入哈希、递归或堆分配。线程共享不可变执行程序，但各自拥有工作区。
内存说明：标量/梯度工作区随执行节点数减少，执行程序容器容量不必减少；
本次不声称编译峰值内存降低。

## 接口

Direct CLI 新参数：

```text
--execution-sharing structural   # direct 默认：开启结构共享
--execution-sharing none         # 原 pointer-DAG 执行程序，对照模式
```

同时用于所请求的硬树查询、LSE 查询和距离审计；三者复用一次编译的程序。
`--boundary-query recursive` 的递归求值本身不受影响；若同时请求编译审计或 LSE，
这些编译程序仍使用该开关。flat/local/volume/Alpha-tet 默认路径不变。
`--distance-field lse` 仍需显式启用；本次没有把平滑函数设为默认输出。

C++ API 保持旧默认不变，显式开启如下：

```cpp
HS::CompiledTree hard(tree.raw, tree.leaves, true);
HS::SmoothTree smooth(tree.raw, tree.leaves, {0.0, 0.001, true});
auto hard_work = hard.makeWorkspace();
auto smooth_work = smooth.makeGradientWorkspace();
const auto sample = smooth.evaluateWithGradient(point, smooth_work);
```

执行指标新增：

| 指标 | 含义 |
| --- | --- |
| `execution_structural_sharing` | 请求的结构共享开关 |
| `compiled_shared_nodes` | 硬树执行程序共享掉的节点；无编译程序时为 0 |
| `smooth_execution_nodes` / `smooth_execution_edges` | LSE 实际执行节点/引用数 |
| `smooth_shared_nodes` | LSE 源节点数减执行节点数 |
| `smooth_exp_calls_per_query` / `smooth_log_calls_per_query` | 当前内核一次成功完整求值的指数/对数调用次数 |
| `smooth_reuses_compiled_program` | LSE 与硬树/审计是否持有同一不可变程序实例 |

原 `tree_nodes` 等仍表示导出的源树规模，`query_node_visits` 表示所选硬求值器实际访问数。
原函数定义 JSON 和严格自由空间 manifest 无需迁移；旧消费者仍可用相同树复现函数，
但只有启用新编译选项才能享受执行共享。

## 测试与 benchmark

新增成对测试涵盖不同指针的同构子表达式、不同 child 顺序、不同算子/元数、
重复项、并列、signed zero、极端有限值、非法 DAG、共享程序的多线程工作区。
函数值、梯度逐分量按 double 位模式对照；自动 beta 和误差界不变。
全部 255 个非空 2×2×2 占据模式也比较共享前后的值和 signed zero。

新的成对 benchmark 使用同一个可执行文件，none/structural 顺序逐轮轮换。
核心时间包含执行共享的编译成本；完整验证不省略；I/O 和额外查询计时独立。
每轮比较树、几何、函数参数、逐点值/梯度 CSV 的字节一致性；距离审计除
`audit_ms` 和 `work_used` 外完全一致。单独查询进程使用相同查询坐标。

```bash
python3 scripts/benchmark_execution_sharing.py output/execution_sharing_comparison
# 默认 9 场景，1 mm 平滑误差预算；1 预热 + 3 构建，3 次独立查询进程。
# 更多预算：追加 --errors .001 .005 .01
```

历史 `benchmark_boundary_direct.py`、`benchmark_boundary_direct_query.py`、
`benchmark_smooth_distance.py` 已显式锁定 `none`，以保留各自历史实验口径。

## 修改文件

- `include/rokae_demo/halfspace_compiled_tree.hpp`：验证后执行层结构共享与统计。
- `include/rokae_demo/halfspace_smooth_tree.hpp`：共享选项、执行计数；LSE 内核不变。
- `include/rokae_demo/boundary_smooth_distance.hpp`、`src/boundary_smooth_distance.cpp`：审计使用同一开关。
- `src/boundary_direct_pipeline.cpp`：可选参数和执行指标。
- 对应 C++/CLI 测试、`tests/execution_sharing_benchmark_test.py`、`CMakeLists.txt`。
- `scripts/benchmark_execution_sharing.py`、历史 benchmark 的显式对照选项及文档。

## 实测结果

2026-09-09，本机 Release。9 场景 × 3 预算 × 2 模式，共 54 组记录；
每组 1 次预热 + 3 次构建取中位数、3 次独立查询进程取中位数。
查询进程内部仍使用 7 组 × 5 轮中位数。构建和查询都轮换 none/structural 顺序。

完整 [成对测量表](../output/execution_sharing_comparison/report.md)、
[统计](../output/execution_sharing_comparison/summary.json)、
[可执行文件及输入指纹](../output/execution_sharing_comparison/manifest.json)。
本次命令追加了 `--errors .001 .005 .01`。

相机数据的同轮测量：

| 平滑预算 | 模式 | 核心构建 ms | 构建+完整验证 ms | 硬查询 μs | LSE 值 μs | 值+梯度 μs |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1 mm | none | 1.137 | 326.865 | 1.673 | 14.625 | 16.463 |
| 1 mm | structural | 1.278 | 285.197 | 1.299 | 12.066 | 13.644 |
| 5 mm | none | 1.192 | 389.701 | 1.672 | 18.396 | 21.071 |
| 5 mm | structural | 1.257 | 322.731 | 1.384 | 15.213 | 17.011 |
| 10 mm | none | 1.115 | 383.351 | 1.660 | 18.586 | 21.579 |
| 10 mm | structural | 1.315 | 322.343 | 1.377 | 15.369 | 17.737 |

1 mm 档：值查询减少约 **17.5%**，值+梯度减少约 **17.1%**，完整构建+验证减少约 **12.7%**。
核心构建没有变快，而是增加约 **0.141 ms**；硬/LSE 编译分别由 0.126/0.142 ms
变为 0.176/0.196 ms。额外结构识别成本计入核心时间，没有隐藏。
仅按这组核心及查询中位数估算，每帧约 51 次值+梯度查询即可抵消增加的核心构建成本；
这是摊销估算，不是实时期限或单次延迟保证。

相机梯度工作区由 16768 B 降到 13728 B；执行程序容量估计两种模式均为 18384 B。
另外 8 个简单场景的执行节点数均未减少，查询未显示稳定收益，还需承担结构匹配开销。
因此不能把相机场景的提速百分比推广到任意点云。

验证结果：

- 30 项 CTest 全部通过，包括相机端到端成对 benchmark。
- 7 项核心 ASan/UBSan 全部通过；ASan leak detection 关闭。
- 54 组记录全部通过完整验证，27 对输出树/几何/函数参数/值与梯度 CSV 完全一致。
- 1 mm 相机优化结果与上一版 `output/smooth_distance_benchmark/camera/error_0.001`
  的原树、边界 OFF、严格自由空间 manifest、函数定义及逐点 CSV 也逐字节一致。
- 数值误差、符号、有限差分和射线诊断保持不变；仅工作量计数与耗时变化。

本阶段完成；后续二元 LSE 专用内核见上方独立报告。
本轮后续收尾已统一 direct LSE 路径的执行程序：流水线在平滑模式下只构造一份
`CompiledTree`，`SmoothTree` 持有其只读共享指针，距离审计直接复用同一份指令、
子节点索引和系数快照；硬树工作区仍独立分配。这样不改变任何树/几何/函数值，
并去掉审计阶段再次序列化、结构匹配和复制平面的开销。旧的 `SmoothTree(tree, leaves, options)`
构造函数仍保留，其他调用路径行为不变。
