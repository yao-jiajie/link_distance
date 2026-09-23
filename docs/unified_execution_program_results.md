# 统一硬树、LSE 与距离审计执行程序

## 范围

这是点云到函数优化的第三项。boundary-direct 的 LSE 路径现在每个几何快照只构造
一份不可变 `CompiledTree`，硬 MIN/MAX 查询、`SmoothTree` 和距离审计共享其指令、
有序子引用及平面系数。各次独立调用和线程仍须各自持有工作区。

原表达树、结构共享规则、LSE 公式、beta、误差界、二元内核、查询运算顺序和全部
验证范围均未改变。导出的树、几何、函数参数、逐点值和梯度保持一致。

## 实现

`SmoothTree` 内部由按值保存 `CompiledTree` 改为持有
`std::shared_ptr<const CompiledTree>`。现有从表达树或序列化树构造的接口仍保留，
它们会自行创建只读程序；流水线使用新增的共享程序构造接口。

预编译程序记录是否启用了结构共享。共享构造时，`SmoothTreeOptions` 的配置必须与
程序匹配，否则立即拒绝，避免静默采用不同的执行规模。空共享指针同样拒绝。
距离审计不再接收另一份 `TreeRepresentation` 并重新编译；它明确审计
`SmoothTree::compiledProgram()` 所代表的同一几何快照。
旧的带 `TreeRepresentation` 参数重载仍保留原行为：独立编译传入的树用于核对，
不会忽略它。只有不带该参数的新接口复用程序。

新增运行指标：

| 指标 | 含义 |
| --- | --- |
| `execution_program_instances` | 当前流水线持有的编译程序实例数；LSE 路径为 1 |
| `smooth_reuses_compiled_program` | `SmoothTree` 是否复用流水线的同一实例 |
| `evaluator_compile_ms` | 唯一执行程序编译，以及所请求的硬查询工作区准备时间 |
| `smooth_compile_ms` | LSE 参数/误差权重和工作区准备；不再包含第二次程序编译 |

## 实测

2026-09-09，本机 Release，相机数据、1 mm 平滑预算、结构共享开启。
旧记录来自二元内核阶段的 3 次中位数；新记录为 1 次预热、5 次构建中位数。
两次测量并非同一进程内的严格成对试验，因此只把直接对应的准备阶段作为证据，
不根据易受负载影响的完整验证时间推断收益。

| 指标 | 独立硬/LSE 程序 | 共享程序 |
| --- | ---: | ---: |
| 唯一程序编译 `evaluator_compile_ms` | 0.171 ms | 0.203 ms |
| LSE 准备 `smooth_compile_ms` | 0.188 ms | 0.020 ms |
| 上述两段合计 | 0.359 ms | 0.223 ms |
| LSE 准备降幅 | - | 89.4% |
| 同时存在的执行程序实例峰值 | 3 | 1 |

旧路径同时保留硬查询程序和 LSE 程序，距离审计期间再创建一份硬程序；新路径只有
一份。按相机记录中每份约 18.4 KB 的程序容量估算，执行程序存储减少约 36.7 KB。
标量和梯度工作区未合并，仍按调用者及线程分别持有。

最终采样记录位于
[output/unified_program_reuse_benchmark](../output/unified_program_reuse_benchmark/report.md)。
其中 generic/binary 两种内核的源节点均为 524、执行节点均为 429，所有原树、几何、
函数值、梯度和除耗时外的距离审计结果一致；与旧版相机 binary 输出也做了同样比较。
`smooth_program_bytes` 包含共享程序，不能再与 `compiled_program_bytes` 相加；
这些容量估算均不包含分配器和共享指针控制块开销。

## 验证

- Release 全工程构建无警告，31 项 CTest 全部通过。
- 7 项核心测试用 ASan/UBSan 重新编译并通过；ASan leak detection 关闭。
- 单元测试验证共享对象身份、复制/生命周期、空程序拒绝、共享配置不匹配拒绝、
  程序快照，以及旧审计接口继续检查传入树的行为。
- 端到端测试要求 `execution_program_instances == smooth_reuses_compiled_program == 1`。
- 相机最终运行保留完整 all-strata 符号证书、距离探针、有限差分和零面射线诊断。

后续第五项已加入[多帧几何、程序与工作区复用](multi_frame_reuse_results.md)。只有
体素尺寸、占据键和全部相关配置精确匹配时才复用；几何、叶系数或构建参数变化仍会
创建新程序并重新认证，不能用旧快照替代新场景的构建和验证。
