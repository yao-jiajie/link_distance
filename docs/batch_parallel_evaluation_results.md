# 批量查询与并行评估

## 范围

这是点云到函数优化的第六项。`CompiledTree`、`SmoothTree` 新增保持输入顺序的批量
接口，把不同查询点分配给独立 worker；每个点内部仍调用原有的硬 MIN/MAX、二元或通用
LSE 以及解析梯度内核。本项没有改变树、系数、结构共享、LSE 参数、误差界或距离含义，
也没有使用近似数学函数。

本项选择点级并行，没有加入 SIMD。这样每点的叶算式、子项顺序、并列选择、指数/对数
调用和 signed-zero 行为不变，可直接与原标量结果逐位核对。

## API 与执行模型

公共类型 `EvaluationBatchWorkspace<T>` 为每个 worker 保存一份节点工作区。它只能通过
对应程序创建、可移动但不可复制，避免无意复制大型缓冲。三个入口为：

```cpp
auto hard_work = program.makeBatchWorkspace(workers);
std::vector<double> hard(points.size());
program.evaluateBatch(points, hard, hard_work);

auto value_work = smooth.makeValueBatchWorkspace(workers);
auto gradient_work = smooth.makeGradientBatchWorkspace(workers);
std::vector<double> values(points.size());
std::vector<HS::ValueGradient> gradients(points.size());
smooth.evaluateBatch(points, values, value_work);
smooth.evaluateWithGradientBatch(points, gradients, gradient_work);
```

`evaluateWithGradientTrustedBatch()` 供已经在批次边界检查所有输入、并会检查所有输出的
核心链路使用；它省去逐 DAG 节点的重复有限性/范围检查，但保持相同的 LSE 运算顺序。
当 `beta*gap > 750` 时，binary64 的 `exp(-beta*gap)` 必为精确 `+0`，该入口直接使用
零并跳过相应的 `exp`；零和 MAX 修正量以及二元 MIN 修正量同样复用逐位相同的结果。
一般调用方应继续使用上面的 checked 接口。

输出必须由调用者预先调整到与输入相同的大小。worker 数在创建工作区时显式给定；值为
1 时走相同批量接口但不创建子线程，可用于隔离 API 开销。大于 1 时按连续、数量相差至多
一个点的区间分块，主线程处理第一个区间，其余区间由工作区持有的常驻线程处理。输出索引因此
始终与输入一致，各 worker 只写自己的区间和工作区。

单个工作区只能顺序复用；同时执行的多个批次必须分别持有工作区和输出数组。工作区
只保存会被覆盖的临时缓冲，因此节点数相同的其他程序也可以复用它，节点数不匹配会拒绝。

线程池和异常槽在工作区创建时建立，后续批次只分派任务；工作区销毁时才回收线程。API 不会
自动使用 `hardware_concurrency()`，避免库内部擅自占满调用方的线程资源。多个外部任务并发
时，调用方应自行减少 worker 数。若任一查询抛出异常，所有已启动线程会先 join，再向调用方
重抛；其他区间可能已经写入，因此错误批次的输出应视为部分结果。

工作区按 worker 数线性增长。相机 429 节点程序的实测容器容量估算为：

| workers | 硬值/LSE 值工作区 | 值梯度工作区 |
| ---: | ---: | ---: |
| 2 | 6,944 B | 27,536 B |
| 4 | 13,856 B | 55,040 B |
| 8 | 27,680 B | 110,048 B |

## 同进程实测

2026-09-09，本机 Release、8 个物理核，相机 429 节点结构共享程序，`h=0.04 m`、
1 mm LSE 标量误差预算、二元内核。每个配置先预热，再对 scalar loop、1-worker batch
和 parallel batch 做 7 轮轮换顺序测量；表中为每点微秒中位数，线程创建包含在 parallel
列中。完整数据、协议和二进制哈希见
[最终报告](../output/halfspace_batch_benchmark_final/report.md)。

| 点数 | workers | 硬值 scalar / parallel us | 加速 | LSE 值 scalar / parallel us | 加速 | 值梯度 scalar / parallel us | 加速 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 256 | 8 | 1.252 / 0.700 | 1.79x | 11.730 / 2.298 | 5.10x | 13.699 / 3.845 | 3.56x |
| 4,096 | 8 | 1.388 / 0.354 | 3.92x | 12.873 / 3.216 | 4.00x | 16.300 / 4.037 | 4.04x |
| 16,384 | 2 | 1.154 / 0.566 | 2.04x | 10.012 / 5.310 | 1.89x | 11.541 / 5.947 | 1.94x |
| 16,384 | 4 | 1.303 / 0.318 | 4.10x | 11.845 / 3.039 | 3.90x | 12.953 / 3.464 | 3.74x |
| 16,384 | 8 | 1.539 / 0.355 | 4.34x | 13.461 / 1.788 | 7.53x | 15.470 / 1.961 | 7.89x |

`1-worker batch` 与手写 scalar loop 在多数配置接近，说明批量接口本身没有固定的大额开销；
具体数值见完整报告。256 点硬值只有 `1.79x`，线程启动与调度已经限制收益。不同进程的
CPU 频率和系统负载也会改变加速比，因此表格是本机历史结果，不是当前常驻线程池或 worker
数对应的吞吐保证。

## 验证

单元测试对 `1/4` workers 的硬值、LSE 值和值梯度批次逐项比较原标量接口，覆盖有序随机
点、signed zero、最小非正规数、空批次、重复使用、错误输出大小、其他程序的工作区、
move 后空工作区、非有限点和 worker 数零。所有 `double` 以及梯度三个分量均按位一致。
另覆盖点数少于 worker 数、worker 异常后恢复，以及关闭结构共享时的 generic/binary 内核。

同进程 benchmark 在计时前及每轮三模式执行后，对相机查询集做硬值、LSE 值和全部梯度
分量的逐位比较，核对耗时不计入测量。JSON 保存各模式逐轮耗时，smoke test 检查报告
中位数与原始记录一致，同时检查 429 节点结构、工作区容量关系和所有计时字段为正。

Release 全部 35 项 CTest 通过。两项批量核心测试、相机批量 benchmark 以及 direct tree、
帧缓存、平滑距离三项回归另以 ASan/UBSan Debug 通过；关闭 leak detection，未做 TSan 检查。

本项不把完整验证并行化，也不引入全局线程池或自动线程数。线程池归调用方创建的工作区
独占；若调用方已有实时调度器，应避免与其过度订阅，并为并发批次分别创建工作区。
