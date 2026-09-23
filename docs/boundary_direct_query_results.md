# Boundary-direct：保持几何不变的 DAG 查询优化

## 本轮改变什么

保持 occupied voxels、真实 boundary planes、canonical cells、切分规则、Φ、常量折叠及所有 MIN/MAX 节点不变。
不做树化简、全局 Boolean 优化、BSP-field、R-function 或任何几何近似。

旧 `HS::evaluate` 递归访问每个引用，遇到共享子树会重复求值。相机树虽然只有 524 个存储节点，
一次查询却访问 51402 次节点。新 `HS::CompiledTree` 按已有序列化 DAG 的后序顺序，
把同一共享节点在**本次查询中只计算一次**。下一查询会覆盖全部工作区，不跨查询缓存结果。

实现文件：[halfspace_compiled_tree.hpp](../include/rokae_demo/halfspace_compiled_tree.hpp)。

后续已添加可选的 [LSE 距离代理与解析梯度](smooth_distance_results.md)。
下文仍是原硬 DAG 的查询性能；不能把它当作含指数/对数的平滑函数耗时。
后续 direct 默认还启用了 [执行层结构共享](execution_sharing_results.md)。
本页历史 benchmark 显式指定 `--execution-sharing none`，不混入新优化的效果。

- 编译时复制平面系数，准备连续指令/子节点索引；节点类型依旧只有 LEAF/MIN/MAX。
- 检查节点 ID、postorder、根、叶子索引、非法 operator、空分支、环/前向引用及不可达节点。
- 查询不递归、不查哈希表、不分配内存；工作量为 `O(stored nodes + child references)`。
- 保持原来的叶子算式、MIN/MAX 孩子顺序及初始正/负无穷值，不重排运算，保留 signed-zero 行为。
- 执行程序是树/系数快照；几何帧变化后必须重新编译。程序可共享，每个线程需独立 workspace。
- 只支持有限查询坐标和有限平面系数；对有限输入中间运算溢出的行为有回归，但不把浮点溢出解释为几何保证。

另外，direct 的 raw 与 simplified 本来就指向同一表达式。符号验证现在仅在这两个指针相等时
复用该 DAG 的验证值，而不是序列化、计算两遍。遍历的全部体/面/边/顶点 strata 不变；
若 raw/simplified 不同，仍然分别检查两者。修改位于 `src/boundary_local_tree.cpp`，不是删除验证。

## 接口和兼容性

Direct CLI 新参数：

| 参数 | 默认 | 用途 |
| --- | --- | --- |
| `--boundary-query dag` | 是 | 编译后求值 |
| `--boundary-query recursive` | 否 | 原递归求值器，对照使用 |
| `--boundary-query-cross-check true` | 否 | 每个世界坐标验证点额外比较两种求值器的完整标量与 signed zero |

默认 DAG 模式继续执行原始点、体素中心/角点、patch 中心及两侧、随机点验证，完整符号证书不能关闭。
cross-check 是额外审计，不是几何验证开关；开启后会额外支付递归求值成本。
`--boundary-max-probe-ops` 和 `--boundary-max-query-ops` 按所选 evaluator 的工作估计限额；
cross-check 验证预算同时计入两种 evaluator。原构建的展开引用上限仍保留，未绕过。

```bash
build/build_halfspace_tree data/realistic_sparse_camera.xyz output/direct_dag_new \
  --pipeline octree --tree-source boundary --boundary-expression direct \
  --voxel-size .04 --boundary-query dag --query-benchmark true

# 额外逐点标量审计，成本不代表正常 DAG 查询
build/build_halfspace_tree data/realistic_sparse_camera.xyz output/direct_audit_new \
  --pipeline octree --tree-source boundary --boundary-expression direct \
  --voxel-size .04 --boundary-query-cross-check true
```

外部 C++ 使用方式（`r` 是既有 `TreeRepresentation`）：

```cpp
#include <rokae_demo/halfspace_compiled_tree.hpp>
namespace HS = rokae_demo::halfspace;

const HS::CompiledTree program(r.raw, r.leaves); // 每帧一次
auto workspace = program.makeWorkspace();      // 每个查询线程各一份
const double q = program.evaluate(point, workspace);
const bool strictly_free = q < 0;              // 仅用于已验证的 boundary-direct 树
```

也接受现有 `SerializedTree` 与 leaves。既有导出的 tree JSON 不变，`strict_free_space.json`
增加 `query_backend` 提示。**仅加载相同 JSON 不会使旧递归消费者自动加速**，消费者需要使用新执行接口。
现有通用 `HS::evaluate` 未修改，flat/local、volume/rect、Alpha-tet 的默认求值器均未替换。
顶层 CLI 默认 pipeline 也没有改变，只是已经选择 direct 后默认使用 DAG 查询。

## 性能口径

- `geometry_core_ms`：原有体素化、边界提取、canonical、构树。
- `evaluator_compile_ms`：DAG 编译和 workspace 准备；无编译的递归模式为 0。
- `core_runtime_ms`：包含上述编译/准备，不隐藏预处理成本。
- `query_node_visits`：本次实际访问节点数；`expanded_node_references` 仍报告原树展开规模。
- `query_work_estimate`：DAG 为存储节点+边，递归为展开节点，供预算保护，不是精确 CPU 指令数。
- `compiled_program_bytes` / `query_workspace_bytes`：拥有的容器容量及对象估算，不是 RSS 或峰值内存。
- `total_ms`：核心构建+全部强制验证，排除 I/O 和可选 query benchmark。

新 benchmark：`scripts/benchmark_boundary_direct_query.py`。
同一可执行文件中 serial rotating triples 对比 recursive/direct、DAG/direct 和旧 flat。
两个 direct 模式的 tree、coefficients、split、occupancy、boundary、符号证书逐字节比较；
验证点数相等；三个模式的计时随机点和边界点数组相等。
每个场景另跑一次不计入测量记录的 scalar/signed-zero 审计。
旧 `benchmark_boundary_direct.py` 已显式固定 `--boundary-query recursive`，避免历史 benchmark 悄悄换求值器。

```bash
python3 scripts/benchmark_boundary_direct_query.py output/direct_query_comparison_new \
  --repeats 5 --warmups 1 --query-repeats 3
```

## 测试内容

- 新增 `halfspace_compiled_tree_test.cpp`：普通和高度共享 DAG 的标量一致性、signed zero、极小/极大有限输入、
  非法序列化输入、错误工作区、系数快照、多线程独立 workspace、指数展开但线性存储的 49-node DAG。
- 扩充 direct 原生测试：全部 255 种非空 2×2×2 占据模式的所有 strata 代表点对比递归与 DAG，
  并覆盖原有负坐标、0.03/0.04 m、随机体素、角点 nextafter 邻点。
- 扩充 CLI 测试：默认 DAG、recursive 审计、双向 cross-check、原样导出和参数拒绝。
- 新增 `boundary_direct_query_benchmark_test.py`：三模式对照、标量审计、时间口径和输出一致性。

本轮只优化求值执行和同指针验证复用；不宣称逻辑树规模下降、全流程达到 1 ms、获得距离场或 CBF 控制保证。

## 实测结果（2026-09-09）

完整 [配对报告](../output/boundary_direct_query_comparison/report.md)、
[逐场景统计](../output/boundary_direct_query_comparison/summary.json)、
[运行方法及输入/二进制指纹](../output/boundary_direct_query_comparison/manifest.json)。
Release，同一二进制，1 次预热 + 5 次测量，query 做 3 个配对进程，各 7×5 轮中位数。
正式性能运行期间未并行构建或运行回归测试。

### 相机样例

1655 原始点，341 occupied voxels，`voxel_size=0.04 m`；
三个模式使用相同输入和计时点集。

| 指标 | direct / recursive | direct / DAG | flat / 原求值器 |
| --- | ---: | ---: | ---: |
| 存储节点 | 524 | 524 | 721 |
| 每点访问节点 | 51402 | 524 | 721 |
| 单点查询 µs | 217.761 | **1.696** | 3.875 |
| 边界点查询 µs | 208.655 | **1.691** | 4.010 |
| 核心构建 ms | 0.8923 | 1.0687 | 1.4468 |
| 其中执行程序准备 ms | 0 | 0.1316 | 0 |
| 构建+完整验证 ms | 1600.100 | **33.097** | 273.457 |
| 人工零面/边/顶点 | 0/0/0 | 0/0/0 | 1181/1725/645 |
| 严格 q=0 边界证书 | 通过 | 通过 | 不通过，依赖零值 guard |

相同 direct 树的查询约快 **128 倍**，构建+验证约快 **48 倍**。
新增执行快照准备有约 **0.13 ms** 成本，已计入核心构建，不应忽略。
程序拥有存储约 18,344 B（17.9 KiB），每个线程工作区 4,192 B（4.1 KiB）；这不是进程峰值。

树几何没有缩小：展开 plane refs 仍为 27,408，61 个有向半空间、524 个节点全部保留。
每点不再递归重复访问，才是查询加速来源。两个 direct 模式均使用本轮同指针符号验证复用，
所以表中的 direct 配对主要隔离 DAG 数值求值收益；第一版历史 205 µs/1520 ms 来自另一轮运行，未冒充本次对照。

### 各场景查询与验证

| 场景 | recursive µs | DAG µs | DAG 核心 ms | DAG 构建+验证 ms |
| --- | ---: | ---: | ---: | ---: |
| cuboid | 0.050 | 0.031 | 0.0866 | 0.512 |
| L-shape | 0.063 | 0.036 | 0.0563 | 0.350 |
| U-shape | 0.081 | 0.047 | 0.0770 | 0.456 |
| stair | 0.080 | 0.045 | 0.0753 | 0.450 |
| cavity | 0.112 | 0.061 | 0.0842 | 0.598 |
| disconnected | 0.143 | 0.064 | 0.0763 | 0.516 |
| aligned gap | 0.062 | 0.035 | 0.0507 | 0.298 |
| edge-touch | 0.094 | 0.050 | 0.0641 | 0.370 |
| vertex-touch | 0.159 | 0.073 | 0.0638 | 0.429 |
| camera | 217.761 | 1.696 | 1.0687 | 33.097 |

所有 10 个场景的两个 direct 模式均通过严格符号验证，人工零面/边/顶点全部 0，raw points 无漏包。
相机样例证书覆盖 12,673 个 strata，额外标量审计逐点核对 **6,678** 个世界坐标点，数值及 signed-zero 不一致数为 **0**。
direct 原始/简化树、切分、半空间系数、边界网格和符号证书均与旧版本相机产物逐字节核对一致。

最终 Release **25/25** CTest 通过；ASan+UBSan 五项核心测试 **5/5** 通过（`detect_leaks=0`，不宣称泄漏或 TSAN 检查）。
没有关闭几何验证换取计时结果。

### 当前剩余瓶颈

相机样例强制验证中，几何/符号证书中位耗时约 **19.7 ms**，世界坐标 probes 约 **12.3 ms**；
正常查询不再是之前的主要瓶颈。完整构建+验证仍约 **33 ms**（且不含 I/O），
不能以约 1.07 ms 的核心构建时间宣称端到端 1 ms。本轮到此止步，未继续改动验证算法或扩展机器人控制层。
