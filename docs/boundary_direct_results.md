# Boundary-direct 原型

本文保留第一版**递归求值器**的实现与历史测量。后续已增加不改变表达树的 DAG 求值优化，
direct CLI 当前默认 `--boundary-query dag`；最新性能与接口见
[DAG 查询优化](boundary_direct_query_results.md)。下文的 205 µs 等数字不是当前默认查询性能。

## 实现范围

新增实验入口 `--pipeline octree --tree-source boundary --boundary-expression direct`。
这里的 `octree` 只是兼容现有命令行路由的名称；实际只建立 sorted occupied voxel keys，**不建立 Octree**。
默认仍为旧流程，不替换 flat/local、volume/rect、Alpha-tet。

```text
XYZ 点云（m）
  → 固定分辨率 occupied voxels
  → exposed faces / coplanar rectangular patches
  → 真实 boundary support plane registry
  → compressed canonical cells（包含所有无界外部区间）
  → 递归真实平面切分 + Φ + 构造期常量折叠
  → 原有 LEAF / MIN / MAX 表达式
  → 完整符号证书 + 世界坐标查询校验
```

没有 free rectangles、free prisms、bbox exterior free branches、flat 参考树或后续全局化简。
不调用 Nef、convex decomposition，不加入 R-function、smoothing、BSP-field 或全局 Boolean 优化。
这仍是用户指定的轴平面递归切分，并不意味着完全没有空间划分。

## 切分规则

每轴收集真实暴露边界平面的整数网格坐标 `c[0..m-1]`，去重排序。
canonical 区间是 `(-∞,c0),(c0,c1),…,(c_last,+∞)`。
组合后的每个开三维 cell 全 free 或全 occupied；依据原始 occupied keys 的计数与整数体积相等验证，不能把 mixed cell 当 occupied。
两端无界 cells 均为 free。这是全空间正确性的必要条件，但并没有额外构造六条 bbox 外部分支。

递归区域若全 free/occupied，返回构造期标签。mixed 区域只尝试其内部的真实 support plane。
用三维 prefix sum 计算子区间中 occupied canonical cell 数，最小化：

```text
(mixed_child_count, max(left_cell_count, right_cell_count), axis, cut_rank)
```

优先得到纯子区域，再平衡 **canonical cell 数**，不是物理体积。最后按 X/Y/Z 与坐标顺序打破平局。
这是确定性的初版启发式，不声称最优树。导出的 `direct_splits.json` 记录每次切分范围、support ID 和子区域占据计数。

## Φ 的实际实现

定义 `s=x_axis-c`，`a` 为负侧子表达，`b` 为正侧子表达。
两个孩子都不是常量时，直接创建三个 MIN 和一个 MAX，保持孩子的 shared_ptr 引用：

```text
MAX(MIN(a,b), MIN(a,-s), MIN(s,b))
```

FREE/OCCUPIED 是类型标签，代数意义分别为负/正无穷，**不是大数，也不是最终常量节点**。

| a | b | 构造结果 |
| --- | --- | --- |
| FREE | FREE | FREE |
| OCCUPIED | OCCUPIED | OCCUPIED |
| FREE | OCCUPIED | s |
| OCCUPIED | FREE | -s |
| FREE | b | MIN(s,b) |
| OCCUPIED | b | MAX(b,-s) |
| a | FREE | MIN(a,-s) |
| a | OCCUPIED | MAX(a,s) |

两个表达式指针相同也直接返回该表达式。正/负平面叶子按 `(support ID, sign)` 缓存，平面定义始终保留真实 patch 来源。
没有新增最终节点类型，也没有更改 TreeRepresentation 或原有查询器。
`tree.raw` 与 `tree.simplified` 指向同一个构造结果；后者只是兼容导出名称，`simplifier_applied=false`。

## 正确性与验证边界

当 `s<0` 时 Φ 的符号等于 a；`s>0` 时等于 b。
在 `s=0` 上：a、b 同负则 Φ 负，同正则 Φ 正，其余组合为零。
对子区域闭包递归应用，并考察一个 stratum 周围所有相邻开 cells：

- 全部 free：q<0，包括人工内部面、边、顶点。
- 全部 occupied：q>0，包括占据区域的人工内部接缝。
- free 与 occupied 并存：q=0，正好属于闭体素并集的真实拓扑边界。

验证不是只靠随机点：复用已有符号诊断内核，但剥离了对 flat 构造的依赖。
对真实平面 arrangement 的所有 open volumes、relative-open faces/edges、vertices（含无界部分），
传播离散 `-1/0/+1` 符号，同时核对原始体素占据、压缩 cell 均匀性、边界 patches 和叶子系数来源。
因此 `q=0 iff x∈∂V` 的证书针对这些真实轴平面定义的实数 MIN/MAX 场，不只是采样结论。

这里的 V 是选定分辨率的**闭占据体素并集**，不是未观测真实物体，不会将稀疏表面点自动填充为实体。
物理坐标采用项目一致的 `double(grid_index)*voxel_size` 边界。额外世界坐标测试覆盖 raw points、体素中心/角点、patch 中心及两侧、随机点；原生回归还包含角点的 `nextafter` 邻点与负坐标/非二进制分辨率。
不宣称任意非法输入、溢出浮点查询或改变系数后仍有证书；没有距离场、可微性或机器人 CBF 可行性保证。

严格模式查询只需 `HS::evaluate(root,leaves,p)<0`，没有零值 occupancy guard。
导出的共享 tree JSON 仍沿用通用 `inside iff root<=0` 元数据；此处描述的是非占据域的闭包，严格自由空间必须遵循 `strict_free_space.json` 的 `root_lt_zero`。
`occupancy_reference.json` 仅为调试资料，不是部署时查询依赖。

## 工作量限制与失败行为

| 参数 | 默认 | 范围 |
| --- | ---: | --- |
| `--boundary-max-cells` | 16384 | 包括无界区间的 canonical cells |
| `--boundary-max-direct-nodes` | 16384 | 构造期间分配的表达式节点 |
| `--boundary-max-expanded-nodes` | 262144 | 递归展开节点引用数 |
| `--boundary-max-split-checks` | 1048576 | 所有候选切分检查 |
| `--boundary-max-direct-depth` | 128 | 切分和表达式深度，最多 256 |
| `--boundary-max-strata` | 200000 | 完整符号验证 strata |
| `--boundary-max-sign-ops` | 200000000 | 符号证书工作量 |
| `--boundary-max-probe-ops` | 1000000000 | 原有递归查询器的世界坐标验证工作量 |
| `--boundary-max-query-ops` | 8000000000 | 可选 query benchmark 工作量 |

超限直接报错，不回退为 flat 或其他算法；不能把 skipped certificate 当成功。
符号证书失败会导出 `boundary_sign_failure.json`（含前 24 个见证、坐标表和 stratum 编码）；世界坐标失败输出首个错误点。
只有全部验证和导出成功后才写 `strict_free_space.json`。失败不会删除已有旧输出，务必检查退出码并为新数据使用新目录。
构建仍是批量原型，不是事务式在线发布服务。

## 修改文件

- 新增 `include/rokae_demo/boundary_direct_tree.hpp`、`src/boundary_direct_tree.cpp`：canonical cells、prefix sum、切分、Φ、验证。
- 新增 `src/boundary_direct_pipeline.cpp`：参数、原始点验证、独立导出、计时与查询 benchmark。
- 增量修改 `octree_boundary_tree.hpp/.cpp`：抽出几何平面 registry，旧 sweep 调用同一实现。
- 增量修改 `boundary_local_tree.hpp/.cpp`：抽出可独立使用的全 strata 符号证书；旧 flat/local 接口保留。
- `src/octree_boundary_tree_pipeline.cpp`：仅增加 direct 分发。
- `tests/boundary_direct_tree_test.cpp`、`tests/boundary_direct_pipeline_test.py`、`tests/boundary_direct_benchmark_test.py`：原生、独立导出及配对 benchmark 回归。
- `scripts/benchmark_boundary_direct.py`、`CMakeLists.txt`、`README.md` 和本文档。

## 复现

```bash
cmake --build build -j2
ctest --test-dir build --output-on-failure -j2
build/build_halfspace_tree data/realistic_sparse_camera.xyz output/direct_new \
  --pipeline octree --tree-source boundary --boundary-expression direct \
  --voxel-size 0.04 --validation-samples 1000 --query-benchmark true
python3 scripts/benchmark_boundary_direct.py output/direct_comparison_new \
  --repeats 5 --warmups 1 --query-repeats 3
```

`core_runtime_ms` 包含 voxel occupancy、边界提取/patch 合并、registry、canonical 和构树；
不含输入输出、所有验证、可选查询 benchmark。`total_ms` 包含构建及完整验证，不含 I/O/查询 benchmark。
配对结果使用相同 occupied geometry、随机查询集及边界查询集；flat 为 voxel occupancy + Y sweep，仍带零值 guard。
两条路径的验证参考和 probe 工作量不同，因此完整验证耗时单独列出，不把它归因于单一构造步骤。
查询继续使用原有递归 evaluator，未添加 DAG memoization；存储节点数与递归展开引用数必须同时报告。

## 本机实测结果（2026-09-09）

Release 构建，1 次预热 + 5 次正式配对运行取中位数，固定 `h=0.04 m`。
查询做 3 次独立配对进程，各进程 7×5 轮中位数，再取跨进程中位数。
没有同时运行构建/测试来干扰计时；这些仍是本机测量，不是最坏情况实时保证。
完整记录：[配对报告](../output/boundary_direct_comparison/report.md)、
[全部统计](../output/boundary_direct_comparison/summary.json)、
[输入/二进制指纹及方法](../output/boundary_direct_comparison/manifest.json)。

下表均为 direct；`nodes` 是唯一存储节点，`plane refs` 是展开叶子引用数。
全部场景在所有 strata 上都满足 free<0、occupied>0、true boundary=0，
人工 free 面/边/顶点零值计数 **全部 0/0/0**，occupied 内部错误零值也全部为 0。

| 场景 | 验证 strata | nodes | plane refs | 核心构建 ms | 单点查询 µs |
| --- | ---: | ---: | ---: | ---: | ---: |
| cuboid | 125 | 11 | 6 | 0.0612 | 0.038 |
| L-shape | 245 | 15 | 8 | 0.0486 | 0.057 |
| U-shape | 315 | 17 | 9 | 0.0655 | 0.068 |
| stair | 405 | 19 | 10 | 0.0521 | 0.070 |
| cavity | 729 | 23 | 12 | 0.0836 | 0.092 |
| disconnected obstacles | 567 | 26 | 18 | 0.0821 | 0.146 |
| aligned gap | 225 | 15 | 8 | 0.0418 | 0.053 |
| edge-touch | 245 | 20 | 12 | 0.0516 | 0.090 |
| vertex-touch | 343 | 24 | 16 | 0.0543 | 0.105 |
| realistic_sparse_camera | 12673 | 524 | 27408 | 0.9139 | 205.469 |

L-shape 的旧 flat 人工零面/边/顶点数为 `4/5/2`，direct 为 `0/0/0`。
独立导出测试也在旧 flat 的这些自由区零值位置求 direct 值，要求严格小于零；真实边界仍为零。

### 相机点云：正确性完成，查询性能明显退步

输入 1655 点，341 个 occupied voxels，782 个 exposed unit faces，318 个 boundary patches。
direct 建立 1800 个 canonical cells，304 次递归切分，2107 次候选检查；
80 次完整 Φ、224 次构造期折叠，切分深度 19，表达式深度 25，61 个唯一有向半空间。

| 指标 | flat（voxel + Y sweep） | direct |
| --- | ---: | ---: |
| 核心构建 ms | 1.4395 | 0.9139 |
| 存储 nodes | 721 | 524 |
| 递归展开 nodes | 721 | 51402 |
| 展开 plane refs | 618 | 27408 |
| 随机点 strict query µs | 3.589 | 205.469 |
| 边界点 strict query µs | 3.854 | 196.615 |
| 构建 + 完整验证 ms（不含 I/O） | 265.422 | 1520.324 |
| 人工零面/边/顶点 | 1181/1725/645 | 0/0/0 |
| root-only 严格符号证书 | 不通过，需 guard | 通过，不需 guard |

**不能用 0.914 ms 宣称全流程耗时，也不能用 524 nodes 宣称查询复杂度下降。**
原有 evaluator 是无缓存的递归遍历，Φ 多次引用相同子树；存储共享并不会自动消除重复求值。
相机点云单点查询约为 flat 的 **57 倍**；完整逐点验证也被重复递归访问拖慢。
本轮没有为了改善成绩改用缓存式 DAG evaluator、全局化简或新的几何算法。
因此结论是：在给定体素几何上实现了所要求的精确零值集合；尚不适合作为已验证的实时优化替代品。

### 回归记录

- Release 全部 **23/23** CTest 通过，包括旧后端、导出和 benchmark 回归。
- ASan + UBSan 的 direct、local、flat boundary、sparse occupancy 四项核心测试 **4/4** 通过（`detect_leaks=0`，未声称泄漏检查）。
- 原生测试穷举全部 **255** 种非空 2×2×2 占据模式；另测 12 个随机 4×4×4 模式、负坐标、0.03/0.04 m 分辨率、角点邻接浮点值及工作量上限。
- 独立 Python 导出测试包含 cuboid/L/U/stair/cavity/disconnected/gap/edge-touch/vertex-touch/negative，按原始占据盒和真实 patch 计算期望符号，不使用 direct 的 canonical 标签作 oracle。
- benchmark 的 10 个场景所有 direct 严格符号检查通过，原始点无漏包，边界 mesh 与 flat 字节一致；计时随机点和边界点数组一致。
- cuboid、L 和 camera 的旧 flat tree 导出与此前 voxel-removal benchmark 的对应文件逐字节核对一致。

没有遇到几何符号失败反例。这个结论有全 strata 符号证书支持，但不延伸为未观测真实障碍物保证或 CBF 控制保证。
