# Phase 4：多场景质量与速度评估协议

本阶段只添加 benchmark、数据生成器和只读查询适配器，不修改 Alpha-tet / Octree
构建算法、ConvexCluster、支撑面约简或 MIN/MAX 求值器。几何后端不调用新的凸分解。

## 数据与矩阵

默认 42 组输入/参数配置，每组比较 `alpha-tet` 与 `octree-rect`：

- 5 个合成场景：0.3 m 立方体表面、L 型、U 型、6 mm 薄壁、椭球表面。
- 每场景 256 / 1024 点。大点集的前缀是小点集，避免把不同采样位置误当成密度变化。
  这是点数两档，不宣称不同场景具有相同的单位面积点密度。
- 噪声每坐标独立均匀分布于 `[-a,a]`，a=0 / 0.003 m；最大位移为 sqrt(3)*a。
  无离群点，不把均匀扰动称为精确相机噪声模型。干净/带噪版本共享干净样本。
- 体素尺寸 0.02 / 0.04 m，即 5×2×2×2=40 组。
- 仓库中 1655 点多视角 RGB-D **模拟**样例，另测两种体素尺寸，共 42 组。
  理想椭球参数从现有 metadata 读取，不把它当成实测扫描。
- 可重复添加 `--measured-input scan.xyz` 接入用户实测数据，坐标必须先转为 m。
  没有提供真实几何时，仅评估原始点、占据体素与输出树，不伪造 ground truth。

合成盒组合仅采样并集的外表面，剔除盒间隐藏接口。盒面按面积加权；椭球为均匀球面
方向缩放，非均匀椭球表面积采样。表面、噪声、查询均使用固定且独立的种子。
输入清单保存生成参数、来源、理想几何、XYZ 内容哈希；不写入用户原始文件。

## 构建计时

每组每后端预热 1 次、正式测量 3 次，轮换后端顺序，串行运行，不并行编译/测试。
记录每次命令、日志、构建统计和进程 wall，输出 min / median / p95 / max。
仅 3 次测量时 p95 就是该组最大值，不能用它推断可靠尾延迟或硬实时上界。

Alpha-tet 使用 voxel corner 输入，alpha=3h、offset=h、exact plane reduction；
实际 alpha/offset、offset retries 也保留，不能把恢复后的 offset 当作请求值。
Octree 使用默认 rect，不开启 Phase 3 可视化边界阶段。

支持 `--include-baselines` 同时加入 `octree-none`、`octree-exact`、`octree-boundary`。
每个模式反复运行后检查凸块、raw/simplified 树哈希不变。Phase 3 的文件还需与 rect
逐字节一致。跨 pruning 只要求判定一致，不要求标量值相同。

- `core_runtime_ms`：后端自身的核心构建计时，不含验证及 I/O。
- `validation_ms`：各后端原有完整验证；原始点检查均保留，随机样本均为 1000。
  Octree 还检查所有体素角点/中心、未剪枝树，工作量与 Alpha-tet 不同。
- `total_ms`：构建＋后端验证，不含 I/O。
- `process_wall_ms`：后端进程启动至退出，包含其文件 I/O，不包含另行运行的查询工具。

不能把核心构建中位数称为端到端机器人避障延迟。

## 统一查询与几何评估

每组建立同一份米制查询文件，所有后端读同一文件：

1. 1024 个均匀随机查询。域包含 occupied union、理想物体以及实际导出的 wrap 顶点，
   外加一个体素边距；不假设 wrap 一定在某个未经证明的 offset 扩张范围内。
2. 全部原始点。
3. 每个 occupied voxel 的 8 个角点和中心，完整枚举，不因剪枝减少。
4. 对有理想几何的场景，独立生成 512 个理想表面点、512 个实体内部点。
5. L/U 型在预先定义的理想凹槽核心区域采样 512 点，用来观察槽被填充的比例。

`halfspace_tree_benchmark` 从现有 JSON 重建相同 Expression，然后调用已有的
`halfspace::evaluate`；不替换成 Python 求值器、AABB 查询器或另一种加速结构。
全部查询在 raw/simplified 树上的标量值和零阈值判定都应一致。

分类统一按 `F(x)<=0`，原始点正值会将该项标为失败；这额外检查了 Alpha-tet 原先
带数值容差的验证器以外，导出的零阈值树是否覆盖原始点。
Octree 还必须通过完整体素证书与查询一致性；Alpha-tet 的体素角点/中心及随机漏包数
如实输出，但采样为零也不能升级为 `V ⊆ Alpha output` 的严格证明。
各组另外报告最大正值，以及 `<=1e-10 m` / `>1e-10 m` 的正值计数，用于识别边界
残差量级。这只是诊断分桶，不把阈值用于更改零阈值分类或消除失败记录。

质量指标区分三个参照，不混为一个“精度”：

- 原始点：全部观测点是否被覆盖，是当前的硬验收条件。
- occupied voxel union V：Octree 有构造证书保证 C=V；Alpha 仅增加公共采样检查。
- 理想实体 / 理想凹槽：用于诊断未观测内部的覆盖及过度填充，不是输入点保守性的
  等价条件。稀疏表面点全部被覆盖，并不证明未知实体内部也被覆盖。

体积：Octree 用互不重叠 AABB 体积和；Alpha-tet 用已由后端验证的定向 wrap 三角形
体积和，坐标先平移后用浮点精细求和。这是网格体积测量，不是新的树等价严格证明。
`output_volume_over_voxels_minus_one` / `output_volume_over_ideal_solid_minus_one`
是有符号相对体积差；仅凭该比值不证明集合包含关系，负值也不截成零。

公共均匀查询分别估计 C\V 的增加体积、V\C 的遗漏体积，并给出 Wilson 95% 区间。
这是有限随机样本的不确定性范围，不是严格几何误差上界或凸包误差预算。
凹槽探测只代表指定区域，不声称验证了全局拓扑。

## 树求值计时

用公共随机查询的前 256 点，预热后执行 5 轮查询、重复计时 7 次，报告每次查询
微秒数的 min / median / max 与原始计时序列。计时包含已有递归求值器及求和防优化
开销，不含 JSON/XYZ 解析、raw 树对照和输出。没有布尔短路，不计 CBF 梯度、QP 求解、
传感器通信或控制线程调度；也不证明控制闭环安全性。

## 输出、错误处理与复现

必须使用新目录；已有目录拒绝覆盖。每组写出完成记录，失败保留日志并继续其他配置，
最终含失败时返回非零。中断后的 `matrix_complete=false` 不能当成全矩阵成功。
二进制 SHA256 与系统平台记录在 summary，方便区分不同时期的测量。

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --target build_halfspace_tree halfspace_tree_benchmark -j2
ctest --test-dir build --output-on-failure

python3 scripts/benchmark_octree_phase4.py output/octree_phase4_benchmark
```

可缩小矩阵，例如：

```bash
python3 scripts/benchmark_octree_phase4.py output/octree_phase4_small \
  --scenes l u thin_wall --densities 256 --noise-levels 0 0.003 \
  --voxel-sizes 0.04 --skip-camera --include-baselines
```

总览在 `summary.json`、`summary.csv`、`report.md`；各 case 内保存查询 XYZ / 分组范围，
各 backend 内保存 geometry、tree、trials、query_evaluation 和 quality。
当前实现不新增近似合并、非 AABB 凸近似或真正空间自适应分辨率。
