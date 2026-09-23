# 不同点云输入测试

在原相机输入之外增加 19 组可复现点云，组成 20 场景测试集。全部 XYZ 已保存于
[data/point_cloud_cases](../data/point_cloud_cases/manifest.json)，可直接作为命令行输入。
生成参数、点数、边界范围、文件哈希和输入之间的关系记录在 manifest 中；
[点云预览图](../data/point_cloud_cases/preview.svg)按同一投影比例显示各组输入，每组最多显示 768 个点。

## 输入范围与含义

所有坐标单位为米。除复制保留的 `camera` 基线使用原多视角深度相机模型外，
其余点云是直接生成的几何表面样本，没有传感器、可见性或真实采集数据的含义。
当前占据集合仍由落入体素的点定义，程序不额外填充模型内部。

- 立方体、L/U 形、台阶、双壁和分离物体按长方体并集的外露表面采样，剔除隐藏接合面。
- 圆柱半径 14 cm、高 32 cm，包含侧面和端面；环面主半径 16 cm、管半径 6 cm，
  按角参数均匀采样，不声称表面积均匀。
- 空腔样本包含外部 36 cm 箱体和内部 18 cm 箱体的表面，中心保留空闲探针。
  薄壁物理厚度为 6 mm；基准体素为 4 cm，结果仍受该离散尺度影响。
- 稀疏、普通和密集椭球分别为同一采样序列的前 256、2048、8192 点；
  方向均匀采样后按椭球轴缩放，不声称表面积均匀。
- 噪声版本在同一 2048 点上加入独立笛卡尔高斯噪声，2/8 mm 是每个坐标轴的标准差，
  并非三维位移 RMS 或相机深度噪声。实测位移 RMS 另存于 manifest。
- 局部视野版本取 `x >= 0` 的 1024 个表面点；重复版本将同一组点复制三份并打乱；
  离群版本添加六个明确列出的独立远端点，测试它们仍被占据集合保留。

## 不含验证的点云到函数耗时

2026-09-10，本机 Release，4 cm 体素、1 mm LSE 标量误差预算、结构共享、二元 LSE、
packed 符号传播、自动验证缓存。预热一轮后测量七轮，每轮轮换场景顺序、串行启动新进程，
每个进程重新构建函数；正式计时期间没有并行运行本次构建或 CTest。

**下表取 `core_runtime_ms` 中位数，排除验证和文件读写。** 该计时从内存原始点云开始，
包含体素化、边界与树构造、共享执行程序编译、LSE 参数和函数值/梯度工作区准备，
在所有验证开始前结束。它不包含第一次查询、缓存命中复用或多线程查询耗时。

| 输入（点击下载/打开 XYZ） | 点数 | 占据体素 | 执行节点 | 到函数 ms |
| --- | ---: | ---: | ---: | ---: |
| [原相机椭球（基线）](../data/point_cloud_cases/camera.xyz) | 1655 | 341 | 429 | 1.017011 |
| [立方体表面](../data/point_cloud_cases/cube_surface.xyz) | 2048 | 324 | 23 | 0.510080 |
| [L 形表面](../data/point_cloud_cases/l_surface.xyz) | 2048 | 162 | 37 | 0.428970 |
| [U 形表面](../data/point_cloud_cases/u_surface.xyz) | 2048 | 222 | 51 | 0.473570 |
| [6 mm 薄壁表面](../data/point_cloud_cases/thin_wall_surface.xyz) | 2048 | 141 | 30 | 0.405171 |
| [三级台阶表面](../data/point_cloud_cases/staircase_surface.xyz) | 2048 | 361 | 87 | 0.623107 |
| [两个分离物体](../data/point_cloud_cases/two_objects_surface.xyz) | 2048 | 302 | 85 | 0.606415 |
| [双壁窄通道](../data/point_cloud_cases/narrow_corridor_surface.xyz) | 2048 | 497 | 152 | 0.871351 |
| [圆柱表面](../data/point_cloud_cases/cylinder_surface.xyz) | 2048 | 317 | 210 | 0.752842 |
| [环面（中孔）](../data/point_cloud_cases/torus_surface.xyz) | 2048 | 259 | 226 | 0.777868 |
| [箱体外表面与内腔表面](../data/point_cloud_cases/hollow_box_surfaces.xyz) | 2048 | 588 | 382 | 1.131389 |
| [开放平面片](../data/point_cloud_cases/plane_patch.xyz) | 1024 | 80 | 11 | 0.195680 |
| [稀疏椭球](../data/point_cloud_cases/ellipsoid_sparse.xyz) | 256 | 169 | 473 | 0.717609 |
| [无噪声椭球](../data/point_cloud_cases/ellipsoid_surface.xyz) | 2048 | 339 | 410 | 1.064290 |
| [密集椭球](../data/point_cloud_cases/ellipsoid_dense.xyz) | 8192 | 378 | 330 | 1.665706 |
| [椭球＋2 mm 各轴噪声](../data/point_cloud_cases/ellipsoid_noise_2mm.xyz) | 2048 | 341 | 424 | 1.020193 |
| [椭球＋8 mm 各轴噪声](../data/point_cloud_cases/ellipsoid_noise_8mm.xyz) | 2048 | 416 | 497 | 1.157479 |
| [椭球局部半面](../data/point_cloud_cases/ellipsoid_partial.xyz) | 1024 | 174 | 213 | 0.575716 |
| [椭球三份重复并打乱](../data/point_cloud_cases/ellipsoid_duplicates.xyz) | 6144 | 339 | 410 | 1.641761 |
| [椭球＋6 个离群点](../data/point_cloud_cases/ellipsoid_outliers.xyz) | 2054 | 345 | 542 | 1.185258 |

本轮所有输入的中位数范围为 0.196–1.666 ms。点数并不能独立决定构造成本：例如稀疏
椭球的占据分布更零散，执行程序有 473 个节点，密集椭球为 330 个；后者体素化的输入点更多。
相机本轮为约 1.02 ms，属于本轮机器状态下的实测值，不用不同轮次的计时差推断实现变化。
逐轮时间、最小值、p95、单独的验证时间、执行命令与二进制哈希见
[正式 benchmark 报告](../output/point_cloud_input_benchmark_final/report.md)。七个样本的最近秩 p95 等于最大值。

## 验证内容

计时排除验证，但测试运行仍执行原有全部验证，默认工作量上限保持不变。

- 所有输入通过原始点包含性、外露面/patch 覆盖和全体、面、边、顶点 strata 严格符号证书。
- 所有输入通过 LSE 保守次序与误差界、解析梯度有限差分和零面射线审计。
- Python 独立解释导出的 MIN/MAX 与 LSE 树，对每个场景从距离 CSV 均匀选取 24 行，
  重算硬值、平滑值、梯度和有限矩形边界的带符号参考距离，并核对梯度范数及平滑误差界。
- 八个显式探针覆盖 L/U 凹口、物体间隙、通道、环面中孔、箱体空腔、开放平面和未观测半面；
  同时检查探针不在占据集合中，硬值与平滑值均判为空闲。
- 同一输入的八次执行（含预热）导出一致；重复且乱序的输入与原输入导出相同的占据键、
  边界、支撑、分割和树。密度序列的占据键严格包含，六个离群点分别增加六个占据体素。
- 生成器检查解析表面方程、噪声缩放关系、重复点计数、XYZ 浮点往返、文件哈希和重生成一致性，
  并拒绝覆盖已有输出目录。

这些检查针对输入所生成的体素集合及其函数；不把有限采样误认为原解析实体的完整重建，
也不把 LSE 标量误差预算当成全局欧氏距离误差保证。

新增两个 CTest：`point_cloud_input_fixtures` 和 `point_cloud_inputs`。后者覆盖全部 20 组输入，
使用 100 个随机验证样本；正式基准使用 1000 个，两者均保留全 strata 证书以及默认距离
审计的梯度和射线预算。Release 全部 **39/39 CTest** 通过。本次新增输入、测试和报告工具，
没有修改 C++ 求值或构造实现。

## 使用与复现

直接运行任意已保存的输入：

```bash
build/build_halfspace_tree data/point_cloud_cases/torus_surface.xyz output/torus_input_new \
  --pipeline octree --tree-source boundary --boundary-expression direct \
  --boundary-query dag --voxel-size .04 --distance-field lse --lse-error .001
```

重建、运行完整测试和正式输入基准：

```bash
cmake -S . -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
ctest --test-dir build --output-on-failure
python3 scripts/benchmark_point_cloud_inputs.py output/point_cloud_inputs_new \
  --repeats 7 --warmups 1 --validation-samples 1000
```

生成一套新的相同输入，或只测选定场景：

```bash
python3 scripts/generate_point_cloud_cases.py output/point_cloud_cases_new
python3 scripts/benchmark_point_cloud_inputs.py output/selected_cloud_results_new \
  --input-dir output/point_cloud_cases_new \
  --cases cylinder_surface torus_surface ellipsoid_dense --repeats 7
```

生成器与基准的输出目录均需尚不存在；已有相机基线和此前 benchmark 产物保留。
