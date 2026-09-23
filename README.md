# Link Distance Node

本项目的核心程序是 `link_distance_node`：输入机械臂球分解 URDF、环境 XYZ 点云和关节角，输出每个关节连杆组相对环境的平滑距离代理值，以及该值对全部关节变量的一阶解析梯度。它目前是独立 C++ 可执行程序，不依赖 ROS；连续帧通过 stdin stream 接口输入。

> 输出是 `smooth_clearance_proxy_not_euclidean_sdf`，不是严格的欧氏符号距离。正值表示平滑模型下有间隙，零附近表示接触带，负值表示平滑模型下发生碰撞。

## 核心链路

```text
球分解 URDF + q ──→ 正运动学、球心 Jacobian ─────────────┐
                                                        ├─→ 各球间隙 ─→ 每连杆 smooth-min ─→ 距离与 ∂d/∂q
XYZ 点云 ─→ 占据体素 ─→ direct halfspace DAG ─→ LSE 场 ─┘
```

节点会完成：

1. 解析 URDF 中的球形 `<collision>`；
2. 将点云体素化，并删除内部公共面、合并共面外表面；
3. 构建环境 direct MIN/MAX DAG；
4. 对环境树执行保守 LSE 平滑，并进行严格符号认证；
5. 计算每个球的环境间隙与空间梯度；
6. 按关节参考坐标系合并固定子连杆；
7. 对同组球执行 smooth-min，输出每组距离和关节梯度。

只有环境场通过严格符号认证后，节点才会写出结果。

## 构建

要求：

- C++17 编译器；
- CMake 3.16 或更新版本；
- Eigen3、Boost、GMP、MPFR；
- CGAL 5.6.3。仓库的 `_deps/cgal` 中提供了固定版本。

Ubuntu 可安装基础依赖：

```bash
sudo apt install build-essential cmake libeigen3-dev libboost-dev libgmp-dev libmpfr-dev
```

配置并构建核心节点：

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCGAL_DIR="$PWD/_deps/cgal/lib/cmake/CGAL" \
  -DBUILD_TESTING=ON

cmake --build build -j"$(nproc)" --target link_distance_node
```

构建全部目标并运行测试：

```bash
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

当前全量测试为 46 项。

## 快速运行

```bash
mkdir -p output

build/link_distance_node \
  resources/vamp/panda/panda_spherized.urdf \
  data/realistic_sparse_camera.xyz \
  output/panda_link_distances.json \
  --voxel-size 0.02 \
  --q 0,-0.4,0,-2.0,0,1.6,0.8 \
  --env-lse-error 0.001 \
  --link-lse-error 0.001 \
  --margin 0.005 \
  --field-from-base 0,0.3,0,0,0,0 \
  --query-workers 6
```

成功后输出：

```text
Wrote 7 link distance results to "output/panda_link_distances.json"
```

Panda 的可动关节顺序由程序从 URDF 读取。固定连接的 `panda_link8`、手和手指球会合并到最近的上游可动关节组，因此最终输出 7 个连杆组，而不是按 URDF link 标签输出 11 组。

## 输入约定

### 球分解 URDF

节点读取每个 link 的球形 collision：

```xml
<collision>
  <origin xyz="..." rpy="..."/>
  <geometry>
    <sphere radius="..."/>
  </geometry>
</collision>
```

支持 fixed、revolute、continuous 和 prismatic joint，支持一个 link 上的多个球。固定子树会合并到最近的上游可动关节。

不支持 box、cylinder、mesh 等非球形 collision，也不支持尚未展开的可动 mimic joint。仓库中的 `resources/vamp/panda/panda_spherized.urdf` 可以直接作为 Panda 球分解模型输入。

### XYZ 环境点云

每个有效行必须恰好包含三列有限浮点数：

```text
x y z
```

允许空行和 `#` 注释。坐标单位统一为米；毫米数据必须先除以 1000。多个非连通障碍物可以放在同一个 XYZ 文件中，节点按占据体素并集处理，不要求点云连通。

### 关节配置

`--q` 的数量必须和 URDF 可动关节数一致。实际顺序写入结果 JSON 的 `joint_order`，不要假定它与 link 名称排序一致。

### 坐标系

点云位于 field frame，URDF 运动学位于 robot base frame。

`--field-from-base x,y,z,roll,pitch,yaw` 给出 base 到 field 的变换，平移单位为米、角度单位为弧度。旋转顺序为：

```text
Rz(yaw) · Ry(pitch) · Rx(roll)
```

默认是单位变换。

## 命令行参数

基本形式：

```text
link_distance_node robot.urdf environment.xyz output.json [options]
```

必需参数：

| 参数 | 含义 |
|---|---|
| `--voxel-size M` | 环境体素边长，单位 m，必须大于 0 |
| `--q q1,q2,...` | URDF 可动关节配置 |
| `--env-lse-error M` | 环境树标量 LSE 平滑上界 |
| `--env-lse-beta 1/M` | 直接指定环境 LSE 逆温度 |
| `--link-lse-error M` | 连杆球集合 smooth-min 上界 |
| `--link-lse-alpha 1/M` | 直接指定连杆 smooth-min 逆温度 |

`--env-lse-error` 与 `--env-lse-beta` 必须且只能选一个；`--link-lse-error` 与 `--link-lse-alpha` 同理。

常用可选参数：

| 参数 | 默认值 | 含义 |
|---|---:|---|
| `--margin M` | `0` | 从每个球间隙中减去的安全裕量 |
| `--field-from-base x,y,z,r,p,y` | 单位变换 | base 到 field 的刚体变换 |
| `--query-workers N` | `min(6, hardware concurrency)` | 球心环境场批查询线程数 |
| `--profile-query true\|false` | `false` | 输出查询分阶段时间，会增加计时开销 |
| `--stream true\|false` | `false` | 启用多帧 stdin 接口 |
| `--reuse-environment true\|false` | `true` | stream 模式下复用完全相同的占据环境 |

## 距离和梯度定义

第 `i` 个球的间隙代理为：

```text
s_i(q) = D_env(c_i(q)) - r_i - margin
```

其中 `c_i(q)` 是 field frame 中的球心，`D_env` 是环境 direct DAG 的保守 LSE 距离代理，`r_i` 是球半径。

包含 `n` 个球的连杆组使用：

```text
d_link(q) = -log(sum_i exp(-alpha · s_i(q))) / alpha
```

解析梯度为：

```text
∂d_link/∂q = sum_i w_i · J_i(q)^T · ∇D_env(c_i(q))
```

`gradient_dq[k]` 对应 `joint_order[k]`。连杆 smooth-min 相对硬最小值的标量差由 `log(n)/alpha` 控制，环境树的 LSE 上界单独写入 JSON。

这些上界描述代数树平滑误差，不是相对真实、未观测物体表面的欧氏距离误差。

## 输出 JSON

顶层关键字段：

| 字段 | 含义 |
|---|---|
| `value_kind` | `smooth_clearance_proxy_not_euclidean_sdf` |
| `joint_order` | 梯度和配置采用的关节顺序 |
| `configuration` | 本次查询关节配置 |
| `voxel_size_m` | 实际体素尺寸 |
| `environment_beta_per_m` | 环境 LSE 参数 |
| `environment_scalar_smoothing_bound_m` | 环境代数平滑上界 |
| `link_alpha_per_m` | 连杆 smooth-min 参数 |
| `maximum_link_smoothing_bound_m` | 最大连杆平滑上界 |
| `environment_build_ms` | 环境构建和认证耗时 |
| `query_ms` | 一次机械臂距离和梯度查询耗时 |
| `voxel_occupancy` | 占据体素及合并后外边界面 |
| `links` | 各关节连杆组的距离与关节梯度 |
| `spheres` | 各碰撞球的球心、半径、环境值和空间梯度 |

每个 `links[]` 元素主要包含：

```json
{
  "name": "panda_link3",
  "joint": "panda_joint3",
  "joint_configuration_index": 2,
  "source_links": ["panda_link3"],
  "sphere_count": 8,
  "distance_proxy_m": 0.12,
  "gradient_dq": [0.0, 0.01, -0.03, 0.0, 0.0, 0.0, 0.0],
  "hard_min_proxy_m": 0.13,
  "nearest_sphere": 17,
  "link_smoothing_bound_m": 0.002
}
```

以上数值仅为字段示例。`source_links` 给出合并到该关节组的 URDF links；`nearest_sphere` 是硬最小间隙球的全局编号。静态 base/root 球仍出现在 `spheres` 诊断中，但不形成可动关节距离组。

## 多帧点云接口

环境逐帧变化时，使用 stream 模式，使 URDF、球模型和查询工作区常驻：

```bash
printf 'data/frame_001.xyz\toutput/frame_001.json\ndata/frame_002.xyz\toutput/frame_002.json\n' |
  build/link_distance_node \
    resources/vamp/panda/panda_spherized.urdf - - \
    --stream true \
    --voxel-size 0.02 \
    --q 0,-0.4,0,-2.0,0,1.6,0.8 \
    --env-lse-error 0.001 \
    --link-lse-error 0.001 \
    --margin 0.005
```

stdin 每行格式：

```text
cloud.xyz<TAB>output.json[<TAB>q1,q2,...]
```

第三列可覆盖该帧默认 `--q`。stdout 每帧输出一行状态 JSON，包括 `cache_hit`、`evaluator_reused`、`occupancy_ms`、`boundary_ms`、`direct_ms`、`environment_ms`、`query_ms` 和 `core_ms`。`core_ms` 排除点云读取和结果 JSON 写入；cache miss 时包含完整环境认证。

缓存只复用体素尺寸、环境选项和占据体素键完全相同的帧。占据发生变化时执行完整环境重建和严格认证；当前没有局部增量更新路径。`--reuse-environment false` 可在 stream 模式下强制每帧完整重建，用于 A/B 测试。

## 查询分析与性能

开启查询阶段分析：

```bash
build/link_distance_node \
  resources/vamp/panda/panda_spherized.urdf \
  data/realistic_sparse_camera.xyz \
  output/profiled.json \
  --voxel-size 0.02 \
  --q 0,-0.4,0,-2.0,0,1.6,0.8 \
  --env-lse-error 0.001 \
  --link-lse-error 0.001 \
  --field-from-base 0,0.3,0,0,0,0 \
  --profile-query true
```

结果会增加 `query_phase_ms.kinematics`、`query_phase_ms.field` 和 `query_phase_ms.link`。`--profile-query true` 会增加三个计时点，不应用于最终性能验收。

在 `data/panda_complex_scene.xyz`、`--voxel-size 0.02`、6 个 worker 上，排除认证、文件 I/O、JSON 和可视化的“点数组到七连杆结果”全量单帧基准中位数约为：

```text
5.755 ms ≈ 173.8 Hz
```

该基准包含占据体素构建、direct DAG、场编译和一次七连杆距离/梯度查询，不等同于节点冷启动端到端墙钟时间。复现实验：

```bash
cmake --build build -j"$(nproc)" --target link_distance_frame_benchmark

build/link_distance_frame_benchmark \
  resources/vamp/panda/panda_spherized.urdf \
  data/panda_complex_scene.xyz \
  0.02 101
```

详细结果见 [动态单帧性能记录](docs/link_distance_dynamic_frame_results.md)。

## 验证球参考坐标系

独立从 URDF XML 重建固定关节链和正运动学，检查每个球所属关节参考系：

```bash
python3 scripts/validate_sphere_joint_frames.py \
  output/panda_link_distances.json \
  --urdf resources/vamp/panda/panda_spherized.urdf
```

脚本会核对 `center_urdf_link_m`、`center_group_link_m`、`joint_group` 和 field-frame 球心。

## 3D 可视化

生成自包含交互式 HTML：

```bash
python3 scripts/visualize_link_distances_3d.py \
  output/panda_link_distances.json \
  output/panda_link_distances_3d.html \
  --cloud data/realistic_sparse_camera.xyz \
  --clearance-scale 0.1
```

可显示原始点云、机械臂碰撞球、各连杆最近球、坐标轴，以及删除内部公共面并合并共面矩形后的占据体素外表面。

正交投影 SVG：

```bash
python3 scripts/visualize_link_distances.py \
  output/panda_link_distances.json \
  output/panda_link_distances.svg \
  --cloud data/realistic_sparse_camera.xyz
```

两个可视化脚本都转换到 robot base frame 显示。

## 高级安全上限

以下参数是拒绝失控工作的上限，不是期望工作量：

- `--boundary-max-cells`
- `--boundary-max-direct-nodes`
- `--boundary-max-expanded-nodes`
- `--boundary-max-split-checks`
- `--boundary-max-direct-depth`，最大 256
- `--boundary-max-strata`
- `--boundary-max-sign-ops`

实际采用的上限写入 `environment_direct_limits`。超过上限时节点报错退出，不会静默切换几何表示。

## 当前限制

- 只使用 URDF sphere collision；
- 不计算机械臂自碰撞；
- 不检查关节上下限；
- 环境表示为体素并集，而不是原始点集的连续曲面；
- 距离是两层 LSE 平滑的保守代理，不是欧氏 SDF；
- 新环境帧全量重建，没有局部增量更新；
- JSON 包含球和合并体素面诊断，复杂场景下文件可能较大。

## 相关文件

核心实现：

- `src/link_distance_node.cpp`：命令行、单帧/多帧接口和 JSON；
- `src/link_distance_evaluator.cpp`：球查询、连杆 smooth-min 和关节梯度；
- `src/spherical_robot_model.cpp`：球分解 URDF、关节分组和正运动学；
- `src/boundary_direct_tree.cpp`：占据体素到 direct halfspace DAG；
- `src/boundary_direct_frame_cache.cpp`：环境快照与相同帧复用；
- `include/rokae_demo/halfspace_smooth_tree.hpp`：环境 LSE 值和空间梯度。

样例和结果：

- [Panda 球分解 URDF](resources/vamp/panda/panda_spherized.urdf)
- [多个非连通不规则点云](docs/panda_disconnected_point_cloud_examples.md)
- [动态单帧性能](docs/link_distance_dynamic_frame_results.md)
- [查询内核优化记录](docs/link_distance_predecoded_kernel_results.md)
- [stream 测试记录](docs/link_distance_stream_results.md)

`build_halfspace_tree`、Alpha Wrap、Octree 和其他 benchmark 仍保留为底层研究、验证和对照工具，但不是本 README 的主入口。
