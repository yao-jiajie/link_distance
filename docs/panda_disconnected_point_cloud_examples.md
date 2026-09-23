# Panda 多个非连通不规则点云示例

`link_distance_node` 的环境接口接收一个 XYZ 文件。多个障碍物不需要分别调用节点：
将所有分量的点放在同一个 XYZ 中即可；体素化、公共面删除和距离场构建会保留空间上
分开的占据分量。

## 六形状复杂场景（推荐压力输入）

在三分量场景的 C 形支架、波浪环面和弯管之外，增加带 1 mm 坐标扰动的扭曲波纹带、
七槽立柱和折叠非凸顶棚。机械臂配置与 base 坐标系均不变。场景包含 9216 个点；
`h=0.02 m` 下得到 3989 个占据体素和 6869 个合并边界面，明显高于三分量场景
同分辨率下的 2228 个体素和 4088 个合并面。

```bash
python3 scripts/generate_panda_complex_scene.py

build/link_distance_node \
  resources/vamp/panda/panda_spherized.urdf \
  data/panda_complex_scene.xyz \
  output/panda_complex_link_distances.json \
  --voxel-size 0.02 \
  --q 0,-0.4,0,-2.0,0,1.6,0.8 \
  --env-lse-error 0.001 \
  --link-lse-error 0.001 \
  --margin 0.005 \
  --field-from-base 0,0,0,0,0,0 \
  --query-workers 6 \
  --boundary-max-cells 262144 \
  --boundary-max-expanded-nodes 32000000 \
  --boundary-max-strata 2000000 \
  --boundary-max-sign-ops 40000000000

python3 scripts/visualize_link_distances_3d.py \
  output/panda_complex_link_distances.json \
  output/panda_complex_link_distances_3d.html \
  --cloud data/panda_complex_scene.xyz \
  --clearance-scale 0.15
```

当前固定姿态下，七组连杆的 smooth clearance proxy 最小值约为 `174.36 mm`，
59 个球的 clearance 均为正；球面到原始输入点的最小直接距离约为 `251 mm`。
较大的 direct-tree 参数仍是明确有界的拒绝上限，用于容纳这个压力场景，并不要求
算法一定执行到该上限。解析定义、随机种子、采样噪声和边界记录在
[`panda_complex_scene.json`](../data/panda_complex_scene.json)，交互式输出为
[`panda_complex_link_distances_3d.html`](../output/panda_complex_link_distances_3d.html)。

## 三分量场景

推荐示例包含三个不与机械臂相撞的不规则表面点云：

1. 左后方的非凸 C 形支架；
2. 右后方的三瓣波浪环面；
3. 机械臂后方的变半径弯管。

机械臂使用固定配置 `0,-0.4,0,-2.0,0,1.6,0.8`，环境直接放在 base 坐标系，
没有修改机械臂位姿。确定性重新生成输入并运行节点：

```bash
python3 scripts/generate_panda_disconnected_scene.py

build/link_distance_node \
  resources/vamp/panda/panda_spherized.urdf \
  data/panda_disconnected_irregular_scene.xyz \
  output/panda_three_disconnected_link_distances.json \
  --voxel-size 0.04 \
  --q 0,-0.4,0,-2.0,0,1.6,0.8 \
  --env-lse-error 0.001 \
  --link-lse-error 0.001 \
  --margin 0.005 \
  --field-from-base 0,0,0,0,0,0 \
  --query-workers 4

python3 scripts/visualize_link_distances_3d.py \
  output/panda_three_disconnected_link_distances.json \
  output/panda_three_disconnected_3d.html \
  --cloud data/panda_disconnected_irregular_scene.xyz \
  --clearance-scale 0.15
```

默认 4608 个点在 `h=0.04 m` 下形成 801 个占据体素。按共享面定义的 6 邻接检查
得到三个连通分量，大小分别为 404、202、195 个体素。当前固定姿态下七组连杆的
smooth clearance proxy 最小值为 `164.705 mm`，七组结果均为正。

输入参数和三个解析表面定义记录在
[`panda_disconnected_irregular_scene.json`](../data/panda_disconnected_irregular_scene.json)。
交互式结果为
[`panda_three_disconnected_3d.html`](../output/panda_three_disconnected_3d.html)，静态结果为
[`panda_three_disconnected.svg`](../output/panda_three_disconnected.svg)。

## 两分量场景

已有的 C 形支架＋波浪环面输入也可以直接用于机械臂节点。保持相同关节角，通过
`field_from_base` 只平移环境：

```bash
build/link_distance_node \
  resources/vamp/panda/panda_spherized.urdf \
  data/disconnected_irregular_surface.xyz \
  output/panda_disconnected_irregular_link_distances.json \
  --voxel-size 0.04 \
  --q 0,-0.4,0,-2.0,0,1.6,0.8 \
  --env-lse-error 0.001 \
  --link-lse-error 0.001 \
  --margin 0.005 \
  --field-from-base 0,0.6,0,0,0,0 \
  --query-workers 4

python3 scripts/visualize_link_distances_3d.py \
  output/panda_disconnected_irregular_link_distances.json \
  output/panda_disconnected_irregular_3d.html \
  --cloud data/disconnected_irregular_surface.xyz \
  --clearance-scale 0.15
```

这个场景的七组连杆距离也全部为正，最小值为 `214.88 mm`。更详细的点云自身构建
说明见[两个不连通不规则点云示例](disconnected_irregular_point_cloud_example.md)。
