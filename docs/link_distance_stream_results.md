# 常驻多帧接口与相同占据复用

`link_distance_node --stream true` 从 stdin 逐行读取点云路径、输出路径及可选关节角；
URDF 只加载一次，`BoundaryDirectFrameCache(1)`、最近一次场对应的
`LinkDistanceEvaluator` 和查询工作区在进程内常驻。每帧仍从文件读取新点云并重新体素化。
若体素尺寸、排序去重后的占据键及构建选项相同，直接复用已认证的场；原始点可以在
原体素内部移动。占据键变化时，环境边界、direct 树和执行场统一全量重建，并在发布前执行完整的
严格符号认证。

## 计时口径

stdout 每帧一条 JSON。`core_ms` 从点云读取完成后到距离/梯度查询结束，
不含 URDF 加载、点云文件读取和结果 JSON 写入。缓存命中时，该时间仍含
原始点对闭合体素并集的包含性检查。`core_no_validation_ms` 只在命中时提供，
从 `core_ms` 减去这项检查耗时；重建帧的其他严格认证不作这种减法，字段为 `null`。
`--reuse-environment false` 在每帧准备场前清空缓存，用于相同进程的完整重建对照。

## 复杂 Panda 场景结果

点云为 `data/panda_complex_scene.xyz`（9216 点，3989 个占据体素），
体素尺寸 0.02 m；第二帧把每个原始点向其所在体素中心移动 25%，占据键不变。
4 组正反顺序交替运行，每进程 3 帧，首帧预热/建场不计入以下稳态统计，
每种模式共 8 个计时帧。运行：

```bash
python3 scripts/benchmark_link_distance_stream.py \
  build/link_distance_node resources/vamp/panda/panda_spherized.urdf \
  data/panda_complex_scene.xyz --voxel-size 0.02 --pairs 4 --frames 3
```

| 同一常驻节点的模式 | 稳态单帧 `core_ms`，组间 Median |
| --- | ---: |
| 关闭复用，每帧完整重建及认证 | 1591.60 ms |
| 开启复用，命中且保留原始点检查 | 4.58 ms |
| 开启复用，命中后扣除原始点检查 | 2.72 ms |

完整重建的 1591.60 ms **包含**严格认证，不能与无验证数字直接比较。
作为无验证的独立参考，原有 `link_distance_frame_benchmark` 在同一场景下
4 次独立运行、每次 11 个计时帧，其完整重建单帧 Median 的跨次中位数为
11.58 ms；新命中路径扣除检查后的 2.72 ms 包含体素化、缓存查找和七连杆
距离/梯度查询。两者不是同一进程内的配对比较，约 4.3 倍只是参考量级。

测试覆盖点坐标改变但占据不变的命中、同场不同关节角的命中、占据变化时的重建，
以及关闭复用的对照。缓存版与重建版的结果 JSON 除计时字段外相同，
`link_distance_stream`、`link_distance_node_smoke`、`boundary_direct_frame_cache`
测试均通过。
