# Abseil 哈希容器试验

## 结论

只在 direct 半空间树的子表达式驻留表使用 `absl::flat_hash_map` 有收益；
体素邻居表改用 `absl::flat_hash_set` 则变慢。因此两个替换均保留为独立编译开关，
默认关闭，以免没有 Abseil 开发包的现有构建失效。推荐只启用前者。

## 构建

需要 Abseil 开发包，例如 Ubuntu 的 `libabsl-dev`。在已有项目依赖配置的基础上：

```bash
cmake -S . -B build-absl -DBUILD_TESTING=ON \
  -DROKAE_USE_ABSEIL_DIRECT_HASH=ON \
  -DROKAE_USE_ABSEIL_BOUNDARY_HASH=OFF
cmake --build build-absl -j4 --target link_distance_node link_distance_frame_benchmark
```

`ROKAE_USE_ABSEIL_BOUNDARY_HASH=ON` 也可用于复现实验，但不推荐用于当前数据。

## 单帧测试

输入 `data/panda_complex_scene.xyz`（9216 点）、Panda URDF、体素尺寸 0.02 m，
机器人模型常驻；单帧包含环境重建和一次七连杆距离/梯度查询，不含验证、文件 I/O。
每个进程预热 1 帧，测 11 帧；标准版和 Abseil 版交替运行 8 组。下表的差值是
每组 Median 或 P95 的配对差值再取中位数，不是两个独立总体中位数的差。

| 方案 | 完整帧 Median 的跨组中位数 | 相对标准版的配对 Median 差值 | 配对 P95 差值 |
| --- | ---: | ---: | ---: |
| `std::unordered_map` | 11.65 ms | 基准 | 基准 |
| 仅 direct `absl::flat_hash_map` | 10.85 ms | -0.78 ms（约 -6.7%） | -0.97 ms |

8 组中 7 组完整帧 Median 更快；direct-tree 阶段的配对 Median 差值为 -0.60 ms。
若同时替换体素邻居表，6 组配对中其提取阶段反而慢约 0.40 ms，
抵消了部分构树收益，因此不把它作为推荐配置。结果只针对上述场景和本机
Abseil 20210324 版本；其他输入规模和版本仍需复测。

两版的 `boundary_direct_tree`、`octree_boundary`、`boundary_frontend` 测试均通过。
对同一 Panda 输入运行 `link_distance_node` 后，去除 `environment_build_ms` 和
`query_ms` 两个计时字段，JSON 逐字节一致；七个连杆的距离与关节梯度没有变化。
