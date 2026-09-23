#!/usr/bin/env python3
import json
import math
import pathlib
import subprocess
import sys
import tempfile


def run(node, urdf, cloud, output, extra=(), stdin=None, q="0.3,0.1"):
    command = [
        str(node), str(urdf), str(cloud), str(output),
        "--voxel-size", "0.1", "--q", q,
        "--env-lse-error", "0.001", "--link-lse-error", "0.001",
        "--margin", "0.005", *extra,
    ]
    completed = subprocess.run(command, input=stdin, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    return completed


def without_timing(path):
    result = json.loads(path.read_text())
    result.pop("environment_build_ms")
    result.pop("query_ms")
    return result


def main():
    node, urdf, cloud = map(pathlib.Path, sys.argv[1:])
    with tempfile.TemporaryDirectory() as folder:
        folder = pathlib.Path(folder)
        points = [tuple(map(float, row.split())) for row in cloud.read_text().splitlines()]
        moved = folder / "same_voxels.xyz"
        moved.write_text("".join(
            " ".join(f"{value + 0.25 * ((math.floor(value / 0.1) + 0.5) * 0.1 - value):.17g}"
                     for value in point) + "\n"
            for point in points
        ))
        changed = folder / "changed_voxels.xyz"
        changed.write_text(moved.read_text() + "0.2 0.2 0.2\n")
        outputs = [folder / f"frame_{index}.json" for index in range(4)]
        rows = [
            f"{cloud}\t{outputs[0]}",
            f"{moved}\t{outputs[1]}",
            f"{moved}\t{outputs[2]}\t0.25,0.12",
            f"{changed}\t{outputs[3]}",
        ]
        streamed = run(node, urdf, "-", "-", ("--stream", "true"),
                       "\n".join(rows) + "\n")
        summaries = [json.loads(row) for row in streamed.stdout.splitlines()]
        assert [item["cache_hit"] for item in summaries] == [False, True, True, False]
        assert [item["evaluator_reused"] for item in summaries] == [False, True, True, False]
        assert [item["frame"] for item in summaries] == list(range(4))
        assert all(item["core_ms"] >= item["query_ms"] >= 0 for item in summaries)
        assert summaries[0]["core_no_validation_ms"] is None
        assert all(item["core_no_validation_ms"] is not None for item in summaries[1:3])
        assert all(item["occupied_voxels"] > 0 for item in summaries)
        assert json.loads(outputs[2].read_text())["configuration"] == [0.25, 0.12]

        standalone = folder / "standalone.json"
        run(node, urdf, moved, standalone)
        assert without_timing(outputs[1]) == without_timing(standalone)
        run(node, urdf, moved, standalone, q="0.25,0.12")
        assert without_timing(outputs[2]) == without_timing(standalone)
        run(node, urdf, changed, standalone)
        assert without_timing(outputs[3]) == without_timing(standalone)

        no_reuse = [folder / f"no_reuse_{index}.json" for index in range(2)]
        rows = [f"{cloud}\t{no_reuse[0]}", f"{moved}\t{no_reuse[1]}"]
        rebuilt = run(node, urdf, "-", "-",
                      ("--stream", "true", "--reuse-environment", "false"),
                      "\n".join(rows) + "\n")
        assert [json.loads(row)["cache_hit"] for row in rebuilt.stdout.splitlines()] == [False, False]
        assert without_timing(no_reuse[1]) == without_timing(outputs[1])


if __name__ == "__main__":
    main()
