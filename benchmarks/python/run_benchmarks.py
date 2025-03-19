# Copyright (C) 2025 Severin Strobl <git@severin-strobl.de>
#
# SPDX-License-Identifier: MIT
#
# Exact calculation of the overlap volume of spheres and mesh elements.
# http://dx.doi.org/10.1016/j.jcp.2016.02.003

import pathlib
import subprocess
import sys
import tempfile

import pyperf


def main() -> None:
    benchmarks: list[pyperf.Benchmark] = []

    for benchmark_script in [p for p in (pathlib.Path(__file__).parent).glob("benchmark_*.py")]:
        name = benchmark_script.stem
        print(f"running benchmark {name}...")

        with tempfile.NamedTemporaryFile() as json_log_file:
            json_log_file.close()

            subprocess.run(
                [sys.executable, benchmark_script, "-o", json_log_file.name],
                stdout=subprocess.DEVNULL,
                cwd=benchmark_script.parent,
                check=True,
            )

            with pathlib.Path(json_log_file.name).open() as f:
                benchmarks += pyperf.BenchmarkSuite.load(f).get_benchmarks()

    benchmark_suite = pyperf.BenchmarkSuite(benchmarks)
    benchmark_suite.dump(
        str((pathlib.Path(__file__).parent) / "overlap_python.json"),
        compact=False,
        replace=True,
    )


if __name__ == "__main__":
    main()
