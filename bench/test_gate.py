import json
import os
import tempfile
import unittest
from pathlib import Path

import gate


def write_estimates(path: Path, point: float, lower: float) -> None:
    path.mkdir(parents=True, exist_ok=True)
    mean = {"point_estimate": point, "confidence_interval": {"lower_bound": lower, "upper_bound": point}}
    (path / "estimates.json").write_text(json.dumps({"mean": mean}))


def write_bench(root: Path, name: str, change_point: float, change_lower: float) -> None:
    bench = root / name
    for run in ("base", "new"):
        write_estimates(bench / run, 1000.0, 1000.0)
        (bench / run / "benchmark.json").write_text(json.dumps({"full_id": name}))
    write_estimates(bench / "change", change_point, change_lower)
    later = os.stat(bench / "base" / "estimates.json").st_mtime_ns + 1_000_000
    os.utime(bench / "change" / "estimates.json", ns=(later, later))


class GateTest(unittest.TestCase):
    def run_gate(self, change_point: float, change_lower: float) -> int:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_bench(root, "parse", change_point, change_lower)
            return gate.main(["--criterion-dir", str(root)])

    def test_noisy_slowdown_with_overlapping_interval_passes(self) -> None:
        self.assertEqual(self.run_gate(0.40, 0.05), 0)

    def test_confident_slowdown_fails(self) -> None:
        self.assertEqual(self.run_gate(0.40, 0.25), 1)

    def test_memory_growth_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_bench(root, "parse", 0.0, -0.01)
            base = root / "base.json"
            head = root / "head.json"
            base.write_text(json.dumps({"files": {"a.cif": {"peak_bytes": 100, "seconds": 1}}, "sweep_peak_bytes": 100}))
            head.write_text(json.dumps({"files": {"a.cif": {"peak_bytes": 110, "seconds": 1}}, "sweep_peak_bytes": 100}))
            arguments = ["--criterion-dir", str(root), "--base-memory", str(base), "--head-memory", str(head)]
            self.assertEqual(gate.main(arguments), 1)


if __name__ == "__main__":
    unittest.main()
