import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

DEFAULT_CRITERION_DIR = Path("target/criterion")
DEFAULT_BASELINE = "base"
DEFAULT_TIME_THRESHOLD_PERCENT = 15.0
DEFAULT_MEMORY_THRESHOLD_PERCENT = 5.0
CURRENT_RUN = "new"
SWEEP_ENTRY = "sweep_peak_bytes"
OK = "ok"
REGRESSED = "REGRESSED"
MISSING_IN_BASE = "missing in base"
MISSING_IN_HEAD = "missing in head"
TIME_UNITS = (("s", 1e9), ("ms", 1e6), ("us", 1e3), ("ns", 1.0))
BYTES_PER_MIB = 1024 * 1024


@dataclass(frozen=True)
class Row:
    name: str
    base: str
    head: str
    change: float | None
    status: str


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fail when criterion benchmarks or mmcif-mem peak heap bytes regress. "
            "Time changes come from <criterion-dir>/**/change/estimates.json written by "
            "`cargo bench -- --baseline NAME`; memory compares two `mmcif-mem` JSON reports. "
            "Exit status: 0 pass, 1 regression, 2 usage error."
        )
    )
    parser.add_argument("--criterion-dir", type=Path, default=DEFAULT_CRITERION_DIR)
    parser.add_argument("--baseline", default=DEFAULT_BASELINE, help="criterion baseline name")
    parser.add_argument(
        "--time-threshold",
        type=float,
        default=DEFAULT_TIME_THRESHOLD_PERCENT,
        metavar="PERCENT",
        help="fail when the 95%% confidence lower bound of the slowdown exceeds this (default %(default)s)",
    )
    parser.add_argument(
        "--memory-threshold",
        type=float,
        default=DEFAULT_MEMORY_THRESHOLD_PERCENT,
        metavar="PERCENT",
        help="largest allowed peak heap growth (default %(default)s)",
    )
    parser.add_argument("--base-memory", type=Path, help="mmcif-mem report of the base commit")
    parser.add_argument("--head-memory", type=Path, help="mmcif-mem report of the head commit")
    return parser.parse_args(argv)


def load_json(path: Path) -> dict:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def mean_estimate(path: Path) -> float:
    return load_json(path)["mean"]["point_estimate"]


def change_lower_bound(path: Path) -> float:
    return load_json(path)["mean"]["confidence_interval"]["lower_bound"]


def bench_dirs(criterion_dir: Path, baseline: str) -> list[Path]:
    runs = {CURRENT_RUN, baseline}
    found = {
        path.parent.parent
        for path in criterion_dir.rglob("benchmark.json")
        if path.parent.name in runs
    }
    return sorted(found)


def bench_id(bench_dir: Path, baseline: str) -> str:
    for run in (CURRENT_RUN, baseline):
        path = bench_dir / run / "benchmark.json"
        if path.exists():
            return load_json(path)["full_id"]
    return bench_dir.name


def fresh_change(bench_dir: Path, baseline: str) -> Path | None:
    change = bench_dir / "change" / "estimates.json"
    base = bench_dir / baseline / "estimates.json"
    if change.exists() and change.stat().st_mtime_ns > base.stat().st_mtime_ns:
        return change
    return None


def format_time(nanoseconds: float) -> str:
    for unit, scale in TIME_UNITS:
        if nanoseconds >= scale:
            return f"{nanoseconds / scale:.4g} {unit}"
    return f"{nanoseconds:.4g} ns"


def format_bytes(count: int) -> str:
    return f"{count / BYTES_PER_MIB:.2f} MiB"


def format_change(change: float | None) -> str:
    return "-" if change is None else f"{change * 100:+.1f}%"


def status_for(change: float, threshold: float) -> str:
    return REGRESSED if change > threshold else OK


def time_row(bench_dir: Path, baseline: str, threshold: float) -> Row:
    name = bench_id(bench_dir, baseline)
    base_estimates = bench_dir / baseline / "estimates.json"
    head_estimates = bench_dir / CURRENT_RUN / "estimates.json"
    if not base_estimates.exists():
        return Row(name, "-", format_time(mean_estimate(head_estimates)), None, MISSING_IN_BASE)
    base_time = format_time(mean_estimate(base_estimates))
    change = fresh_change(bench_dir, baseline)
    if change is None:
        return Row(name, base_time, "-", None, MISSING_IN_HEAD)
    relative = mean_estimate(change)
    head_time = format_time(mean_estimate(head_estimates))
    return Row(name, base_time, head_time, relative, status_for(change_lower_bound(change), threshold))


def time_rows(criterion_dir: Path, baseline: str, threshold: float) -> list[Row]:
    rows = [time_row(path, baseline, threshold) for path in bench_dirs(criterion_dir, baseline)]
    return sorted(rows, key=lambda row: row.name)


def memory_entries(path: Path | None) -> dict[str, int]:
    if path is None or not path.exists():
        return {}
    report = load_json(path)
    entries = {name: value["peak_bytes"] for name, value in report.get("files", {}).items()}
    if SWEEP_ENTRY in report:
        entries[SWEEP_ENTRY] = report[SWEEP_ENTRY]
    return entries


def relative_growth(base: int, head: int) -> float:
    if base == 0:
        return 0.0 if head == 0 else math.inf
    return (head - base) / base


def memory_row(name: str, base: int | None, head: int | None, threshold: float) -> Row:
    if base is None:
        return Row(name, "-", format_bytes(head), None, MISSING_IN_BASE)
    if head is None:
        return Row(name, format_bytes(base), "-", None, MISSING_IN_HEAD)
    change = relative_growth(base, head)
    return Row(name, format_bytes(base), format_bytes(head), change, status_for(change, threshold))


def memory_rows(base_path: Path | None, head_path: Path | None, threshold: float) -> list[Row]:
    base = memory_entries(base_path)
    head = memory_entries(head_path)
    names = sorted(base.keys() | head.keys())
    return [memory_row(name, base.get(name), head.get(name), threshold) for name in names]


def print_table(title: str, rows: list[Row]) -> None:
    header = (title, "base", "head", "change", "status")
    cells = [header] + [
        (row.name, row.base, row.head, format_change(row.change), row.status) for row in rows
    ]
    widths = [max(len(cell[column]) for cell in cells) for column in range(len(header))]
    for cell in cells:
        name, base, head, change, status = cell
        print(
            f"{name:<{widths[0]}}  {base:>{widths[1]}}  {head:>{widths[2]}}"
            f"  {change:>{widths[3]}}  {status}"
        )
    print()


def regressions(rows: list[Row]) -> list[str]:
    return [f"{row.name} {format_change(row.change)}" for row in rows if row.status == REGRESSED]


def count_status(rows: list[Row], status: str) -> int:
    return sum(1 for row in rows if row.status == status)


def missing_summary(kind: str, rows: list[Row]) -> str | None:
    only_head = count_status(rows, MISSING_IN_BASE)
    only_base = count_status(rows, MISSING_IN_HEAD)
    if only_head == only_base == 0:
        return None
    return f"{kind}: {only_head} only in head, {only_base} only in base (not gated)"


def report_missing_files(paths: list[Path | None]) -> None:
    for path in paths:
        if path is not None and not path.exists():
            print(f"note: memory report {path} not found; its entries count as missing")


def pass_summary(time: list[Row], memory: list[Row], args: argparse.Namespace) -> str:
    parts = [f"{count_status(time, OK)} benchmarks within +{args.time_threshold:g}%"]
    if memory:
        parts.append(f"{count_status(memory, OK)} memory entries within +{args.memory_threshold:g}%")
    return ", ".join(parts)


def verdict(time: list[Row], memory: list[Row], args: argparse.Namespace) -> int:
    slow = regressions(time)
    grown = regressions(memory)
    notes = [missing_summary("benchmarks", time), missing_summary("memory entries", memory)]
    for note in filter(None, notes):
        print(f"note: {note}")
    if slow:
        print(f"time regressions confidently above +{args.time_threshold:g}%: {', '.join(slow)}")
    if grown:
        print(f"memory regressions above +{args.memory_threshold:g}%: {', '.join(grown)}")
    if slow or grown:
        print("verdict: FAIL")
        return 1
    print(f"verdict: PASS ({pass_summary(time, memory, args)})")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.criterion_dir.is_dir():
        print(f"error: criterion directory {args.criterion_dir} not found", file=sys.stderr)
        return 2
    time = time_rows(args.criterion_dir, args.baseline, args.time_threshold / 100)
    print_table("benchmark", time)
    memory = []
    if args.base_memory is not None or args.head_memory is not None:
        report_missing_files([args.base_memory, args.head_memory])
        memory = memory_rows(args.base_memory, args.head_memory, args.memory_threshold / 100)
        print_table("peak heap", memory)
    return verdict(time, memory, args)


if __name__ == "__main__":
    sys.exit(main())
