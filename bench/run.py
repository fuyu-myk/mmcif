import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import IO, NoReturn

REPO_ROOT = Path(__file__).resolve().parent.parent
GATE_SCRIPT = REPO_ROOT / "bench" / "gate.py"
MEMORY_TOOL = Path("tools") / "mmcif-mem"
BASELINE = "base"
DEFAULT_BASE_REF = "main"
DEFAULT_WARM_UP_SECONDS = 1.0
DEFAULT_MEASUREMENT_SECONDS = 3.0
DEFAULT_TIME_THRESHOLD_PERCENT = 15.0
DEFAULT_MEMORY_THRESHOLD_PERCENT = 5.0
DESCRIPTION = """\
Benchmark a base commit against the working tree and gate regressions.

steps:
  1. check out REF in a temporary git worktree that shares this CARGO_TARGET_DIR
  2. base:  cargo bench --workspace --bench '*' -- --save-baseline base
  3. base:  mmcif-mem --corpus, built from this tree's tools/mmcif-mem
  4. head:  cargo bench --workspace --bench '*' -- --baseline-lenient base
  5. head:  mmcif-mem --corpus
  6. bench/gate.py fails when a mean time grows more than --time-threshold percent
     or a peak heap grows more than --memory-threshold percent

Results land in <target>/bench, wiped first: criterion/, memory-base.json and
memory-head.json. The corpus comes from MMCIF_CORPUS_DIR or corpus/data; fetch it
with harness/.venv/bin/python corpus/fetch.py.

manual comparison without a worktree:
  cargo bench --workspace --bench '*' -- --save-baseline base
  cargo bench --workspace --bench '*' -- --baseline base     (after changing code)
  python3 bench/gate.py
"""


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base", default=DEFAULT_BASE_REF, metavar="REF", help="commit to compare against (default %(default)s)")
    parser.add_argument("--quick", action="store_true", help="pass --quick to criterion instead of fixed times")
    parser.add_argument("--warm-up-time", type=float, default=DEFAULT_WARM_UP_SECONDS, metavar="SECONDS")
    parser.add_argument("--measurement-time", type=float, default=DEFAULT_MEASUREMENT_SECONDS, metavar="SECONDS")
    parser.add_argument("--time-threshold", type=float, default=DEFAULT_TIME_THRESHOLD_PERCENT, metavar="PERCENT")
    parser.add_argument("--memory-threshold", type=float, default=DEFAULT_MEMORY_THRESHOLD_PERCENT, metavar="PERCENT")
    return parser.parse_args(argv)


def fail(message: str) -> NoReturn:
    print(f"error: {message}", file=sys.stderr)
    sys.exit(2)


def say(message: str) -> None:
    print(message, flush=True)


def run(command: list[str], cwd: Path, env: dict[str, str], stdout: IO[str] | None = None) -> None:
    subprocess.run(command, cwd=cwd, env=env, stdout=stdout, check=True)


def corpus_dir() -> Path:
    path = Path(os.environ.get("MMCIF_CORPUS_DIR", REPO_ROOT / "corpus" / "data")).resolve()
    if not path.is_dir():
        fail(f"corpus directory {path} is missing: run harness/.venv/bin/python corpus/fetch.py or set MMCIF_CORPUS_DIR")
    return path


def target_dir() -> Path:
    return Path(os.environ.get("CARGO_TARGET_DIR", REPO_ROOT / "target")).resolve()


def bench_env(target: Path, criterion: Path, corpus: Path) -> dict[str, str]:
    return {**os.environ, "CARGO_TARGET_DIR": str(target), "CRITERION_HOME": str(criterion), "MMCIF_CORPUS_DIR": str(corpus)}


def criterion_timing(args: argparse.Namespace) -> list[str]:
    if args.quick:
        return ["--quick"]
    return ["--warm-up-time", f"{args.warm_up_time:g}", "--measurement-time", f"{args.measurement_time:g}"]


def has_bench_targets(tree: Path, env: dict[str, str]) -> bool:
    command = ["cargo", "metadata", "--no-deps", "--format-version", "1"]
    completed = subprocess.run(command, cwd=tree, env=env, stdout=subprocess.PIPE, check=True, text=True)
    packages = json.loads(completed.stdout)["packages"]
    return any("bench" in target["kind"] for package in packages for target in package["targets"])


def run_benches(tree: Path, env: dict[str, str], criterion_args: list[str]) -> None:
    if not has_bench_targets(tree, env):
        say(f"{tree} has no bench targets; its benchmarks count as missing")
        return
    run(["cargo", "bench", "--workspace", "--bench", "*", "--", "--noplot", *criterion_args], tree, env)


def measure_memory(tree: Path, env: dict[str, str], report: Path) -> None:
    with open(report, "w", encoding="utf-8") as output:
        run(["cargo", "run", "--release", "-p", "mmcif-mem", "--", "--corpus"], tree, env, stdout=output)


def overlay_memory_tool(tree: Path) -> None:
    destination = tree / MEMORY_TOOL
    shutil.rmtree(destination, ignore_errors=True)
    shutil.copytree(REPO_ROOT / MEMORY_TOOL, destination)


def measure_base_memory(tree: Path, env: dict[str, str], report: Path) -> None:
    overlay_memory_tool(tree)
    try:
        measure_memory(tree, env, report)
    except subprocess.CalledProcessError as error:
        report.unlink(missing_ok=True)
        say(f"warning: mmcif-mem failed on the base ({error}); its memory entries count as missing")


def bench_base(ref: str, env: dict[str, str], timing: list[str], report: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="mmcif-bench-") as scratch:
        tree = Path(scratch) / "base"
        run(["git", "worktree", "add", "--detach", str(tree), ref], REPO_ROOT, env)
        try:
            say(f"\n==> criterion benches of base {ref}")
            run_benches(tree, env, ["--save-baseline", BASELINE, *timing])
            say(f"\n==> mmcif-mem of base {ref}")
            measure_base_memory(tree, env, report)
        finally:
            run(["git", "worktree", "remove", "--force", str(tree)], REPO_ROOT, env)


def bench_head(env: dict[str, str], timing: list[str], report: Path) -> None:
    say("\n==> criterion benches of the working tree")
    run_benches(REPO_ROOT, env, ["--baseline-lenient", BASELINE, *timing])
    say("\n==> mmcif-mem of the working tree")
    measure_memory(REPO_ROOT, env, report)


def gate(args: argparse.Namespace, out: Path) -> int:
    say("\n==> regression gate")
    command = [
        sys.executable, str(GATE_SCRIPT),
        "--criterion-dir", str(out / "criterion"),
        "--baseline", BASELINE,
        "--time-threshold", f"{args.time_threshold:g}",
        "--memory-threshold", f"{args.memory_threshold:g}",
        "--base-memory", str(out / "memory-base.json"),
        "--head-memory", str(out / "memory-head.json"),
    ]
    return subprocess.run(command, cwd=REPO_ROOT).returncode


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    target = target_dir()
    out = target / "bench"
    env = bench_env(target, out / "criterion", corpus_dir())
    timing = criterion_timing(args)
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    try:
        bench_base(args.base, env, timing, out / "memory-base.json")
        bench_head(env, timing, out / "memory-head.json")
    except subprocess.CalledProcessError as error:
        fail(f"{' '.join(error.cmd)} exited with status {error.returncode}")
    return gate(args, out)


if __name__ == "__main__":
    sys.exit(main())
