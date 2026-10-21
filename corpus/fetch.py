import argparse
import sys
import time
import tomllib
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CANDIDATES_PATH = REPO_ROOT / "corpus" / "candidates.toml"
DATA_DIR = REPO_ROOT / "corpus" / "data"
USER_AGENT = "mmcif-corpus-fetch/1.0 (differential test corpus)"
REQUEST_TIMEOUT_SECONDS = 60
DEFAULT_DELAY_SECONDS = 0.5
MAX_ATTEMPTS = 3

SOURCES = (
    ("cif", "https://files.rcsb.org/download/{upper}.cif"),
    ("bcif", "https://models.rcsb.org/{lower}.bcif"),
)


def load_candidate_ids(path: Path) -> list[str]:
    with open(path, "rb") as handle:
        classes = tomllib.load(handle)
    ids = {identifier.lower() for entry in classes.values() for identifier in entry["ids"]}
    return sorted(ids)


def download(url: str, destination: Path) -> int:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
                payload = response.read()
            partial = destination.with_suffix(destination.suffix + ".part")
            partial.write_bytes(payload)
            partial.replace(destination)
            return len(payload)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise
            last_error = error
        except (urllib.error.URLError, TimeoutError) as error:
            last_error = error
        time.sleep(attempt)
    raise last_error


def fetch_one(identifier: str, extension: str, template: str, delay: float) -> tuple[str, int]:
    destination = DATA_DIR / f"{identifier.lower()}.{extension}"
    if destination.exists():
        return "present", 0
    url = template.format(upper=identifier.upper(), lower=identifier.lower())
    size = download(url, destination)
    time.sleep(delay)
    return "downloaded", size


def fetch_all(ids: list[str], delay: float) -> tuple[list[tuple[str, str, str]], int]:
    failures = []
    downloaded_bytes = 0
    for identifier in ids:
        for extension, template in SOURCES:
            try:
                status, size = fetch_one(identifier, extension, template, delay)
            except (urllib.error.URLError, TimeoutError) as error:
                failures.append((identifier, extension, str(error)))
                print(f"FAIL {identifier}.{extension}: {error}", file=sys.stderr)
                continue
            downloaded_bytes += size
            if status == "downloaded":
                print(f"fetched {identifier}.{extension} ({size} bytes)")
    return failures, downloaded_bytes


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ids", nargs="*", help="fetch only these ids instead of all candidates")
    parser.add_argument("--delay", type=float, default=DEFAULT_DELAY_SECONDS)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    ids = [i.lower() for i in arguments.ids] if arguments.ids else load_candidate_ids(CANDIDATES_PATH)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    failures, downloaded_bytes = fetch_all(ids, arguments.delay)
    print(f"{len(ids)} ids, {downloaded_bytes} bytes downloaded, {len(failures)} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
