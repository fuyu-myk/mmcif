import argparse
import gzip
import hashlib
import http.client
import json
import sys
import time
import tomllib
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CANDIDATES_PATH = REPO_ROOT / "corpus" / "candidates.toml"
DATA_DIR = REPO_ROOT / "corpus" / "data"
USER_AGENT = "mmcif-corpus-fetch/1.0 (differential test corpus)"
REQUEST_TIMEOUT_SECONDS = 60
DEFAULT_DELAY_SECONDS = 0.5
MAX_ATTEMPTS = 3
HASH_CHUNK_BYTES = 1 << 20
TRANSFER_ERRORS = (OSError, http.client.HTTPException, EOFError)

SOURCES = (
    ("cif", "https://files.rcsb.org/download/{upper}.cif.gz", True),
    ("bcif", "https://models.rcsb.org/{lower}.bcif", False),
)


def load_candidate_ids(path: Path) -> list[str]:
    with open(path, "rb") as handle:
        classes = tomllib.load(handle)
    ids = {identifier.lower() for entry in classes.values() for identifier in entry["ids"]}
    return sorted(ids)


def candidate_ids(corpus: dict) -> list[str]:
    pooled = [identifier for pool in corpus["pools"].values() for identifier in pool["ids"]]
    return sorted({identifier.lower() for identifier in corpus["baseline"] + pooled})


def corpus_ids(corpus: dict) -> list[str]:
    return sorted(set(candidate_ids(corpus)) | set(corpus.get("pins", {})))


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(HASH_CHUNK_BYTES):
            digest.update(chunk)
    return digest.hexdigest()


def changed_files(pins: dict[str, dict[str, str]]) -> list[str]:
    changed = []
    for identifier, pin in sorted(pins.items()):
        for extension, _, _ in SOURCES:
            path = DATA_DIR / f"{identifier}.{extension}"
            if path.exists() and sha256_of(path) != pin[extension]:
                changed.append(path.name)
    return changed


def download(url: str, destination: Path, compressed: bool) -> int:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
                payload = response.read()
            if compressed:
                payload = gzip.decompress(payload)
            partial = destination.with_suffix(destination.suffix + ".part")
            partial.write_bytes(payload)
            partial.replace(destination)
            return len(payload)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise
            last_error = error
        except TRANSFER_ERRORS as error:
            last_error = error
        time.sleep(attempt)
    raise last_error


def fetch_one(identifier: str, extension: str, template: str, compressed: bool, delay: float) -> tuple[str, int]:
    destination = DATA_DIR / f"{identifier.lower()}.{extension}"
    if destination.exists():
        return "present", 0
    url = template.format(upper=identifier.upper(), lower=identifier.lower())
    size = download(url, destination, compressed)
    time.sleep(delay)
    return "downloaded", size


def fetch_entry(identifier: str, delay: float) -> tuple[list[tuple[str, str, str]], int]:
    failures = []
    downloaded_bytes = 0
    for extension, template, compressed in SOURCES:
        try:
            status, size = fetch_one(identifier, extension, template, compressed, delay)
        except TRANSFER_ERRORS as error:
            failures.append((identifier, extension, str(error)))
            print(f"FAIL {identifier}.{extension}: {error}", file=sys.stderr)
            continue
        downloaded_bytes += size
        if status == "downloaded":
            print(f"fetched {identifier}.{extension} ({size} bytes)")
    return failures, downloaded_bytes


def fetch_all(ids: list[str], delay: float, workers: int) -> tuple[list[tuple[str, str, str]], int]:
    pool = ThreadPoolExecutor(workers)
    try:
        results = list(pool.map(lambda identifier: fetch_entry(identifier, delay), ids))
    except KeyboardInterrupt:
        pool.shutdown(wait=False, cancel_futures=True)
        raise
    pool.shutdown()
    failures = [failure for entry_failures, _ in results for failure in entry_failures]
    return failures, sum(size for _, size in results)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ids", nargs="*", help="fetch only these ids instead of all candidates")
    parser.add_argument("--corpus", type=Path, help="fetch every id in this corpus file, such as corpus/full.json")
    parser.add_argument("--delay", type=float, default=DEFAULT_DELAY_SECONDS)
    parser.add_argument("--workers", type=int, default=1)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    pins = {}
    if arguments.ids:
        ids = sorted({i.lower() for i in arguments.ids})
    elif arguments.corpus:
        corpus = json.loads(arguments.corpus.read_text())
        ids = corpus_ids(corpus)
        pins = corpus.get("pins", {})
    else:
        ids = load_candidate_ids(CANDIDATES_PATH)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    failures, downloaded_bytes = fetch_all(ids, arguments.delay, arguments.workers)
    changed = changed_files(pins)
    for name in changed:
        print(f"CHANGED {name} differs from its pin in {arguments.corpus}")
    if changed:
        print("run corpus/select_entries.py classify to pin the new files, and commit the result")
    print(f"{len(ids)} ids, {downloaded_bytes} bytes downloaded, {len(failures)} failures, {len(changed)} changed")
    return 1 if failures or changed else 0


if __name__ == "__main__":
    sys.exit(main())
