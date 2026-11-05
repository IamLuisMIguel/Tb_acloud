"""Build a reproducible candidate archive from an explicit public-file allowlist."""

import argparse
import gzip
import hashlib
import io
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PUBLIC_FILES = [
    "instructions.md", "task.md", "benchmark/inputs.json", "benchmark/requirements.txt",
    "benchmark/results.schema.json", "benchmark/results_format.md", "benchmark/validate_results.py",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    protected = [(ROOT / name).resolve() for name in PUBLIC_FILES]
    protected += list((ROOT / "author_only").glob("*.py"))
    protected += list((ROOT / "author_only").glob("*.json"))
    if args.output.resolve() in protected:
        parser.error("Output must not overwrite a benchmark source or reference file")
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w", format=tarfile.USTAR_FORMAT) as archive:
        for name in PUBLIC_FILES:
            data = (ROOT / name).read_bytes()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mode = 0o644
            info.mtime = 0
            info.uid = info.gid = 0
            archive.addfile(info, io.BytesIO(data))
    with args.output.open("wb") as destination:
        with gzip.GzipFile(filename="", mode="wb", fileobj=destination, mtime=0) as compressed:
            compressed.write(buffer.getvalue())
    with tarfile.open(args.output, "r:gz") as archive:
        assert archive.getnames() == PUBLIC_FILES
    print("Files:", ", ".join(PUBLIC_FILES))
    print("SHA-256:", hashlib.sha256(args.output.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
