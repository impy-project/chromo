#!/usr/bin/env python3
# save as fix_utf8.py ;  chmod +x fix_utf8.py
import argparse
import sys
from pathlib import Path

from charset_normalizer import CharsetMatch, from_path


def main():
    ap = argparse.ArgumentParser(
        description="Convert every text file to UTF-8 in-place"
    )
    ap.add_argument("root", help="Top-level directory to process")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    if not root.is_dir():
        sys.exit(f"❌ {root} is not a directory")

    for file in root.rglob("*"):
        if not file.is_file():
            continue
        match: CharsetMatch | None = from_path(file).best()
        if match and match.encoding.lower() != "utf_8":
            file.write_bytes(match.output())
            print(  # noqa: T201
                f"✔ {file.relative_to(root)}  ({match.encoding} → UTF-8)"
            )


if __name__ == "__main__":
    main()
