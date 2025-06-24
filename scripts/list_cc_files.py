import sys
from pathlib import Path

src_dir = Path(sys.argv[1])
print("\n".join(str(p) for p in src_dir.rglob("*.cc")))  # noqa: T201
