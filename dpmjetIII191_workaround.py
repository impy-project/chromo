from pathlib import Path
import re
import sys
import shutil

this_program = Path(sys.argv[0]).name

PATTERN = r"INCLUDE '\(([A-Z]+)\)'"

if not sys.argv[1:]:  # for testing
    cmake_binary_dir = "build"
    inc_dir = "src/fortran/dpmjetIII-19.1/include"
    src = []
    dirs = [
        "src/fortran/dpmjetIII-19.1/src/phojet",
        "src/fortran/dpmjetIII-19.1/src/pythia",
        "src/fortran/dpmjetIII-19.1/src/dpmjet",
        "src/fortran/dpmjetIII-19.1/common",
    ]
    for d in dirs:
        for fn in Path(d).rglob("*.f"):
            src.append(fn)
else:
    with Path(sys.argv[3]).open() as f:
        src = f.read().split()
    cmake_binary_dir, inc_dir = sys.argv[1:3]

src = [Path(x) for x in src]
inc_dir = Path(inc_dir)
cmake_binary_dir = Path(cmake_binary_dir)

dst_include = cmake_binary_dir / "include"
dst_include.mkdir(parents=True, exist_ok=True)
for fn in inc_dir.rglob("(*)"):
    shutil.copy(fn, dst_include / fn.name[1:-1])

needs_mod = []
unchanged = []
for fn in src:
    with fn.open() as f:
        c = f.read()
    m = re.search(PATTERN, c, re.MULTILINE)
    if m:
        needs_mod.append(fn)
    else:
        unchanged.append(fn)

# ensure filenames are not colliding
assert len(set([x.name for x in needs_mod])) == len(needs_mod)

modded = []
for fn in needs_mod:
    fn2 = cmake_binary_dir / fn.name
    modded.append(fn2)
    shutil.copy(fn, fn2)

for fn in modded:
    with fn.open() as f:
        c = f.read()
    c = re.sub(PATTERN, r"INCLUDE '\1'", c)
    with fn.open("w") as f:
        f.write(c)
    sys.stderr.write(f"-- {this_program}: generated {fn}\n")

out = list(map(str, unchanged))
out += list(map(str, modded))


with open(sys.argv[3], "w") as f:
    f.write(";".join(out).replace("\\", "/"))
