import pathlib

try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # For older Python with 'tomli' installed

pyproject = pathlib.Path(__file__).resolve().parents[1] / "pyproject.toml"

with pyproject.open("rb") as f:
    data = tomllib.load(f)

version = data["project"]["version"]
print(version)  # noqa: T201
