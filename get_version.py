import pathlib
import importlib

tomllib = importlib.import_module(
    "tomllib" if importlib.util.find_spec("tomllib") else "tomli"
)

pyproject = pathlib.Path(__file__).resolve().parents[0] / "pyproject.toml"

with pyproject.open("rb") as f:
    data = tomllib.load(f)

version = data["project"]["version"]
print(version)
