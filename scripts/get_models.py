import argparse
import pathlib

try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # For older Python with 'tomli' installed


def get_models(enabled):
    pyproject = pathlib.Path(__file__).resolve().parents[1] / "pyproject.toml"

    with pyproject.open("rb") as f:
        data = tomllib.load(f)

    enabled_models = data["tool"]["chromo"]["enabled-models"]
    disabled_models = data["tool"]["chromo"]["disabled-models"]
    if enabled:
        return enabled_models
    return disabled_models


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get enabled and disabled models from pyproject.toml"
    )
    parser.add_argument("--enabled", action="store_true", help="Print enabled models")
    parser.add_argument("--disabled", action="store_true", help="Print disabled models")

    args = parser.parse_args()

    if args.enabled:
        print("\n".join(get_models(True)))  # noqa: T201
    if args.disabled:
        print("\n".join(get_models(False)))  # noqa: T201
    if not args.enabled and not args.disabled:
        print(  # noqa: T201
            "Please specify --enabled or --disabled to see the respective models."
        )
