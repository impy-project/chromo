#!/usr/bin/env python3
"""
Helper script to generate f2py wrapper code for meson build system.
This replaces the f2py_add_module functionality from CMake.
"""

import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional


def setup_logging(module_name: str, output_dir: str):
    """Set up logging to file and console for f2py generation."""
    output_path = Path(output_dir)
    output_path.mkdir(
        parents=True, exist_ok=True
    )  # Create directory if it doesn't exist
    log_file = output_path / f"{module_name}-f2py_log.txt"

    # Clear any existing handlers
    logging.getLogger().handlers.clear()

    # Configure logging to write to both file and console
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="w"),
            # logging.StreamHandler(sys.stdout)
        ],
    )

    return logging.getLogger(__name__)


def generate_f2py_module(
    module_name: str,
    functions: list,
    sources: list,
    include: str,
    flags: str,
    output_dir: str,
    logger: Optional[logging.Logger] = None,
):
    """
    Generate f2py wrapper code for a Fortran module.

    Args:
        module_name: Name of the Python module to generate
        functions: List of function names to expose
        sources: List of source files to compile List of source files for interface generation
        include: String of include directive
        output_dir: Directory to place generated files
    """

    # Set up logging
    logger = logger or setup_logging(module_name, output_dir)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    include = include.split(" ") if include else []

    logger.info(f"Generating f2py wrapper for module: {module_name}")
    logger.info(f"Functions: {functions}")
    logger.info(f"Sources: {sources}")
    logger.info(f"Output dir: {output_path}")

    # Generate signature file
    pyf_file = output_path / f"{module_name}.pyf"

    # Add interface sources
    source_list = [
        str(src) for src in sources if src and not src.endswith(".c")
    ]  # Skip empty and C sources
    flags = [f.strip() for f in flags.split(" ") if f]  # Split and clean flags
    logger.info(f"Source files: {source_list}")
    logger.info(f"Include directories: {include}")
    logger.info(f"Flags: {flags}")
    # Run gfortran preprocessor to evaluate all includes and flags
    # Concatenate output into a single temporary file to be processed by f2py
    preprocessed_sources = output_path / f"{module_name}pp.f"
    logger.info(
        f"Running gfortran preprocessor: {' '.join(['gfortran', '-cpp', '-E', *flags, *include, *source_list])}"
    )
    cmd = ["gfortran", "-cpp", "-E", *flags, *include, *source_list]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, cwd=output_path
        )
        if result.stdout:
            with open(preprocessed_sources, "w") as f:
                f.write(result.stdout)
            logger.info(f"Preprocessed sources written to: {preprocessed_sources}")
        if result.stderr:
            logger.info(f"Preprocessing stderr: {result.stderr}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error preprocessing file: {e}")  # noqa: TRY400
        # logger.error(f"stdout: {e.stdout}")
        logger.error(f"stderr: {e.stderr}")  # noqa: TRY400
        return False

    cmd = [
        sys.executable,
        "-m",
        "numpy.f2py",
        "-h",
        str(pyf_file),
        "-m",
        module_name,
        " ".join(include),
        f"only: {' '.join(functions)} :",
        "--overwrite-signature",
        "--lower",
        str(preprocessed_sources.absolute()),
    ]
    logger.info(f"Running f2py signature generation: {' '.join(cmd)}")
    # Run f2py to generate signature from the parent directory
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, cwd=output_path.parent
        )
        logger.info(f"Generated signature file: {pyf_file}")
        if result.stdout:
            logger.info(f"f2py stdout: {result.stdout}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Error generating signature file: {e}")  # noqa: TRY400
        logger.error(f"stdout: {e.stdout}")  # noqa: TRY400
        logger.error(f"stderr: {e.stderr}")  # noqa: TRY400
        return False

    # Generate wrapper source files using f2py
    # Create expected output file paths
    wrapper_c = output_path / f"{module_name}module.c"
    wrapper_f = output_path / f"{module_name}-f2pywrappers.f"

    # Use f2py to generate the wrapper files without compiling
    cmd = [sys.executable, "-m", "numpy.f2py", str(pyf_file.absolute()), "--lower"]

    logger.info(f"Running f2py wrapper generation: {' '.join(cmd)}")

    try:
        # Run from the output directory
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, cwd=output_path
        )

        if result.stdout:
            logger.info(f"f2py stdout: {result.stdout}")

        # Check if the files were generated
        if not wrapper_c.exists():
            logger.warning(f"Expected wrapper C file not found: {wrapper_c}")
            # Look for any generated .c files
            c_files = list(output_path.glob("*.c"))
            if c_files:
                actual_c = c_files[0]
                logger.info(f"Found C file: {actual_c}, renaming to {wrapper_c}")
                actual_c.rename(wrapper_c)

        if not wrapper_f.exists():
            # Look for any generated wrapper .f files
            f_files = list(output_path.glob("*-f2pywrappers*.f"))
            if f_files:
                actual_f = f_files[0]
                logger.info(f"Found F file: {actual_f}, renaming to {wrapper_f}")
                actual_f.rename(wrapper_f)

        logger.info(f"Generated wrapper files: {wrapper_c}, {wrapper_f}")

        # Output the generated file paths for meson
        generated_files = {
            "module_c": str(wrapper_c) if wrapper_c.exists() else None,
            "wrapper_f": str(wrapper_f) if wrapper_f.exists() else None,
            "pyf": str(pyf_file),
        }
        logger.info("GENERATED_FILES: " + json.dumps(generated_files))

    except subprocess.CalledProcessError as e:
        logger.error(f"Error generating wrapper files: {e}")  # noqa: TRY400
        logger.error(f"stdout: {e.stdout}")  # noqa: TRY400
        logger.error(f"stderr: {e.stderr}")  # noqa: TRY400
        return False

    return True


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(  # noqa: T201
            "Usage: generate_f2py.py <module_name> <functions> <logging_dir>"
            " <includes> <output_dir> <sources>"
        )
        sys.exit(1)

    module_name = sys.argv[1]
    functions = (
        [f.strip() for f in sys.argv[2].split(",") if f.strip()] if sys.argv[2] else []
    )
    logging_dir = sys.argv[3]
    logger = setup_logging(module_name, logging_dir)
    includes = sys.argv[4]
    flags = sys.argv[5]
    output_dir = sys.argv[6]
    # Ther remaining arguments are the source files
    if len(sys.argv) < 7:
        logger.error("At least one source file must be provided.")
        sys.exit(1)
    sources = sys.argv[7:]
    # Check that all files in sources exist
    for src in sources:
        if not Path(src).exists():
            logger.error(f"Source file does not exist: {src}")
            assert Path(src).exists()
            sys.exit(1)
    # Set up logging for the main script

    assert output_dir is not None, "output_dir must be provided"
    assert module_name, "module_name must be provided"
    assert functions, "functions must be provided"
    assert sources, "sources must be provided"

    logger.info("Script arguments:")
    logger.info(f"  module_name: {module_name}")
    logger.info(f"  functions: {functions}")
    logger.info(f"  sources: {sources}")
    logger.info(f"  include_dirs: {includes}")
    logger.info(f"  output_dir: {output_dir}")
    logger.info(f"  flags: {flags}")

    success = generate_f2py_module(
        module_name=module_name,
        functions=functions,
        sources=sources,
        include=includes,
        output_dir=output_dir,
        flags=flags,
        logger=logger,
    )

    sys.exit(0 if success else 1)
