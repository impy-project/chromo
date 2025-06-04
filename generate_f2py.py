#!/usr/bin/env python3
"""
Helper script to generate f2py wrapper code for meson build system.
This replaces the f2py_add_module functionality from CMake.
"""

import subprocess
import sys
import json
from pathlib import Path


def fix_pyf_structure(pyf_file: Path, module_name: str):
    """
    Fix the structure of a .pyf file to wrap all content in a python module block.
    This is needed because f2py sometimes generates individual subroutine blocks
    instead of proper python module structure.
    """
    print(f"Fixing .pyf file structure: {pyf_file}")

    # Read the original content
    with open(pyf_file, "r") as f:
        content = f.read()

    # Check if it already has a python module block
    if "python module " in content:
        print("File already has python module structure")
        return

    # Split content into lines
    lines = content.strip().split("\n")

    # Find the header comments and keep them
    header_lines = []
    content_lines = []
    in_header = True

    for line in lines:
        if in_header and (
            line.startswith("!") or line.startswith("//") or line.strip() == ""
        ):
            header_lines.append(line)
        else:
            in_header = False
            content_lines.append(line)

    # Create the new content with proper module structure
    new_content = []

    # Add header
    new_content.extend(header_lines)
    new_content.append("")

    # Add python module block
    new_content.append(f"python module {module_name} ! in")
    new_content.append(f"    interface  ! in :{module_name}")

    # Add all the original content indented
    for line in content_lines:
        if line.strip():
            new_content.append("        " + line)
        else:
            new_content.append("")

    # Close the blocks
    new_content.append("    end interface")
    new_content.append(f"end python module {module_name}")
    new_content.append("")

    # Write the fixed content back
    with open(pyf_file, "w") as f:
        f.write("\n".join(new_content))

    print(f"Fixed .pyf file structure for module: {module_name}")


def generate_f2py_module(
    module_name: str,
    functions: list,
    sources: list,
    interface_sources: list = None,
    include_dirs: list = None,
    output_dir: str = ".",
):
    """
    Generate f2py wrapper code for a Fortran module.

    Args:
        module_name: Name of the Python module to generate
        functions: List of function names to expose
        sources: List of source files to compile
        interface_sources: List of source files for interface generation
        include_dirs: List of include directories
        output_dir: Directory to place generated files
    """

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # If interface_sources not provided, use sources
    if interface_sources is None:
        interface_sources = sources

    print(f"Generating f2py wrapper for module: {module_name}")
    print(f"Functions: {functions}")
    print(f"Sources: {sources}")
    print(f"Output dir: {output_path}")

    # Generate signature file
    pyf_file = output_path / f"{module_name}.pyf"

    cmd = [
        sys.executable,
        "-m",
        "numpy.f2py",
        "-h",
        str(pyf_file),
        "--overwrite-signature",
    ]

    if functions:
        cmd.extend(["only:"] + functions + [":"])

    # Add include directories
    if include_dirs:
        for inc_dir in include_dirs:
            cmd.extend(["-I", str(inc_dir)])

    # Add interface sources
    for src in interface_sources:
        if src:  # Skip empty sources
            cmd.append(str(src))

    print(f"Running f2py signature generation: {' '.join(cmd)}")

    # Run f2py to generate signature from the parent directory
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, cwd=output_path.parent
        )
        print(f"Generated signature file: {pyf_file}")
        if result.stdout:
            print(f"f2py stdout: {result.stdout}")

        # Fix the .pyf file structure - wrap everything in a python module block
        fix_pyf_structure(pyf_file, module_name)

    except subprocess.CalledProcessError as e:
        print(f"Error generating signature file: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return False

    # Generate wrapper source files using f2py
    # Create expected output file paths
    wrapper_c = output_path / f"{module_name}module.c"
    wrapper_f = output_path / f"{module_name}-f2pywrappers.f"

    # Use f2py to generate the wrapper files without compiling
    cmd = [sys.executable, "-m", "numpy.f2py", str(pyf_file.absolute()), "--lower"]

    print(f"Running f2py wrapper generation: {' '.join(cmd)}")

    try:
        # Run from the output directory
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, cwd=output_path
        )

        if result.stdout:
            print(f"f2py stdout: {result.stdout}")

        # Check if the files were generated
        if not wrapper_c.exists():
            print(f"Warning: Expected wrapper C file not found: {wrapper_c}")
            # Look for any generated .c files
            c_files = list(output_path.glob("*.c"))
            if c_files:
                actual_c = c_files[0]
                print(f"Found C file: {actual_c}, renaming to {wrapper_c}")
                actual_c.rename(wrapper_c)

        if not wrapper_f.exists():
            # Look for any generated wrapper .f files
            f_files = list(output_path.glob("*-f2pywrappers*.f"))
            if f_files:
                actual_f = f_files[0]
                print(f"Found F file: {actual_f}, renaming to {wrapper_f}")
                actual_f.rename(wrapper_f)

        print(f"Generated wrapper files: {wrapper_c}, {wrapper_f}")

        # Output the generated file paths for meson
        generated_files = {
            "module_c": str(wrapper_c) if wrapper_c.exists() else None,
            "wrapper_f": str(wrapper_f) if wrapper_f.exists() else None,
            "pyf": str(pyf_file),
        }
        print("GENERATED_FILES:", json.dumps(generated_files))

    except subprocess.CalledProcessError as e:
        print(f"Error generating wrapper files: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return False

    return True


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(
            "Usage: generate_f2py.py <module_name> <functions> <sources> [interface_sources] [include_dirs] [output_dir]"
        )
        sys.exit(1)

    module_name = sys.argv[1]
    functions = (
        [f.strip() for f in sys.argv[2].split(",") if f.strip()] if sys.argv[2] else []
    )
    sources = [
        s.strip() for s in sys.argv[3].split(",") if s.strip() and s.strip() != ""
    ]
    interface_sources = (
        [s.strip() for s in sys.argv[4].split(",") if s.strip() and s.strip() != ""]
        if len(sys.argv) > 4 and sys.argv[4]
        else None
    )
    include_dirs = (
        [d.strip() for d in sys.argv[5].split(",") if d.strip() and d.strip() != ""]
        if len(sys.argv) > 5 and sys.argv[5]
        else None
    )
    output_dir = sys.argv[6] if len(sys.argv) > 6 else "."

    print("Script arguments:")
    print(f"  module_name: {module_name}")
    print(f"  functions: {functions}")
    print(f"  sources: {sources}")
    print(f"  interface_sources: {interface_sources}")
    print(f"  include_dirs: {include_dirs}")
    print(f"  output_dir: {output_dir}")

    success = generate_f2py_module(
        module_name=module_name,
        functions=functions,
        sources=sources,
        interface_sources=interface_sources,
        include_dirs=include_dirs,
        output_dir=output_dir,
    )

    sys.exit(0 if success else 1)
