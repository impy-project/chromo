#!/usr/bin/env python3
"""
Windows Build Validation Test Script

This script helps validate the Windows build fixes by checking:
1. Command line lengths in build.ninja
2. Source list files are generated
3. Extension modules are built
4. Modules can be imported

Usage:
    python test_windows_build.py [--phase N]

Examples:
    python test_windows_build.py           # Check current build
    python test_windows_build.py --phase 2 # Run Phase 2 validation
"""

import argparse
import subprocess
import sys
from pathlib import Path


def check_command_line_lengths(build_dir="build"):
    """Check that custom_target commands are under Windows limit."""
    build_ninja = Path(build_dir) / "build.ninja"

    if not build_ninja.exists():
        print(f"❌ {build_ninja} not found. Run build first.")
        return False

    print("\n🔍 Checking command line lengths in build.ninja...")

    max_length = 0
    max_line = ""

    with open(build_ninja, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if 'generate_f2py.py' in line:
                length = len(line)
                if length > max_length:
                    max_length = length
                    max_line = line[:100] + "..." if len(line) > 100 else line

    WINDOWS_LIMIT = 8191
    SAFE_THRESHOLD = 8000

    print(f"   Max command length: {max_length} characters")

    if max_length > WINDOWS_LIMIT:
        print(f"   ❌ FAIL: Exceeds Windows limit of {WINDOWS_LIMIT} chars!")
        print(f"   Command: {max_line}")
        return False
    elif max_length > SAFE_THRESHOLD:
        print(f"   ⚠️  WARNING: Close to limit ({SAFE_THRESHOLD} threshold)")
        return True
    else:
        print(f"   ✅ PASS: Well under Windows limit")
        return True


def check_source_list_files(build_dir="build"):
    """Check that source list files were generated."""
    print("\n🔍 Checking for source list files...")

    build_path = Path(build_dir)
    source_lists = list(build_path.glob("*_sources.txt"))

    if not source_lists:
        print("   ⚠️  No source list files found (expected for models with few files)")
        return True

    print(f"   Found {len(source_lists)} source list files:")
    for sl in sorted(source_lists):
        line_count = len(sl.read_text().strip().split('\n'))
        print(f"   ✅ {sl.name}: {line_count} source files")

    # Check DPMJET files specifically (should have 700+ lines)
    dpmjet_lists = [sl for sl in source_lists if 'dpmjet' in sl.name.lower()]
    for dl in dpmjet_lists:
        lines = len(dl.read_text().strip().split('\n'))
        if '191' in dl.name or '193' in dl.name:
            if lines > 600:
                print(f"   ✅ {dl.name}: {lines} files (DPMJET large model)")
            else:
                print(f"   ⚠️  {dl.name}: Only {lines} files (expected 700+)")

    return True


def check_extension_modules(build_dir="build"):
    """Check that extension modules were built."""
    print("\n🔍 Checking for built extension modules...")

    build_path = Path(build_dir)

    # Find the cpXX directory
    cp_dirs = list(build_path.glob("cp*"))
    if not cp_dirs:
        print("   ❌ No cp* directory found in build/")
        return False

    models_dir = cp_dirs[0] / "chromo" / "models"
    if not models_dir.exists():
        print(f"   ❌ Models directory not found: {models_dir}")
        return False

    # Find .pyd files (Windows) or .so files (Unix)
    extensions = list(models_dir.glob("*.pyd")) + list(models_dir.glob("*.so"))

    if not extensions:
        print(f"   ❌ No extension modules found in {models_dir}")
        return False

    print(f"   Found {len(extensions)} extension modules:")
    for ext in sorted(extensions):
        size_kb = ext.stat().st_size / 1024
        print(f"   ✅ {ext.name} ({size_kb:.1f} KB)")

    return True


def check_imports():
    """Try to import the chromo models."""
    print("\n🔍 Testing module imports...")

    try:
        import chromo.models as models
        print("   ✅ chromo.models package imported")
    except ImportError as e:
        print(f"   ❌ Failed to import chromo.models: {e}")
        return False

    # Try to import specific models
    model_names = [
        '_sib21', '_sib23d', '_qgs2_04', '_sophia', '_pythia6',
        '_eposlhc', '_dpmjet191', '_dpmjet193', '_dpmjet307', '_urqmd34'
    ]

    imported = []
    not_built = []
    failed = []

    for model in model_names:
        try:
            mod = getattr(models, model, None)
            if mod is None:
                not_built.append(model)
            else:
                imported.append(model)
        except Exception as e:
            failed.append((model, str(e)))

    if imported:
        print(f"   ✅ Successfully imported {len(imported)} models:")
        for m in imported:
            print(f"      • {m}")

    if not_built:
        print(f"   ⚠️  {len(not_built)} models not built (disabled in config):")
        for m in not_built:
            print(f"      • {m}")

    if failed:
        print(f"   ❌ {len(failed)} models failed to import:")
        for m, err in failed:
            print(f"      • {m}: {err}")
        return False

    return True


def check_response_file_usage(build_dir="build"):
    """Check that commands use --source-file-list."""
    print("\n🔍 Checking response file usage...")

    build_ninja = Path(build_dir) / "build.ninja"

    if not build_ninja.exists():
        print(f"   ❌ {build_ninja} not found")
        return False

    with open(build_ninja, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    f2py_commands = [line for line in content.split('\n') if 'generate_f2py.py' in line]

    if not f2py_commands:
        print("   ❌ No generate_f2py.py commands found")
        return False

    response_file_count = sum('--source-file-list' in cmd for cmd in f2py_commands)

    print(f"   Total generate_f2py.py commands: {len(f2py_commands)}")
    print(f"   Using --source-file-list: {response_file_count}")

    if response_file_count == len(f2py_commands):
        print("   ✅ All commands use response files")
        return True
    elif response_file_count > 0:
        print("   ✅ Some commands use response files (expected)")
        return True
    else:
        print("   ⚠️  No commands use response files (unexpected)")
        return False


def run_phase_validation(phase):
    """Run validation for a specific phase."""
    phase_tests = {
        2: ["test_sibyll.py", "-k", "Sibyll21"],
        3: ["test_sibyll.py"],
        4: ["test_qgsjet.py"],
        5: ["-k", "sophia or pythia6"],
        6: ["test_epos.py"],
        7: ["test_dpmjetIII.py"],
        8: ["-k", "urqmd"],
    }

    if phase not in phase_tests:
        print(f"❌ Unknown phase: {phase}")
        return False

    print(f"\n🧪 Running Phase {phase} tests...")

    test_args = ["python", "-m", "pytest", "-n", "0", "-vv"] + phase_tests[phase]

    print(f"   Command: {' '.join(test_args)}")

    result = subprocess.run(test_args)

    if result.returncode == 0:
        print(f"\n✅ Phase {phase} tests passed!")
        return True
    else:
        print(f"\n❌ Phase {phase} tests failed!")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Validate Windows build fixes for chromo"
    )
    parser.add_argument(
        "--phase",
        type=int,
        choices=[2, 3, 4, 5, 6, 7, 8],
        help="Run validation for specific phase"
    )
    parser.add_argument(
        "--build-dir",
        default="build",
        help="Build directory (default: build)"
    )

    args = parser.parse_args()

    print("=" * 70)
    print("Windows Build Validation for chromo")
    print("=" * 70)

    all_passed = True

    # Always run these checks
    all_passed &= check_command_line_lengths(args.build_dir)
    all_passed &= check_source_list_files(args.build_dir)
    all_passed &= check_response_file_usage(args.build_dir)
    all_passed &= check_extension_modules(args.build_dir)
    all_passed &= check_imports()

    # Run phase-specific tests if requested
    if args.phase:
        all_passed &= run_phase_validation(args.phase)

    print("\n" + "=" * 70)
    if all_passed:
        print("✅ All validation checks passed!")
        print("=" * 70)
        return 0
    else:
        print("❌ Some validation checks failed")
        print("=" * 70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
