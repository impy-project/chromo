# AGENTS.md

## Purpose
This document provides guidelines for AI coding agents working on this project. It outlines best practices, project-specific conventions, and important context to ensure high-quality, consistent contributions.

---

## Project Overview
This repository is a scientific computing project involving compiled model modules (Fortran/C/C++), Python wrappers, and extensive testing. The build system uses Meson, and compiled shared objects (`.so` files) are placed in the `build/` directory. Python tests are located in `tests/`.

---

## Key Conventions

- **Testing:**  
  - All tests are in the `tests/` directory and use `pytest`.
  - Tests may dynamically discover and import compiled modules from the build output.

- **Build Artifacts:**  
  - Compiled modules are found in `build/cp*/`.
  - Python modules are imported by extracting the base name from the `.so` file.

- **Symbol Checking:**  
  - Some tests (e.g., `test_exposed_symbols.py`) check for the presence of exported symbols in compiled modules.
  - For certain models (e.g., "phojet"), ignore symbols with specific prefixes as defined in the test.

- **Meson Build Integration:**  
  - Model function lists may be parsed from `meson.build`.
  - When updating model functions, ensure consistency between the build configuration and test expectations.

---

## Agent Workflow

1. **Context Gathering:**  
   - Always gather relevant context before making changes. Use semantic search or file search to locate definitions, usages, and related files.

2. **Editing Files:**  
   - Use concise, minimal edits. Represent unchanged code with comments (`# ...existing code...`).
   - Follow Python and project-specific style conventions.

3. **Testing:**  
   - After making changes, run the test suite (`pytest`) to ensure nothing is broken.
   - Address any test failures or errors before considering the task complete.

4. **Build System Awareness:**  
   - Be aware of the build system (Meson) and the location of build artifacts.
   - When working with compiled modules, ensure the correct paths and import logic are used.

5. **Documentation:**  
   - Update or create documentation as needed when adding new features or changing behavior.

---

## Special Notes

- **Dynamic Imports:**  
  - When importing modules from `.so` files, always update `sys.path` as needed and extract the module name correctly.

- **Error Handling:**  
  - Assert meaningful error messages, especially when expected files or symbols are missing.

---

## Example: Adding a New Model

1. Update `meson.build` with the new model and its function list.
2. Ensure the model is built and the `.so` file appears in `build/cp*/`.
3. Update or add tests in `tests/` to check for the new model's symbols.
4. Run all tests and verify success.

---
