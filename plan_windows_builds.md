
## Phase 6: Add EPOS Models

**Goal:** Test multi-file models with moderate complexity.

### Configuration Change:
Add to `enabled-models`:
```toml
enabled-models = [
    # ... previous models ...
    "eposlhc",
    "eposlhcr",
]
```

### Commands:
```bash
rm -rf build/
pip install --no-build-isolation -v -e .[test]
python -m pytest -n 32 -vv -k "epos"
```

### Success Criteria:
- [ ] Both EPOS modules built
- [ ] EPOS tests pass

---

## Phase 7: Add DPMJET Models 🔥 CRITICAL TEST

**Goal:** Validate the command line length fix with 700+ file models.

This is the critical phase that tests the response file implementation!

### Phase 7a: DPMJET 307 (Simple)
Add to `enabled-models`:
```toml
enabled-models = [
    # ... previous models ...
    "dpmjet_phojet307",
]
```

```bash
rm -rf build/
pip install --no-build-isolation -v -e .[test]
python -m pytest -n 0 -vv -k "307"
```

### Phase 7b: DPMJET 191 (700+ files)
Add to `enabled-models`:
```toml
    "dpmjet_phojet191",
```

```bash
rm -rf build/
pip install --no-build-isolation -v -e .[test]

# Check command line lengths are under limit
grep "CUSTOM_COMMAND" build/build.ninja | awk '{print length, $0}' | sort -rn | head -5

# Verify all commands are under 8000 characters (should be ~200)
python -m pytest -n 0 -vv -k "191"
```

### Phase 7c: DPMJET 193
Add to `enabled-models`:
```toml
    "dpmjet_phojet193",
```

```bash
rm -rf build/
pip install --no-build-isolation -v -e .[test]
python -m pytest -n 0 -vv -k "193"
```

### Success Criteria:
- [ ] All DPMJET modules built without command line length errors
- [ ] Command lines in build.ninja are under 8,000 characters (ideally ~200)
- [ ] Source list files created: `build/dpmjet_phojet_*_sources.txt`
- [ ] DPMJET tests pass

### Critical Validation:
```bash
# Verify the response file approach is working
cat build/dpmjet_phojet_191_sources.txt | wc -l  # Should show ~700+ lines
cat build/dpmjet_phojet_193_sources.txt | wc -l  # Should show ~700+ lines

# Check command line lengths
grep "generate_f2py.py" build/build.ninja | grep "dpmjet" | awk '{print length}'
# All should be around 200 characters, not 42,000!
```

---

## Phase 8: Add UrQMD (Optional)

**Goal:** Test UrQMD (may fail due to non-standard Fortran).

### Configuration Change:
Add to `enabled-models`:
```toml
    "urqmd34",
```

### Commands:
```bash
rm -rf build/
pip install --no-build-isolation -v -e .[test]
python -m pytest -n 0 -vv -k "urqmd"
```

### Success Criteria:
- [ ] UrQMD builds successfully
- [ ] UrQMD tests pass

### If UrQMD Fails:
This is acceptable. Document the failure and move `urqmd34` to `disabled-models` in `pyproject.toml` with a comment:
```toml
disabled-models = [
    "urqmd34",  # Non-standard Fortran, incompatible with gfortran on Windows
]
```

---

## Final Validation: Full Build

### Enable All Models (Except PYTHIA8)
Edit `pyproject.toml` to enable all models except PYTHIA8:
```toml
enabled-models = [
    "sib21", "sib23c01", "sib23d", "sib23d_star", "sib23e", "sib23e_star",
    "qgs01", "qgs2_03", "qgs2_04", "qgs3",
    "sophia", "pythia6",
    "eposlhc", "eposlhcr",
    "dpmjet_phojet191", "dpmjet_phojet193", "dpmjet_phojet307",
    "urqmd34",  # Omit if Phase 8 failed
]
disabled-models = [
    "pythia8",  # Skipped on Windows
]
```

### Full Test Suite:
```bash
rm -rf build/
pip install --no-build-isolation -v -e .[test]

# Verify all extension modules
ls -lh build/cp*/chromo/models/*.pyd

# Run full test suite with Windows-appropriate arguments
python scripts/download_data.py
python -m pytest -vv -k "not (Pythia8 or test_decay_handler)" -n 2
```

### Success Criteria:
- [ ] All enabled models build successfully
- [ ] All extension modules load
- [ ] Test suite passes (except known xfails)
- [ ] No command line length errors in build process

---

## 🐛 Troubleshooting

### Build fails with "command line too long"
- **Cause:** The response file fix didn't apply correctly
- **Check:** Look at `build/build.ninja` for custom_target commands
- **Fix:** Ensure `--source-file-list` is in the command and source list files exist

### Import fails: "DLL load failed"
- **Cause:** Missing MinGW runtime DLLs
- **Check:** Ensure MinGW bin directory is in PATH
- **Fix:** `export PATH="/c/path/to/mingw64/bin:$PATH"`

### gfortran preprocessor errors
- **Cause:** Include paths not found
- **Check:** Verify submodules are initialized
- **Fix:** `git submodule update --init --recursive`

### Tests fail with "model already initialized"
- **Cause:** Fortran models using global state
- **Expected:** Normal behavior, tests run in subprocesses
- **Fix:** Run tests with `-n 0` for serial execution when debugging

### PYTHIA8 attempts to build on Windows
- **Cause:** Platform check not working
- **Check:** Line 356 in `meson.build` should have `and host_machine.system() != 'windows'`
- **Verify:** `python -c "import platform; print(platform.system())"` returns "Windows"

---

## 📊 Validation Checklist

After completing all phases:

### Build System
- [ ] Command line lengths under 8,000 characters for all models
- [ ] Source list files generated for all models
- [ ] PYTHIA8 skipped on Windows
- [ ] All enabled models have `.pyd` files in build directory

### Functionality
- [ ] All extension modules import successfully
- [ ] Test suite passes with expected filters
- [ ] Cross-section calculations work
- [ ] Event generation works

### CI/CD
- [ ] Windows CI enabled in `.github/workflows/test.yml`
- [ ] Windows wheels enabled in `.github/workflows/release.yml`
- [ ] CI tests pass with `-k "not (Pythia8 or test_decay_handler)" -n 2`

### Documentation
- [ ] `CLAUDE.md` updated with Windows build info
- [ ] Known limitations documented
- [ ] Response file implementation explained

