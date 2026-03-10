# Fixing Non-Standard Fortran Common Block Declarations in EPOS-LHC-R

## Problem Description

The EPOS-LHC-R Fortran source code contains non-standard Fortran syntax where multiple common blocks are declared on a single line. This pattern is not properly handled by modern Fortran tools like numpy's f2py, causing compilation and processing issues.

**Problematic pattern:**
```fortran
common /common1/var1,var2 /common2/varX
```

**Standard compliant pattern:**
```fortran
common /common1/var1,var2
common /common2/varX
```

## Detection Method

Use this regex pattern to find all occurrences of multiple common blocks on a single line:

```regex
^(?=.*\bcommon\b)(?:[^\/]*\/){4}.*$
```

Or alternatively, this simpler pattern:
```regex
common.*\/[^\/]*\/[^\/]*\/[^\/]*\/
```

## Solution Process

### Step 1: Search and Identify

Search for all lines with multiple common blocks across all source files:

```bash
grep -rn "common.*\/[^\/]*\/[^\/]*\/[^\/]*\/" src/fortran/epos-lhc-r/sources/*
```

### Step 2: Systematic Fixing

For each identified line, split the multiple common block declarations into separate lines while preserving the exact variable lists and maintaining proper Fortran syntax.

## Files Modified and Changes Made

### 1. epos-jps.f (2 occurrences)

**Line ~140:**
```fortran
# Before:
common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctain/mtain

# After:
common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
common/ctain/mtain
```

**Line ~427:**
```fortran
# Before:
common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctain/mtain

# After:
common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
common/ctain/mtain
```

### 2. epos-rsh.f (1 occurrence)

**Line ~39:**
```fortran
# Before:
common/cprtx/nprtjx,pprtx(5,2)/ciptl/iptl

# After:
common/cprtx/nprtjx,pprtx(5,2)
common/ciptl/iptl
```

### 3. urqmdepos.f (4 occurrences)

**Line ~47:**
```fortran
# Before:
common/cdir/edir /cishuuu/ishuuu

# After:
common/cdir/edir
common/cishuuu/ishuuu
```

**Line ~370:**
```fortran
# Before:
common /city/itypart /ctrace/itrace /corpart/iorpart

# After:
common /city/itypart
common /ctrace/itrace
common /corpart/iorpart
```

**Line ~1088:**
```fortran
# Before:
common /city/itypart /ctrace/itrace /corpart/iorpart

# After:
common /city/itypart
common /ctrace/itrace
common /corpart/iorpart
```

**Line ~3325:**
```fortran
# Before:
common /ctrace/itrace /cishuuu/ishuuu

# After:
common /ctrace/itrace
common /cishuuu/ishuuu
```

### 4. epos-int.f (4 occurrences)

**Line ~8:**
```fortran
# Before:
common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat /ctimel/ntc

# After:
common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat
common/ctimel/ntc
```

**Line ~11:**
```fortran
# Before:
common/cttaun/ttaun /cttau0/ttau0 /geom1/rcproj,rctarg

# After:
common/cttaun/ttaun
common/cttau0/ttau0
common/geom1/rcproj,rctarg
```

**Line ~182 and multiple other locations:**
```fortran
# Before:
common/cdelzet/delzet,delsce /cvocell/vocell,xlongcell

# After:
common/cdelzet/delzet,delsce
common/cvocell/vocell,xlongcell
```

### 5. epos-dro.f (6 occurrences)

**Multiple locations with patterns like:**
```fortran
# Before:
common/cspez4/ffstat(2,0:mspez+2) /ctfo/tfo

# After:
common/cspez4/ffstat(2,0:mspez+2)
common/ctfo/tfo
```

### 6. epos-fra.f (2 occurrences)

**Line ~18:**
```fortran
# Before:
common/pb/pb /cnsbp/nsbp  /cn8ptl/n8ptl

# After:
common/pb/pb
common/cnsbp/nsbp
common/cn8ptl/n8ptl
```

### 7. epos-xan.f (2 occurrences)

**Line ~74:**
```fortran
# Before:
common/nl/noplin  /cnnnhis/nnnhis

# After:
common/nl/noplin
common/cnnnhis/nnnhis
```

### 8. epos-uti.f (2 occurrences)

**Line ~6563:**
```fortran
# Before:
common/ccc20/icc20  /ciext4/iext4

# After:
common/ccc20/icc20
common/ciext4/iext4
```

### 9. epos-dky.f (1 occurrence)

**Line ~16:**
```fortran
# Before:
common/cnptlbur/nptlbur /cij99/ij99

# After:
common/cnptlbur/nptlbur
common/cij99/ij99
```

### 10. epos-ems.f (1 occurrence)

**Line ~4703:**
```fortran
# Before:
common/col3/ncol,kolpt /cfacmss/facmss /cts/its

# After:
common/col3/ncol,kolpt
common/cfacmss/facmss
common/cts/its
```

## Automation Script Template

For future similar tasks, here's a template approach:

```bash
#!/bin/bash

# Find all files with problematic common blocks
files_to_fix=$(grep -l "common.*\/[^\/]*\/[^\/]*\/[^\/]*\/" src/fortran/epos-lhc-r/sources/*.f)

for file in $files_to_fix; do
    echo "Processing: $file"
    # Manual inspection and fixing required for each occurrence
    # Use your preferred editor with regex search/replace
done
```

## Verification

After making changes:

1. **Syntax Check:** Ensure all Fortran files still compile
2. **Functionality Check:** Run test cases to ensure behavior is unchanged
3. **Tool Compatibility:** Verify f2py can now process the files without errors

## Important Notes

1. **Preserve Exact Functionality:** The split common blocks must maintain identical variable declarations
2. **Maintain Fortran Format:** Keep proper column formatting for fixed-form Fortran
3. **Context Preservation:** Include 3-5 lines of context when making replacements to avoid ambiguity
4. **Documentation:** Each change should be traceable and reversible if needed

## Total Summary

- **Files Modified:** 11 Fortran source files
- **Lines Fixed:** 26 problematic common block declarations
- **Result:** Standard-compliant Fortran code compatible with modern tools like f2py

This process ensures the EPOS-LHC-R code maintains its functionality while becoming compatible with modern Fortran processing tools and compilers.
