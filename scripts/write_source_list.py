#!/usr/bin/env python3
"""
Helper script to write source file list to a text file.
Used by meson.build to avoid Windows command line length limits.

Usage: python write_source_list.py <source1> <source2> ... <sourceN>
Outputs: one source file path per line to stdout
"""

import sys

if __name__ == "__main__":
    # Skip script name, output one source file path per line
    for source in sys.argv[1:]:
        print(source)
