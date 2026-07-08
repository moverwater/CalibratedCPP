#!/usr/bin/env python3
"""
Extract the first N% of sites (by alignment column) from a NEXUS DNA matrix,
keeping all taxa, and write the result out as a new, well-formed NEXUS file.

Usage:
    python truncate_alignment.py [input.nex] [output.nex] [fraction]

Defaults:
    input.nex  = primates.nex
    output.nex = primates_first10pct.nex
    fraction   = 0.1 (first 10% of sites)
"""
import re
import sys


def parse_nexus_matrix(text):
    """Return (ntax, nchar, format_line, taxa) where taxa is a list of (name, sequence)."""
    dims = re.search(r"Dimensions\s+ntax=(\d+)\s+nchar=(\d+)\s*;", text, re.IGNORECASE)
    if not dims:
        raise ValueError("Could not find a Dimensions line (ntax=.../nchar=...)")
    ntax, nchar = int(dims.group(1)), int(dims.group(2))

    fmt = re.search(r"(Format[^\n]*;)", text, re.IGNORECASE)
    format_line = fmt.group(1).strip() if fmt else "Format datatype=nucleotide gap=-;"

    matrix_start = re.search(r"Matrix\s*", text, re.IGNORECASE)
    if not matrix_start:
        raise ValueError("Could not find the Matrix keyword")

    lines = text[matrix_start.end():].splitlines()
    taxa = []
    for line in lines:
        line = line.strip()
        if not line or line in (";", "End;") or line.lower().startswith("end"):
            continue
        if line == ";":
            break
        name, seq = line.split(None, 1)
        taxa.append((name, seq.strip().rstrip(";").strip()))
        if len(taxa) == ntax:
            break

    if len(taxa) != ntax:
        raise ValueError(f"Expected {ntax} taxa, found {len(taxa)}")
    for name, seq in taxa:
        if len(seq) != nchar:
            raise ValueError(f"Taxon {name}: expected {nchar} characters, got {len(seq)}")

    return ntax, nchar, format_line, taxa


def write_nexus(path, ntax, nchar, format_line, taxa):
    name_width = max(len(name) for name, _ in taxa) + 3
    with open(path, "w", newline="\n") as fh:
        fh.write("#NEXUS\n\n")
        fh.write("BEGIN DATA;\n")
        fh.write(f"\tDimensions ntax={ntax} nchar={nchar};\n")
        fh.write(f"\t{format_line}\n")
        fh.write("\tMatrix\n")
        for name, seq in taxa:
            fh.write(f"\t{name.ljust(name_width)}{seq}\n")
        fh.write(";\n")
        fh.write("END;\n")


def main():
    input_path = sys.argv[1] if len(sys.argv) > 1 else "primates.nex"
    output_path = sys.argv[2] if len(sys.argv) > 2 else "primates_first10pct.nex"
    fraction = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1

    with open(input_path) as fh:
        text = fh.read()

    ntax, nchar, format_line, taxa = parse_nexus_matrix(text)
    new_nchar = int(nchar * fraction)
    truncated = [(name, seq[:new_nchar]) for name, seq in taxa]

    write_nexus(output_path, ntax, new_nchar, format_line, truncated)
    print(f"Wrote {output_path}: {ntax} taxa, {new_nchar} of {nchar} sites "
          f"({fraction:.0%})")


if __name__ == "__main__":
    main()
