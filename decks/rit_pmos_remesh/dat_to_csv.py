#!/usr/bin/env python3
"""Convert Genius TCAD .dat files to CSV.

Usage:
    python3 dat_to_csv.py idvg_vd_m0p1.dat
    python3 dat_to_csv.py idvd_prebias.dat idvg_vd_m0p1.dat
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


VARIABLE_RE = re.compile(r"^#\s*(\d+)\s+(.+?)\s*$")


def parse_headers(path: Path) -> list[str]:
    numbered: list[tuple[int, str]] = []
    with path.open() as f:
        for line in f:
            match = VARIABLE_RE.match(line)
            if match:
                numbered.append((int(match.group(1)), match.group(2).strip()))

    numbered.sort()
    return [label for _index, label in numbered]


def data_rows(path: Path):
    with path.open() as f:
        for line_number, line in enumerate(f, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                yield line_number, [float(value) for value in stripped.split()]
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: could not parse numeric row") from exc


def convert(dat_path: Path) -> Path:
    csv_path = dat_path.with_suffix(".csv")
    headers = parse_headers(dat_path)
    rows = list(data_rows(dat_path))

    if not rows:
        raise ValueError(f"{dat_path}: no numeric data rows found")

    width = len(rows[0][1])
    for line_number, row in rows:
        if len(row) != width:
            raise ValueError(
                f"{dat_path}:{line_number}: row has {len(row)} columns, expected {width}"
            )

    if not headers:
        headers = [f"col_{i}" for i in range(1, width + 1)]
    elif len(headers) != width:
        raise ValueError(
            f"{dat_path}: header has {len(headers)} columns, data has {width} columns"
        )

    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for _line_number, row in rows:
            writer.writerow(row)

    print(f"{dat_path.name} -> {csv_path.name} ({len(rows)} rows, {width} columns)")
    return csv_path


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("Usage: python3 dat_to_csv.py <file.dat> [file2.dat ...]", file=sys.stderr)
        return 2

    for arg in argv[1:]:
        convert(Path(arg))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
