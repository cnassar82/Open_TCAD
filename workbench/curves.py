"""Read numeric CSV, TSV, and Genius gnuplot tables."""

import csv
import math
import re


def read_table(path):
    columns = []
    data = []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        if path.suffix.lower() in {".csv", ".tsv"}:
            rows = csv.reader(stream, delimiter="," if path.suffix.lower() == ".csv" else "\t")
        else:
            rows = []
            for line in stream:
                if line.lstrip().startswith("#"):
                    variable = re.match(r"\s*#\s*\t(\d+)\t(.+)", line)
                    if variable:
                        columns.append(variable[2].strip())
                elif line.strip():
                    rows.append(line.split())
        for row in rows:
            if not any(cell.strip() for cell in row) or row[0].lstrip().startswith("#"):
                continue
            try:
                values = tuple(float(cell) if cell.strip() else math.nan for cell in row)
            except ValueError:
                if not columns and not data:
                    columns = [cell.strip() for cell in row]
                    continue
                raise ValueError(f"Non-numeric table row in {path.name}.") from None
            if not columns:
                columns = [f"Column {i}" for i in range(1, len(values) + 1)]
            if len(values) != len(columns):
                raise ValueError(f"Inconsistent column count in {path.name}.")
            data.append(values)
    if len(columns) < 2 or not data:
        raise ValueError(f"{path.name} needs at least two numeric columns and one row.")
    return columns, list(zip(*data))
