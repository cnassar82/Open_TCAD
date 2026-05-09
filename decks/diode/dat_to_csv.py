"""
Convert Genius TCAD .dat files to .csv using pandas.

Usage:
    python dat_to_csv.py forward_iv.dat
    python dat_to_csv.py forward_iv.dat reverse_iv.dat
    python dat_to_csv.py *.dat
"""

import sys
import re
import pandas as pd
from pathlib import Path


def parse_headers(path):
    headers = []
    with open(path) as f:
        for line in f:
            m = re.match(r'^#\s+\d+\s+(.+)$', line)
            if m:
                headers.append(m.group(1).strip())
    return headers


def convert(dat_path):
    dat_path = Path(dat_path)
    csv_path = dat_path.with_suffix('.csv')

    headers = parse_headers(dat_path)
    df = pd.read_csv(dat_path, comment='#', sep=r'\s+', header=None, names=headers)
    df.to_csv(csv_path, index=False)

    print(f"{dat_path.name} -> {csv_path.name}  ({len(df)} rows, {len(df.columns)} columns)")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python dat_to_csv.py <file.dat> [file2.dat ...]")
        sys.exit(1)

    for arg in sys.argv[1:]:
        convert(arg)
