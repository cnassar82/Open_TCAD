#!/usr/bin/env python3
"""Check this engine's active-dopant codes at the SUPREM/Genius boundary."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'suprem_viewer'))
from make_genius_contacts import SolutionRow,solution_name_to_index,acceptor_donor
# Carrier densities, velocities, grain size and tracer must never become doping.
codes=[5,23,3,21,2,20,4,22,7,8,9,24,25]
values=[8e17,2e16,9e17,3e16,8e18,4e16,9e18,5e16,1e20,2e20,-3e20,4e20,5e20]
accept,donor=acceptor_donor(SolutionRow(0,3,values),solution_name_to_index(codes))
assert accept==2e16 and donor==12e16,(accept,donor)
accept,donor=acceptor_donor(SolutionRow(0,3,[2e16,3e16]),solution_name_to_index([5,3]))
assert accept==2e16 and donor==3e16
print('PASS active-dopant export and non-dopant exclusion')
