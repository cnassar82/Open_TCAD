#!/usr/bin/env python3
"""Check quadrature exactness, affine energy, and the physical nullspaces."""
import argparse,subprocess,tempfile,os
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--src',type=Path,default=root/'src');args=p.parse_args();src=args.src.resolve()
with tempfile.TemporaryDirectory() as d:
 exe=Path(d)/'test';cmd=['gcc','-std=gnu89','-O1','-g','-fcommon','-ffunction-sections','-fdata-sections','-fsanitize=address,undefined','-fno-pie','-no-pie','-Wl,--gc-sections','-I',str(src/'include'),str(root/'tests/fe_triangle7.c'),str(src/'finel/tri7.c'),str(src/'oxide/triox.c'),'-lm','-o',str(exe)];subprocess.run(cmd,check=True)
 lines=subprocess.check_output([str(exe)],text=True,env={**os.environ,'ASAN_OPTIONS':'detect_leaks=0:halt_on_error=1','UBSAN_OPTIONS':'halt_on_error=1'}).splitlines()
 for name,part,expected in [('elastic',lines[:14],11),('oxide',lines[14:],17)]:
  a=np.array([list(map(float,line.split())) for line in part]);assert np.max(abs(a-a.T))<=64*np.finfo(float).eps*np.max(abs(a))
  eig=np.linalg.eigvalsh(a);tol=max(abs(eig))*1e-10;assert np.count_nonzero(eig>tol)==expected,(name,eig);assert min(eig)>-tol
  if float(os.environ.get('SUPREM_TEST_BULK','1'))>=1e6:
   pressure_modes=3 if name=='oxide' and os.environ.get('SUPREM4_FE_BBAR_P1') else 5
   assert np.count_nonzero(eig>float(os.environ['SUPREM_TEST_BULK'])**.5)==pressure_modes,(name,eig)
   print(f'PASS {name}: {pressure_modes} volumetric modes, no added zero modes')
  print(f'PASS {name}: rank={expected}, physical nullity={len(eig)-expected}')
 print('PASS degree-five moments, affine fields, obtuse and translated thin element energies; ASan/UBSan')
