"""Compare the engine predicate with exact rational signs, including wide exponents."""
from fractions import Fraction
from pathlib import Path
import random
import struct
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


class OrientationTest(unittest.TestCase):
    def test_exact_stored_float_signs(self):
        rng=random.Random(140026)
        def value():
            bits=rng.getrandbits(32)
            # Avoid infinities and NaNs; retain subnormal and extreme exponents.
            if (bits>>23)&255 == 255:bits ^= 1<<23
            return struct.unpack('f',struct.pack('I',bits))[0]
        cases=[[0.,0.,1.,0.,0.,1.], [0.,0.,1.,1.,2.,2.],
               [-2.8819578709955665e-25,-3.1026742362882942e-5,
                -4.440892465941686e-23,-3.3673830330371857e-5,
                -.001,-.00087734760018065572]]
        cases += [[value() for _ in range(6)] for _ in range(10000)]
        cases=[[struct.unpack('f',struct.pack('f',v))[0] for v in row] for row in cases]
        expected=[]
        for row in cases:
            ax,ay,bx,by,cx,cy=map(Fraction,row)
            det=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax)
            expected.append(int(det>0)-int(det<0))
        with tempfile.TemporaryDirectory() as directory:
            work=Path(directory); driver=work/'driver.c';binary=work/'orient'
            driver.write_text('''#include <stdio.h>
extern int orient2d_float_sign(float *,float *,float *);
int main(void){float a[2],b[2],c[2];
while(scanf("%f %f %f %f %f %f",a,a+1,b,b+1,c,c+1)==6)
printf("%d\\n",orient2d_float_sign(a,b,c));return 0;}
''')
            subprocess.run(['gcc','-std=gnu89','-O1','-ffunction-sections','-I',str(ROOT/'src/include'),
                            str(driver),str(ROOT/'src/geom/misc.c'),'-Wl,--gc-sections','-lm','-o',str(binary)],check=True)
            data='\n'.join(' '.join(v.hex() for v in row) for row in cases)+'\n'
            result=subprocess.run([str(binary)],input=data,text=True,capture_output=True,check=True)
            actual=list(map(int,result.stdout.split()))
            self.assertEqual(actual,expected)


if __name__=='__main__':unittest.main()
