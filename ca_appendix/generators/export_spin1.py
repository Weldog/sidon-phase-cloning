import sys; sys.argv = ['x', '0']
src = open('cert_spin1_O.py').read().split("cert = dict(")[0]
exec(src)
import json
from fractions import Fraction as Fr
ivt = open('/workspace/r9/cert_spin1_O_ivt.py').read()
part = ivt.split("SQ2 = iv.sqrt(2)")[0].split("cert = json.load(open('cert_spin1_O_seed0.json'))")[1]
cert = json.load(open('cert_spin1_O_seed0.json'))
assert [[str(a), str(b)] for a, b in fq] == cert['fiducial_rational']
exec(part)
t = Fr(1, 10 ** 13); DEN = 10 ** 15
def shifted(sgn):
    return [(a + sgn * t * Fr(round(z.real * DEN), DEN), b + sgn * t * Fr(round(z.imag * DEN), DEN)) for (a, b), z in zip(fq, w)]
cl = lambda z: [float(z.real), float(z.imag)]
data = dict(fiducial=[[str(a), str(b)] for a, b in fq], v_top=[cl(z) for z in Vx[:, -1]], mu=float(lh[-1]) + 1e-12,
            B4=[[cl(B4[i, j]) for j in range(6)] for i in range(6)],
            f_minus=[[str(a), str(b)] for a, b in shifted(-1)], f_plus=[[str(a), str(b)] for a, b in shifted(1)])
json.dump(data, open('spin1_data.json', 'w'), indent=0)
print('saved', float(gap.a))
