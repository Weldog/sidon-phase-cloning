import sys, json, math, time, os, importlib.util
from fractions import Fraction as Fr
folder, budget = sys.argv[1], float(sys.argv[2])
spec = importlib.util.spec_from_file_location("chk", os.path.join(folder, "check_torus_lp.py"))
chk = importlib.util.module_from_spec(spec); spec.loader.exec_module(chk)
data = json.load(open(os.path.join(folder, "cert.json")))
state = f"/home/claude/ca/state_{os.path.basename(folder)}.txt"; out = f"/home/claude/ca/lines_{os.path.basename(folder)}.txt"
start = int(open(state).read()) if os.path.exists(state) else 0
t0 = time.time(); i = start
with open(out, "a") as fo:
    while i < len(data) and time.time() - t0 < budget:
        case = data[i]; bd, ws, na, nsh, ok = chk.check(case)
        fo.write('d=%d T=%d k=%d  bound<=%.10f  max_unbal_single<=%.6f  active_types=%d shifts=%d  %s\n' % (
            case['d'], case['T'], case['k'], math.ceil(bd * 10**10) / 10**10, math.ceil(ws * 10**6) / 10**6, na, nsh, 'PASS' if ok else 'FAIL'))
        fo.flush(); i += 1
open(state, "w").write(str(i)); print(f"{os.path.basename(folder)}: {i}/{len(data)} cases done")
