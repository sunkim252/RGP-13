"""OpenFOAM 질량·에너지 전역 보존오차 측정 (오프라인 후처리).

    python3 budget_conservation.py <case>

E_t = integral(rho*h - p + 0.5*rho|U|^2) dV      순유입 = -sum_patch phi*(h + 0.5|U|^2)
벽이 h zeroGradient + noSlip 이면 벽 열유속·점성일 0 이라 인렛/출구만으로 닫혀야 한다.
필요: reconstructPar 된 시간들 + `foamPostProcess -func "writeCellVolumes writeCellCentres"`

★함정 2건 (2026-08-06 실측에서 오차 204% 를 만든 원인):
  1. 오차를 *순유입*(작은 차)으로 정규화하면 100% 처럼 보인다. 처리량/2 로 나눌 것.
  2. `zeroGradient` 패치는 파일에 `value` 블록이 없다 -> 출구가 통째로 누락된다.
     구조격자라 x 최대 열의 인접 셀값으로 대체한다 (아래 OUT).

실측 결과(tip prism 케이스): 오차 = 처리량의 8~11%, 정상상태에서 일정.
에너지오차/질량오차 = -2.81e5 J/kg ~= h_ox -> 에너지 비보존은 질량 비보존의 종속량.
"""
import glob, os, re, sys
def rd(path):
    t = open(path).read()
    m = re.search(r"internalField\s+(uniform|nonuniform)", t)
    if not m: return None, {}
    if m.group(1) == "uniform":
        v = re.search(r"uniform\s+(\([^)]*\)|[-\d.eE+]+)", t[m.start():]).group(1)
        vals = [[float(x) for x in v.strip("()").split()]]
    else:
        i = t.index("(", m.end()); j = t.index("\n)", i); b = t[i+1:j]
        if b.lstrip().startswith("("):
            f = [float(x) for x in b.replace("("," ").replace(")"," ").split()]
            vals = [f[k:k+3] for k in range(0, len(f), 3)]
        else:
            vals = [[float(x)] for x in b.split()]
    # 경계 패치 값
    bf = {}
    k = t.find("boundaryField")
    if k > 0:
        for pm in re.finditer(r"\n    (\w+)\s*\n    \{(.*?)\n    \}", t[k:], re.S):
            name, body = pm.group(1), pm.group(2)
            vm = re.search(r"value\s+(uniform|nonuniform)([^;]*);", body, re.S)
            if not vm: continue
            if vm.group(1) == "uniform":
                s = vm.group(2).strip()
                bf[name] = [[float(x) for x in s.strip("()").split()]]
            else:
                s = vm.group(2)
                a = s.index("("); z = s.rindex(")")
                bb = s[a+1:z]
                if bb.lstrip().startswith("("):
                    f = [float(x) for x in bb.replace("("," ").replace(")"," ").split()]
                    bf[name] = [f[q:q+3] for q in range(0, len(f), 3)]
                else:
                    bf[name] = [[float(x)] for x in bb.split()]
    return vals, bf

D = sys.argv[1]
ts = sorted((float(os.path.basename(p)), os.path.basename(p))
            for p in glob.glob(f"{D}/[0-9]*") if os.path.isdir(p)
            and os.path.exists(f"{D}/{os.path.basename(p)}/rho"))
vc = glob.glob(f"{D}/*/Vc"); V = [v[0] for v in rd(vc[0])[0]]
cx = glob.glob(f"{D}/*/Ccx")
CX = [v[0] for v in rd(cx[0])[0]] if cx else None
# zeroGradient 패치(=출구)는 파일에 value 가 없다 -> 인접 셀값으로 대체.
# 구조격자이므로 x 최대 열이 출구 인접 셀.
OUT = None
if CX:
    xm = max(CX); dxm = 1e-6
    OUT = [i for i in range(len(CX)) if CX[i] > xm - dxm]
    print("출구 인접 셀 %d개 (x>%.5f mm)" % (len(OUT), 1e3*(xm-dxm)))
def bcast(v):
    return v*len(V) if (v and len(v) == 1) else v
print("%-9s %11s %11s %9s | %12s %12s %9s" %
      ("t[s]", "dM/dt", "net_m_in", "질량오차%", "dEt/dt[W]", "net_E_in[W]", "E오차%"))
prev = None
for tv, t in ts:
    rho = bcast(rd(f"{D}/{t}/rho")[0]); h, hb = rd(f"{D}/{t}/h")
    p = bcast(rd(f"{D}/{t}/p")[0]); U, Ub = rd(f"{D}/{t}/U")
    phi, phib = rd(f"{D}/{t}/phi")
    h = bcast(h); U = bcast(U)
    if not rho or len(rho) < len(V): continue
    M = sum(rho[i][0]*V[i] for i in range(len(V)))
    Et = sum((rho[i][0]*h[i][0] - p[i][0]
              + 0.5*rho[i][0]*sum(c*c for c in U[i]))*V[i] for i in range(len(V)))
    nm = nE = thru_m = thru_E = 0.0
    for pat, pf in phib.items():
        if pat in hb and pat in Ub:
            hh, uu = hb[pat], Ub[pat]
        elif OUT and len(OUT) == len(pf):        # zeroGradient -> 인접 셀값
            hh = [h[i] for i in OUT]; uu = [U[i] for i in OUT]
        else:
            continue
        if len(hh) == 1: hh = hh*len(pf)
        if len(uu) == 1: uu = uu*len(pf)
        for k in range(min(len(pf), len(hh), len(uu))):
            f = pf[k][0]
            ht = hh[k][0] + 0.5*sum(c*c for c in uu[k])
            nm -= f; nE -= f*ht
            thru_m += abs(f); thru_E += abs(f*ht)
    if prev:
        dM = (M - prev[1])/(tv - prev[0]); dE = (Et - prev[2])/(tv - prev[0])
        print("%-9.5g %11.3e %11.3e %8.2f%% | %12.4e %12.4e %8.2f%%"
              % (tv, dM, nm, 100*abs(dM-nm)/max(thru_m/2, 1e-30),
                 dE, nE, 100*abs(dE-nE)/max(thru_E/2, 1e-30)))
    prev = (tv, M, Et)
