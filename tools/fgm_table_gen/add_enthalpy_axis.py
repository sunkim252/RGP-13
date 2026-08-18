"""Non-adiabatic FPV (method b): add an ENTHALPY-DEFECT 4th axis to a 3-D FGM
table, so the manifold can be queried at inlet/heat-loss enthalpies different
from the 800 K flamelet boundary -- in particular cryogenic LOx + warm kerosene.

Coordinate (the clean part): TOTAL enthalpy h is a CONSERVED scalar (like Z), so
the adiabatic reference is a mixing line  h_ad(Z) = (1-Z) h_ox + Z h_fuel  with
h_ox, h_fuel the pure-stream total enthalpies at the 800 K flamelet boundary.
The 4th axis is the enthalpy DEFECT  dh = h - h_ad(Z)  (dh=0 -> adiabatic 800 K
manifold = the original 3-D table EXACTLY; dh<0 -> colder inlet / heat loss).

Build (frozen-composition, no extra flamelet solves -> avoids the stiff cryo
flamelet Newton stall). The defect is the SENSIBLE-enthalpy change referenced to
the flamelet's OWN adiabatic state, so the formation-enthalpy / differential-
diffusion (Takahashi/Chung) / dual-gas offset CANCELS (h(Tad,Y) need NOT equal
the mixing line):
  T(Z,gZ,c,dh):  solve  h(T,Y) - h(Tad,Y) = dh   at frozen Y(Z,gZ,c)   [RK EOS]
  Y(Z,gZ,c,dh) = Y(Z,gZ,c)                                    (frozen)
  omega_C(...,dh) = omega_C(Z,gZ,c) * ignition_gate(T(...,dh))  (inert when cold)
Solver side: dh = h_transported - ((1-Z)hOx + Z hFuel); at a stream boundary the
composition IS the pure stream so h(Tad,Y)=h_mix there -> the two references
coincide, and in the interior the heat-loss departure is measured consistently.

Solver reads: nH, fourthAxis=enthalpy, hOx, hFuel, the dh axis. It transports
total enthalpy h, forms dh = h - ((1-Z)hOx + Z hFuel), and does the 4-D lookup.

** CAVEAT: NASA-7 thermo floors at ~300 K; below that (cold O2) the h<->T
inversion is EXTRAPOLATED. EOS density is physical, cold-stream sensible h is
approximate -> swap a NIST/REFPROP h(T) for the pure cold streams for production. **

Usage: python add_enthalpy_axis.py <src_3d.npz> <out_prefix> [n_h] [T_ox_cold] [T_fuel_cold]
"""
import sys, time
import numpy as np, cantera as ct

SRC   = sys.argv[1]
OUT   = sys.argv[2]
N_H   = int(sys.argv[3])   if len(sys.argv) > 3 else 11
T_OXC = float(sys.argv[4]) if len(sys.argv) > 4 else 100.0   # cold oxidizer (LOx)
T_FUC = float(sys.argv[5]) if len(sys.argv) > 5 else 300.0   # cold fuel (kerosene)
T_BOUND = 800.0
FUEL = {'NC10H22': 0.74, 'PHC3H7': 0.15, 'CYC9H18': 0.11}
OXID = {'O2': 1.0}
import os as _os
MECH = _os.environ.get('RGP_MECH', 'data/wang2011_srk_v32.yaml')
T_IGN, DT_IGN = 900.0, 150.0     # source ignition gate (reaction frozen below ~crossover)
def log(m): print(m, flush=True)

d = np.load(SRC, allow_pickle=True)
Z_axis = d['Z_axis']; g_axis = d['gZ_axis']; C_axis = d['C_axis']
P = float(d['P']); species = [str(s) for s in d['species']]
T3 = d['T']; om3 = d['omega_C']
has_hrr = 'hrr' in d.files
hrr3 = d['hrr'] if has_hrr else None
Yk = {sp: d[f'Y_{sp}'] for sp in species}
nZ, nGz, nC = T3.shape
log(f"loaded 3-D table {SRC}: {nZ}x{nGz}x{nC}, P={P/1e5:.1f}bar, {len(species)} species")

g = ct.Solution(MECH)
# Ideal-gas twin for the h(T) inversion: the SRK cubic root solver fails to
# converge for some (T,Y) with the real composition (e.g. T~1630 K) and aborts
# the build. The enthalpy used for the T inversion is dominated by formation +
# sensible terms (both identical in ideal/real); the real-fluid departure is
# small and largely cancels in the sensible DIFFERENCE h(T,Y)-h(Tad,Y). So fall
# back to the ideal twin wherever the cubic fails -> robust, no crashes.
gi = ct.Solution('data/wang2011_ideal_v32.yaml')
def h_at(sol, T, Ydict):
    sol.TPY = T, P, Ydict; return sol.enthalpy_mass
# boundary total enthalpies at the flamelet boundary T (defect reference)
g.TPX = T_BOUND, P, OXID;  h_ox = g.enthalpy_mass
g.TPX = T_BOUND, P, FUEL;  h_fuel = g.enthalpy_mass
# cold-stream defects set the axis span
g.TPX = T_OXC, P, OXID;  dh_ox = g.enthalpy_mass - h_ox
g.TPX = T_FUC, P, FUEL;  dh_fu = g.enthalpy_mass - h_fuel
# 저온 유입(LOx 100 K / 케로신 300 K)이 요구하는 결손.
# 주의: 축 좌표는 혼합선이 아니라 dhRef(Z,c) 기준이므로, 실제 요구량은
#   dh_need(Z,c) = cold(Z) - dhRef(Z,c)
# 이고 dhRef 가 양인 곳(차등확산으로 살찐 희박·과농측)에서 더 깊어진다.
# w14 실측: cold 최저 -1.460, dhRef 최대 +0.795 -> 요구 -2.26 MJ/kg 인데
# 옛 공식 1.05*min(dh_ox,dh_fu) = -1.533 이라 0.6 MJ/kg 모자랐다.
dh_min = 1.05 * min(dh_ox, dh_fu)         # most negative defect (+5% margin)
if len(sys.argv) > 7:
    dh_min = float(sys.argv[7])           # 명시 하한 (요구량 기반)
# The axis was one-sided (deficit only), which silently clamps every lookup that
# needs MORE enthalpy than the dh=0 slice. Measured on the 134-member family:
# the flamelet defect vs the mixing line is positive for 100% of members at
# Z=0.02-0.08 and Z=0.7 (median +386 kJ/kg at Z=0.04) -- preferential diffusion
# enriches the lean and rich edges -- and the low-chi members sit ~340 kJ/kg
# above the dh=0 slice. Those lookups clamped and came back 300-460 K too cold.
DH_MAX = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0   # J/kg, 0 = 옛 동작
h_axis = (np.linspace(dh_min, 0.0, N_H) if DH_MAX <= 0 else
          np.concatenate([np.linspace(dh_min, 0.0, N_H)[:-1],
                          np.linspace(0.0, DH_MAX, max(2, N_H//2 + 1))]))
log(f"h_ox(800K)={h_ox/1e6:.3f} h_fuel(800K)={h_fuel/1e6:.3f} MJ/kg")
log(f"cold defects: O2@{T_OXC:g}K dh={dh_ox/1e6:.3f}  fuel@{T_FUC:g}K dh={dh_fu/1e6:.3f} MJ/kg")
N_H = len(h_axis)                          # 양의 구간을 붙이면 길이가 늘어난다
log(f"enthalpy-defect axis: {N_H} pts, "
    f"[{h_axis[0]/1e6:.3f}, {h_axis[-1]/1e6:+.3f}] MJ/kg")

# ---- build the 4-D arrays (frozen composition, ROBUST h->T inversion) ----
# Cantera's HPY Newton lands on spurious roots below the 300 K NASA-poly floor
# (non-monotonic extrapolated h(T)). Instead build a MONOTONE-ENFORCED h(T) grid
# per composition and invert by interpolation, clamped to [T_FLOOR, T_ad]. This
# never returns garbage; over-cold targets (e.g. a burnt state pushed to the
# coldest defect -- a corner the CFD never visits) simply clamp to T_FLOOR.
T_FLOOR = 80.0          # coldest physical (LOx ~90 K); covers O2 @100 K inlet
NT = 32                 # h(T) grid points per composition
T_CUT = 600.0           # hard reaction cutoff: omega = 0 below this T
T4  = np.zeros((nZ, nGz, nC, N_H))
om4 = np.zeros((nZ, nGz, nC, N_H))
# h(T_ad, Y) per node: the manifold's OWN adiabatic enthalpy. With differential
# diffusion this is NOT the mixing line (measured -0.67 MJ/kg at Z_st for
# SRK+HP-Chung, +1.40 for mixture-averaged), yet the solver forms its defect as
# h - mixing line. Ship the difference as dhRef so dh=0 <=> adiabatic there too.
hAd3 = np.zeros((nZ, nGz, nC))
hrr4 = np.zeros((nZ, nGz, nC, N_H)) if has_hrr else None
t0 = time.time(); npt = 0
for iZ in range(nZ):
    for iG in range(nGz):
        for iC in range(nC):
            Yvec = np.array([Yk[sp][iZ, iG, iC] for sp in species])
            s = Yvec.sum()
            Tad = T3[iZ, iG, iC]
            if s <= 1e-6:                      # empty/edge -> replicate adiabatic
                hAd3[iZ, iG, iC] = np.nan          # -> dhRef 0 (mixing line)
                T4[iZ, iG, iC, :] = Tad
                om4[iZ, iG, iC, :] = om3[iZ, iG, iC]
                if has_hrr: hrr4[iZ, iG, iC, :] = hrr3[iZ, iG, iC]
                continue
            Yvec /= s
            # ** set composition BY NAME ** -- the table's species list is sorted
            # alphabetically, NOT in Cantera's order, so g.TPY=T,P,Yvec (array)
            # SCRAMBLES the composition. A name->value dict is order-independent.
            Ydict = dict(zip(species, Yvec))
            # monotone h(T) from T_FLOOR up to the local adiabatic T.
            # With a two-sided dh axis the grid must reach ABOVE Tad, otherwise
            # np.interp clamps every positive-dh slice back to Tad and the whole
            # extension is a no-op (measured: w15 came out bit-identical to w14).
            # cp >= ~800 J/kg/K here, so DH_MAX/800 bounds the extra rise.
            Ttop = max(Tad, T_FLOOR + 50.0)
            if DH_MAX > 0:
                Ttop = min(Ttop + DH_MAX/800.0, 4800.0)
            Tgrid = np.linspace(T_FLOOR, Ttop, NT)
            hgrid = np.empty(NT)
            for j, Tg in enumerate(Tgrid):
                try:
                    hgrid[j] = h_at(g, Tg, Ydict)        # SRK (real-fluid)
                except Exception:
                    hgrid[j] = h_at(gi, Tg, Ydict)       # ideal fallback (cubic failed)
                npt += 1
            hgrid = np.maximum.accumulate(hgrid)   # enforce dh/dT>0 (clip extrap wiggle)
            # SENSIBLE-enthalpy defect referenced to the flamelet's OWN adiabatic
            # state h(Tad,Y) = hgrid[-1]:  solve  h(T,Y) - h(Tad,Y) = dh.  The
            # formation-enthalpy (and the differential-diffusion / dual-gas) offset
            # CANCELS in the difference, so dh=0 -> T=Tad exactly and a cold defect
            # cools T by the right sensible amount regardless of that offset. (An
            # absolute h_ad=(1-Z)hOx+ZhFuel reference fails because h(Tad,Y) != h_ad
            # for real-fluid flamelets.) Matches the solver's dh = h - h_mix(Z).
            # 기준은 h(Tad) 지 격자 끝이 아니다 — 양방향 축이면 끝이 Tad 위에 있다
            try:
                hAd = h_at(g, Tad, Ydict)
            except Exception:
                hAd = h_at(gi, Tad, Ydict)
            hAd = float(np.clip(hAd, hgrid[0], hgrid[-1]))   # 단조화 격자 안으로
            hAd3[iZ, iG, iC] = hAd
            targets = hAd + h_axis
            Tloc = np.interp(targets, hgrid, Tgrid)  # clamps to [T_FLOOR, Tad]
            Tloc[np.abs(h_axis) < 1e-9] = Tad        # dh=0 slice exact
            T4[iZ, iG, iC, :] = Tloc
            gate = 0.5 * (1.0 + np.tanh((Tloc - T_IGN) / DT_IGN))
            gate[Tloc < T_CUT] = 0.0                 # hard inert cutoff when cold
            om4[iZ, iG, iC, :] = om3[iZ, iG, iC] * gate
            if has_hrr: hrr4[iZ, iG, iC, :] = hrr3[iZ, iG, iC] * gate
    if iZ % 10 == 0:
        log(f"  iZ={iZ}/{nZ}  ({time.time()-t0:.0f}s, {npt} h(T) evals)")
# Y frozen: broadcast along the new axis
Y4 = {sp: np.repeat(Yk[sp][..., None], N_H, axis=3) for sp in species}
log(f"4-D build done: {npt} h(T) evals in {time.time()-t0:.0f}s, shape {T4.shape}")

# sanity: dh=0 slice must equal the original 3-D T exactly
# 양방향 축에서는 dh=0 이 마지막이 아니라 중간이다. 인덱스를 찾아 쓴다.
i_ad = int(np.argmin(np.abs(h_axis)))
err = np.abs(T4[..., i_ad] - T3).max()
log(f"dh=0 slice (idx {i_ad}/{N_H-1}) vs original 3-D Tmax-diff = {err:.3e} K (should be ~0)")

# ---- write npz + OpenFOAM dict (enthalpy-axis 4-D format) ----
def fmt(a): return "\n".join(f"    {v:.8e}" for v in np.asarray(a).reshape(-1))

# dhRef(Z,gZ,c) = h(T_ad,Y) - ((1-Z) hOx + Z hFuel), 4th축으로 복제해 저장
h_mix = (1.0 - Z_axis)[:, None, None]*h_ox + Z_axis[:, None, None]*h_fuel
dhRef3 = hAd3 - h_mix
dhRef3[~np.isfinite(dhRef3)] = 0.0
dhRef4 = np.repeat(dhRef3[..., None], N_H, axis=3)
log(f"dhRef: [{dhRef3.min()/1e6:.3f}, {dhRef3.max()/1e6:.3f}] MJ/kg "
    f"(|dhRef| > 0.05 MJ/kg at {100*np.mean(np.abs(dhRef3) > 5e4):.1f}% of nodes)")

np.savez_compressed(OUT + '.npz', dhRef=dhRef4,
    Z_axis=Z_axis, gZ_axis=g_axis, C_axis=C_axis, h_axis=h_axis,
    species=np.array(species, dtype=object), P=P, h_ox=h_ox, h_fuel=h_fuel,
    fourthAxis='enthalpy', T=T4, omega_C=om4,
    **({'hrr': hrr4} if has_hrr else {}),
    **{f'Y_{sp}': Y4[sp] for sp in species})
log(f"[write] {OUT}.npz")

nH = N_H
blocks = [f"sourcePV\n(\n{fmt(om4)}\n);", f"T\n(\n{fmt(T4)}\n);",
          f"dhRef\n(\n{fmt(dhRef4)}\n);"]
blocks.append("species ( " + " ".join(species) + " );")
for sp in species:
    blocks.append(f"Y_{sp}\n(\n{fmt(Y4[sp])}\n);")
body = "\n\n".join(blocks)
header = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| FGM 4-D NON-ADIABATIC (Z~, gZ, c, dh=enthalpy-defect) -- add_enthalpy_axis.py
| P_ref = {P:.6g} Pa | h_ox(800K)={h_ox:.6g} h_fuel(800K)={h_fuel:.6g} J/kg
| dh = h_total - ((1-Z) hOx + Z hFuel) - dhRef(Z,gZ,c);  dh=0 -> adiabatic
| dhRef = h(T_ad,Y)|manifold - mixing line (0 for a unity-Lewis manifold)
| flat C-order: idx = ((iZ*nGz + iGz)*nC + iC)*nH + iH
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      fgmProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

nZ      {nZ};
nGz     {nGz};
nC      {nC};
nH      {nH};
fourthAxis  enthalpy;
hOx     {h_ox:.8e};
hFuel   {h_fuel:.8e};

Le
{{
    Z   0.63;
    C   0.60;
    h   0.72;
}}

Z
(
{fmt(Z_axis)}
);

gZ
(
{fmt(g_axis)}
);

C
(
{fmt(C_axis)}
);

enthalpy
(
{fmt(h_axis)}
);

{body}

// ************************************************************************* //
"""
open(OUT, 'w').write(header)
log(f"[write] {OUT}  (4-D non-adiabatic, nH={nH})")
