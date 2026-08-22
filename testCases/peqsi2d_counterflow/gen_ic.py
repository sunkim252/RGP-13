"""Initial fields for the 2-D counterflow-geometry FGM case.

Z runs linearly 0 -> 1 across x at fixed progress variable, so every cell
sits on a manifold node in (gZ, c, dh) and only Z varies.  That gives the
strong composition and enthalpy gradients the energy audit needs, with
the state still exactly on the table at t = 0.
"""
import os, re, sys

T = "/home/sunkim/openfoam/RGP-13/RGP-13-realFluid/tools/fgm_table_gen/data/fpvwork/peqsi_dp2b_lo30_coef_rg.dict"
D, C0 = sys.argv[1], float(sys.argv[2])
NX, NY = 100, 20
nZ, nG, nC, nH = 51, 11, 41, 16
iG, iH = 0, 10
node = lambda z, g, c, h: ((z*nG + g)*nC + c)*nH + h

scal, axes, species, blocks = {}, {}, [], {}
with open(T) as f:
    for line in f:
        m = re.match(r"^(\w+)\s+([-0-9.eE+]+);", line)
        if m:
            scal[m.group(1)] = float(m.group(2)); continue
        m = re.match(r"^(\w+)\s+\(", line)
        if not m:
            continue
        k = m.group(1)
        if k == "species":
            species = line.split("(", 1)[1].split(")")[0].split()
        elif k in ("Z", "C", "enthalpy", "Cnorm"):
            axes[k] = [float(x) for x in line.split("(", 1)[1].split(")")[0].split()]

iC = min(range(nC), key=lambda i: abs(axes["C"][i] - C0))
# every cell sits on a Z node: pick NX nodes spread over the axis
zidx = [round(i*(nZ - 1)/(NX - 1)) for i in range(NX)]
want = {"T", "dhRef"} | {"Y_" + s for s in species}
with open(T) as f:
    for line in f:
        m = re.match(r"^(\w+)\s+\(", line)
        if m and m.group(1) in want:
            blocks[m.group(1)] = line.split("(", 1)[1].split(")")[0].split()
            want.discard(m.group(1))
        if not want:
            break

pRef = scal["pressure"]
col = {}
for iZ in set(zidx):
    n = node(iZ, iG, iC, iH)
    Zc = axes["Z"][iZ]
    hMix = (1.0 - Zc)*scal["hOx"] + Zc*scal["hFuel"]
    col[iZ] = dict(
        Z=Zc, T=float(blocks["T"][n]),
        h=hMix + float(blocks["dhRef"][n]),
        Yc=axes["C"][iC]*axes["Cnorm"][iZ],
        Y={s: float(blocks["Y_" + s][n]) for s in species},
    )

cells = [col[zidx[i % NX]] for i in range(NX*NY)]
print(f"  c = {axes['C'][iC]:.4f}   Z {cells[0]['Z']:.4f} .. {cells[NX-1]['Z']:.4f}")
print(f"  T {min(c['T'] for c in cells):.1f} .. {max(c['T'] for c in cells):.1f} K"
      f"   h {min(c['h'] for c in cells):.4g} .. {max(c['h'] for c in cells):.4g} J/kg")

BC = ("boundaryField {{ left {{ type {b}; }} right {{ type {b}; }} "
      "top {{ type {b}; }} bottom {{ type {b}; }} empties {{ type empty; }} }}\n")
os.makedirs(os.path.join(D, "0"), exist_ok=True)
def wr(obj, dim, vals, cls="volScalarField", bc="zeroGradient"):
    body = ("uniform " + vals if isinstance(vals, str)
            else "nonuniform List<scalar>\n%d\n(\n%s\n)" % (len(vals),
                 "\n".join(repr(v) for v in vals)))
    open(os.path.join(D, "0", obj), "w").write(
        f"FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n"
        f"dimensions {dim};\ninternalField   {body};\n" + BC.format(b=bc))

wr("Z", "[0 0 0 0 0 0 0]", [c["Z"] for c in cells])
wr("Yc", "[0 0 0 0 0 0 0]", [c["Yc"] for c in cells])
wr("h", "[0 2 -2 0 0 0 0]", [c["h"] for c in cells])
wr("T", "[0 0 0 1 0 0 0]", [c["T"] for c in cells])
wr("Zvar", "[0 0 0 0 0 0 0]", "0")
wr("p", "[1 -1 -2 0 0 0 0]", repr(pRef))
open(os.path.join(D, "0", "U"), "w").write(
    "FoamFile { version 2.0; format ascii; class volVectorField; object U; }\n"
    "dimensions [0 1 -1 0 0 0 0];\ninternalField   uniform (0 0 0);\n"
    + BC.format(b="fixedValue; value uniform (0 0 0)"))
for s in species:
    wr(s, "[0 0 0 0 0 0 0]", [c["Y"][s] for c in cells])
wr("Ydefault", "[0 0 0 0 0 0 0]", "0")
print(f"  wrote {len(species)+8} fields, {NX*NY} cells -> {D}/0")
