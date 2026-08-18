# Production config — recommended stack (2026-08-16 critical review)

Evidence: 3D injector A/B (L3/JA/JB/JH1, cirius 96-128 ranks, 300-step legs from
checkpoint 2000), PEP 8000-step long-run audit, Mayer case-3 benchmark
(mayer2d M0/M1/M2 + wedge-band + 3D sector). See Obsidian
"RGP-13 Adversarial Code Review" §2026-08-16 for the full review.

## fvSolution PIMPLE — dissipation stack

```
// thermoPerCorrector false: hoist the manifold+EOS refresh from per-pressure-
// corrector (6x/step) to per-outer (3x/step). Measured -31% s/step (A/B
// 2026-08-17), mild contErr increase (+36% mean), no spike in 60-step leg.
// Validate on a 500+ step run before making it the campaign default.
thermoPerCorrector false;

// consistency pair (coefficient-free, keep on)
rAUfDensityWeighted true;
psisFreezeStep  true;
pepAdvectUpwind true;

// LAD — keep for 3D hygiene (removal accumulates interface wiggles: boosted
// faces +32%/300 steps at LAD=0). Sweep ladder JA/JH1/JH2/JB + Mayer M3 shows
// monotone tradeoff; 0.15/0.25 halves plateau erosion at ~baseline hygiene.
// Do NOT set to 0 in production.
LADCoeff        0.15;   // sweep verdict 2026-08-16 (JH1/JH2/M3 ladder)
LADUCoeff       0.25;   // halves Mayer plateau erosion, 3D hygiene ~baseline
LADrhoCoeff     0;
LADbulkCoeff    0.5;

// JST sensor-gated face dissipation — validated harmless-or-better everywhere:
// inert on healthy fields (Mayer: 0 boosted faces), -82% spurious dp in 1D,
// contErr -9% in 3D. Always-on recommended.
jstFaceDissipation true;
jstCoeff        16;
jstBoostMax     8;
```

## constant/momentumTransport — REQUIRED FIX

Production 3M/5Mv2-family cases currently run `simulationType laminar`
(no SGS at all). Mayer benchmark proved resolved-fluctuation-free fields get
zero turbulent mixing — the current 3D mixing is numerics + LAD only.

```
simulationType  LES;
LES
{
    model           WALE;      // or SIGMA via coded FO (see project note)
    turbulence      on;
    delta           cubeRootVol;
}
```
Requires 0/nut (calculated/zeroGradient walls) and LES thermophysicalTransport:
```
LES { model eddyDiffusivity; Prt 0.85; }
```

## system/fvConstraints — dead entry

`limitT` (limitTemperature) binds to `thermo.he()`; fgmFluid never solves he,
so the 50-4000 K clamp HAS NEVER ACTED (runtime warning confirms). Effective T
bounding comes from the manifold table axis range. Remove the entry or leave a
comment; do not rely on it.

`limitP min 1e6`, `pMinPa 3.5e6`, `limitU max 800`: retained — measured margins
(swirl low-pressure core 43.7 bar, |U|max 221 m/s) are comfortable.

## Transport — steam kappa HOLD

The 2026-08-12 Chung B6/B7 coefficient fix removed the negative-kappa pathology
but WORSENED NIST agreement for steam (kappa MAPE 40 -> 55 % at 250 bar, dense
branch saturates the undocumented clamp `min(max(kappa,1e-4),1.0)` at
chungTransportI.H:383). Re-derive against the original Chung 1988 tables before
any steam-diluted (Paper B) run; long-term fix is an IAPWS correlation for the
H2O component. O2/kerosene/N2/H2 are unaffected (polar terms vanish).

## Tables — FPV formulation

Production tables still rely on the error-cancellation pair (transportYc off +
legacy dhRef step). Regenerate with `transportYc` + `--eq-dh-file` two-pass
before any quantitative reacting run.
