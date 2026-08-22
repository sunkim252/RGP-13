# peqsi2d_counterflow — 2-D reacting case with gradients

The gap every earlier reacting test left open.  `peqsi1d_react` is
uniform, so the two halves of the S_Y pairing see identical states and
cancel trivially; `peqsi1d_reactCP` has gradients but is open, so the
energy budget needs flux accounting.  This one has strong gradients AND
is closed, so the internal energy must hold exactly and the reference
needs no benchmark data.

## Setup

100 x 20 cells, 2 mm x 0.4 mm, 2000 cells, ~45 s serial.  Geometry and
boundary pattern follow `testCases/counterflowDiffusion_FGM`; the physics
does not -- that case is `multicomponentFluid` on `hePsiThermo` with a
different table, so only the layout carries over.

Z runs linearly 0 -> 1 across x at fixed progress variable c = 0.80, on
the burning branch.  Every cell sits on a manifold NODE in (gZ, c, dh)
and only Z varies, so the initial state is exactly on the table and
carries no interpolation error of its own.  The result is a strong
composition and enthalpy gradient with reaction everywhere:

    T   285.7 .. 2829.8 K
    h   -1.967e+06 .. +1.725e+05 J/kg

Closed on all sides.  y is uniform, so top/bottom are cyclic -- that is
physically identical to walls here and gives WENO a periodic direction to
extend into.  With 20 cells and non-cyclic ends the WENO bank build
overruns and segfaults in `WENOBase::buildIntBasFlat`; the schemes are
kept identical to the other PEQSI cases rather than downgraded.

## What it establishes

    peqsiCompSource   max |d(int rho e)/int rho e|
    no                9.96e-08
    yes               2.65e-04

The base scheme conserves internal energy on a smooth gradient.  The
composition source breaks it by a factor of 2700.

That corrects a reading taken from a 1-D case with a DISCONTINUOUS
initial enthalpy (half at c = 0.80, half at c = 0.10), where the source
appeared to account for only half the leak.  The other half was the
discontinuity, which the base scheme handles badly on its own.  On a
smooth gradient the base level is 1e-7 and the source is the whole
signal.

## Why

The pairing itself is exact.  `iL*S_Y` is added to `Lp` and `Lrh` with
the same value, so within an RK stage `rho e = rho h - p` gains `X - X`.

The loss is downstream, in the acoustic update.  With
`h = (rhoStar hStar - (coef hN - 1) dp)/rho` and `p = pStar + dp`,

    d(rho e) = -coef * hN * dp

per step, and unlike the mass identity `int(coef dp) = dt int(sComp)`
that term is weighted by h and does not telescope.  It is proportional to
`dp`, which is why it sat at 1e-7 before: nothing was driving `dp`
hard.  S_Y raises `dp` by orders of magnitude and the pre-existing defect
comes with it.

So the source exposes a property of the fractional step rather than
introducing one -- but the practical consequence is the source's, and
`peqsiCompSource` stays off by default until it is resolved.  Fixing it
means revisiting the Eq. 11-13 update pair, which is the scheme's core.

## Run

    ./Allrun

Regenerate the initial fields with `python3 gen_ic.py . 0.80`.
