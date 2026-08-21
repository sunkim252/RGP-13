# peqsi1d_fgm — multicomponent pressure-equilibrium contact (FGM stage 2a gate)

The wiring test for the manifold path.  A smooth periodic Z pulse
(0 -> 1 -> 0, tanh width 0.05) advects at uniform u = 1 m/s in uniform
p = 52.5 bar, with h set exactly on the mixing line
(1-Z) hOx + Z hFuel.  Nothing should happen: composition rides along,
pressure stays flat, mass and energy stay put.  That is the whole test.

Passing numbers (100 cells, 2000 steps):

    mass rel     -8.9e-14
    rho*h rel    -4.8e-11
    dp            0.036 Pa      (case A, single species: 0.04 Pa)

Reads the verification bake `peqsi_dpm2_lo30_coef.dict` as
constant/fgmProperties (30 species + PEQSI coefficient blocks + W) and
the 106-species Wang-2011 thermo from the rd0110 case.

The FGM coefficient diagnostic is expected to look BAD here and that is
not a failure: this contact sits at 143-283 K on the dh axis origin,
where the verification table's self-reproduction error is largest
(|T_tbl - T| ~ 22 K, dBeta ~ 97%).  Those numbers are reported, never
consumed -- the solver runs on the runtime SRK evaluation -- so the
wiring, conservation and stability verdict is independent of table
quality.  Switching the solver onto the tabulated coefficients waits
for the production bake.

Run: needs libWENOEXT.so in controlDict libs (the WENOUpwindFit01
bounded scheme transports Z).
