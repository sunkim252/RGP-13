/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2024-2025 RGP-13
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "elyHanleyTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::elyHanleyTransport<Thermo>::elyHanleyTransport
(
    const word& name,
    const dictionary& dict
)
:
    Thermo(name, dict),
    TcM_(this->Tc()),
    // rfSpecie stores Vc in cm^3/mol; the rest of this class works in
    // SI (m^3/kmol). Convert here so VcM_ / VcO_ ratios stay dimensionless.
    VcM_(this->Vc()*1e-3),
    MM_(this->W()),
    omegaM_(this->omega()),
    ZcM_(0.0),                          // computed below
    TcO_(190.564),
    VcO_(0.0992),                       // m^3/kmol (methane)
    MO_(16.043),
    omegaO_(0.011),
    ZcO_(0.286),
    // Chung fallback auxiliaries (per-species values; overwritten by
    // SRKchungTakaMixture::calcMixture through updateTRANS for mixtures).
    sigmaM_(0.809*pow(this->Vc(), 1.0/3)),
    epsilonkM_(this->Tc()/1.2593),
    miuiM_(this->miui()),
    kappaiM_(this->kappai()),
    Ymd_(2),
    Xmd_(2),
    Tcmd_(2),
    Pcmd_(2),
    Mmd_(2),
    sigmd_(2)
{
    // Compute mixture Zc = pc Vc / (R Tc).
    // pc is not stored on this class directly; retrieve through rfSpecie Pc().
    const scalar Ru = 8.31446;
    const scalar Tc_ = this->Tc();
    const scalar Pc_ = this->Pc();
    const scalar Vc_ = this->Vc();  // m^3/kmol here?
    // Vc in rfSpecie is stored in cm^3/mol == 1e-3 m^3/kmol? Check:
    //   For O2 we have Vc = 73.529 (cm^3/mol) which equals 7.3529e-2 m^3/kmol.
    // So convert to m^3/kmol assuming stored unit is cm^3/mol:
    const scalar Vc_SI = Vc_*1e-3;          // -> m^3/kmol
    if (Tc_ > 1e-6 && Pc_ > 1e-6)
    {
        // Zc = p_c V_c / (R T_c), with V_c per kmol and R_u = 8314.46
        ZcM_ = Pc_*Vc_SI/(Ru*1e3*Tc_);
    }
    else
    {
        ZcM_ = 0.29;
    }

    // Initialise diffusivity bookkeeping symmetrically to chungTransport
    forAll(Ymd_, i)
    {
        Ymd_[i] = this->Y();
    }
    forAll(Xmd_, i)
    {
        Xmd_[i] = this->Y()/this->W();
    }
    List<scalar> TcTemp(2);
    forAll(TcTemp, i) TcTemp[i] = this->Tc();
    forAll(Tcmd_, i) Tcmd_[i] = TcTemp;

    List<scalar> PcTemp(2);
    forAll(PcTemp, i) PcTemp[i] = this->Pc();
    forAll(Pcmd_, i) Pcmd_[i] = PcTemp;

    List<scalar> MTemp(2);
    forAll(MTemp, i) MTemp[i] = this->W();
    forAll(Mmd_, i) Mmd_[i] = MTemp;

    List<scalar> sigTemp(2);
    forAll(sigTemp, i) sigTemp[i] = this->sigmvi();
    forAll(sigmd_, i) sigmd_[i] = sigTemp;
}


template<class Thermo>
Foam::elyHanleyTransport<Thermo>::elyHanleyTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t),
    TcM_(this->Tc()),
    VcM_(this->Vc()*1e-3),      // cm^3/mol -> m^3/kmol
    MM_(this->W()),
    omegaM_(this->omega()),
    ZcM_(0.0),
    TcO_(190.564), VcO_(0.0992), MO_(16.043), omegaO_(0.011), ZcO_(0.286),
    sigmaM_(0.809*pow(this->Vc(), 1.0/3)),
    epsilonkM_(this->Tc()/1.2593),
    miuiM_(this->miui()),
    kappaiM_(this->kappai()),
    Ymd_(2), Xmd_(2), Tcmd_(2), Pcmd_(2), Mmd_(2), sigmd_(2)
{
    const scalar Ru = 8.31446;
    const scalar Vc_SI = this->Vc()*1e-3;
    if (this->Tc() > 1e-6 && this->Pc() > 1e-6)
    {
        ZcM_ = this->Pc()*Vc_SI/(Ru*1e3*this->Tc());
    }
    else
    {
        ZcM_ = 0.29;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

template<class Thermo>
void Foam::elyHanleyTransport<Thermo>::write(Ostream& os) const
{
    Thermo::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const elyHanleyTransport<Thermo>& tr
)
{
    tr.write(os);
    return os;
}
