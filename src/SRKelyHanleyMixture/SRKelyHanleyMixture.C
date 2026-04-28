/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2024-2026 RGP-13
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "SRKelyHanleyMixture.H"
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::SRKelyHanleyMixture<ThermoType>::calculateRealGas
(
    const List<scalar>& X,
    scalar& bM,
    scalar& coef1,
    scalar& coef2,
    scalar& coef3,
    scalar& cM,
    scalar& MM,
    scalar& VcM,
    scalar& TcM,
    scalar& omegaM,
    scalar& ZcM
) const
{
    // ---------------- SRK EoS mixing rules (identical to chung) ----------
    forAll(BM_, i)
    {
        bM += X[i]*BM_[i];
    }

    forAll(CM_, i)
    {
        cM += X[i]*CM_[i];
    }

    forAll(COEF1_, i)
    {
        const scalar Xi2 = X[i]*X[i];
        coef1 += Xi2*COEF1_[i][i];
        coef2 += Xi2*COEF2_[i][i];
        coef3 += Xi2*COEF3_[i][i];
    }
    for (label i = 0; i < COEF1_.size(); i++)
    {
        for (label j = i + 1; j < COEF1_.size(); j++)
        {
            const scalar twoXiXj = 2.0*X[i]*X[j];
            coef1 += twoXiXj*COEF1_[i][j];
            coef2 += twoXiXj*COEF2_[i][j];
            coef3 += twoXiXj*COEF3_[i][j];
        }
    }

    // ---------------- Ely-Hanley mixture critical state ------------------
    // Pseudo-pure mole-weighted averaging — see notes in
    // RGP-13_ECS_deviations.ipynb Item 11.
    //   T_c,m = sum_i X_i T_c,i
    //   V_c,m = sum_i X_i V_c,i
    //   omega_m = sum_i X_i omega_i
    //   M_m   = sum_i X_i W_i           (mole-weighted molar mass)
    //   P_c,m = sum_i X_i P_c,i         (auxiliary, only for Z_c,m)
    //   Z_c,m = P_c,m * V_c,m^(SI) / (R T_c,m)
    // V_c,i is stored in cm^3/mol by rfSpecie; convert to m^3/kmol for Z.
    scalar PcM = 0;
    forAll(X, i)
    {
        TcM    += X[i]*ListTc_[i];
        VcM    += X[i]*ListVc_[i];      // cm^3/mol
        omegaM += X[i]*ListOmega_[i];
        MM     += X[i]*ListW_[i];
        PcM    += X[i]*ListPc_[i];
    }

    if (TcM == 0) { TcM = 1e-30; }
    if (VcM == 0) { VcM = 1e-30; }
    if (MM  == 0) { MM  = 1e-30; }

    const scalar Ru   = Foam::constant::thermodynamic::RR;   // J/(kmol K)
    const scalar VcMSI = VcM*1e-3;                           // -> m^3/kmol
    if (PcM > 1e-6 && TcM > 1e-6)
    {
        ZcM = PcM*VcMSI/(Ru*TcM);
    }
    else
    {
        // Fall back to Pitzer if any species lacks P_c
        ZcM = 0.291 - 0.08*omegaM;
        if (ZcM <= 0) ZcM = 0.29;
    }
}


template<class ThermoType>
const typename Foam::SRKelyHanleyMixture<ThermoType>::thermoMixtureType&
Foam::SRKelyHanleyMixture<ThermoType>::calcMixture
(
    const scalarFieldListSlice& Y
) const
{
    // Mass-fraction weighted base thermo
    mixture_ = Y[0]*this->specieThermos()[0];
    for (label n = 1; n < Y.size(); n++)
    {
        mixture_ += Y[n]*this->specieThermos()[n];
    }

    const label nSpecies = Y.size();
    List<scalar> X(nSpecies);
    List<scalar> Yl(nSpecies);
    scalar sumXb = 0.0;

    forAll(X, i)
    {
        sumXb += Y[i]/ListW_[i];
    }
    if (sumXb == 0) { sumXb = 1e-30; }

    forAll(X, i)
    {
        X[i] = (Y[i]/ListW_[i])/sumXb;
        Yl[i] = Y[i];
        if (X[i]  <= 0) X[i]  = 0;
        if (Yl[i] <= 0) Yl[i] = 0;
    }

    // SRK + Ely-Hanley mixture aggregates
    scalar bM = 0, coef1 = 0, coef2 = 0, coef3 = 0, cM = 0;
    scalar MM = 0, VcM = 0, TcM = 0, omegaM = 0, ZcM = 0;

    calculateRealGas
    (
        X, bM, coef1, coef2, coef3, cM,
        MM, VcM, TcM, omegaM, ZcM
    );

    // Push SRK + Peneloux into the EoS slot
    mixture_.updateEoS(bM, coef1, coef2, coef3, cM);

    // Composition-dependent Tc, Pc matrices for Fuller + Takahashi
    // (mole-weighted pair averages; identical to SRKchungTakaMixture)
    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(X, i)
    {
        X[i] += 1e-40;
        sumXcorrected += X[i];
    }
    forAll(X, i)
    {
        X[i] /= sumXcorrected;
        WmixCorrect += X[i]*ListW_[i];
    }
    forAll(Yl, i)
    {
        Yl[i] = X[i]*ListW_[i]/WmixCorrect;
    }

    List<List<scalar>> TCMD(nSpecies);
    List<List<scalar>> PCMD(nSpecies);
    forAll(TCMD, i)
    {
        TCMD[i].setSize(nSpecies);
        PCMD[i].setSize(nSpecies);
    }
    for (label i = 0; i < nSpecies; i++)
    {
        TCMD[i][i] = ListTc_[i];
        if (TCMD[i][i] == 0) TCMD[i][i] = 1e-40;
        PCMD[i][i] = ListPc_[i];
        if (PCMD[i][i] == 0) PCMD[i][i] = 1e-40;

        for (label j = i + 1; j < nSpecies; j++)
        {
            const scalar XiPlusXj = X[i] + X[j];
            const scalar invSum = (XiPlusXj == 0) ? 0.0 : 1.0/XiPlusXj;

            scalar tcVal = (X[i]*ListTc_[i] + X[j]*ListTc_[j])*invSum;
            if (tcVal == 0) tcVal = 1e-40;
            scalar pcVal = (X[i]*ListPc_[i] + X[j]*ListPc_[j])*invSum;
            if (pcVal == 0) pcVal = 1e-40;

            TCMD[i][j] = tcVal; TCMD[j][i] = tcVal;
            PCMD[i][j] = pcVal; PCMD[j][i] = pcVal;
        }
    }

    // Drive elyHanleyTransport via its native 5-parameter signature
    mixture_.updateTRANS
    (
        TcM, VcM*1e-3, MM, omegaM, ZcM,   // VcM -> m^3/kmol for elyHanley
        Yl, X, TCMD, PCMD, MMD_, SIGMD_
    );

    return mixture_;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::SRKelyHanleyMixture<ThermoType>::SRKelyHanleyMixture
(
    const dictionary& dict
)
:
    multicomponentMixture<ThermoType>(dict),
    mixture_("mixture", this->specieThermos()[0]),
    numberOfSpecies_(this->specieThermos().size()),
    ListW_(numberOfSpecies_),
    ListTc_(numberOfSpecies_),
    ListPc_(numberOfSpecies_),
    ListVc_(numberOfSpecies_),
    ListOmega_(numberOfSpecies_),
    BM_(numberOfSpecies_),
    COEF1_(numberOfSpecies_),
    COEF2_(numberOfSpecies_),
    COEF3_(numberOfSpecies_),
    CM_(numberOfSpecies_, scalar(0)),
    KIJ_(numberOfSpecies_),
    MMD_(numberOfSpecies_),
    SIGMD_(numberOfSpecies_)
{
    const scalar RR = Foam::constant::thermodynamic::RR;

    // Pure-species pulls
    forAll(BM_, i)
    {
        ListW_[i]     = this->specieThermos()[i].W();
        ListTc_[i]    = this->specieThermos()[i].Tc();
        ListPc_[i]    = this->specieThermos()[i].Pc();
        ListVc_[i]    = this->specieThermos()[i].Vc();
        ListOmega_[i] = this->specieThermos()[i].omega();
        CM_[i]        = this->specieThermos()[i].c();
        BM_[i]        = 0.08664*RR
                       *this->specieThermos()[i].Tc()
                       /this->specieThermos()[i].Pc();
    }

    // k_ij parsing (identical to SRKchungTakaMixture)
    forAll(KIJ_, i)
    {
        KIJ_[i].setSize(numberOfSpecies_, scalar(0));
    }
    if (dict.found("binaryInteraction"))
    {
        const dictionary& bd = dict.subDict("binaryInteraction");
        forAllConstIter(dictionary, bd, iter)
        {
            if (!iter().isDict()) continue;
            const word& pair = iter().keyword();
            if (pair == "default") continue;
            const std::size_t sep = pair.find('_');
            if (sep == std::string::npos) continue;
            const word a = pair.substr(0, sep);
            const word b = pair.substr(sep + 1);
            label ia = -1, ib = -1;
            forAll(this->specieThermos(), s)
            {
                if (this->specieThermos()[s].name() == a) ia = s;
                if (this->specieThermos()[s].name() == b) ib = s;
            }
            if (ia < 0 || ib < 0) continue;
            const scalar kij = iter().dict().lookupOrDefault<scalar>("kij", 0);
            KIJ_[ia][ib] = kij;
            KIJ_[ib][ia] = kij;
        }
    }

    // Pre-compute SRK + diffusivity pair matrices
    List<scalar> nCOEF1(numberOfSpecies_);
    List<scalar> nCOEF2(numberOfSpecies_);
    List<scalar> nCOEF3(numberOfSpecies_);
    List<scalar> nMMD(numberOfSpecies_);
    List<scalar> nSIGMD(numberOfSpecies_);

    forAll(COEF1_, i)
    {
        const auto& si = this->specieThermos()[i];
        const scalar miOmega = si.omega();
        const scalar miTc    = si.Tc();
        const scalar miPc    = si.Pc();
        const scalar miW     = si.W();
        const scalar miSigmv = si.sigmvi();

        const scalar mFi = 0.48508 + 1.5517*miOmega - 0.15613*sqr(miOmega);
        const scalar Ai  = 0.42747*sqr(RR*miTc)/miPc;

        forAll(nCOEF1, j)
        {
            const auto& sj = this->specieThermos()[j];
            const scalar mjOmega = sj.omega();
            const scalar mjTc    = sj.Tc();
            const scalar mjPc    = sj.Pc();
            const scalar mjW     = sj.W();
            const scalar mjSigmv = sj.sigmvi();

            const scalar mFj = 0.48508 + 1.5517*mjOmega - 0.15613*sqr(mjOmega);
            const scalar Aj  = 0.42747*sqr(RR*mjTc)/mjPc;

            const scalar oneMinusKij = (i == j) ? 1.0 : (1.0 - KIJ_[i][j]);
            const scalar sqAij = sqrt(Ai*Aj);

            nCOEF1[j] = sqAij*(1 + mFi)*(1 + mFj)*oneMinusKij;
            nCOEF2[j] =
                sqAij
               *(
                    (1.0 + mFj)*mFi/sqrt(miTc)
                  + (1.0 + mFi)*mFj/sqrt(mjTc)
                )*oneMinusKij;
            nCOEF3[j] = sqAij*mFi*mFj/sqrt(miTc*mjTc)*oneMinusKij;

            nMMD[j]  = 1/miW + 1/mjW;
            nSIGMD[j] = pow(miSigmv, 1.0/3) + pow(mjSigmv, 1.0/3);
        }

        COEF1_[i] = nCOEF1;
        COEF2_[i] = nCOEF2;
        COEF3_[i] = nCOEF3;
        MMD_[i]   = nMMD;
        SIGMD_[i] = nSIGMD;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const typename Foam::SRKelyHanleyMixture<ThermoType>::thermoMixtureType&
Foam::SRKelyHanleyMixture<ThermoType>::thermoMixture
(
    const scalarFieldListSlice& Y
) const
{
    return calcMixture(Y);
}


template<class ThermoType>
const typename Foam::SRKelyHanleyMixture<ThermoType>::transportMixtureType&
Foam::SRKelyHanleyMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    return thermoMixture(Y);
}


template<class ThermoType>
const typename Foam::SRKelyHanleyMixture<ThermoType>::transportMixtureType&
Foam::SRKelyHanleyMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice&,
    const thermoMixtureType& mixture
) const
{
    return mixture;
}


// ************************************************************************* //
