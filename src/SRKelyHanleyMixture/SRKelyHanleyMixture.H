/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2024-2026 RGP-13
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::SRKelyHanleyMixture

Description
    Multicomponent mixture class with real-gas mixing rules tailored to
    the Ely-Hanley extended-corresponding-states transport backend:

        EoS:        SRK (with optional Peneloux volume translation, k_ij)
        Thermo:     janafThermo
        Transport:  Ely-Hanley ECS via methane reference, plus Fuller +
                    Takahashi for binary diffusion coefficients

    Mixture-averaged critical state (T_{c,m}, V_{c,m}, omega_m, M_m,
    Z_{c,m}) is supplied by mole-fraction-weighted (pseudo-pure)
    averaging of pure-species critical data. Z_{c,m} is computed from
    the mole-weighted (P_{c,m}, V_{c,m}, T_{c,m}) so the Enskog
    correction X_eta, X_lam inside elyHanleyTransport sees a consistent
    Z_{c,m}/Z_{c,o} ratio (not the Pitzer fallback used by the previous
    Chung-signature adapter).

    SRK pre-computed pair-interaction matrices (BM, COEF1..3, KIJ, CM)
    are kept identical to SRKchungTakaMixture, so EoS results are
    bit-for-bit interchangeable. Only the transport-side pre-computed
    matrices differ:
        - SRKchungTakaMixture stores SIGMA3M, EPSILONKM0, OMEGAM0, MM0,
          MIUIM0, KAPPAIM (Chung Lennard-Jones combining rules)
        - SRKelyHanleyMixture stores ListTc, ListVc, ListOmega — that
          is all the ECS pipeline needs.

SourceFiles
    SRKelyHanleyMixture.C

\*---------------------------------------------------------------------------*/

#ifndef SRKelyHanleyMixture_H
#define SRKelyHanleyMixture_H

#include "multicomponentMixture.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class SRKelyHanleyMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class SRKelyHanleyMixture
:
    public multicomponentMixture<ThermoType>
{
public:

    // Public Typedefs

        typedef ThermoType thermoMixtureType;
        typedef ThermoType transportMixtureType;


private:

    // Private Data

        //- Per-cell scratch mixture
        mutable thermoMixtureType mixture_;

        label numberOfSpecies_;


    // Pure-species data (constant in time)

        List<scalar> ListW_;     // molecular weights [kg/kmol]
        List<scalar> ListTc_;    // critical temperature [K]
        List<scalar> ListPc_;    // critical pressure [Pa]
        List<scalar> ListVc_;    // critical volume [cm^3/mol]
                                 // (rfSpecie native unit; the
                                 // elyHanleyTransport adapter converts
                                 // to m^3/kmol internally)
        List<scalar> ListOmega_; // acentric factor [-]


    // SRK equation-of-state pre-computed quantities (identical to
    // SRKchungTakaMixture)

        List<scalar> BM_;
        List<List<scalar>> COEF1_;
        List<List<scalar>> COEF2_;
        List<List<scalar>> COEF3_;

        //- Per-species Peneloux volume shift [m^3/kmol]
        List<scalar> CM_;

        //- Binary interaction matrix k_ij (1 - k_ij) applied to SRK
        //  cross terms. Diagonal = 0.
        List<List<scalar>> KIJ_;


    // Pre-computed matrices for Fuller / Takahashi binary diffusion
    // (identical to SRKchungTakaMixture)

        List<List<scalar>> MMD_;
        List<List<scalar>> SIGMD_;


    // Private Member Functions

        //- Compute SRK mixture parameters (a_m via three coefficients,
        //  b_m, Peneloux c_m) and Ely-Hanley mixture critical state
        //  (T_c,m, V_c,m, M_m, omega_m, Z_c,m) from mole fractions.
        void calculateRealGas
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
        ) const;

        //- Build the per-cell mixture object
        const thermoMixtureType& calcMixture
        (
            const scalarFieldListSlice& Y
        ) const;


public:

    // Constructors

        SRKelyHanleyMixture(const dictionary&);

        SRKelyHanleyMixture(const SRKelyHanleyMixture<ThermoType>&) = delete;


    //- Destructor
    virtual ~SRKelyHanleyMixture()
    {}


    // Member Functions

        static word typeName()
        {
            return "SRKelyHanleyMixture<" + ThermoType::typeName() + '>';
        }

        inline const label& numberOfSpecies() const
        {
            return numberOfSpecies_;
        }

        const thermoMixtureType& thermoMixture
        (
            const scalarFieldListSlice&
        ) const;

        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&
        ) const;

        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&,
            const thermoMixtureType&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SRKelyHanleyMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
