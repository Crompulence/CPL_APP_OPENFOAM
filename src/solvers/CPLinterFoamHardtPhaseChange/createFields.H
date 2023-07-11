#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

Info<< "Reading field mDot\n" << endl;
    volScalarField mDot
    (
        IOobject
        (
            "mDot",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField mDot0
    (
        IOobject
        (
            "mDot0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mDot,
        mDot.boundaryField().types()
    );

#include "createPhi.H"

IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

dimensionedScalar nu
(
	transportProperties.lookup("nu")
);

scalar CONDENSATION
(
 readScalar(transportProperties.lookup("CONDENSATION"))
 );


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

dictionary phase1 = transportProperties.subDict(mixture.phase1Name());
    dictionary phase2 = transportProperties.subDict(mixture.phase2Name());
    dictionary pcProperties = transportProperties.subDict("PhaseChangeProperties");   

    dimensionedScalar k1 (phase1.lookup("kappa"));
    dimensionedScalar cp1 (phase1.lookup("Cp"));
    dimensionedScalar hf1 (phase1.lookup("Hf"));

    dimensionedScalar k2  (phase2.lookup("kappa"));
    dimensionedScalar cp2 (phase2.lookup("Cp"));
    dimensionedScalar hf2 (phase2.lookup("Hf"));

    dimensionedScalar Tsat (pcProperties.lookup("Tsat"));
    dimensionedScalar R (pcProperties.lookup("R"));
    dimensionedScalar sigmaEvap (pcProperties.lookup("sigmaEvap"));
    dimensionedScalar sigmaCond (pcProperties.lookup("sigmaCond"));
    dimensionedScalar DmDot (pcProperties.lookup("DmDot"));


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

Info << "Initializing rhoCp" << nl << endl;
    volScalarField rhoCp
    (
        IOobject
        (
            "rhoCp",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1*cp1 + alpha2*rho2*cp2,
        alpha1.boundaryField().types()
    );
    rhoCp.oldTime();


   //Energy flux
    surfaceScalarField rhoPhiCp
    (
        IOobject
        (
            "rhoPhiCp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rhoCp)*phi
    );



// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

volScalarField continuitySource
    (
        IOobject
        (
            "continuitySource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("continuitySource",dimensionSet(1,-3,-1,0,0,0,0),0.0)
    );

    volScalarField energySource
    (
        IOobject
        (
            "energySource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("energySource",dimensionSet(1,-1,-3,-1,0,0,0),0.0)
    );

dimensionedScalar evapRate("evapRate",dimensionSet(1,0,-1,0,0,0,0), 0.0);
dimensionedScalar totMdot("totMdot",dimensionSet(1,0,-1,0,0,0,0), 0.0);
dimensionedScalar mDotL("mDotL",dimensionSet(1,0,-1,0,0,0,0), 0.0);
dimensionedScalar mDotV("mDotV",dimensionSet(1,0,-1,0,0,0,0), 0.0);

#include "createMRF.H"
#include "createFvOptions.H"
