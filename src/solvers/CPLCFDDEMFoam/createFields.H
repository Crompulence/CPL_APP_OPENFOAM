    Info<< "Reading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField beta
    (
        IOobject
        (
            "beta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalar(1) - alpha
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Ub\n" << endl;
    volVectorField Ub
    (
        IOobject
        (
            "Ub",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Setting up place-holder for Ua\n" << endl;
    volVectorField Ua
    (
        IOobject
        (
            "Ua",
            runTime.timeName(),
            mesh,
            // IOobject::READ_IF_PRESENT,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
         ),
        Ub - Ub
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha*Ua + beta*Ub
    );

    Info<< "Reading field F\n" << endl;
    volVectorField F
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    // default max alpha if not defined in dictionary
    scalar maxPossibleAlpha = 0.99;
    maxPossibleAlpha = transportProperties.lookupOrDefault<scalar>("maxPossibleAlpha", 0.99);
    if (!transportProperties.found("maxPossibleAlpha"))
    {
        Info<< "*** Max possible alpha set to 0.99 by default" << endl;
    }
    else
    {
        Info<< "*** Max possible alpha set to " << maxPossibleAlpha << endl;
    }

    dimensionedScalar rhob
    (
        transportProperties.lookup("rhob")
    );

    dimensionedScalar nub
    (
        transportProperties.lookup("nub")
    );


#   include "createPhia.H"
#   include "createPhib.H"

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alpha)*phia
      + fvc::interpolate(beta)*phib
    );

    volVectorField dPdr
    (
        IOobject
        (
            "dPdr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(p)
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

    surfaceScalarField betaf("betaf", fvc::interpolate(beta));
    surfaceScalarField betaPhib("betaPhib", betaf*phib);

