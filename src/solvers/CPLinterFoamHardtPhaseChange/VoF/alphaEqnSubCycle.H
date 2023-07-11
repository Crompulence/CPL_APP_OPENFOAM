if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    surfaceScalarField rhoPhiCpSum("rhoPhiCpSum", 0.0*rhoPhiCp);

    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    }

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        rhoPhiCpSum += (runTime.deltaT()/totalDeltaT)*rhoPhiCp;
    }

    rhoPhi = rhoPhiSum;
    rhoPhiCp = rhoPhiCpSum;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;
rhoCp == alpha1*rho1*cp1 + alpha2*rho2*cp2;

Info << "Total Evaporation Rate (before smoothing) = " << evapRate  << endl;
Info << "Total Evaporation Rate (after smoothing) = " << totMdot  << endl;
//Info << "N = " << N  << endl;
//Info << "Nv = " << NV  << endl;
//Info << "Nl = " << NL << nl << endl;
