 if (runTime.outputTime())
{
    #include "CourantNo.H"
    #include "alphaEqn.H"

    //Update dPdr to write out
    dPdr = fvc::grad(p);

    volVectorField Ur
    (
        IOobject
        (
            "Ur",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Ua - Ub
    );

    runTime.write();

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = "   << runTime.elapsedClockTime() << " s" << endl;
}
