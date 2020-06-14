/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverPiso

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CPLSocketFOAM.H"
/*
#include "singlePhaseTransportModel.H"

#include "OFversion.H"
#if defined(version30)
    #include "turbulentTransportModel.H"
    #include "pisoControl.H"
#else
    #include "turbulenceModel.H"
#endif
#if defined(versionv1606plus) || defined(version40)
    #include "fvOptions.H"
#else
    #include "fvIOoptionList.H"
#endif
#include "fixedFluxPressureFvPatchScalarField.H"
#include "cfdemCloud.H"

#if defined(anisotropicRotation)
    #include "cfdemCloudRotation.H"
#endif
#if defined(superquadrics_flag)
    #include "cfdemCloudRotationSuperquadric.H"
#endif
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    /*
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #if defined(version30)
        pisoControl piso(mesh);
        #include "createTimeControls.H"
    #endif
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    create cfdemCloud
    #include "readGravitationalAcceleration.H"
    #include "checkImCoupleM.H"
    #if defined(anisotropicRotation)
        cfdemCloudRotation particleCloud(mesh);
    #elif defined(superquadrics_flag)
        cfdemCloudRotationSuperquadric particleCloud(mesh);
    #else
        cfdemCloud particleCloud(mesh);
    #endif
    #include "checkModelType.H"
    */

    // This Command turns off solver output 
    solverPerformance::debug=0;

    // Create a CPL object (not used if uncoupled) and intialises MPI
    CPLSocketFOAM CPL;
    MPI_Init(&argc, &argv);
    CPL.initComms(argc, argv);

    // Assume coupled simulation
    bool coupled=true;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "createFields.H"
    #include "readPISO.H"
    #include "initContinuityErrs.H"
    
    // Declare required MD fields
    volScalarField dragCoef = alpha*dimensionedScalar("dum", dimensionSet(1, -3, -1, 0, 0), 0.0);

    // Initialise CFD
    CPL.initCFD(runTime, mesh);
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Main time loop
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.run())
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        /*
        #if defined(version30)
            #include "readTimeControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"
        #else
            #include "readPISOControls.H"
            #include "CourantNo.H"
        #endif
        
        // do particle stuff
        particleCloud.clockM().start(1,"Global");
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(beta,Ua,Ub);

        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).impParticleForces());
        }
    
        dragCoef = particleCloud.momCoupleM(particleCloud.registryM().getProperty("implicitCouple_index")).impMomSource();
        dragCoef.correctBoundaryConditions();
        */

        // Pack and send CFD fields
        CPL.pack(Ub, p, nub, mesh, CPL.VEL | CPL.GRADPRESSURE | CPL.DIVSTRESS);
        CPL.send();

        // Recieve and unpack MD fields
        CPL.recv();
        CPL.unpackPorousVelForceCoeff(Ua, F, dragCoef, beta, maxPossibleAlpha, mesh);
        alpha = scalar(1) - beta;
        surfaceScalarField betaf = fvc::interpolate(beta);
        betaPhib = betaf*phib;

        /*
        //Force Checks
        #include "forceCheckIm.H"

        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");
        */

        // if(particleCloud.solveFlow())
        {
            // Pressure-velocity PISO corrector
            {
                // Momentum predictor
                fvVectorMatrix UbEqn
                (
                    fvm::ddt(beta,Ub) - fvm::Sp(fvc::ddt(beta),Ub)
                  + fvm::div(betaPhib,Ub) - fvm::Sp(fvc::div(betaPhib),Ub)
                  // + turbulence->divDevReff(U)
                  // + particleCloud.divVoidfractionTau(U, voidfraction)
                  - fvm::laplacian(nub*beta, Ub)
                  - fvc::div(nub*beta*dev2(fvc::grad(Ub)().T()))
                  ==
                  - fvm::Sp(dragCoef/rhob,Ub)
                  // + fvOptions(Ub)
                );

                UbEqn.relax();
                // fvOptions.constrain(UbEqn);

                /*
                #if defined(version30)
                    if (piso.momentumPredictor())
                #else
                    if (momentumPredictor)
                #endif
                {
                    if (modelType=="B" || modelType=="Bfull")
                        solve(UbEqn == - fvc::grad(p) + dragCoef/rhob*Ua);
                    else
                        solve(UbEqn == - beta*fvc::grad(p) + dragCoef/rhob*Ua);

                    fvOptions.correct(Ub);
                }
                */

                solve(UbEqn == - beta*fvc::grad(p) + dragCoef/rhob*Ua);


                // --- PISO loop
                /*
                #if defined(version30)
                    while (piso.correct())
                #else
                    for (int corr=0; corr<nCorr; corr++)
                #endif
                */
                for (int corr=0; corr<nCorr; corr++)
                {
                    volScalarField rUbA = 1.0/UbEqn.A();

                    surfaceScalarField rUbAf = fvc::interpolate(rUbA);
                    volScalarField rUbAbeta = rUbA*beta;
                    surfaceScalarField rUbAbetaf = fvc::interpolate(rUbAbeta);

                    /*
                    surfaceScalarField rUbAf("(1|A(Ub))", fvc::interpolate(rUbA));
                    volScalarField rUbAbeta("(voidfraction2|A(Ub))",rUbA*beta);
                    surfaceScalarField rUbAbetaf("(voidfraction2|A(Ub)F)", fvc::interpolate(rUbAbeta));
                    */

                    Ub = rUbA*UbEqn.H();

                    /*
                    #ifdef version23
                        betaPhib = ( fvc::interpolate(Ub) & mesh.Sf() )
                            + rUbAbetaf*fvc::ddtCorr(Ub, phib);
                    #else
                        betaPhib = ( fvc::interpolate(Ub) & mesh.Sf() )
                            + fvc::ddtPhiCorr(rUbAbeta, Ub, phib);
                    #endif
                    */        
                    betaPhib = (fvc::interpolate(Ub) & mesh.Sf())
                             + rUbAbetaf*fvc::ddtCorr(Ub, phib);

                    surfaceScalarField phia = (fvc::interpolate(Ua) & mesh.Sf());
                    betaPhib += rUbAf*(fvc::interpolate(dragCoef/rhob)*phia);

                    /*
                    if (modelType=="A")
                        rUbAbeta = volScalarField("(voidfraction2|A(Ub))",rUbA*beta*beta);

                    // Update the fixedFluxPressure BCs to ensure flux consistency
                    #include "fixedFluxPressureHandling.H"
                    */

                    // Non-orthogonal pressure corrector loop
                    /*
                    #if defined(version30)
                        while (piso.correctNonOrthogonal())
                    #else
                        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                    #endif
                    */
                    rUbAbeta = rUbA*beta*beta;
                    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                    {
                        // Pressure corrector
                        fvScalarMatrix pEqn
                        (
                            // fvm::laplacian(rUbAbeta, p) == fvc::div(betaf*betaPhib) + particleCloud.ddtVoidfraction()
                            fvm::laplacian(rUbAbeta, p) == fvc::div(betaf*betaPhib) + fvc::ddt(beta)
                        );
                        pEqn.setReference(pRefCell, pRefValue);

                        /*
                        #if defined(version30)
                            pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
                            if (piso.finalNonOrthogonalIter())
                            {
                                phib = betaPhib - pEqn.flux()/betaf;
                            }
                        #else
                            if( corr == nCorr-1 && nonOrth == nNonOrthCorr )
                                #if defined(versionExt32)
                                    pEqn.solve(mesh.solutionDict().solver("pFinal"));
                                #else
                                    pEqn.solve(mesh.solver("pFinal"));
                                #endif
                            else
                                pEqn.solve();

                            if (nonOrth == nNonOrthCorr)
                            {
                                phib = betaPhib - pEqn.flux()/betaf;
                            }
                        #endif
                        */
                            
                        if (corr == nCorr-1 && nonOrth == nNonOrthCorr)
                            pEqn.solve(mesh.solver("pFinal"));    
                        else
                            pEqn.solve();

                        if (nonOrth == nNonOrthCorr)
                        {
                            phib = betaPhib - pEqn.flux()/betaf;
                        }


                    } // end non-orthogonal corrector loop

                    betaPhib = betaf*phib;
                    Ub -= beta*rUbA*fvc::grad(p) - dragCoef/rhob*Ua*rUbA;

                    /*
                    #include "continuityErrorPhiPU.H"

                    if (modelType=="B" || modelType=="Bfull")
                        Ub -= rUbA*fvc::grad(p) - dragCoef/rhob*Ua*rUbA;
                    else
                        Ub -= beta*rUbA*fvc::grad(p) - dragCoef/rhob*Ua*rUbA;
                    
                    Ub.correctBoundaryConditions();
                    fvOptions.correct(Ub);
                    */
                    Ub.correctBoundaryConditions();
                } // end piso loop
            }

            /*
            laminarTransport.correct();
            turbulence->correct();
            */
        }// end solveFlow
        /*
        else
        {
            Info << "skipping flow solution." << endl;
        }
        */

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        /*
        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
        */
    }

    Info<< "End\n" << endl;

    CPL.finalize();
    if (! Pstream::parRun())  MPI_Finalize();
    return 0;
}


// ************************************************************************* //
