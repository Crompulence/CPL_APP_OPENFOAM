/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    lammpsFoam

Description
    Solver for a system of an incompressible fluid phase with one
    phase dispersed, i.e. particles in a liquid. The dispersed phase
    is modeled with DEM (LAMMPS)

If you are interested in a partial description of the numerical methodology 
used in bubbleFoam/twoPhaseEulerFoam, and in particular on the phase intensive 
formulation of the momentum equation, you can read the following paper:

    P. J. Oliveira, R. I. Issa, Numerical aspects of an algorithm for the
    Eulerian simulation, of two-phase flows, International Journal of Numerical 
    Methods in Fluids, 2003; 43:1177–1198 (DOI: 10.1002/fld.508)


The details of the implementation in OpenFOAM(r) are described in an internal 
report of OpenCFD (NOT AVAILABLE ANYWHERE - GREAT SCIENCE!) and are summed
up in Henrik Rusche PhD thesis.

The basic ideas behind the solution algorithm are the following

    > do not include the pressure gradient in the momentum predictor, 
      which is not solved (only a Jacobi iteration is done to find a 
      guess of the velocity before solving for the pressure equation 
      based on the mixture)
    > move the gravity and the explicit part of the drag term to the 
      pressure equation (this approach is known in the literature as 
      semi-implicit coupling)
    > solve the pressure equation to enforce mass conservation for the mixture
    > do not correct the velocity field directly, instead correct the flux
    > obtain the velocity correction from a reconstruction of the flux correction

The last two steps are the key point, if you want to have a stable solution when sharp 
interfaces with large density gradients are present.

Description from 
https://www.cfd-online.com/Forums/openfoam-solving/58178-twophaseeulerfoam-documentation-2.html

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Kmesh.H"
#include "UOprocess.H"
#include "fft.H"
#include "singlePhaseTransportModel.H"
#include "PhaseIncompressibleTurbulenceModel.H"
#include "nearWallDist.H"
#include "wallFvPatch.H"
#include "Switch.H"
#include "CPLSocketFOAM.H"

//#include "enhancedCloud.H"
//#include "chPressureGrad.H"

// #define RANDOM_TURB
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    Info << "start main = " << endl;
    CPLSocketFOAM CPL;
    CPL.initComms(argc, argv);

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "createFields.H"

    scalar t0 = runTime.elapsedCpuTime();
    //#include "createParticles.H"
    // Also update/create MD-related fields.
    beta = scalar(1) - alpha;
    volScalarField dragCoef = alpha*dimensionedScalar("dum", dimensionSet(1, -3, -1, 0, 0), 0.0);
    // Lift (place holder only. Changed in fluid loop)
    volVectorField liftCoeff = Cl*beta*rhob*(Ub ^ fvc::curl(U));
    #include "initContinuityErrs.H"

	// MPI_Init is called somewhere in the PStream library
    CPL.initCFD(runTime, mesh);
    scalarList splitTime(5,0.0);

    Info<< "\nStarting time loop\n" << endl;
    #include "liftDragCoeffs.H"

    splitTime[1] += runTime.elapsedCpuTime() - t0;

    while (runTime.run())
    {
        t0 = runTime.elapsedCpuTime();
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        //Packup and send velocity, gradident of pressure and divergence of stress
        CPL.pack(Ub, p, nu, mesh, CPL.VEL | CPL.GRADPRESSURE | CPL.DIVSTRESS);
        CPL.send();

        //Recieve and unpack particle velocity, force, 
        //sum of force weightings and porosity
        CPL.recv();
        CPL.unpackPorousVelForceCoeff(Ua, F, dragCoef, beta, mesh);
        alpha = scalar(1) - beta;

        #include "readPISO.H"
        #include "CourantNo.H"
        #include "alphaEqn.H"

        fvVectorMatrix UbEqn(Ub, Ub.dimensions()*dimVol/dimTime);
        betaf = fvc::interpolate(beta);
        betaPhib = betaf*phib;

        // See H. Xiao and J. Sun / Commun. Comput. Phys., 9 (2011), pp. 297-323 
        // For explaination of various terms omega and A from Cloud
        // \sum DragCoef*U_b where U_b is fluid velocity 
        UbEqn =
        (
            fvm::ddt(beta, Ub)
          + fvm::div(betaPhib, Ub, "div(phib,Ub)")
          - fvm::Sp(fvc::ddt(beta) + fvc::div(betaPhib), Ub)
          // This term is different in Anderson & Jackson or Kafui et al: e*du VS. d(e*u)
          - fvm::laplacian(nub*beta, Ub)
         ==
          - beta*fvm::Sp(dragCoef/rhob, Ub)   // Implicit drag transfered to p-equation
        );
        UbEqn.relax();

        // --- PISO loop
        volScalarField rUbA = 1.0/UbEqn.A()*beta;

        for (int corr = 0; corr < nCorr; corr++)
        {
            surfaceScalarField alphaf = fvc::interpolate(alpha);
            surfaceScalarField betaf = scalar(1) - alphaf;
            surfaceScalarField rUbAf = fvc::interpolate(rUbA);
            Ub = rUbA*UbEqn.H()/beta;
            surfaceScalarField phiDragb = fvc::interpolate(rUbA/rhob) 
                                         *(fvc::interpolate(F) & mesh.Sf())
                                         + rUbAf*(g & mesh.Sf());   //This term here applies gravity!!

            forAll(p.boundaryField(), patchi)
            {
                if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
                {
                    phiDragb.boundaryField()[patchi] = 0.0;
                }
            }

            Ua.correctBoundaryConditions();

            phia = (fvc::interpolate(Ua) & mesh.Sf());
            phib = (fvc::interpolate(Ub) & mesh.Sf())
                  + rUbAf*fvc::ddtCorr(Ub, phib)
                  + phiDragb;

            phi = alphaf*phia + betaf*phib;

            surfaceScalarField Dp("(rho*(1|A(U)))", betaf*rUbAf/rhob);

            //Assume a uniform grid so non-othogonal corrections    
            fvScalarMatrix pEqn
            (
                fvm::laplacian(Dp, p) == fvc::div(phi)
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            surfaceScalarField SfGradp = pEqn.flux()/Dp;

            phib -= rUbAf*SfGradp/rhob;
            phi = alphaf*phia + betaf*phib;

            p.relax();
            SfGradp = pEqn.flux()/Dp;
            Ub += (fvc::reconstruct(phiDragb - rUbAf*SfGradp/rhob));
            Ub.correctBoundaryConditions();
            U = alpha*Ua + beta*Ub;
        }

        Ub.correctBoundaryConditions();

        {
            DDtUa =
                fvc::ddt(Ua)
              + fvc::div(phia, Ua)
              - fvc::div(phia)*Ua;

            DDtUb =
                fvc::ddt(Ub)
              + fvc::div(phib, Ub)
              - fvc::div(phib)*Ub;
        }


        //Define dPdr to write out
        dPdr = fvc::grad(p);
        
        splitTime[0] += runTime.elapsedCpuTime() - t0;
        t0 = runTime.elapsedCpuTime();
        splitTime[1] += runTime.elapsedCpuTime() - t0;
        t0 = runTime.elapsedCpuTime();

        #include "liftDragCoeffs.H"
        #include "write.H"

        splitTime[2] += runTime.elapsedCpuTime() - t0;
        t0 = runTime.elapsedCpuTime();

        #include "writeCPUTime.H"

    }

    Info<< "End\n" << endl;

    if (! Pstream::parRun())  MPI_Finalize();
    return(0);
}


// ************************************************************************* //
