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
report of OpenCFD (NOT AVAILABLE ANYWHERE) but are apparently summed
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
#include "CPLSocketFOAM.H"

int main(int argc, char *argv[])
{

    //This Command turns off solver output 
    solverPerformance::debug=0;

    //Check if coupled based on cpl/COUPLER.in input file
    bool coupled;
    if (file_exists("./cpl/COUPLER.in")) {
        Info<< "Assuming coupled run as cpl/COUPLER.in exists\n" << endl;
        coupled=true;
    } else {
        Info<< "Assuming uncoupled run as cpl/COUPLER.in does not exist\n" << endl;
        coupled=false;
    }

    // Create a CPL object (not used if uncoupled) and intialise MPI
    CPLSocketFOAM CPL;
    MPI_Init(&argc, &argv);
    if (coupled)
        CPL.initComms(argc, argv);

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "createFields.H"
    #include "readPISO.H"
    #include "initContinuityErrs.H"

    // Also update/create MD-related fields.
    beta = scalar(1) - alpha;
    volScalarField dragCoef = alpha*dimensionedScalar("dum", dimensionSet(1, -3, -1, 0, 0), 0.0);

	// MPI_Init is called somewhere in the PStream library
    if (coupled)
        CPL.initCFD(runTime, mesh);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime++;
        
        if (runTime.outputTime())
            Info<< "Time = " << runTime.timeName() << endl;

        if (coupled){
            //Packup and send velocity, gradident of pressure and divergence of stress
            CPL.pack(Ub, p, nub, mesh, CPL.VEL | CPL.GRADPRESSURE | CPL.DIVSTRESS);
            CPL.send();

            //Recieve and unpack particle velocity, force, 
            //sum of force weightings and porosity
            CPL.recv();
            CPL.unpackPorousVelForceCoeff(Ua, F, dragCoef, beta, maxPossibleAlpha, mesh);
            alpha = scalar(1) - beta;
        }

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


        // E.S. A full solve seems to be needed here in place of relax to give correct answer
        // for plain Couette flow solver
        solve(UbEqn == -fvc::grad(p));

        // Only a Jacobi iteration is done to find a guess of the velocity 
        // before solving for the pressure equation based on the mixture
        //UbEqn.relax();

        // --- PISO loop
        volScalarField rUbA = 1.0/UbEqn.A()*beta;

        // Iterate over number of nCorr specified by PISO input
        for (int corr = 0; corr < nCorr; corr++)
        {
            surfaceScalarField alphaf = fvc::interpolate(alpha);
            surfaceScalarField betaf = scalar(1) - alphaf;
            surfaceScalarField rUbAf = fvc::interpolate(rUbA);
            Ub = rUbA*UbEqn.H()/beta;

            // The gravity and explicit part of drag are moved to the
            // pressure equation (this is known as semi-implicit coupling)
            surfaceScalarField phiDragb = fvc::interpolate(rUbA/rhob) 
                                         *(fvc::interpolate(F) & mesh.Sf())
                                         + rUbAf*(g & mesh.Sf());  
            forAll(p.boundaryField(), patchi)
            {
                if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
                    phiDragb.boundaryField()[patchi] = 0.0;
            }
            Ua.correctBoundaryConditions();

            // Solve the pressure equation to enforce mass conservation for the mixture
            phia = (fvc::interpolate(Ua) & mesh.Sf());
            phib = (fvc::interpolate(Ub) & mesh.Sf())
                  + rUbAf*fvc::ddtCorr(Ub, phib)
                  + phiDragb;
            phi = alphaf*phia + betaf*phib;
            surfaceScalarField Dp("(rhob*(1|A(U)))", betaf*rUbAf/rhob);
            fvScalarMatrix pEqn(fvm::laplacian(Dp, p) == fvc::div(phi));
            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            // Do not correct the velocity field directly, instead correct the flux
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

        #include "write.H"

    }

    Info<< "End\n" << endl;

    CPL.finalize();
    if (! Pstream::parRun())  MPI_Finalize();
    return(0);
}


// ************************************************************************* //
