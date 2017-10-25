/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "CPLSocketFOAM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    CPLSocketFOAM CPL;
    CPL.initComms(argc, argv);
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"
    
	// MPI_Init is called somewhere in the PStream library
    CPL.initCFD(runTime, mesh);
    //Foam::dimensionedScalar mu(CPL.CPLDensity*nu);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //dimensionedScalar rho("rho",  dimensionSet(1, -3, 0, 0, 0, 0, 0), 1.0);
    scalar rho=1000.0; //Water density

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {

        //Pack data to send
        CPL.pack(U, p, nu, mesh, CPL.VEL | CPL.GRADPRESSURE | CPL.DIVSTRESS);
        CPL.send();

        //Recieve data from particle code
        CPL.recv();
        //CPL.unpackPorousForce(F, eps, mesh);
        CPL.unpackPorousVelForceCoeff(Up, F, dragCoef, eps, mesh);

        //Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Get momentum divided by eps
        F = CPL.divideFieldsVectorbyScalar(F, eps, mesh);
        F.correctBoundaryConditions();

        //Get Stress times grad eps term
	    //Foam::volSymmTensorField sigma(nu*2*dev(symm(fvc::grad(U))));
        //Foam::volVectorField sigmagradeps(CPL.divideFieldsVectorbyScalar(sigma*fvc::grad(eps), eps, mesh));

        //Main part of the NS equation with no pressure solve
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
          - F/rho               //Explicit Force
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());

            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn.H();
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                (fvc::interpolate(HbyA) & mesh.Sf())
               + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
