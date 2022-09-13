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
#include "PstreamGlobals.H"
#include "mpi.h"
#include "cpl.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    //Define variables
    bool cflag;
    int CFD_realm = 1;
    MPI_Comm CFD_COMM, CART_COMM;
    CPL::ndArray<double> send_array, recv_array;

    int flag = 0;
    int ierr = MPI_Initialized(&flag);

    Foam::Info << "MPI_Initialized(&flag) " << flag << Foam::endl;
    if (flag == 0)
		MPI_Init(&argc, &argv);

    //Initialise CPL library
    CPL::init(CFD_realm, CFD_COMM); 

    //If you want to use shared MPI_COMM_WORLD, this line is set
    Info<< "\nSetting CPLRealmComm\n" << endl;
	Foam::PstreamGlobals::CPLRealmComm = CFD_COMM;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // Initial communication to initialize domains
    int npxyz[3] = {1, 1, 1}; int periods[3] = {1, 1, 1};
    MPI_Cart_create(CFD_COMM, 3, npxyz, periods, 1, &CART_COMM);

    double xyzL[3] = {1.0, 1.0, 1.0}; double xyz_orig[3] = {0.0, 0.0, 0.0};
    int ncxyz[3] = {32, 32, 32}; 
    CPL::setup_cfd(CART_COMM, xyzL, xyz_orig, ncxyz);
    CPL::get_arrays(&recv_array, 4, &send_array, 1);

    Info<< "\nStarting time loop\n" << endl;
    int time=0;
    while (runTime.loop())
    {

        cflag = CPL::recv(&recv_array);
        std::cout << "CPLTestFoam " << time << " " << recv_array(0,0,0,0) << std::endl;
        send_array(0,0,0,0) = 2.*time;
        cflag = CPL::send(&send_array);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
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

        time += 1;
    }
    Info<< "End\n" << endl;
	CPL::finalize(); // uncomment if CPL is installed

    return 0;
}


// ************************************************************************* //
