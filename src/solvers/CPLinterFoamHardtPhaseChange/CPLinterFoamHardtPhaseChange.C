/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    interFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "CPLSocketFOAM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  Info << "Hello, this is CPLinterFoamHardtPhaseChange from ~/demo/CPL..." << endl;
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    //Check if coupled based on cpl/COUPLER.in input file
    bool coupled;
    if (file_exists("./cpl/COUPLER.in")) {
        Info<< "Assuming coupled run as cpl/COUPLER.in exists\n" << endl;
        coupled=true;
    } else {
        Info<< "Assuming uncoupled run as cpl/COUPLER.in does not exist\n" << endl;
        coupled=false;
    }

    // Create a CPL object (not used if uncoupled) and intialises MPI
    CPLSocketFOAM CPL;

    if (coupled)
        CPL.initComms(argc, argv);
    
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    // MPI_Init is called somewhere in the PStream library
    if (coupled)
      CPL.initCFD(runTime, mesh);
    
    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Initial communication to initialize domains
    if (coupled){
      CPL.pack(U, p, nu, mesh, CPL.VEL);
      CPL.send();
      CPL.recvVelocityPressure();
      CPL.unpackVelocityVOF(U, alpha1, alpha2, 
			    rho1, rho2, 
			    T, mesh);
    }
    
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

	    //Not sure whether it should stay here or after next mixture.correct()
	    if (coupled){
	      CPL.pack(U, p, nu, mesh, CPL.VEL);
	      //CPL.pack(U, p, nu, mesh, CPL.STRESS);
	      CPL.send();
	      CPL.recvVelocityPressure();
	      CPL.unpackVelocityVOF(U, alpha1, alpha2, 
				    rho1, rho2, 
				    T, mesh);
	    }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "calcMassSource.H"

            #include "UEqn.H"

            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    if (coupled)
      CPL.finalize();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
