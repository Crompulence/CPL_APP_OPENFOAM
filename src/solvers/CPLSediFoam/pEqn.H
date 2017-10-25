        // --- PISO loop
        volScalarField rUbA = 1.0/UbEqn.A()*beta;

        for (int corr = 0; corr < nCorr; corr++)
        {
            surfaceScalarField alphaf = fvc::interpolate(alpha);
            surfaceScalarField betaf = scalar(1) - alphaf;
            surfaceScalarField rUbAf = fvc::interpolate(rUbA);
            Ub = rUbA*UbEqn.H()/beta;
            surfaceScalarField phiDragb =
                fvc::interpolate(rUbA/rhob)*(fvc::interpolate(F) & mesh.Sf())
              + rUbAf*(g & mesh.Sf());

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

            for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(Dp, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve();
                // Not completely clear about how the nonOrth correction
                // should be modified. Need more thinking.
                if (nonOrth == nNonOrthCorr)
                {
                    // Info<< "Dp: " << Dp << endl;
                    surfaceScalarField SfGradp = pEqn.flux()/Dp;

                    // phia -= rUaAf*SfGradp/rhoa;
                    phib -= rUbAf*SfGradp/rhob;
                    phi = alphaf*phia + betaf*phib;

                    p.relax();
                    SfGradp = pEqn.flux()/Dp;
                    Ub += (fvc::reconstruct(phiDragb - rUbAf*SfGradp/rhob));
                    Ub.correctBoundaryConditions();
                    U = alpha*Ua + beta*Ub;
                }
            }
        }

        gradP.adjust(rUbA);

        Ub.correctBoundaryConditions();

        volTensorField gradUb("gradUb",fvc::grad(Ub));
        volScalarField nuEff("nuEff",continuousPhaseTurbulence->nuEff());
        volScalarField k("k",continuousPhaseTurbulence->k());

        B = ((2.0/3.0)*I)*k - nuEff*(twoSymm(gradUb));
