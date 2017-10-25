        fvVectorMatrix UbEqn(Ub, Ub.dimensions()*dimVol/dimTime);
        betaf = fvc::interpolate(beta);
        betaPhib = betaf*phib;

        #ifdef RANDOM_TURB
            #include "calcDNSForce.H"
        #endif

        {
            UbEqn =
            (
                fvm::ddt(beta, Ub)
              + fvm::div(betaPhib, Ub, "div(phib,Ub)")
              - fvm::Sp(fvc::ddt(beta) + fvc::div(betaPhib), Ub)
              + (Cvm*rhob*alpha*beta/rhob)*
                (
                    fvm::ddt(Ub)
                  + fvm::div(phib, Ub, "div(phib,Ub)")
                  - fvm::Sp(fvc::div(phib), Ub)
                )

                // divDevReff(U) = - laplacian(beta*U)
              + continuousPhaseTurbulence->divDevReff(Ub)
              + continuousPhaseTurbulence->nuEff()*(fvc::grad(beta) & fvc::grad(Ub))
             ==
            //  g                                 // Buoyancy term transfered to p-equation
            // See H. Xiao and J. Sun / Commun. Comput. Phys., 9 (2011), pp. 297-323 
            // For explaination of various terms omega and A from Cloud
            // DragCoef*U_b where U_b is fluid velocity 
              - beta*fvm::Sp(dragCoef/rhob, Ub)   // Implicit drag transfered to p-equation
            //+ alpha/rhob*dragCoef*Ua             // Explicit drag transfered to p-equation
	        + beta*alpha/rhob*(liftCoeff + Cvm*rhob*DDtUa)
            //  + beta*alpha/rhob*liftCoeff
            //  + beta*alpha*Cvm*rhob*DDtUa
              + fvc::average(beta)*gradP.flowDirection()*gradP.value() // Adding pressure gradient
            #ifdef RANDOM_TURB
              + fvc::average(beta)*turbulenceForce // Adding DNS force
            #endif
            );

            if (addIBMForce)
            {
                UbEqn -= fvm::Sp(-ibmIndicatorPtr()/ibmRelaxTime, Ub);
            }

            UbEqn.relax();

        }
