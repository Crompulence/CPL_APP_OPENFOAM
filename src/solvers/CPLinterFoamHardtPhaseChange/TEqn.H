{

  volScalarField kappa = alpha1*k1 + alpha2*k2;

  fvScalarMatrix TEqn
   (
       fvm::ddt(rhoCp,T)
     + fvm::div(rhoPhiCp,T)
     - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
     - fvm::laplacian(kappa,T)
     - fvm::Sp(energySource,T)
//     - fvc::DDt(phi,p_rgh)
     ==
       - energySource*Tsat
   );


   TEqn.relax();
   TEqn.solve();

   Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() <<endl;

   //update continuity source with new T
   #include "calcMassSource.H"

}
