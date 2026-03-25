/*=========================================================================================================*/
/*---------------------------------------------Energy equation---------------------------------------------*/
/*=========================================================================================================*/
fvScalarMatrix rhoEEqn
(
//fvm::ddt(rhoE)  == -fvc::div(phiEp) + fvc::div(phiEpv)
fvm::ddt(rhoE) + fvm::div(jacobianRhoE,rhoE) - fvc::div(jacobianRhoE,rhoE)- fvm::div(jacobianRhoEv,rhoE) + fvc::div(jacobianRhoEv,rhoE)  == -fvc::div(phiEp) + fvc::div(phiEpv)
 );

// --- Insert equation in rhoUEEqn
rhoUEEqn.insertEquation(4, rhoEEqn);

// --- Just to fill the vector 6
fvScalarMatrix zero
(
  fvm::ddt(zeroVec6) == fvc::ddt(zeroVec6)
);

rhoUEEqn.insertEquation(5,zero);
