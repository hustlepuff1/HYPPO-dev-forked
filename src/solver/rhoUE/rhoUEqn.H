/*=========================================================================================================*/
/*--------------------------------------------Momentum equation--------------------------------------------*/
/*=========================================================================================================*/
fvVectorMatrix rhoUEqn
(
 //fvm::ddt(rhoU)  == -fvc::div(phiUp) + fvc::div(phiUpv)
 fvm::ddt(rhoU) + fvm::div(jacobianRhoU,rhoU) - fvc::div(jacobianRhoU,rhoU) - fvm::div(jacobianRhoUv,rhoU) + fvc::div(jacobianRhoUv,rhoU) == -fvc::div(phiUp) + fvc::div(phiUpv)
 );

// --- Insert equation in rhoUEEqn
rhoUEEqn.insertEquation(1, rhoUEqn);
