/*=========================================================================================================*/
/*---------------------------------------------Density equation--------------------------------------------*/
/*=========================================================================================================*/
fvScalarMatrix rhoEqn
(
//fvm::ddt(rho)  == -fvc::div(phi)
fvm::ddt(rho) + fvm::div(jacobianRho,rho) - fvc::div(jacobianRho,rho)  == -fvc::div(phi)
 );

// --- Insert equation in rhoUEEqn
rhoUEEqn.insertEquation(0, rhoEqn);
