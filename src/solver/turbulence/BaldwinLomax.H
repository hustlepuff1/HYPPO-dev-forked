/*==============================================================================================================*/
/*----------------------------------------------Baldwin-Lomax model---------------------------------------------*/
/*==============================================================================================================*/
/*
For more details please refer to "Thin Layer Approximation and Algebraic Model
for Separated Turbulent Flows" B.S. Baldwin and H. Lomax, NASA Ames Research Center
*/

// --- Computes enthalpy far from the wall
turbulence.farEnthalpy(mesh,H);

// --- Computes the vorticity
magVorticity = turbulence.magVorticity(U);

// --- Computes uTau
turbulence.uTauBoundary(uTau,mu,U,rho,patches);


turbulence.boundaryWall(mesh,uTau,Cf,patches,rhoW,muW);

// --- Manage a parallel computation
turbulence.parControl();

// --- Computes the distance to the wall
turbulence.wallDistance(mesh,uTau,C,rhoW,muW);

// --- Computes y+
yPlus = turbulence.yPlus(uTau,mesh,rhoW,muW);

// --- Computes F
F = turbulence.F(magVorticity,yPlus,mesh,Aplus);

// --- Computes the turbulence using Baldwin-Lomax model
turbulence.baldwinLomax(mesh,yPlus,muTinner,muTouter,Ckleb,shockDistance,muT,rho,Aplus,H,Cf,C,patches,magVorticity);

