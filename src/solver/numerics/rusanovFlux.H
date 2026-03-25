/*=========================================================================================================*/
/*-----------------------------------------Rusanov flux computation----------------------------------------*/
/*=========================================================================================================*/

// --- Surface flux computation
// --- Density flux
surfaceScalarField f_Rho = 0.5*( (contrrhoU_own + contrrhoU_nei) - BetaSf * (rho_nei - rho_own) );

surfaceScalarField f_RhoN2 = 0.5*( (contrrhoUN2_own + contrrhoUN2_nei) - BetaSf * (rho_neiN2 - rho_ownN2) );
surfaceScalarField f_RhoN = 0.5*( (contrrhoUN_own + contrrhoUN_nei) - BetaSf * (rho_neiN - rho_ownN) );
surfaceScalarField f_RhoNO = 0.5*( (contrrhoUNO_own + contrrhoUNO_nei) - BetaSf * (rho_neiNO - rho_ownNO) );
surfaceScalarField f_RhoO = 0.5*( (contrrhoUO_own + contrrhoUO_nei) - BetaSf * (rho_neiO - rho_ownO) );
surfaceScalarField f_RhoO2 = 0.5*( (contrrhoUO2_own + contrrhoUO2_nei) - BetaSf * (rho_neiO2 - rho_ownO2) );

// --- Momentum flux
// --- Temporary variables
surfaceVectorField fv_own = U_own*contrrhoU_own + p_own*normalVector;
surfaceVectorField fv_nei = U_nei*contrrhoU_nei + p_nei*normalVector;
surfaceVectorField f_RhoU = 0.5*( (fv_own + fv_nei) - BetaSf * (rhoU_nei - rhoU_own) );

// --- Energie flux
// --- Temporary variables
surfaceScalarField  fs_own = contrU_own * (rhoE_own + p_own);
surfaceScalarField  fs_nei = contrU_nei * (rhoE_nei + p_nei);
surfaceScalarField f_RhoE = 0.5*(fs_own + fs_nei - BetaSf * (rhoE_nei - rhoE_own) );

// --- Set time parameters :: amaxSf is useful to calculate the Courant number
surfaceScalarField amaxSf("amaxSf", BetaSf * mesh.magSf() );

// --- Set jaconbian values
surfaceScalarField jacobianRho("jacobianRho", (contrU_own*(rho_own + hrho_own) - contrrhoU_own ) * mesh.magSf() / hrho_own );

surfaceScalarField contrrhoU_jac = ( (rhoU_own & normalVector) + hrhoU_own );
surfaceVectorField fv_own_jac = U_own*contrrhoU_jac + p_own*normalVector;

surfaceScalarField jacobianRhoU("jacobianRhoU", ( (fv_own_jac - fv_own) & mesh.Sf())/ hrhoU_own );

surfaceScalarField fs_own_jac = contrU_own * ( (rhoE_own + hrhoE_own) + p_own);

surfaceScalarField jacobianRhoE("jacobianRhoE", (fs_own_jac - fs_own) * mesh.magSf() / hrhoE_own);



// --- Set jaconbian values explicit part
surfaceScalarField explicitJacobianRho("explicitJacobianRho", (contrU_own*(rho_own + hrho_own) - contrrhoU_own ) * rho_own * mesh.magSf() / hrho_own );

contrrhoU_jac = ( (rhoU_own & normalVector) + hrhoU_own );
fv_own_jac = U_own*contrrhoU_jac + p_own*normalVector;
surfaceVectorField explicitJacobianRhoU("explicitJacobianRhoU", ( (fv_own_jac - fv_own) & mesh.Sf())*rhoU_own/ hrhoU_own );


fs_own_jac = contrU_own * ( (rhoE_own + hrhoE_own) + p_own);
surfaceScalarField explicitJacobianRhoE("explicitJacobianRhoE", (fs_own_jac - fs_own) * mesh.magSf() * rhoE_own / hrhoE_own);

// --- Set rho species jacobians
surfaceScalarField jacobianRhoN2("jacobianRhoN2", (contrU_own*(rho_ownN2 + hrho_own) - contrrhoUN2_own ) * mesh.magSf() / hrho_own );
surfaceScalarField jacobianRhoN("jacobianRhoN", (contrU_own*(rho_ownN + hrho_own) - contrrhoUN_own ) * mesh.magSf() / hrho_own );
surfaceScalarField jacobianRhoNO("jacobianRhoNO", (contrU_own*(rho_ownNO + hrho_own) - contrrhoUNO_own ) * mesh.magSf() / hrho_own );
surfaceScalarField jacobianRhoO("jacobianRhoO", (contrU_own*(rho_ownO + hrho_own) - contrrhoUO_own ) * mesh.magSf() / hrho_own );
surfaceScalarField jacobianRhoO2("jacobianRhoO2", (contrU_own*(rho_ownO2 + hrho_own) - contrrhoUO2_own ) * mesh.magSf() / hrho_own );
