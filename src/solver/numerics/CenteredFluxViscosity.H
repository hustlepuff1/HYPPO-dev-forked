/*=========================================================================================================*/
/*----------------------------------------------Viscous flux-----------------------------------------------*/
/*=========================================================================================================*/

// --- Surface flux computation
// --- Density flux
//surfaceScalarField fv1 = 0.0; this flux does not exist.
surfaceVectorField fv2_own = contrTau_own;
surfaceVectorField fv2_nei = contrTau_nei;
surfaceVectorField f_RhoU_visc = 0.5*((fv2_own + fv2_nei));

surfaceScalarField fv3_own = contrTauVq_own;
surfaceScalarField fv3_nei = contrTauVq_nei;
surfaceScalarField f_RhoE_visc = 0.5*((fv3_nei + fv3_own));

// --- Set jacobian values

surfaceVectorField fv2_own_jac = contrTau_own + contrTau_own_jac;
surfaceScalarField jacobianRhoUv("jacobianRhoUv", ( (fv2_own_jac - fv2_own) & mesh.Sf())/ hrhoU_own );

surfaceScalarField fv3_own_jac = contrTauVq_own + contrTauVq_own_jac;
surfaceScalarField jacobianRhoEv("jacobianRhoEv", ( (fv3_own_jac - fv3_own) * mesh.magSf())/ hrhoE_own);
