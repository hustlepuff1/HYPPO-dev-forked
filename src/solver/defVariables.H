//Initialization
// Get face-to-cell addressing: face area point from owner to neighbour

// --- Default définitions :: to get the good dimensions
surfaceScalarField rho_own =
  fvc::interpolate(rho, pos, "reconstruct(rho)");
surfaceScalarField rho_nei =
  fvc::interpolate(rho, neg, "reconstruct(rho)");

/**************************************************/

surfaceVectorField rhoU_own =
  fvc::interpolate(rhoU, pos, "reconstruct(U)");
surfaceVectorField rhoU_nei =
  fvc::interpolate(rhoU, neg, "reconstruct(U)");

surfaceTensorField tauMC_own =
  fvc::interpolate(tauMC,pos, "reconstruct(tauMC)");
surfaceTensorField tauMC_nei =
  fvc::interpolate(tauMC,neg, "reconstruct(tauMC)");

surfaceVectorField tauVq_own =
  fvc::interpolate(tauVq,pos, "reconstruct(tauVq)");
surfaceVectorField tauVq_nei =
  fvc::interpolate(tauVq,neg, "reconstruct(tauVq)");

surfaceScalarField rhoE_own =
  fvc::interpolate(rhoE, pos, "reconstruct(U)");
surfaceScalarField rhoE_nei =
  fvc::interpolate(rhoE, neg, "reconstruct(U)");

surfaceScalarField E_own =
  fvc::interpolate(E, pos, "reconstruct(U)");
surfaceScalarField E_nei =
  fvc::interpolate(E, neg, "reconstruct(U)");

surfaceScalarField rPsi_own =
  fvc::interpolate(rPsi, pos, "reconstruct(T)");
surfaceScalarField rPsi_nei =
  fvc::interpolate(rPsi, neg, "reconstruct(T)");

// --- Total energy
surfaceScalarField e_own =
  fvc::interpolate(e, pos, "reconstruct(T)");
surfaceScalarField e_nei =
  fvc::interpolate(e, neg, "reconstruct(T)");

// --- To build jacobian
surfaceScalarField hrho_own =
  fvc::interpolate(hrho, pos, "reconstruct(hro)");
surfaceScalarField hrhoU_own =
  fvc::interpolate(hrhoU, pos, "reconstruct(hrhoU)");
surfaceScalarField hrhoE_own =
  fvc::interpolate(hrhoE, pos, "reconstruct(hrhoE)");
surfaceTensorField htauMC_own =
  fvc::interpolate(htauMC, pos, "reconstruct(htauMC)");
surfaceVectorField htauVq_own =
  fvc::interpolate(htauVq, pos, "reconstruct(htauVq)");

//--- Update speed of sound values to calculate the Jacobian's eigen values

if (model != "readTable") c = sqrt(rPsi* Cp/Cv);

Mach=mag(U)/c;

// --- Speed of sound's directed interpolation onto faces
surfaceScalarField cSf_own = fvc::interpolate(c, pos, "reconstruct c");
surfaceScalarField cSf_nei = fvc::interpolate(c, neg, "reconstruct c");
// --- Update values

if (sim.order(runTime) == 1 )
  {
    forAll (owner, faceI)
      {
	const label own = owner[faceI];
	const label nei = neighbour[faceI];

	rho_own[faceI] = rho[own];
	rho_nei[faceI] = rho[nei];

	rhoU_own[faceI] = rhoU[own];
	rhoU_nei[faceI] = rhoU[nei];

	rhoE_own[faceI] = rhoE[own];
	rhoE_nei[faceI] = rhoE[nei];

	E_own[faceI] = E[own];
	E_nei[faceI] = E[nei];

	rPsi_own[faceI] = rPsi[own];
	rPsi_nei[faceI] = rPsi[nei];

	e_own[faceI] = e[own];
	e_nei[faceI] = e[nei];

	cSf_own[faceI] = c[own];
	cSf_nei[faceI] = c[nei];

	hrho_own[faceI] = hrho[own];
	hrhoU_own[faceI] = hrhoU[own];
	hrhoE_own[faceI] = hrhoE[own];
      };
  }
surfaceVectorField U_own = rhoU_own/rho_own;
surfaceVectorField U_nei = rhoU_nei/rho_nei;
surfaceScalarField p_own = rho_own*rPsi_own;
surfaceScalarField p_nei = rho_nei*rPsi_nei;


// --- Contravariant vectors :: inner product between two vectors <==> dot_product
// --- Enables to project the surfaceVectorFields on the faces in order to build RusanovFlux

// --- rhoU
const surfaceScalarField contrrhoU_own = (rhoU_own & normalVector);
const surfaceScalarField contrrhoU_nei = (rhoU_nei & normalVector);
// --- U
const surfaceScalarField contrU_own = (U_own & normalVector);
const surfaceScalarField contrU_nei = (U_nei & normalVector);
//tau
const surfaceVectorField contrTau_own = (tauMC_own & normalVector);
const surfaceVectorField contrTau_nei = (tauMC_nei & normalVector);
//tauVq
const surfaceScalarField contrTauVq_own = (tauVq_own & normalVector);
const surfaceScalarField contrTauVq_nei = (tauVq_nei & normalVector);
//htauMC
const surfaceVectorField contrTau_own_jac = (htauMC_own & normalVector);
//htauVq
const surfaceScalarField contrTauVq_own_jac = (htauVq_own & normalVector);

// --- Jacobian's eigen values
// --- max eigen value in each cell
surfaceScalarField ap = contrU_own + cSf_own;
surfaceScalarField am = contrU_nei + cSf_nei;


// --- The biggest of the eigen values between 2 cells
surfaceScalarField BetaSf("BetaSf", max(mag(am), mag(ap)));
