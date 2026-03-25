/**************************************************/
surfaceScalarField rho_ownN2 =
  fvc::interpolate(rhoN2, pos, "reconstruct(rhoN2)");
surfaceScalarField rho_neiN2 =
  fvc::interpolate(rhoN2, neg, "reconstruct(rhoN2)");

surfaceScalarField rho_ownO2 =
  fvc::interpolate(rhoO2, pos, "reconstruct(rhoO2)");
surfaceScalarField rho_neiO2 =
  fvc::interpolate(rhoO2, neg, "reconstruct(rhoO2)");

surfaceScalarField rho_ownN =
  fvc::interpolate(rhoN, pos, "reconstruct(rhoN)");
surfaceScalarField rho_neiN =
  fvc::interpolate(rhoN, neg, "reconstruct(rhoN)");

surfaceScalarField rho_ownNO =
  fvc::interpolate(rhoNO, pos, "reconstruct(rhoNO)");
surfaceScalarField rho_neiNO =
  fvc::interpolate(rhoNO, neg, "reconstruct(rhoNO)");

surfaceScalarField rho_ownO =
  fvc::interpolate(rhoO, pos, "reconstruct(rhoO)");
surfaceScalarField rho_neiO =
  fvc::interpolate(rhoO, neg, "reconstruct(rhoO)");
/**************************************************/

surfaceVectorField rhoUN2_own =
  fvc::interpolate(rhoUN2, pos, "reconstruct(U)");
surfaceVectorField rhoUN2_nei =
  fvc::interpolate(rhoUN2, neg, "reconstruct(U)");

surfaceVectorField rhoUN_own =
  fvc::interpolate(rhoUN, pos, "reconstruct(U)");
surfaceVectorField rhoUN_nei =
  fvc::interpolate(rhoUN, neg, "reconstruct(U)");

surfaceVectorField rhoUNO_own =
  fvc::interpolate(rhoUNO, pos, "reconstruct(U)");
surfaceVectorField rhoUNO_nei =
  fvc::interpolate(rhoUNO, neg, "reconstruct(U)");

surfaceVectorField rhoUO2_own =
  fvc::interpolate(rhoUO2, pos, "reconstruct(U)");
surfaceVectorField rhoUO2_nei =
  fvc::interpolate(rhoUO2, neg, "reconstruct(U)");

surfaceVectorField rhoUO_own =
  fvc::interpolate(rhoUO, pos, "reconstruct(U)");
surfaceVectorField rhoUO_nei =
  fvc::interpolate(rhoUO, neg, "reconstruct(U)");

// --- Update values

if (sim.order(runTime) == 1 )
  {
    forAll (owner, faceI)
      {
	const label own = owner[faceI];
	const label nei = neighbour[faceI];

	rho_ownN2[faceI] = rhoN2[own];
	rho_neiN2[faceI] = rhoN2[nei];

  rho_ownN[faceI] = rhoN[own];
  rho_neiN[faceI] = rhoN[nei];

  rho_ownNO[faceI] = rhoNO[own];
  rho_neiNO[faceI] = rhoNO[nei];

  rho_ownO[faceI] = rhoO[own];
  rho_neiO[faceI] = rhoO[nei];

  rho_ownO2[faceI] = rhoO2[own];
  rho_neiO2[faceI] = rhoO2[nei];

	rhoUN2_own[faceI] = rhoUN2[own];
	rhoUN2_nei[faceI] = rhoUN2[nei];

  rhoUN_own[faceI] = rhoUN[own];
  rhoUN_nei[faceI] = rhoUN[nei];

  rhoUNO_own[faceI] = rhoUNO[own];
  rhoUNO_nei[faceI] = rhoUNO[nei];

  rhoUO2_own[faceI] = rhoUO2[own];
  rhoUO2_nei[faceI] = rhoUO2[nei];

  rhoUO_own[faceI] = rhoUO[own];
  rhoUO_nei[faceI] = rhoUO[nei];
      };
  }


// --- rhoU
const surfaceScalarField contrrhoUN2_own = (rhoUN2_own & normalVector);
const surfaceScalarField contrrhoUN2_nei = (rhoUN2_nei & normalVector);

const surfaceScalarField contrrhoUN_own = (rhoUN_own & normalVector);
const surfaceScalarField contrrhoUN_nei = (rhoUN_nei & normalVector);

const surfaceScalarField contrrhoUNO_own = (rhoUNO_own & normalVector);
const surfaceScalarField contrrhoUNO_nei = (rhoUNO_nei & normalVector);

const surfaceScalarField contrrhoUO_own = (rhoUO_own & normalVector);
const surfaceScalarField contrrhoUO_nei = (rhoUO_nei & normalVector);

const surfaceScalarField contrrhoUO2_own = (rhoUO2_own & normalVector);
const surfaceScalarField contrrhoUO2_nei = (rhoUO2_nei & normalVector);
