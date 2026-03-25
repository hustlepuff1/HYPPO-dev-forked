/*=========================================================================================================*/
/*------------------------------------------Functions definitions------------------------------------------*/
/*=========================================================================================================*/

#include "turbulence/turbulence.H"

// --- Constructor
Turbulence::Turbulence() : matchTol(1e-10), kappa(0.4), K(0.0168), Ccp(1.6)
{
	cellFaceDist = 0;
	yTemp        = 0;
	diffDist     = 0;
	ymax         = 0;

	gatheredStuff.resize(Pstream::nProcs());
};


// --- Compute the magnitude of the vorticity
volScalarField Turbulence::magVorticity(volVectorField& field) const
{
	return mag(fvc::curl(field));
};


// --- Compute enthalpy far from wall
void Turbulence::farEnthalpy(fvMesh& mesh, volScalarField& H)
{
    forAll(wallDist(mesh).y(),cellI)
    {
      if(wallDist(mesh).y()[cellI] == Ymax)
      {
        cellInfinity = cellI;
      }
    }
	Htot = H[cellInfinity];
};


// --- Computes F (cf. Baldwin-Lomax model)
volScalarField Turbulence::F(volScalarField& magVorticity, volScalarField& yPlus, fvMesh& mesh, volScalarField Aplus) const
{
	return magVorticity * wallDist(mesh).y() * (1 - exp(-yPlus / Aplus));
};


// --- Computes uTau boundary field
void Turbulence::uTauBoundary(volScalarField& uTau, volScalarField& mu, volVectorField& U, volScalarField& rho, const fvPatchList& patches)
{
	forAll(patches, patchi)
	{
		const fvPatch& currPatch = patches[patchi];

		if (isType<wallFvPatch>(currPatch))
		{
		  uTau.boundaryField()[patchi] =sqrt(mu.boundaryField()[patchi]*mag(U.boundaryField()[patchi].snGrad()) /rho.boundaryField()[patchi]);
		}
	}
};


// --- Set boundary wall
void Turbulence::boundaryWall(fvMesh& mesh, volScalarField& uTau, const surfaceVectorField& Cf, const fvPatchList& patches,volScalarField& rhoW, volScalarField& muW)
{
	forAll(patches,patchi)
	{
	  const fvPatch& currPatch = patches[patchi];
	  nFaces_ = mesh.boundaryMesh()[patchi].size();

	  if(isType<wallFvPatch>(currPatch))
	  {
		nFaces2=nFaces_;
		for(int facei = 0; facei<nFaces_; facei++)
		{
		  xyz.append((Cf.boundaryField()[patchi][facei]));
		  uTauPatch.append(uTau.boundaryField()[patchi][facei]);
		  rhoWPatch.append(rhoW.boundaryField()[patchi][facei]);
		  muWPatch.append(muW.boundaryField()[patchi][facei]);
		}
	  }
	}
};


// --- Computes the distance from the wall
void Turbulence::wallDistance(fvMesh& mesh, volScalarField& uTau, const volVectorField& C,volScalarField& rhoW, volScalarField& muW)
{
	forAll(mesh.C(),cellI)
	{
		for(int facei = 0; facei<nFaces2; facei++)
		{
			cellFaceDist = Foam::sqrt(sqr(C[cellI].x()-xyz[facei].x()) + sqr(C[cellI].y()-xyz[facei].y()) + sqr(C[cellI].z()-xyz[facei].z()) );
			yTemp = wallDist(mesh).y()[cellI];
			diffDist = mag(cellFaceDist - yTemp);

			if( diffDist <= matchTol)
			{
				uTau[cellI] = uTauPatch[facei];
				rhoW[cellI] = rhoWPatch[facei];
				muW[cellI] = muWPatch[facei];
			}
		}
	}
};


// --- Computes yPlus
volScalarField Turbulence::yPlus(volScalarField& uTau, fvMesh& mesh,volScalarField& rhoW, volScalarField& muW) const
{
	return wallDist(mesh).y() * uTau * rhoW / muW;
};


// --- Computes muT using Baldwin-Lomax model
void Turbulence::baldwinLomax(fvMesh& mesh, volScalarField& yPlus, volScalarField& muTinner, volScalarField& muTouter, volScalarField& Ckleb, scalar shockDistance, volScalarField& muT, volScalarField& rho, volScalarField& Aplus,  volScalarField& H, const surfaceVectorField& Cf, const volVectorField& C, const fvPatchList& patches, volScalarField& magVorticity)
{
	forAll(patches,patchi)
	{
		const fvPatch& currPatch = patches[patchi];
		nFaces_ = mesh.boundaryMesh()[patchi].size();

		if(isType<wallFvPatch>(currPatch))
		{
			for(int facei = 0; facei<nFaces_; facei++)
			{
				listOfCells.clear();
				Fvalues.clear();
				ycrossoverList.clear();

				forAll(wallDist(mesh).y(),cellI)
				{
					cellFaceDistScalar = Foam::sqrt(sqr(C[cellI].x()-Cf.boundaryField()[patchi][facei].x()) + sqr(C[cellI].y()-Cf.boundaryField()[patchi][facei].y())+sqr(C[cellI].z()-Cf.boundaryField()[patchi][facei].z()));
					yTemp = wallDist(mesh).y()[cellI];
					diffDist = mag(cellFaceDistScalar - yTemp);

					if( diffDist < matchTol)
					{
						listOfCells.append(cellI);
					}
				}

				for(int k=0; k<NumberOfCells;k++)
				{
					cellID = listOfCells[k];
					HCell = H[cellID];

					if((mag(HCell - Htot) <= 0.05*Htot)|| wallDist(mesh).y()[cellID] > shockDistance)
					{
					  Fvalues.append(0);
					}
					else
					{
					  Fvalues.append(magVorticity[cellID] * wallDist(mesh).y()[cellID] * (1 - std::exp(-yPlus[cellID] / Aplus[cellID])));
					}
				}

				Fmax = max(Fvalues);

				for(int j = 0; j<NumberOfCells; j++)
				{
					cellID = listOfCells[j];
					if( Fvalues[j] == Fmax)
					{
					  ymax = wallDist(mesh).y()[cellID];
					}
				}
				
				Fwake = ymax * Fmax ;

				for(int j = 0; j<NumberOfCells; j++)
				{
					cellID = listOfCells[j];
					l = wallDist(mesh).y()[cellID] * kappa * (1 - std::exp(-yPlus[cellID] / Aplus[cellID]));
					muTinner[cellID] = rho[cellID] * l * l *magVorticity[cellID];
					muTouter[cellID] = K* rho[cellID] * Ccp * Fwake * (1/ (1+5.5*pow(Ckleb[cellID] * wallDist(mesh).y()[cellID]/ymax,6)));
				}
				for(int k=0; k<NumberOfCells;k++)
				{
					cellID = listOfCells[k];
					if(muTinner[cellID] > muTouter[cellID])
					{
					  ycrossoverList.append(wallDist(mesh).y()[cellID]);
					}
				}
				
				ycrossover = min(ycrossoverList);

				for(int j=0; j<NumberOfCells;j++)
				{
					cellID = listOfCells[j];
					if(wallDist(mesh).y()[cellID] > ycrossover)
					{
					  muT[cellID] = muTouter[cellID];
					}
					else
					{
					  muT[cellID] = muTinner[cellID];
					}
				}
			}
		}
	}
};


// --- Deal with a parralel computation
void Turbulence::parControl()
{
	for (int i = 0 ; i<Pstream::nProcs() ; i++)
	{
	  gatheredStuff[i] = xyz;
	}

	xyz=gatheredStuff[0];
	Pstream::gatherList(gatheredStuff);
	Pstream::scatterList(gatheredStuff);
};
