/*==============================================================================================================*/
/*-----------------------------------------------Turbulence class-----------------------------------------------*/
/*==============================================================================================================*/

class Turbulence
{
	private:
		// --- Surface face center coordinates
		DynamicList<vector> xyz;

		// --- Values on the wall
		DynamicList<scalar> uTauPatch;
		DynamicList<scalar> rhoWPatch;
		DynamicList<scalar> muWPatch;

		// --- List of cells on the wall
		DynamicList<scalar> listOfCells;

		// --- 
		DynamicList<scalar> Fvalues;

		// --- Y crossover list
		DynamicList<scalar> ycrossoverList;

		// --- Constants used in Baldwin-Lomax model
		const scalar kappa;
		const scalar K;
		const scalar Ccp;

		// --- Tolerance for diffDist
		const scalar matchTol;

		// --- Enthalpy far from the shock
		scalar Htot;

		// --- Distance from the wall
		scalar cellFaceDist;

		// --- Stores the distance to the wall
		scalar yTemp;

		// --- Difference between Foam distance and cartesian distance
		scalar diffDist;

		// --- Max value of y
		scalar ymax;

		// --- Distance from the wall
		scalar cellFaceDistScalar;

		// --- y corresponding to Fmax
		scalar Ymax;

		// --- Current cell's enthalpy
		scalar HCell;

		// --- Max value of F
		scalar Fmax;

		// --- ymax * Fmax
		scalar Fwake;

		// --- y[cellID] * kappa * (1 - std::exp(-yPlus[cellID] / Aplus[cellID]))
		scalar l;

		// --- y corresponding to muTinner=muTouter
		scalar ycrossover;

		// --- Number of faces on the wall
		int nFaces2;

		// --- Cell far from the wall
		int cellInfinity;

		// --- Number of cells
		int NumberOfCells;

		// --- Current cell
		int cellID;

		// --- Number of faces on boundary patches
		label nFaces_;

		// --- Gathered values to parralel computation
        List<List<vector>> gatheredStuff;

	public:
		// --- Constructor
		Turbulence();
	
		// --- Computes enthalpy far from the wall
		void farEnthalpy(fvMesh& mesh, volScalarField& H);

		// --- Computes uTau boundary field
		void uTauBoundary(volScalarField& uTau, volScalarField& mu, volVectorField& U, volScalarField& rho, const fvPatchList& patches);

		// --- Set boundary wall
		void boundaryWall(fvMesh& mesh, volScalarField& uTau, const surfaceVectorField& Cf, const fvPatchList& patches,volScalarField& rhoW, volScalarField& muW);

		// --- Computes the distance from the wall
		void wallDistance(fvMesh& mesh, volScalarField& uTau, const volVectorField& C, volScalarField& rhoW, volScalarField& muW);

		// --- Computes muT using Baldwin-Lomax model
		void baldwinLomax(fvMesh& mesh, volScalarField& yPlus, volScalarField& muTinner, volScalarField& muTouter, volScalarField& CKleb, scalar shockDistance, volScalarField& muT, volScalarField& rho, volScalarField& Aplus, volScalarField& H, const surfaceVectorField& Cf, const volVectorField& C, const fvPatchList& patches, volScalarField& magVorticity);

		// --- Used to deal with a parralel computation
		void parControl();
	
		// --- Computes yPlus
		volScalarField yPlus(volScalarField& uTau, fvMesh& mesh,volScalarField& rhoW, volScalarField& muW) const;

		// --- Computes F (cf. Baldwin-Lomax model)
		volScalarField F(volScalarField& magVorticity, volScalarField& yPlus, fvMesh& mesh, volScalarField Aplus) const;

		// --- Computes the magnitude of the vorticity
		volScalarField magVorticity(volVectorField& field) const;	
};
