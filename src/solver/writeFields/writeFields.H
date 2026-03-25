/*==============================================================================================================*/
/*----------------------------------------------Write fields class----------------------------------------------*/
/*==============================================================================================================*/

// --- Class to manage the wrting of the fields
class Write
{
	public:
		// --- Print values corresponding to the negative energy
		void msgNegativeEnergy (int cellI, volScalarField& rho, volScalarField& rhoE, volScalarField& e, volVectorField& U, volScalarField& T, volScalarField& psi);

		// --- Write fields for sutherland model
		void writeFields (volScalarField& mu);

		// --- Write fields for table reading air 5 model or finite rate
		void writeFields (volScalarField& mu, volScalarField& O2, volScalarField& N2, volScalarField& NO, volScalarField& N, volScalarField& O);

		// --- Write fields for table reading air 11 model 
		void writeFields (volScalarField& mu, volScalarField& O2, volScalarField& N2, volScalarField& NO, volScalarField& N, volScalarField& O, volScalarField& O2Plus, volScalarField& N2Plus, volScalarField& NOPlus, volScalarField& NPlus, volScalarField& OPlus, volScalarField& eMoins);

		// --- Write fields for turbulence model
		void writeFields (volScalarField& muT, volScalarField& yPlus);
};
