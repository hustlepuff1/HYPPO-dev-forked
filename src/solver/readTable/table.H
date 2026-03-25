/*=========================================================================================================*/
/*-------------------------------------------Reading table class-------------------------------------------*/
/*=========================================================================================================*/

#include <list>

class ReadTable
{
	private:
		// --- Possible readable lists
		DynamicList<scalar> rhoList,
							energyList,
							rhoeList,
							temperatureList,
							pressureList,
							soundList,
							enthalpyList,
							CpList,
							muList,
							lambdaList,
							CvList,
							NList,
							NOList,
							OList,
							N2List,
							O2List,
							eMoinsList,
							NPlusList,
							OPlusList,
							N2PlusList,
							O2PlusList,
							NOPlusList;

		// --- Number of rho values
		int nRho;

		// --- Size of the lists
		int sizeTable_;

		// --- Values used to bound
		scalar maxRho_; 
		scalar min_e;
		scalar max_e;

		// --- Coefficients to interpolate
		scalar rhoCoef;
		scalar eCoef;	
	
		// --- Index used to parse lists
		int min_eIndex;
		int max_eIndex;
		int rhoIndex;

		// --- Open variables for table
		std::fstream air5;
		std::fstream air11;

		// --- Used to stop when correct value found in interpolation
		bool stop_;

		// --- Air 11 model
		bool air11_;

		// --- Variables used to compute boundary corrections
		label face;
		label own;

	public:
		// --- If air11 is used
		bool airModel11 () const;
		
		// --- Open file
		bool openFile (word model);
		
		// --- Table Reading
		void reading ();

		// --- Counting of values of rho in the table
		void countRho ();

		// --- Size of the lists
		void sizeTable ();

		// --- Maximum of rho, used to clamp
		void maxRho ();

		// --- Setting the value stop_
		void stop (bool value);

		// --- Coefficient of the interpolation for rho
		void interRho (scalar value);

		// --- Coefficient of the interpolation for the normalized energy
		void inter_e(scalar value);

		// --- Correct boundary conditions of fields
		void boundaryCorrection(volScalarField& field, fvMesh& mesh, const polyBoundaryMesh& boundaryMesh);

		// --- Selection of the list to use for interpolation
		DynamicList<scalar>& fieldToList (volScalarField& field) ;

		scalar interpolation (volScalarField& field);

		// --- Interpolation table function for air 5
		void interpolationField (
			int cellI,
			volScalarField& T, 
			volScalarField& p ,
			volScalarField& mu ,
			volScalarField& Cv ,
			volScalarField& Cp ,
			volScalarField& lambda,
			volScalarField& N2,
			volScalarField& N,
			volScalarField& NO,
			volScalarField& O2,
			volScalarField& O,
			volScalarField& c);
		
		// --- Interpolation table function for air 11
		void interpolationField (
			int cellI,
			volScalarField& T, 
			volScalarField& p ,
			volScalarField& mu ,
			volScalarField& Cv ,
			volScalarField& Cp ,
			volScalarField& lambda,
			volScalarField& N2,
			volScalarField& N,
			volScalarField& NO,
			volScalarField& O2,
			volScalarField& O,
			volScalarField& c,
			volScalarField& N2Plus,
			volScalarField& NPlus,
			volScalarField& NOPlus,
			volScalarField& O2Plus,
			volScalarField& OPlus,
			volScalarField& eMoins);

		// --- Interpolation table function for air 5
		void boundaryCorrectionField (
			fvMesh& mesh, 
			const polyBoundaryMesh& boundaryMesh,
			volScalarField& T, 
			volScalarField& p ,
			volScalarField& mu ,
			volScalarField& Cv ,
			volScalarField& Cp ,
			volScalarField& lambda,
			volScalarField& N2,
			volScalarField& N,
			volScalarField& NO,
			volScalarField& O2,
			volScalarField& O,
			volScalarField& c);
		
		// --- Interpolation table function for air 11
		void boundaryCorrectionField (
			fvMesh& mesh, 
			const polyBoundaryMesh& boundaryMesh,
			volScalarField& T, 
			volScalarField& p ,
			volScalarField& mu ,
			volScalarField& Cv ,
			volScalarField& Cp ,
			volScalarField& lambda,
			volScalarField& N2,
			volScalarField& N,
			volScalarField& NO,
			volScalarField& O2,
			volScalarField& O,
			volScalarField& c,
			volScalarField& N2Plus,
			volScalarField& NPlus,
			volScalarField& NOPlus,
			volScalarField& O2Plus,
			volScalarField& OPlus,
			volScalarField& eMoins);		

};
