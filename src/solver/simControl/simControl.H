/*=========================================================================================================*/
/*-----------------------------------------Simulation control class----------------------------------------*/
/*=========================================================================================================*/

class SimControl 
{
	private:
		// --- Courant number
		scalar coNum_;

		// --- Mean Courant number
		scalar meanCoNum;

		// --- Variables to stock old Courant number
    	scalar maxCo_oldTime;
    	scalar maxCo_oldTime2;

		// --- Minimum Courant number defined
		scalar minCo_;

		// --- Maximum delta T 
		scalar maxDeltaTFact;

		// --- Delta T
		scalar deltaTFact;

		// --- Order of the scheme (0 or 1)
		scalar order_;

		// --- Convergence criterion
		scalar convergenceCriterion;

		// --- Enthalpy formation
		scalar hf;
	
		// --- Variables to manage the loopback
		bool stopMinCo_;
		bool negative_energy_;
		bool restart_;

		// --- Variables to manage iteration number
		int k;
		int k1;
		int k2;
			
	public:		
		// --- Maximum residual, public to avoid issue with foam class
		scalar maxResidual;

		// --- Order reading
		scalar order (Time& runTime) ;

		// --- Courant number computation
		void coNum (scalarField sumPhi, Time& runTime, fvMesh& mesh);

		// --- Methods to manage iteration number
		void switchK2K1 ();
		void switchKK1 ();

		// --- Minimum Courant number reading
		void minCo (Time& runTime);

		// --- Setting of delta T
		void setDeltaT (bool adjustTimeStep, scalar maxCo, scalar maxDeltaT, Time& runTime);

		// --- Setting energy boolean
		void negative_energy (bool value);

		// --- Setting restart boolean
		void restart (bool value);

		// --- Setting stopMinCo boolean
		void stopMinCo (bool value);

		// --- Incrementation of K2
		void incremK2 ();

		// --- Method used to manage loopback in parallel
		void parControl ();

		// --- Checking convergence
		void convergenceCheck (Time& runTime);

		// --- Set the value of hf
		void setHf(scalar value);

		// --- test wrong energy
		bool testEnergy(scalar e, word model);

		// --- Give the value of restart
		bool restart () const;



};
