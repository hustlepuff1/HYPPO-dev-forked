/*=========================================================================================================*/
/*------------------------------------------Functions definitions------------------------------------------*/
/*=========================================================================================================*/

#include "simControl/simControl.H"

// --- Computes and display the max and the mean Courant number
void SimControl::coNum(scalarField sumPhi, Time& runTime, fvMesh& mesh) 
{
	if (mesh.nInternalFaces())
	{      

		//Courant number :: CFL >= (lambda * Sum(CellInterfacesSurfaceAreas) ) / (CellVolume)) * Dt 
		//with lambda==max(max(Jacobian's eigenvalues))

		coNum_ = gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();//*0.5;

		//Mean courant number
		meanCoNum =(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();//*0.5;
	}

	Info<< "Courant Number mean: " << meanCoNum << " max: " << coNum_  << nl << endl;
};

// --- Switch K2 and K1 values
void SimControl::switchK2K1()
{
	if(k2>0)
	{
		k = k1;
	}
};

// --- Switch K and K1 values
void SimControl::switchKK1()
{
	if(k<=9)
	{
		k = k+1;
	}
	else
	{
		k=0;
	}
	k1 = k;
};

// --- Minimum Courant number computation
void SimControl::minCo(Time& runTime)
{
	minCo_ = runTime.controlDict().lookupOrDefault<scalar>("minCo", 0.02);
};

// --- Computation and regulation of delta T
void SimControl::setDeltaT(bool adjustTimeStep, scalar maxCo, scalar maxDeltaT, Time& runTime)
{
	if (adjustTimeStep )
	{
		if(k2==0 && k==0)
		{
			maxCo_oldTime = coNum_;
		}
		else if (maxCo_oldTime2 < maxCo)
		{
			maxCo_oldTime = maxCo_oldTime2;
		}
		else
		{
			maxCo_oldTime = maxCo;
		}
	
	// --- Loopback
	if (negative_energy_)
	{
		if (maxCo_oldTime>minCo_/0.8)
		{
			maxCo_oldTime = maxCo_oldTime*0.8;
		}
		// --- If coNum is not equal to minCo
		else if (!stopMinCo_)
		{
			maxCo_oldTime = minCo_;
			stopMinCo_=true;
		}
		else
		{
			FatalErrorIn
			(
			"setDeltaT"
			)
			<< "Lower Courant number needed"
			<< abort(FatalError);
		}
		k = 0;
		maxCo_oldTime2 = maxCo_oldTime;
	}
	else
	{
		if(k>9 && k2>100)
		{
			maxCo_oldTime = maxCo_oldTime/0.8;
			k = 0;
		}

		maxCo_oldTime2 = maxCo_oldTime;
	}

	maxDeltaTFact = maxCo_oldTime/(coNum_);
	deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
	runTime.setDeltaT(min(deltaTFact*runTime.deltaT().value(),maxDeltaT));

	Info<< "deltaT = " <<  runTime.deltaT().value() << nl << endl;
	}
};

// --- Setting energy boolean
void SimControl::negative_energy(bool value)
{
	negative_energy_=value;
};

// --- Setting stopMinCo boolean
void SimControl::stopMinCo(bool value)
{
	stopMinCo_=value;
};

// --- Setting restart boolean
void SimControl::restart(bool value)
{
	restart_=value;
};

// --- Give restart value
bool SimControl::restart() const
{
	return restart_;
};

// --- Incrementation of K2
void SimControl::incremK2()
{
	k2 = k2 + 1;
	Info << "Time iteration number ::  " << k2 << nl << endl;
};

// --- Manage loopback in parallel
void SimControl::parControl()
{
	reduce(restart_ , sumOp<bool>());
	reduce(negative_energy_ , sumOp<bool>());
};

// --- Read the order of the simulation
scalar SimControl::order(Time& runTime)
{
    return order_ = runTime.controlDict().lookupOrDefault<scalar>("order", 1.0);
};

// --- Checking convergence
void SimControl::convergenceCheck(Time& runTime)
{
	if (maxResidual < convergenceCriterion)
	{
		Info<< "reached convergence criterion: " << convergenceCriterion << endl;
		runTime.writeAndEnd();
		Info<< "latestTime = " << runTime.timeName() << endl;
	}
};

// --- Test wrong energy
bool SimControl::testEnergy(scalar e, word model)
{
	if (model == "none" || model == "sutherland")
	{
		if (e < 0) return true;
	}
	else if (model == "readTable" || model == "finiteRate")
	{
		if (e < -300000) return true;		// --- Enthalpy formation from mutation++
	}
	return false;
};

void SimControl::setHf(scalar value)
{
	hf = value;
};





