/*=========================================================================================================*/
/*------------------------------------------Functions definitions------------------------------------------*/
/*=========================================================================================================*/

#include "readTable/table.H"

// --- Open File
bool ReadTable::openFile(word model)
{
	air5  = std::fstream("constant/air5.txt");

	if (air11 = std::fstream("constant/air11.txt")) air11_=true;

	if(air5 && air11) //Cannot store files in the folder constant at the same time
	{
		FatalErrorIn
		(
		    "openFile"
		)
		<< "You can't use air5 and air11 at the same time"
		<< abort(FatalError);
	}

	if(model=="readTable" && air5 || air11 ) //If the table is found
	{
		Info << "Thermophysical properties from table" << nl << endl;
		return true;
	}
	else //If cannot find the table
	{
	  	Info << "Constant thermophysical properties used" << nl << endl;
		return false;
	}
};


// --- Bool to test if air 11 is used
bool ReadTable::airModel11() const
{
	return air11_;
};


// --- Setting the variable stop
void ReadTable::stop(bool value)
{
	stop_ = value;
};


// --- Table reading
void ReadTable::reading ()
{
	std::string ligne;
	scalar value;

	if (air5) 
	{
		while(getline(air5,ligne))
		{
			air5 >> value;
			if(air5.eof()) break;
			rhoList.append(value);
			air5 >> value;
			rhoeList.append(value);
			air5 >> value;
			energyList.append(value);
			air5 >> value;
			soundList.append(value);
			air5 >> value;
			temperatureList.append(value);
			air5 >> value;
			pressureList.append(value);
			air5 >> value;
			CpList.append(value);
			air5 >> value;
			CvList.append(value);
			air5 >> value;
			muList.append(value);
			air5 >> value;
			lambdaList.append(value);
			air5 >> value;
			NList.append(value);
			air5 >> value;
			OList.append(value);
			air5 >> value;
			NOList.append(value);
			air5 >> value;
			N2List.append(value);
			air5 >> value;
			O2List.append(value);
		}
	}
	else if(air11)
	{
		while(getline(air11,ligne))
		{
			air11 >> value;  
			if(air11.eof()) break;
			rhoList.append(value);
			air11 >> value;
			rhoeList.append(value);
			air11 >> value;
			energyList.append(value);
			air11 >> value;
			soundList.append(value);
			air11 >> value;
			temperatureList.append(value);
			air11 >> value;
			pressureList.append(value);
			air11 >> value;
			CpList.append(value);
			air11 >> value;
			CvList.append(value);
			air11 >> value;
			muList.append(value);
			air11 >> value;
			lambdaList.append(value);
			air11 >> value;
			eMoinsList.append(value);
			air11 >> value;
			NList.append(value);
			air11 >> value;
			NPlusList.append(value);
			air11 >> value;
			OList.append(value);
			air11 >> value;
			OPlusList.append(value);
			air11 >> value;
			NOList.append(value);
			air11 >> value;
			N2List.append(value);
			air11 >> value;
			N2PlusList.append(value);
			air11 >> value;
			O2List.append(value);
			air11 >> value;
			O2PlusList.append(value);
			air11 >> value;
			NOPlusList.append(value);
		}
	}
};


// --- Selection of the list to use for interpolation
DynamicList<scalar>& ReadTable::fieldToList (volScalarField& field) 
{
	if 		(field.name()=="mu") 		 return muList;
	else if (field.name()=="rho") 		 return rhoList;
	else if (field.name()=="rhoe") 		 return rhoeList;
	else if (field.name()=="energy")   	 return energyList;
	else if (field.name()=="enthalpy") 	 return enthalpyList;
	else if (field.name()=="T") 		 return temperatureList;
	else if (field.name()=="p") 		 return pressureList;
	else if (field.name()=="Cp") 		 return CpList;
	else if (field.name()=="Cv") 		 return CvList;
	else if (field.name()=="c") 		 return soundList;
	else if (field.name()=="lambda") 	 return lambdaList;
	else if (field.name()=="N2") 		 return N2List;
	else if (field.name()=="O2") 		 return O2List;
	else if (field.name()=="O") 		 return OList;
	else if (field.name()=="N") 		 return NList;
	else if (field.name()=="NO") 		 return NOList;
	else if (field.name()=="NOPlus") 	 return NOPlusList;
	else if (field.name()=="NPlus") 	 return NPlusList;
	else if (field.name()=="N2Plus") 	 return N2PlusList;
	else if (field.name()=="O2Plus") 	 return O2PlusList;
	else if (field.name()=="OPlus") 	 return OPlusList;
	else if (field.name()=="eMoins") 	 return eMoinsList;
};


// --- Interpolation table function for air 5
void ReadTable::interpolationField(
	int cellI,
	volScalarField& T, 
	volScalarField& p,
	volScalarField& mu,
	volScalarField& Cv,
	volScalarField& Cp,
	volScalarField& lambda,
	volScalarField& N2,
	volScalarField& N,
	volScalarField& NO,
	volScalarField& O2,
	volScalarField& O,
	volScalarField& c)
{
	mu[cellI]     = interpolation(mu);
	Cp[cellI]     = interpolation(Cp);
	lambda[cellI] = interpolation(lambda);
	Cv[cellI]     = interpolation(Cv);
	T[cellI]      = interpolation(T);
	p[cellI]      = interpolation(p);
	O[cellI]      = interpolation(O);
	NO[cellI]     = interpolation(NO);
	N[cellI]      = interpolation(N);
	N2[cellI]     = interpolation(N2);
	O2[cellI]     = interpolation(O2);	
	c[cellI]	  = interpolation(c);
};


// --- Interpolation table function for air 11
void ReadTable::interpolationField(
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
	volScalarField& eMoins)
{
	mu[cellI]     = interpolation(mu);
	Cp[cellI]     = interpolation(Cp);
	lambda[cellI] = interpolation(lambda);
	Cv[cellI]     = interpolation(Cv);
	T[cellI]      = interpolation(T);
	p[cellI]      = interpolation(p);
	O[cellI]      = interpolation(O);
	NO[cellI]     = interpolation(NO);
	N[cellI]      = interpolation(N);
	N2[cellI]     = interpolation(N2);
	O2[cellI]     = interpolation(O2);
	c[cellI]	  = interpolation(c);
	eMoins[cellI] = interpolation(eMoins);
	NOPlus[cellI] = interpolation(NOPlus);
	NPlus[cellI]  = interpolation(NPlus);
	O2Plus[cellI] = interpolation(O2Plus);
	OPlus[cellI]  = interpolation(OPlus);
	N2Plus[cellI] = interpolation(N2Plus);

};


// --- Boundary field update function for air 5
void ReadTable::boundaryCorrectionField (
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
	volScalarField& c)
{
	boundaryCorrection(T,mesh, boundaryMesh);
	boundaryCorrection(p,mesh, boundaryMesh);
	boundaryCorrection(mu,mesh, boundaryMesh);
	boundaryCorrection(Cv,mesh, boundaryMesh);
	boundaryCorrection(Cp,mesh, boundaryMesh);
	boundaryCorrection(lambda,mesh, boundaryMesh);
	boundaryCorrection(N2,mesh, boundaryMesh);
	boundaryCorrection(N,mesh, boundaryMesh);
	boundaryCorrection(O2,mesh, boundaryMesh);
	boundaryCorrection(O,mesh, boundaryMesh);
	boundaryCorrection(NO,mesh, boundaryMesh);
	boundaryCorrection(c,mesh, boundaryMesh);
};


// --- Boundary field update function for air 11
void ReadTable::boundaryCorrectionField (
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
	volScalarField& eMoins)
{
	boundaryCorrection(T,mesh, boundaryMesh);
	boundaryCorrection(p,mesh, boundaryMesh);
	boundaryCorrection(mu,mesh, boundaryMesh);
	boundaryCorrection(Cv,mesh, boundaryMesh);
	boundaryCorrection(Cp,mesh, boundaryMesh);
	boundaryCorrection(lambda,mesh, boundaryMesh);
	boundaryCorrection(N2,mesh, boundaryMesh);
	boundaryCorrection(N,mesh, boundaryMesh);
	boundaryCorrection(O2,mesh, boundaryMesh);
	boundaryCorrection(O,mesh, boundaryMesh);
	boundaryCorrection(NO,mesh, boundaryMesh);
	boundaryCorrection(c,mesh, boundaryMesh);
	boundaryCorrection(NOPlus,mesh, boundaryMesh);
	boundaryCorrection(NPlus,mesh, boundaryMesh);
	boundaryCorrection(OPlus,mesh, boundaryMesh);
	boundaryCorrection(N2Plus,mesh, boundaryMesh);
	boundaryCorrection(O2Plus,mesh, boundaryMesh);
	boundaryCorrection(eMoins,mesh, boundaryMesh);
};


// --- Interpolation function
scalar ReadTable::interpolation(volScalarField& field)
{
	return rhoCoef*(eCoef*fieldToList(field)[max_eIndex] + (1-eCoef)*fieldToList(field)[min_eIndex]) + (1-rhoCoef)*(eCoef*fieldToList(field)[max_eIndex - nRho] + (1-eCoef)*fieldToList(field)[min_eIndex-nRho]);
};


// --- Counting of values of rho in the table
void ReadTable::countRho()
{
	while(rhoList[nRho] == rhoList[nRho+1])
	{
    	nRho++;
  	}	
		nRho++;
};


// --- Size of the lists
void ReadTable::sizeTable()
{
	sizeTable_ = temperatureList.size();
};


// --- Maximum of rho, used to clamp
void ReadTable::maxRho ()
{
	maxRho_ = max(rhoList);
};


// --- Coefficient of the interpolation for rho
void ReadTable::interRho(scalar value)
{
	for (int i = 0 ; i < sizeTable_ ; i+=nRho) //Browse rho values
	{
		if (value > maxRho_)      //Bound high
		{
			rhoIndex = sizeTable_ - nRho;
			i = sizeTable_;
			rhoCoef = 1;
		}
		else if (value < rhoList[i])   //Computation and bound low
		{
  			if( i == 0) //Bound low
  			{
				rhoIndex = nRho;  
				i = sizeTable_;
				rhoCoef = 0;
  			}
  			else //Computation
  			{
				rhoIndex = i; 
				i = sizeTable_;
				rhoCoef = (value - rhoList[rhoIndex-1])/(rhoList[rhoIndex] - rhoList[rhoIndex-1]);
			}
		}	
	}	
};


// --- Coefficient of the interpolation for the normalized energy
void ReadTable::inter_e (scalar value)
{
	for (int j = rhoIndex ; j <= rhoIndex + nRho-1 ; j++)
	{
		min_e = energyList[rhoIndex];
		max_e = energyList[rhoIndex + nRho -1];

		if(value <= min_e) //Bound low
		{
			min_eIndex = rhoIndex;
			max_eIndex = rhoIndex;
			eCoef = 0;
			j = sizeTable_;
		}
		else if(value > max_e) //Bound high
		{
			min_eIndex = rhoIndex + nRho -1;
			max_eIndex = rhoIndex + nRho -1;
			eCoef = 1;
			j = sizeTable_;
		}		
		else if( (value < energyList[j]) && (stop_ ==false)) //Computation
		{
			max_eIndex = j;
			min_eIndex = j-1;
			eCoef = (value-energyList[min_eIndex])/(energyList[max_eIndex] - energyList[min_eIndex] );
			stop_=true;
		}
	}
};


// --- Boundary correction
void ReadTable::boundaryCorrection(volScalarField& field, fvMesh& mesh, const polyBoundaryMesh& boundaryMesh)
{
	forAll(mesh.boundary(), patch)
	{
	  forAll(mesh.boundary()[patch], facei)
	  {
		face = boundaryMesh[patch].start() + facei;
		own = mesh.owner()[face];
		field.boundaryField()[patch][facei] = field[own];
	  }
	}
};
