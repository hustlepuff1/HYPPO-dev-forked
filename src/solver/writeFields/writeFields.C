/*=========================================================================================================*/
/*------------------------------------------Functions definitions------------------------------------------*/
/*=========================================================================================================*/

#include "writeFields/writeFields.H"

// --- Print the values if the energy is wrong
void Write::msgNegativeEnergy(int cellI, volScalarField& rho, volScalarField& rhoE, volScalarField& e, volVectorField& U, volScalarField& T, volScalarField& psi)
{
	Info <<"celli = " << cellI << " ; rho =" << rho[cellI] << endl;
	Info <<"celli = " << cellI << " ; rhoE =" << rhoE[cellI] << endl;
	Info <<"celli = " << cellI << " ; e =" << e[cellI] << endl;
	Info <<"celli = " << cellI << " ; rhoE/rho = " << rhoE[cellI]/rho[cellI] << endl;
	Info <<"celli = " << cellI << " ; 0.5*u2 = " << 0.5*magSqr(U[cellI]) << endl;
	Info <<"celli = " << cellI << " ; T = " << T[cellI] << endl;
	Info <<"celli = " << cellI << " ; psi = " << psi[cellI] << endl;
};


// --- Writes mu 
void Write::writeFields(volScalarField& mu)
{
	mu.write();
};


// --- Writes mu and air5 mass fraction species
void Write::writeFields(volScalarField& mu, volScalarField& O2, volScalarField& N2, volScalarField& NO, volScalarField& N, volScalarField& O)
{
	mu.write();
	O.write();
	O2.write();
	NO.write();
	N.write();
	N2.write();
};


// --- Writes mu and air11 mass fraction species
void Write::writeFields(volScalarField& mu, volScalarField& O2, volScalarField& N2, volScalarField& NO, volScalarField& N, volScalarField& O, volScalarField& O2Plus, volScalarField& N2Plus, volScalarField& NOPlus, volScalarField& NPlus, volScalarField& OPlus, volScalarField& eMoins)
{
	mu.write();
	O.write();
	O2.write();
	NO.write();
	N.write();
	N2.write();
	NOPlus.write();
	OPlus.write();
	eMoins.write();
	N2Plus.write();
	O2Plus.write();
	NPlus.write();
};


// --- Writes muT and y+
void Write::writeFields(volScalarField& muT, volScalarField& yPlus)
{
	muT.write();
	yPlus.write();
};


