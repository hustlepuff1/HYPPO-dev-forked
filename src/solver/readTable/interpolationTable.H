/*=========================================================================================================*/
/*------------------------------------------Interpolation from table---------------------------------------*/
/*=========================================================================================================*/

#include <fstream>
#include <iostream>
#include <string>

// --- Interpolation of the table values
forAll(p,cellI)
{
	// -- Used to stop parsing values when correct value found
 	table.stop(false);
	
	// --- Coefficients for interpolation
	table.interRho(rho[cellI]);
	table.inter_e(e[cellI]);

	// --- Interpolation
	if(table.airModel11())
		table.interpolationField(cellI,T,p,mu,Cv,Cp,lambda,N2,N,NO,O2,O,c,eMoins,N2Plus,O2Plus,NPlus,OPlus,NOPlus);
		
	else
		table.interpolationField(cellI,T,p,mu,Cv,Cp,lambda,N2,N,NO,O2,O,c);
}

// --- Boundary correction
if(table.airModel11())
	table.boundaryCorrectionField(mesh,boundaryMesh,T,p,mu,Cv,Cp,lambda,N2,N,NO,O2,O,c);
else
	table.boundaryCorrectionField(mesh,boundaryMesh,T,p,mu,Cv,Cp,lambda,N2,N,NO,O2,O,c,eMoins,N2Plus,O2Plus,NPlus,OPlus,NOPlus);

