/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    hyppo

Description
    Transient solver for supersonic and hypersonic flow, with implicit
    coupling between density, momentum and energy achieved by fvBlockMatrix.


\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "basicPsiThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"
#include "interpolation.H"
#include "volPointInterpolation.H"
#include "leastSquaresVolPointInterpolation.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "fixedValuePointPatchFields.H"
#include "wallDist.H"
#include "error.H"
#include "nearWallDist.H"
#include "wallFvPatch.H"
#include <iostream>
#include <string>
#include <fstream>

#undef Log
#include "mutation++.h"

#include "readTable/table.C"
#include "simControl/simControl.C"
#include "turbulence/turbulence.C"
#include "writeFields/writeFields.C"
#include "numerics/Diamond.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{

    #include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "finiteRate/createChemistryFields.H"
	#include "turbulence/createTurbulenceField.H"

/*-------------------------------------Constants used for cumputation--------------------------------------*/
    const fvPatchList& patches = mesh.boundary();  					// --- initialize patches
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh(); 	// --- All boundaries
    const unallocLabelList& owner = mesh.owner();					// --- label to current cell
    const unallocLabelList& neighbour = mesh.neighbour();			// --- label to neighbour cells
	const surfaceVectorField& Cf = mesh.Cf();   					// --- Face center coordinates
	const volVectorField& C = mesh.C(); 							// --- cells center coordinates
    const surfaceVectorField normalVector( mesh.Sf()/mesh.magSf()); // --- Normal vector to calculate contravariant vectors
	const pointField& points = mesh.points(); 						// --- nodes coordinates
	const surfaceScalarField& magSf = mesh.magSf();					// --- Faces of surfaces
	const faceList& faces = mesh.faces(); 							// --- Faces of cells
	const volPointInterpolation vpi(mesh);							// --- Point interpolation of mesh

/*---------------------------------------Initialize internal energy----------------------------------------*/
	if (model == "finiteRate")
	{
		Mutation::Mixture& mix = *mixPtr;
		forAll(e, cellI)
		{
		    rhoi[mix.speciesIndex("N2")] = rhoN2[cellI];
			rhoi[mix.speciesIndex("N")] = rhoN[cellI];
			rhoi[mix.speciesIndex("NO")] = rhoNO[cellI];
			rhoi[mix.speciesIndex("O")] = rhoO[cellI];
			rhoi[mix.speciesIndex("O2")] = rhoO2[cellI];
			mix.setState(rhoi,&T[cellI],1);
			e[cellI]=mix.mixtureEnergyMass();
		}
		E = e + 0.5*magSqr(U);
		rhoE = rho*E;
	}
/*-----------------------------------------Open file for readTable-----------------------------------------*/
	if (table.openFile(model))
	{
		// --- Add all data in list
		table.reading();

		// --- Count rhoe values for one rho
		table.countRho();

		// --- Get the total size of the lists
		table.sizeTable();

		// --- Get the maximum value of rho
		table.maxRho();
	}

/*---------------------------------------Initialize internal fields----------------------------------------*/
	if (model == "readTable")
	{
		if (table.airModel11())
		{
			opts.setSpeciesDescriptor("e- N N+ O O+ NO N2 N2+ O2 O2+ NO+"); // --- 11 species air11 mixture
			opts.setMechanism("parkair01");									// --- Mechanism used for air 11 model
		}
		else opts.setMechanism("air5");										// --- Mechanism used for air 5 model
		opts.setStateModel("EquilTP");         								// --- chemical equilibrium
		Mutation::Mixture mix(opts);

		forAll(e,cellI)
		{
			mix.setState(&T[cellI], &p[cellI]);
			e[cellI]=mix.mixtureEnergyMass();
		}

		E = e + 0.5*magSqr(U);
		rhoE = rho*E;
		#include "readTable/interpolationTable.H"
	}

/*------------------------------------------------Time Loop------------------------------------------------*/
	Info<< "\nStarting time loop\n" << endl;
    
	while (runTime.loop())
    {
		Info << "Time = " << runTime.timeName() << nl << endl;

/*------------------------------------Control parameters initialization------------------------------------*/
		// --- Variable to restart an iteration to recalculate the courant number
		sim.restart(true);
		// --- Boolean becoming true if internal energy is negative
		sim.negative_energy(false);	
		// --- Boolean to stop simulation if frozen
		sim.stopMinCo(false);		

/*------------------------------------------Gradient computation-------------------------------------------*/
      	pointVectorField Up = vpi.interpolate(U);
		pointScalarField Tp = vpi.interpolate(T); 

		grad.diamondFace(mesh, C, U, T, faces, points, Up,  Tp,  owner,  neighbour,  magSf, GradU, GradT);
		grad.diamondBoundary(mesh, C, U, T, faces,  points,  Up,  Tp,  owner,  neighbour,  magSf,  GradU,  GradT, boundaryMesh, Cf);

/*-------------------------------------Secondary variables computation-------------------------------------*/
		volTensorField volGradU("volGradU", fvc::reconstruct(GradU * magSf));
		volVectorField volGradT("volGradT", fvc::reconstruct(GradT * magSf));
		volTensorField tauMC("tauMC",(mu)*volGradU + (mu)*dev2(volGradU.T()));
		volTensorField htauMC("htauMC", tauMC * SMALL);
		volVectorField q("q", lambda*volGradT);
		volVectorField tauVq("tauVq" , (tauMC & U) + q);
		volVectorField htauVq("htauVq", tauVq * SMALL);

		H = e + p/rho + 0.5*magSqr(U);

/*------------------------------------------Computes Rusanov flux------------------------------------------*/
		#include "defVariables.H"	
		#include "finiteRate/defVariablesChemistry.H"
		#include "numerics/rusanovFlux.H"
		#include "numerics/CenteredFluxViscosity.H"

/*------------------------------------------------Turbulence-----------------------------------------------*/
	// --- If turbulence in controlDict and required time
	if(turbulenceModel=="baldwinLomax" && (runTime.value()>activateTurbulence))
	{
		#include "turbulence/BaldwinLomax.H"
	}

/*-----------------------------------------------Finite Rate-----------------------------------------------*/
	// --- If chemistry model is finite rate
	if(model=="finiteRate")
	{
		#include "finiteRate/finiteRate.H"
	}

/*--------------------------------------------Table equilibrium--------------------------------------------*/
	//If not, constantProperties are used
	if(model=="readTable")
	{
		#include "readTable/interpolationTable.H"
	}

/*------------------------------------------------Sutherland-----------------------------------------------*/
	if(model=="sutherland")
	{
		mu = muValue*pow((T/Tref) , 1.5)*(Tref + S)/(T+S);  // --- sutherLand law, refer to openfoam doc
	}

	lambda = Cp*(mu/Pr + muT/PrT);

/*---------------------------------------------Flux computation--------------------------------------------*/
	surfaceScalarField phi("phi", f_Rho*mesh.magSf() );
	surfaceVectorField phiUp("phiUp", f_RhoU*mesh.magSf());
	surfaceScalarField phiEp("phiEp", f_RhoE*mesh.magSf());
	surfaceVectorField phiUpv("phiUpv", f_RhoU_visc*mesh.magSf());
	surfaceScalarField phiEpv("phiEpv", f_RhoE_visc*mesh.magSf());

/*----------------------------------------Deal with past computation---------------------------------------*/
	// --- Variables to save the solution values calculated at the previous time iteration
	volScalarField rho_oldTime = rho;
	volVectorField rhoU_oldTime = rhoU;
	volScalarField rhoE_oldTime = rhoE;

	// --- Update k value from the previous time iteration
	sim.switchK2K1();

/*----------------------------------------------Restart loop----------------------------------------------*/
	// --- This loop enables to recalculate the solutions with a lower courant number if e becomes wrong.
	while(sim.restart())
	{
		// --- By default it does the loop just once
		sim.restart(false);

		// ---Initialize the rhoUE block system (matrix, source and reference to Up)
		fvBlockMatrix<vector6> rhoUEEqn(rhoUE);

		// --- Summation of fluxes for each cells
		scalarField sumPhi
		(
		fvc::surfaceSum(mag(amaxSf))().internalField()
		);

		// --- Computes the Courant number
		sim.coNum(sumPhi, runTime, mesh);

		// --- Read time controls
		#include "readTimeControls.H"

		// --- CoMin initialization
		sim.minCo(runTime);

		// --- Set the new delta T
		sim.setDeltaT(adjustTimeStep, maxCo, maxDeltaT, runTime);

		// --- Reset the current solution values if the Courant number has been reset
		rho = rho_oldTime;
		rhoU = rhoU_oldTime;
		rhoE = rhoE_oldTime;

		// --- Assemble and insert density equation
		#include "rhoUE/rhoEqn.H"
		// --- Assemble and insert momentum equation
		#include "rhoUE/rhoUEqn.H"
		// --- Assemble and insert total energy
		#include "rhoUE/rhoEEqn.H"

		// --- Solve the block matrix
		sim.maxResidual = cmptMax(rhoUEEqn.solve().initialResidual());

		// --- Retrieve solution
		rhoUEEqn.retrieveSolution(0, rho.internalField());
		rhoUEEqn.retrieveSolution(1, rhoU.internalField());
		rhoUEEqn.retrieveSolution(4, rhoE.internalField());

		// --- Update the velocity field
		U.dimensionedInternalField() = rhoU.dimensionedInternalField()/rho.dimensionedInternalField();
		U.correctBoundaryConditions();
		rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

		// --- Update the internal energy field
		e = rhoE/rho - 0.5*magSqr(U);
		e.correctBoundaryConditions();

		// --- Update rhoE boundaryConditions
		rhoE.boundaryField() ==   rho.boundaryField()* ( e.boundaryField() + 0.5*magSqr(U.boundaryField()));

/*-------------------------Update pressure, temperature and properties of the gas--------------------------*/
/*------------------------------------------------Ideal gas------------------------------------------------*/
		if(model == "none" || model == "sutherland")
		{
			T = e/Cv;
			T.correctBoundaryConditions();
			r = Cp - Cv;
			rPsi = T*r;

			// --- Update rPsi and psi Fields
			forAll(rPsi,celli)
			{
				psi[celli] = 1.0/rPsi[celli];
			};

			rPsi.correctBoundaryConditions();
			psi.correctBoundaryConditions();
			T.correctBoundaryConditions();
			// --- Update pressure Field
			p.dimensionedInternalField() = rho.dimensionedInternalField()*rPsi.dimensionedInternalField() ;
			p.correctBoundaryConditions();
			// --- Update density boundary conditions
			rho.boundaryField() == psi.boundaryField()* p.boundaryField();
		}

		T.correctBoundaryConditions();
		psi = rho/p;
		psi.correctBoundaryConditions();
		p.correctBoundaryConditions();
		rPsi = 1/psi;
		rPsi.correctBoundaryConditions();
		rho.boundaryField() == psi.boundaryField()* p.boundaryField();

/*----------------------------------------Finite rate update field-----------------------------------------*/
		if(model == "finiteRate")
		{
			Mutation::Mixture& mix = *mixPtr;
			forAll(T,cellI)
			{
				rhoi[mix.speciesIndex("N2")] = rhoN2[cellI];
				rhoi[mix.speciesIndex("N")] = rhoN[cellI];
				rhoi[mix.speciesIndex("NO")] = rhoNO[cellI];
				rhoi[mix.speciesIndex("O")] = rhoO[cellI];
				rhoi[mix.speciesIndex("O2")] = rhoO2[cellI];

				scalar etot = e[cellI] * rho[cellI];

				// --- Update T,p,mu,lambda with e and rhospecies
				mix.setState(rhoi,&etot);
				T[cellI] = mix.T();
				p[cellI] = mix.P();
				mu[cellI]=mix.viscosity();
				lambda[cellI]=mix.heavyThermalConductivity();
			}

			p.correctBoundaryConditions();
			T.correctBoundaryConditions();
			rhoUN2 = rhoN2*U;
			rhoUN2.boundaryField() = rhoN2.boundaryField()*U.boundaryField();
			rhoUN2.correctBoundaryConditions();

			rhoUN= rhoN*U;
			rhoUN.boundaryField() = rhoN.boundaryField()*U.boundaryField();
			rhoUN.correctBoundaryConditions();

			rhoUO = rhoO*U;
			rhoUO.boundaryField() = rhoO.boundaryField()*U.boundaryField();
			rhoUO.correctBoundaryConditions();

			rhoUO2 = rhoO2*U;
			rhoUO2.boundaryField() = rhoO2.boundaryField()*U.boundaryField();
			rhoUO2.correctBoundaryConditions();

			rhoUNO = rhoNO*U;
			rhoUNO.boundaryField() = rhoNO.boundaryField()*U.boundaryField();
			rhoUNO.correctBoundaryConditions();

			rhoN2.correctBoundaryConditions();
			rhoN.correctBoundaryConditions();
			rhoO2.correctBoundaryConditions();
			rhoO.correctBoundaryConditions();
			rhoNO.correctBoundaryConditions();
		}

/*------------------------------Verification of the energy field----------------------------------*/
		forAll (e, cellI)
		{
			// --- Select the correct condtion depending of the model
			if(sim.testEnergy(e[cellI],model))
			{
				// --- Print a message
				fields.msgNegativeEnergy(cellI,rho,rhoE,e,U,T,psi);

				// --- Stay in the loop
				sim.restart(true);

				// --- Used to decrease the time step
				sim.negative_energy(true);

				break;
			}
			else
			{
				// --- Exit the loop
				sim.restart(false);
			}
		};

		// --- Deal with a parralel computation if energy is wrong
		sim.parControl();

/*-----------------------------------------Write Fields---------------------------------------------*/
		runTime.write();
		if(runTime.outputTime())
		{
			// --- If mu is not constant
			if(model == "sutherland")
				fields.writeFields(mu);

			// --- If finite rate chemistry is used
			else if(model == "finiteRate")
				fields.writeFields(mu,rhoN2,rhoO2,rhoNO,rhoN,rhoO);

			// --- If air 5 model is used
			else if(model == "readTable")
			{
				fields.writeFields(mu,N2,O2,NO,N,O);

				// --- If air 11 model is used
				if(table.airModel11())
					fields.writeFields(mu,N2,O2,NO,N,O,N2Plus,O2Plus,NOPlus,NPlus,OPlus,eMoins);
			}

			// --- Values for turbulence
			if(turbulenceModel == "baldwinLomax")
				fields.writeFields(muT,yPlus);
		}

		// --- Manage some iteration values
		sim.switchKK1();
	}

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

		// --- Increase k2
		sim.incremK2();

		// --- Check convergence
		sim.convergenceCheck(runTime);
    }

    Info << "End\n" << endl;

    return 0;
}
