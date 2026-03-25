//update mass concentration
/*
C_N2 = rhoN2/rho;
C_N = rhoN/rho;
C_O2 = rhoO2/rho;
C_O = rhoO/rho;
C_NO = rhoNO/rho;

GradC_N2 = fvc::grad(C_N2);
GradC_N = fvc::grad(C_N);
GradC_NO = fvc::grad(C_NO);
GradC_O2 = fvc::grad(C_O2);
GradC_O = fvc::grad(C_O);
*/

Mutation::Mixture& mix = *mixPtr;

surfaceScalarField phiN2("phiN2", f_RhoN2*mesh.magSf() );
surfaceScalarField phiN("phiN", f_RhoN*mesh.magSf() );
surfaceScalarField phiNO("phiNO", f_RhoNO*mesh.magSf() );
surfaceScalarField phiO("phiO", f_RhoO*mesh.magSf() );
surfaceScalarField phiO2("phiO2", f_RhoO2*mesh.magSf() );

forAll(rho,cellI)
{
    rhoi[mix.speciesIndex("N2")] = rhoN2[cellI];
    rhoi[mix.speciesIndex("N")] = rhoN[cellI];
    rhoi[mix.speciesIndex("NO")] = rhoNO[cellI];
    rhoi[mix.speciesIndex("O")] = rhoO[cellI];
    rhoi[mix.speciesIndex("O2")] = rhoO2[cellI];

    mix.setState(rhoi,&T[cellI],1);
    //mix.forwardRateCoefficients(kf);
    //mix.backwardRateCoefficients(kb);
	mix.netProductionRates(wdot);
	mix.jacobianRho(jac);

	//d_wi/d_rhoi
	jac_N2_N2[cellI]=jac[0*mix.nSpecies()+0];
	jac_N2_N[cellI]=jac[0*mix.nSpecies()+1];
	jac_N2_NO[cellI]=jac[0*mix.nSpecies()+2];
	jac_N2_O[cellI]=jac[0*mix.nSpecies()+3];
	jac_N2_O2[cellI]=jac[0*mix.nSpecies()+4];

	jac_N_N2[cellI]=jac[1*mix.nSpecies()+0];
	jac_N_N[cellI]=jac[1*mix.nSpecies()+1];
	jac_N_NO[cellI]=jac[1*mix.nSpecies()+2];
	jac_N_O[cellI]=jac[1*mix.nSpecies()+3];
	jac_N_O2[cellI]=jac[1*mix.nSpecies()+4];	

	jac_NO_N2[cellI]=jac[2*mix.nSpecies()+0];
	jac_NO_N[cellI]=jac[2*mix.nSpecies()+1];
	jac_NO_NO[cellI]=jac[2*mix.nSpecies()+2];
	jac_NO_O[cellI]=jac[2*mix.nSpecies()+3];
	jac_NO_O2[cellI]=jac[2*mix.nSpecies()+4];

	jac_O_N2[cellI]=jac[3*mix.nSpecies()+0];
	jac_O_N[cellI]=jac[3*mix.nSpecies()+1];
	jac_O_NO[cellI]=jac[3*mix.nSpecies()+2];
	jac_O_O[cellI]=jac[3*mix.nSpecies()+3];
	jac_O_O2[cellI]=jac[3*mix.nSpecies()+4];

	jac_NO_N2[cellI]=jac[4*mix.nSpecies()+0];
	jac_NO_N[cellI]=jac[4*mix.nSpecies()+1];
	jac_NO_NO[cellI]=jac[4*mix.nSpecies()+2];
	jac_NO_O[cellI]=jac[4*mix.nSpecies()+3];
	jac_NO_O2[cellI]=jac[4*mix.nSpecies()+4];
	
    //mix.J(Dm);  //Mutation++ is able to give species diffusion flux.
    //scalar test = 0;
    //Average diffusion coeff
    //mix.averageDiffusionCoeffs(&test);
    //D_average[cellI] = test;

    //std::cout << mix.computeJacobian() << std::endl;
    //Full diffusionMatrix
/*
    m_Dij = mix.diffusionMatrix();

    D_N2 = (m_Dij(mix.speciesIndex("N2"),mix.speciesIndex("N2"))*GradC_N2[cellI] + m_Dij(mix.speciesIndex("N2"),mix.speciesIndex("N"))*GradC_N[cellI] +
    m_Dij(mix.speciesIndex("N2"),mix.speciesIndex("NO"))*GradC_NO[cellI] + m_Dij(mix.speciesIndex("N2"),mix.speciesIndex("O"))*GradC_O[cellI] +
    m_Dij(mix.speciesIndex("N2"),mix.speciesIndex("O2"))*GradC_O2[cellI])*DiffGradC;
    D_N = (m_Dij(mix.speciesIndex("N"),mix.speciesIndex("N"))*GradC_N[cellI] + m_Dij(mix.speciesIndex("N"),mix.speciesIndex("N2"))*GradC_N2[cellI] +
    m_Dij(mix.speciesIndex("N"),mix.speciesIndex("NO"))*GradC_NO[cellI] + m_Dij(mix.speciesIndex("N"),mix.speciesIndex("O"))*GradC_O[cellI] +
    m_Dij(mix.speciesIndex("N"),mix.speciesIndex("O2"))*GradC_O2[cellI])*DiffGradC;
    D_NO = (m_Dij(mix.speciesIndex("NO"),mix.speciesIndex("N"))*GradC_N[cellI] + m_Dij(mix.speciesIndex("NO"),mix.speciesIndex("N2"))*GradC_N2[cellI] +
    m_Dij(mix.speciesIndex("NO"),mix.speciesIndex("NO"))*GradC_NO[cellI] + m_Dij(mix.speciesIndex("NO"),mix.speciesIndex("O"))*GradC_O[cellI] +
    m_Dij(mix.speciesIndex("NO"),mix.speciesIndex("O2"))*GradC_O2[cellI])*DiffGradC;
    D_O2 = (m_Dij(mix.speciesIndex("O2"),mix.speciesIndex("N"))*GradC_N[cellI] + m_Dij(mix.speciesIndex("O2"),mix.speciesIndex("N2"))*GradC_N2[cellI] +
    m_Dij(mix.speciesIndex("O2"),mix.speciesIndex("NO"))*GradC_NO[cellI] + m_Dij(mix.speciesIndex("O2"),mix.speciesIndex("O"))*GradC_O[cellI] +
    m_Dij(mix.speciesIndex("O2"),mix.speciesIndex("O2"))*GradC_O2[cellI])*DiffGradC;
    D_O = (m_Dij(mix.speciesIndex("O"),mix.speciesIndex("N"))*GradC_N[cellI] + m_Dij(mix.speciesIndex("O"),mix.speciesIndex("N2"))*GradC_N2[cellI] +
    m_Dij(mix.speciesIndex("O"),mix.speciesIndex("NO"))*GradC_NO[cellI] + m_Dij(mix.speciesIndex("O"),mix.speciesIndex("O"))*GradC_O[cellI] +
    m_Dij(mix.speciesIndex("O"),mix.speciesIndex("O2"))*GradC_O2[cellI])*DiffGradC;
*/
/*
    D_averageN2[cellI] = rho[cellI]*D_average[cellI]*GradC_N2[cellI];
    D_averageN[cellI] = rho[cellI]*D_average[cellI]*GradC_N[cellI];
    D_averageNO[cellI] = rho[cellI]*D_average[cellI]*GradC_NO[cellI];
    D_averageO[cellI] = rho[cellI]*D_average[cellI]*GradC_O[cellI];
    D_averageO2[cellI] = rho[cellI]*D_average[cellI]*GradC_O2[cellI];
*/
    
	omega_N2[cellI] = wdot[mix.speciesIndex("N2")];
	omega_N[cellI] = wdot[mix.speciesIndex("N")];
	omega_NO[cellI] = wdot[mix.speciesIndex("NO")];
	omega_O[cellI] = wdot[mix.speciesIndex("O")];
	omega_O2[cellI] = wdot[mix.speciesIndex("O2")];
}


fvBlockMatrix<vector6> rhoSEqn(rhoS);
//no diffusion

fvScalarMatrix rhoN2Eqn (fvm::ddt(rhoN2) + fvm::div(jacobianRhoN2,rhoN2) - fvc::div(jacobianRhoN2,rhoN2) + fvc::div(phiN2) == 
omega_N2
/*+ fvm::Sp(jac_N2_N2,rhoN2) 
- jac_N2_N2*rhoN2
- jac_N2_N*rhoN 
- jac_N2_NO*rhoNO 
- jac_N2_O*rhoO 
- jac_N2_O2*rhoO2*/
);

fvScalarMatrix rhoNEqn (fvm::ddt(rhoN) + fvm::div(jacobianRhoN,rhoN) - fvc::div(jacobianRhoN,rhoN) + fvc::div(phiN) == 
omega_N 
/*+ fvm::Sp(jac_N_N,rhoN) 
- jac_N_N*rhoN
- jac_N_N2*rhoN2 
- jac_N_NO*rhoNO 
- jac_N_O*rhoO 
- jac_N_O2*rhoO2*/
);

fvScalarMatrix rhoNOEqn (fvm::ddt(rhoNO) + fvm::div(jacobianRhoNO,rhoNO) - fvc::div(jacobianRhoNO,rhoNO) + fvc::div(phiNO) == 
omega_NO
/*+ fvm::Sp(jac_NO_NO,rhoNO) 
- jac_NO_N*rhoN
- jac_NO_N2*rhoN2 
- jac_NO_NO*rhoNO 
- jac_NO_O*rhoO 
- jac_NO_O2*rhoO2*/
);

fvScalarMatrix rhoOEqn (fvm::ddt(rhoO) + fvm::div(jacobianRhoO,rhoO) - fvc::div(jacobianRhoO,rhoO) + fvc::div(phiO) == 
omega_O
/*+ fvm::Sp(jac_O_O,rhoO) 
- jac_O_N*rhoN
- jac_O_N2*rhoN2 
- jac_O_NO*rhoNO 
- jac_O_O*rhoO 
- jac_O_O2*rhoO2*/
);

fvScalarMatrix rhoO2Eqn (fvm::ddt(rhoO2) + fvm::div(jacobianRhoO2,rhoO2) - fvc::div(jacobianRhoO2,rhoO2) + fvc::div(phiO2) == 
omega_O2
/*+ fvm::Sp(jac_O2_O2,rhoO2) 
- jac_O2_N*rhoN
- jac_O2_N2*rhoN2 
- jac_O2_NO*rhoNO 
- jac_O2_O*rhoO 
- jac_O2_O2*rhoO2*/
);

/*fvScalarMatrix rhoN2Eqn (fvm::ddt(rhoN2) + fvm::div(jacobianRhoN2,rhoN2) - fvc::div(jacobianRhoN2,rhoN2) + fvc::div(phiN2) ==
	M_N2*(Kb_1*pow(rhoN/M_N,2) - Kf_1*rhoN2/M_N2*rhoDivM
		+ Kb_4*((rhoN/M_N) * (rhoNO/M_NO)) - Kf_4*((rhoN2/M_N2) * rhoO/M_O)));

fvScalarMatrix rhoNEqn (fvm::ddt(rhoN) + fvm::div(jacobianRhoN,rhoN) - fvc::div(jacobianRhoN,rhoN) + fvc::div(phiN) ==
	M_N*(-2*(Kb_1*pow(rhoN/M_N,2) - rhoDivM*Kf_1*rhoN2/M_N2) - Kb_4*((rhoN/M_N) * (rhoNO/M_NO) ) + Kf_4*((rhoN2/M_N2) * (rhoO/M_O))
+ Kb_3*((rhoNO/M_NO) * (rhoO/M_O) ) - Kf_3*((rhoN/M_N) * (rhoO2/M_O2)) -1*(Kb_5*((rhoN/M_N) * (rhoO/M_O) ) - Kf_5*((rhoNO/M_NO) * rhoDivM))));

fvScalarMatrix rhoNOEqn (fvm::ddt(rhoNO) + fvm::div(jacobianRhoNO,rhoNO) - fvc::div(jacobianRhoNO,rhoNO) + fvc::div(phiNO) ==
	M_NO*(-1*(Kb_4*((rhoN/M_N) * (rhoNO/M_NO) ) - Kf_4*((rhoN2/M_N2) * (rhoO/M_O)))
	- 1*( Kb_3*((rhoNO/M_NO) * (rhoO/M_O) ) - Kf_3*((rhoN/M_N) * (rhoO2/M_O2)))
	+ Kb_5*((rhoN/M_N) * (rhoO/M_O) ) - Kf_5*((rhoNO/M_NO) * rhoDivM)));

fvScalarMatrix rhoOEqn (fvm::ddt(rhoO) + fvm::div(jacobianRhoO,rhoO) - fvc::div(jacobianRhoO,rhoO) + fvc::div(phiO) ==
	M_O*(-2*(Kb_2*pow(rhoO/M_O,2) - rhoDivM*Kf_2*rhoO2/M_O2)
	+ Kb_4*((rhoN/M_N) * (rhoNO/M_NO) ) - Kf_4*((rhoN2/M_N2) * (rhoO/M_O))
	- Kb_3*((rhoNO/M_NO) * (rhoO/M_O) ) + Kf_3*((rhoN/M_N) * (rhoO2/M_O2)) -1*(Kb_5*((rhoN/M_N) * (rhoO/M_O) ) - Kf_5*((rhoNO/M_NO) * rhoDivM))));

fvScalarMatrix rhoO2Eqn (fvm::ddt(rhoO2) + fvm::div(jacobianRhoO2,rhoO2) - fvc::div(jacobianRhoO2,rhoO2) + fvc::div(phiO2) ==
	M_O2*(Kb_2*pow(rhoO/M_O,2) - rhoDivM*Kf_2*rhoO2/M_O2
	+ Kb_3*((rhoO/M_O) * (rhoNO/M_NO) ) - Kf_3*((rhoO2/M_O2) * rhoN/M_N)));*/


//Full diffusion Matrix
/*
fvScalarMatrix rhoN2Eqn (fvm::ddt(rhoN2) + fvm::div(jacobianRhoN2,rhoN2) - fvc::div(jacobianRhoN2,rhoN2) + fvc::div(phiN2) + fvc::div(rho*D_N2) ==
	M_N2*(Kb_1*pow(rhoN/M_N,2) - Kf_1*rhoN2/M_N2*rhoDivM
		+ Kb_4*((rhoN/M_N) * (rhoNO/M_NO)) - Kf_4*((rhoN2/M_N2) * rhoO/M_O)));

fvScalarMatrix rhoNEqn (fvm::ddt(rhoN) + fvm::div(jacobianRhoN,rhoN) - fvc::div(jacobianRhoN,rhoN) + fvc::div(phiN) + fvc::div(rho*D_N) ==
	M_N*(-2*(Kb_1*pow(rhoN/M_N,2) - rhoDivM*Kf_1*rhoN2/M_N2) - Kb_4*((rhoN/M_N) * (rhoNO/M_NO) ) + Kf_4*((rhoN2/M_N2) * (rhoO/M_O))  // + devant kf4
+ Kb_3*((rhoNO/M_NO) * (rhoO/M_O) ) - Kf_3*((rhoN/M_N) * (rhoO2/M_O2)) -1*(Kb_5*((rhoN/M_N) * (rhoO/M_O) ) - Kf_5*((rhoNO/M_NO) * rhoDivM))));

fvScalarMatrix rhoNOEqn (fvm::ddt(rhoNO) + fvm::div(jacobianRhoNO,rhoNO) - fvc::div(jacobianRhoNO,rhoNO) + fvc::div(phiNO)  + fvc::div(rho*D_NO) ==
	M_NO*(-1*(Kb_4*((rhoN/M_N) * (rhoNO/M_NO) ) - Kf_4*((rhoN2/M_N2) * (rhoO/M_O)))
	- 1*( Kb_3*((rhoNO/M_NO) * (rhoO/M_O) ) - Kf_3*((rhoN/M_N) * (rhoO2/M_O2)))
	+ Kb_5*((rhoN/M_N) * (rhoO/M_O) ) - Kf_5*((rhoNO/M_NO) * rhoDivM)));

fvScalarMatrix rhoOEqn (fvm::ddt(rhoO) + fvm::div(jacobianRhoO,rhoO) - fvc::div(jacobianRhoO,rhoO) + fvc::div(phiO)  + fvc::div(rho*D_O) ==
	M_O*(-2*(Kb_2*pow(rhoO/M_O,2) - rhoDivM*Kf_2*rhoO2/M_O2)
	+ Kb_4*((rhoN/M_N) * (rhoNO/M_NO) ) - Kf_4*((rhoN2/M_N2) * (rhoO/M_O))
	- Kb_3*((rhoNO/M_NO) * (rhoO/M_O) ) + Kf_3*((rhoN/M_N) * (rhoO2/M_O2)) -1*(Kb_5*((rhoN/M_N) * (rhoO/M_O) ) - Kf_5*((rhoNO/M_NO) * rhoDivM))));

fvScalarMatrix rhoO2Eqn (fvm::ddt(rhoO2) + fvm::div(jacobianRhoO2,rhoO2) - fvc::div(jacobianRhoO2,rhoO2) + fvc::div(phiO2)  + fvc::div(rho*D_O2) ==
	M_O2*(Kb_2*pow(rhoO/M_O,2) - rhoDivM*Kf_2*rhoO2/M_O2
	+ Kb_3*((rhoO/M_O) * (rhoNO/M_NO) ) - Kf_3*((rhoO2/M_O2) * rhoN/M_N)));
*/

//Average diffusion
/*
fvScalarMatrix rhoN2Eqn (fvm::ddt(rhoN2) + fvm::div(jacobianRhoN2,rhoN2) - fvc::div(jacobianRhoN2,rhoN2) + fvc::div(phiN2) + fvc::div(D_averageN2) ==
	M_N2*(Kb_1*pow(rhoN/M_N,2) - Kf_1*rhoN2/M_N2*rhoDivM
		+ Kb_4*((rhoN/M_N) * (rhoNO/M_NO)) - Kf_4*((rhoN2/M_N2) * rhoO/M_O)));

fvScalarMatrix rhoNEqn (fvm::ddt(rhoN) + fvm::div(jacobianRhoN,rhoN) - fvc::div(jacobianRhoN,rhoN) + fvc::div(phiN) + fvc::div(D_averageN) ==
	M_N*(-2*(Kb_1*pow(rhoN/M_N,2) - rhoDivM*Kf_1*rhoN2/M_N2) - Kb_4*((rhoN/M_N) * (rhoNO/M_NO) ) + Kf_4*((rhoN2/M_N2) * (rhoO/M_O))  // + devant kf4
+ Kb_3*((rhoNO/M_NO) * (rhoO/M_O) ) - Kf_3*((rhoN/M_N) * (rhoO2/M_O2)) -1*(Kb_5*((rhoN/M_N) * (rhoO/M_O) ) - Kf_5*((rhoNO/M_NO) * rhoDivM))));

fvScalarMatrix rhoNOEqn (fvm::ddt(rhoNO) + fvm::div(jacobianRhoNO,rhoNO) - fvc::div(jacobianRhoNO,rhoNO) + fvc::div(phiNO)  + fvc::div(D_averageNO) ==
	M_NO*(-1*(Kb_4*((rhoN/M_N) * (rhoNO/M_NO) ) - Kf_4*((rhoN2/M_N2) * (rhoO/M_O)))
	- 1*( Kb_3*((rhoNO/M_NO) * (rhoO/M_O) ) - Kf_3*((rhoN/M_N) * (rhoO2/M_O2)))
	+ Kb_5*((rhoN/M_N) * (rhoO/M_O) ) - Kf_5*((rhoNO/M_NO) * rhoDivM)));

fvScalarMatrix rhoOEqn (fvm::ddt(rhoO) + fvm::div(jacobianRhoO,rhoO) - fvc::div(jacobianRhoO,rhoO) + fvc::div(phiO) + fvc::div(D_averageO) ==
	M_O*(-2*(Kb_2*pow(rhoO/M_O,2) - rhoDivM*Kf_2*rhoO2/M_O2)
	+ Kb_4*((rhoN/M_N) * (rhoNO/M_NO) ) - Kf_4*((rhoN2/M_N2) * (rhoO/M_O))
	- Kb_3*((rhoNO/M_NO) * (rhoO/M_O) ) + Kf_3*((rhoN/M_N) * (rhoO2/M_O2)) -1*(Kb_5*((rhoN/M_N) * (rhoO/M_O) ) - Kf_5*((rhoNO/M_NO) * rhoDivM))));

fvScalarMatrix rhoO2Eqn (fvm::ddt(rhoO2) + fvm::div(jacobianRhoO2,rhoO2) - fvc::div(jacobianRhoO2,rhoO2) + fvc::div(phiO2) + fvc::div(D_averageO2) ==
	M_O2*(Kb_2*pow(rhoO/M_O,2) - rhoDivM*Kf_2*rhoO2/M_O2
	+ Kb_3*((rhoO/M_O) * (rhoNO/M_NO) ) - Kf_3*((rhoO2/M_O2) * rhoN/M_N)));
*/


fvScalarMatrix zero(fvm::ddt(zeroVec6) == fvc::ddt(zeroVec6));

rhoSEqn.insertEquation(0,rhoN2Eqn);
rhoSEqn.insertEquation(1,rhoNEqn);
rhoSEqn.insertEquation(2,rhoNOEqn);
rhoSEqn.insertEquation(3,rhoOEqn);
rhoSEqn.insertEquation(4,rhoO2Eqn);
rhoSEqn.insertEquation(5,zero);

/*rhoSEqn.insertEquationCoupling(0,1,fvm::Sp(jac_N2_N,rhoN));
rhoSEqn.insertEquationCoupling(0,2,fvm::Sp(jac_N2_NO,rhoNO));
rhoSEqn.insertEquationCoupling(0,3,fvm::Sp(jac_N2_O,rhoO));
rhoSEqn.insertEquationCoupling(0,4,fvm::Sp(jac_N2_O2,rhoO2));

rhoSEqn.insertEquationCoupling(0,0,fvm::Sp(jac_N_N2,rhoN2));
rhoSEqn.insertEquationCoupling(0,2,fvm::Sp(jac_N_NO,rhoNO));
rhoSEqn.insertEquationCoupling(0,3,fvm::Sp(jac_N_O,rhoO));
rhoSEqn.insertEquationCoupling(0,4,fvm::Sp(jac_N_O2,rhoO2));

rhoSEqn.insertEquationCoupling(2,0,fvm::Sp(jac_NO_N2,rhoN2));
rhoSEqn.insertEquationCoupling(2,1,fvm::Sp(jac_NO_N,rhoN));
rhoSEqn.insertEquationCoupling(2,3,fvm::Sp(jac_NO_O,rhoO));
rhoSEqn.insertEquationCoupling(2,4,fvm::Sp(jac_NO_O2,rhoO2));

rhoSEqn.insertEquationCoupling(3,0,fvm::Sp(jac_O_N2,rhoN2));
rhoSEqn.insertEquationCoupling(3,1,fvm::Sp(jac_O_N,rhoN));
rhoSEqn.insertEquationCoupling(3,2,fvm::Sp(jac_O_NO,rhoNO));
rhoSEqn.insertEquationCoupling(3,4,fvm::Sp(jac_O_O2,rhoO2));

rhoSEqn.insertEquationCoupling(4,0,fvm::Sp(jac_O2_N2,rhoN2));
rhoSEqn.insertEquationCoupling(4,1,fvm::Sp(jac_O2_N,rhoN));
rhoSEqn.insertEquationCoupling(4,2,fvm::Sp(jac_O2_NO,rhoNO));
rhoSEqn.insertEquationCoupling(4,3,fvm::Sp(jac_O2_O,rhoO));*/

//rhoSEqn.updateSourceCoupling();

sim.maxResidual = cmptMax(rhoSEqn.solve().initialResidual());

rhoSEqn.retrieveSolution(0, rhoN2.internalField());
rhoSEqn.retrieveSolution(1, rhoN.internalField());
rhoSEqn.retrieveSolution(2, rhoNO.internalField());
rhoSEqn.retrieveSolution(3, rhoO.internalField());
rhoSEqn.retrieveSolution(4, rhoO2.internalField());

/*Info << omega_N2.average() << endl;
Info << omega_N.average() << endl;
Info << min(omega_N) << endl;
Info << omega_O.average() << endl;
Info << omega_O2.average() << endl;

Info << "sum rho_s = " << rho.average()-(rhoN2.average()+rhoN.average()+rhoNO.average()+rhoO.average()+rhoO2.average()) << endl;

Info << rhoN2.average() << endl;
Info << rhoN.average() << endl;
Info << rhoNO.average() << endl;
Info << rhoO.average() << endl;
Info << rhoO2.average() << endl;
cin.ignore();*/
