/*==================================================================================*/
/*---------------------------Creation of chemistry fields---------------------------*/
/*==================================================================================*/

/*---------------------Mass fraction for equilibrium air5/air11---------------------*/
Info << "Creating chemistry fields \n" << endl;
// --- Oxygen atom
volScalarField O
(
    IOobject
    (
        "O",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("O", dimless, 0.0)
);

// --- Oxygen
volScalarField O2
(
    IOobject
    (
        "O2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("O2", dimless, 0.0)
);

// --- Nitrogen
volScalarField N2
(
    IOobject
    (
        "N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("N2", dimless, 0.0)
);

// --- Nitrogen atom
volScalarField N
(
    IOobject
    (
        "N",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("N", dimless, 0.0)
);

// --- Nitrogen monoxide
volScalarField NO
(
    IOobject
    (
        "NO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("NO", dimless, 0.0)
);

// --- Electron
volScalarField eMoins
(
    IOobject
    (
        "eMoins",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("eMoins", dimless, 0.0)
);

// --- N+
volScalarField NPlus
(
    IOobject
    (
        "NPlus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("NPlus", dimless, 0.0)
);

// --- O+
volScalarField OPlus
(
    IOobject
    (
        "OPlus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("OPlus", dimless, 0.0)
);

// --- N2+
volScalarField N2Plus
(
    IOobject
    (
        "N2Plus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("N2Plus", dimless, 0.0)
);

// --- O2+
volScalarField O2Plus
(
    IOobject
    (
        "O2Plus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("O2Plus", dimless, 0.0)
);

// --- NO+
volScalarField NOPlus
(
    IOobject
    (
        "NOPlus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("NOPlus", dimless, 0.0)
);


/*------------------------Mass fraction for finite rate model------------------------*/
// --- Nitrogen
volScalarField Y_N2
(
    IOobject
    (
        "Y_N2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
	mesh
);

// --- Define rho species boundray type
const volScalarField::GeometricBoundaryField& rhoSpeciesBf = Y_N2.boundaryField();
wordList rhoSpeciesBoundaryType = rhoSpeciesBf.types();

// --- Oxygen
volScalarField Y_O2
(
    IOobject
    (
        "Y_O2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
	mesh
);

// --- Nitrogen atom
volScalarField Y_N
(
    IOobject
    (
        "Y_N",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
	Y_N2*0,
    rhoSpeciesBoundaryType
);

// --- Nitrogen monoxide
volScalarField Y_NO
(
	 IOobject
    (
        "Y_NO",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
	Y_N2*0,
    rhoSpeciesBoundaryType
);

// --- Oxygen atom
volScalarField Y_O
(
    IOobject
    (
        "Y_O",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
	Y_N2*0,
    rhoSpeciesBoundaryType
);

// --- Test if mass fraction sum = 1 (only for finite rate model)
if (model == "finiteRate")
{
    if (Y_O[0]+Y_NO[0]+Y_N[0]+Y_N2[0]+Y_O2[0] != 1){
	    FatalErrorIn
        (
            "Create chemistry fields"
        )
        << "sum of mass fraction is not equal to 1"
        << abort(FatalError);
    }
}

/*-----------------------------Species density-----------------------------*/
// --- Nitrogen density
volScalarField rhoN2
(
    IOobject
    (
        "rhoN2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*Y_N2,
	rhoSpeciesBoundaryType
);

// --- Oxygen density
volScalarField rhoO2
(
    IOobject
    (
        "rhoO2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*Y_O2,
	rhoSpeciesBoundaryType
);

// --- Nitrogen atom density
volScalarField rhoN
(
    IOobject
    (
        "rhoN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhoN2*Y_N,
    rhoSpeciesBoundaryType
);

// --- Nitrogen monoxide density
volScalarField rhoNO
(
	 IOobject
    (
        "rhoNO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhoN2*Y_NO,
    rhoSpeciesBoundaryType
);

// --- Oxygen density
volScalarField rhoO
(
    IOobject
    (
        "rhoO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhoN2*Y_O,
    rhoSpeciesBoundaryType
);

/*-----------------------------Species momentum-----------------------------*/
// --- rhoN2*U
volVectorField rhoUN2
(
    IOobject
    (
        "rhoUN2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhoN2*U,
    UBoundaryTypes
);
rhoUN2.boundaryField() = rhoN2.boundaryField()*U.boundaryField();

// --- rhoN*U
volVectorField rhoUN
(
    IOobject
    (
        "rhoUN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhoN*U,
    UBoundaryTypes
);
rhoUN.boundaryField() = rhoN.boundaryField()*U.boundaryField();

// --- rhoNO*U
volVectorField rhoUNO
(
    IOobject
    (
        "rhoUNO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhoNO*U,
    UBoundaryTypes
);
rhoUNO.boundaryField() = rhoNO.boundaryField()*U.boundaryField();

// --- rhoO2*U
volVectorField rhoUO2
(
    IOobject
    (
        "rhoUO2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhoO2*U,
    UBoundaryTypes
);
rhoUO2.boundaryField() = rhoO2.boundaryField()*U.boundaryField();

// --- rhoO*U
volVectorField rhoUO
(
    IOobject
    (
        "rhoUO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhoO*U,
    UBoundaryTypes
 );
rhoUO.boundaryField() = rhoO.boundaryField()*U.boundaryField();

/*---------------------------Species diffusivity---------------------------*/
volVectorField D_N2
(
    IOobject
    (
        "D_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("D_N2", dimensionSet(0,1,-1,0,0,0,0), vector::zero)  //wrong dimension to correct grad C dimension
);

volScalarField D_average
(
    IOobject
    (
        "D_average",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("D_average", dimensionSet(0,2,-1,0,0,0,0), SMALL)  //dimension is ok here.
);

volScalarField C_N2
(
    IOobject
    (
        "C_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("C_N2", dimless, SMALL)
);

volVectorField GradC_N2
(
    IOobject
    (
        "GradC_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("GradC_N2", dimensionSet(0,-1,0,0,0,0,0), vector::zero)
);

volVectorField D_averageN2
(
    IOobject
    (
        "D_averageN2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("D_averageN2", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)  //wrong dimension to correct grad C dimension
);

volScalarField C_N = C_N2;
volScalarField C_NO = C_N2;
volScalarField C_O2 = C_N2;
volScalarField C_O = C_N2;

volVectorField GradC_N = GradC_N2;
volVectorField GradC_NO = GradC_N2;
volVectorField GradC_O2 = GradC_N2;
volVectorField GradC_O = GradC_N2;

volVectorField D_N = D_N2;
volVectorField D_NO = D_N2;
volVectorField D_O2 = D_N2;
volVectorField D_O = D_N2;

volVectorField D_averageN = D_averageN2;
volVectorField D_averageNO = D_averageN2;
volVectorField D_averageO = D_averageN2;
volVectorField D_averageO2 = D_averageN2;

/*--------------------------------Forward coeffcients--------------------------------*/
// --- Forward reaction rate for N2 = 2N
/*volScalarField Kf_1
(
    IOobject
    (
        "Kf_1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Kf_1", dimensionSet(0,3,-1,0,-1,0,0), SMALL)
);

// --- Forward reaction rate for O2 = 2O
volScalarField Kf_2
(
    IOobject
    (
        "Kf_2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Kf_2", dimensionSet(0,3,-1,0,-1,0,0), SMALL)
);

// --- Forward reaction rate for NO = N + O
volScalarField Kf_3
(
    IOobject
    (
        "Kf_3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Kf_3", dimensionSet(0,3,-1,0,-1,0,0), SMALL)
);

// --- Forward reaction rate for N2 + O = NO + N
volScalarField Kf_4
(
    IOobject
    (
        "Kf_4",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Kf_4", dimensionSet(0,3,-1,0,-1,0,0), SMALL)
);

// --- Forward reaction rate for O2 + N = NO + O
volScalarField Kf_5
(
    IOobject
    (
        "Kf_5",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Kf_5", dimensionSet(0,3,-1,0,-1,0,0), SMALL)
);*/

/*---------------------------------Production rates---------------------------------*/
volScalarField omega_N
(
    IOobject
    (
        "omega_N",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("omega_N", dimensionSet(1,-3,-1,0,0,0,0), SMALL)
);

volScalarField omega_NO
(
    IOobject
    (
        "omega_NO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("omega_NO", dimensionSet(1,-3,-1,0,0,0,0), SMALL)
);


volScalarField omega_O
(
    IOobject
    (
        "omega_O",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("omega_O", dimensionSet(1,-3,-1,0,0,0,0), SMALL)
);


volScalarField omega_O2
(
    IOobject
    (
        "omega_O2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("omega_O2", dimensionSet(1,-3,-1,0,0,0,0), SMALL)
);


volScalarField omega_N2
(
    IOobject
    (
        "omega_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("omega_N2", dimensionSet(1,-3,-1,0,0,0,0), SMALL)
);


/*-------------------------------Backward coeffcients-------------------------------*/
/*
volScalarField Kb_1 = Kf_1;
volScalarField Kb_2 = Kf_1;
volScalarField Kb_3 = Kf_1;
volScalarField Kb_4 = Kf_1;
volScalarField Kb_5 = Kf_1;
*/


/*--------------------------Definition of a mixture object--------------------------*/
Mutation::MixtureOptions opts;
Mutation::Mixture* mixPtr = NULL;
if (model == "finiteRate" || model == "readTable")
{
    opts.setSpeciesDescriptor("N2 N NO O O2"); // --- 5 species air5 mixture
    opts.setThermodynamicDatabase("RRHO");     // --- Thermo database is RRHO
    opts.setMechanism("air5");                 // --- Mechanism used for air5 model
    opts.setStateModel("ChemNonEq1T");         // --- chemical noneq. w/ 1T
    mixPtr = new Mutation::Mixture(opts);
}

/*-----------------------------Variable for finite rate-----------------------------*/
// --- Vector of rho_s
double rhoi [5];

// --- Vector of omega_s
double wdot [5];

// --- Jacobian matrix in row major ordering
double jac [5*5];

// --- Block Coupled Solution Vector
volVector6Field rhoS   
(
    IOobject
    (
        "rhoS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector6("zero", dimless, vector6::zero)
);

/*--------------------------------Species properties--------------------------------*/
// --- Molar mass
dimensionedScalar M_N("M_N",dimensionSet(1,0,0,0,-1,0,0), 0.014);
dimensionedScalar M_NO("M_NO",dimensionSet(1,0,0,0,-1,0,0), 0.030);
dimensionedScalar M_N2("M_N2",dimensionSet(1,0,0,0,-1,0,0), 0.028);
dimensionedScalar M_O("M_O",dimensionSet(1,0,0,0,-1,0,0), 0.016);
dimensionedScalar M_O2("M_O2",dimensionSet(1,0,0,0,-1,0,0), 0.032);

// --- Species gas constant
volScalarField RN = R/0.014;
volScalarField RNO = R/0.03;
volScalarField RN2 = R/0.028;
volScalarField RO = R/0.016;
volScalarField RO2 = R/0.032;

// --- Specific species gas constant
volScalarField r_N = RN/M_N;
volScalarField r_NO = RNO/M_NO;
volScalarField r_N2 = RN2/M_N2;
volScalarField r_O = RO/M_O;
volScalarField r_O2 = RO2/M_O2;

// --- Number to correct dimensions
dimensionedScalar rhoDivM("rhoDivM",dimensionSet(0,-3,0,0,1,0,0), 1);
dimensionedScalar DiffGradC("rhoDivM",dimensionSet(0,1,-1,0,0,0,0), 1); //Correct dimension

//Diffusion matrix (air_diffusion_comparison.cpp -- voir Scoggins sur dropbox)
Eigen::MatrixXd m_Dij(5, 5);

/*------------------------------Species jacobian values------------------------------*/
// --- Jac_X_Y = d(omega_X)/d(rho_Y)

volScalarField jac_N_N2
(
    IOobject
    (
        "jac_N_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_N_N2", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_N_N
(
    IOobject
    (
        "jac_N_N",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_N_N", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_N_NO 
(
    IOobject
    (
        "jac_N_NO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_N_NO", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);
volScalarField jac_N_O 
(
    IOobject
    (
        "jac_N_O",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_N_O", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);
volScalarField jac_N_O2 
(
    IOobject
    (
        "jac_N_O2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_N_O2", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_N2_N2
(
    IOobject
    (
        "jac_N2_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_N2_N2", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_N2_N
(
    IOobject
    (
        "jac_N2_N",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_N2_N", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_N2_NO 
(
    IOobject
    (
        "jac_N2_NO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_N2_NO", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);
volScalarField jac_N2_O 
(
    IOobject
    (
        "jac_N2_O",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_N2_O", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);
volScalarField jac_N2_O2 
(
    IOobject
    (
        "jac_N2_O2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_N2_O2", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_NO_N2
(
    IOobject
    (
        "jac_NO_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_NO_N2", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_NO_N
(
    IOobject
    (
        "jac_NO_N",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_NO_N", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_NO_NO 
(
    IOobject
    (
        "jac_NO_NO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_NO_NO", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);
volScalarField jac_NO_O 
(
    IOobject
    (
        "jac_NO_O",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_NO_O", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);
volScalarField jac_NO_O2 
(
    IOobject
    (
        "jac_O2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_O2", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_O_N2
(
    IOobject
    (
        "jac_O_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_O_N2", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_O_N
(
    IOobject
    (
        "jac_O_N",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_O_N", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_O_NO 
(
    IOobject
    (
        "jac_O_NO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_O_NO", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);
volScalarField jac_O_O 
(
    IOobject
    (
        "jac_O_O",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_O_O", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);
volScalarField jac_O_O2 
(
    IOobject
    (
        "jac_O_O2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_O_O2", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_O2_N2
(
    IOobject
    (
        "jac_O2_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_O2_N2", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_O2_N
(
    IOobject
    (
        "jac_O2_N",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_O2_N", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);

volScalarField jac_O2_NO 
(
    IOobject
    (
        "jac_O2_NO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_O2_NO", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);
volScalarField jac_O2_O 
(
    IOobject
    (
        "jac_O2_O",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_O2_O", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);
volScalarField jac_O2_O2 
(
    IOobject
    (
        "jac_O2_O2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("jac_O2_O2", dimensionSet(0,0,-1,0,0,0,0), SMALL)
);
