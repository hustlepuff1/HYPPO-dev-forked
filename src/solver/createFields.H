/*==================================================================================*/
/*--------------------Creation of files and fields and variables--------------------*/
/*==================================================================================*/

/*--------------------------------Creation of files---------------------------------*/
// --- Definition of properties of air
IOdictionary constantProperties
(
    IOobject
    (
        "constantProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// --- Choice of the model
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// --- Activation of turbulence
IOdictionary turbulenceProperties
(
    IOobject
    (
        "turbulenceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

/*------------------------------------Models------------------------------------*/
// --- Gas/Chemistry model
word model
(
  transportProperties.lookup("model")
);

// --- Turbulence model
word turbulenceModel
(
  turbulenceProperties.lookup("turbulenceModel")
);

// --- Starting time for the turbulence model
scalar activateTurbulence
(
readScalar(turbulenceProperties.lookup("activateTurbulence"))
);

// --- Prediction of the shock distance
scalar shockDistance
(
readScalar(turbulenceProperties.lookup("shockDistance"))
);

/*-----------------------------Reading values in files------------------------------*/
// --- Molar mass
dimensionedScalar MValue
(
    constantProperties.lookup("M")
);

// --- Heat capacity
dimensionedScalar Cpvalue
(
    constantProperties.lookup("Cp")
);

// --- Viscosity
dimensionedScalar muValue
(
    constantProperties.lookup("mu")
);

// --- Prandtl number
dimensionedScalar Pr
(
    constantProperties.lookup("Pr")
);

// --- Turbulent Pr for baldwin lomax
dimensionedScalar PrT = dimensionedScalar("PrT", dimless, 0.9); 

/*-------------------------Definition of secondary fields--------------------------*/
// --- Molar mass
volScalarField M
(
    IOobject
    (
        "M",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("M", dimensionSet(1,0,0,0,-1,0,0),MValue.value())
);

// --- Viscosity
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("mu", dimensionSet(1,-1,-1,0,0,0,0),muValue.value())
);

// --- Heat capacity at constant pressure
volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Cp", dimensionSet(0,2,-2,-1,0,0,0),Cpvalue.value())
);

// --- Heat capacity at constant volume
volScalarField Cv
(
    IOobject
    (
        "Cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Cv", dimensionSet(0,2,-2,-1,0,0,0),Cpvalue.value()/1.4)
);

// --- Specific gas constant
volScalarField r = Cp - Cv;

// --- Gas constant
volScalarField R = M*r;

// --- Turbulent viscosity
volScalarField muT
(
    IOobject
    (
        "muT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("muT", dimensionSet(1,-1,-1,0,0,0,0), SMALL)
);

// --- Enthalpy
volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("H", dimensionSet(0,2,-2,0,0,0,0), SMALL)
);

// --- Volume getter
volScalarField volPos
(
    IOobject
    (
        "volPos",
        runTime.timeName(),
        mesh
     ),
     mesh,
     dimensionedScalar("volPos", dimless, 1.0)
);

// --- Current cell getter
surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

// --- Neighbour cell getter
surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

// --- Small step for rho jacobian
volScalarField hrho
(
    IOobject
    (
        "hrho",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("hrho",  dimMass/dimVolume, 1e-09)
);

// --- Small step for rhoU jacobian
volScalarField hrhoU
(
    IOobject
    (
        "hrhoU",
        runTime.timeName(),
        mesh
     ),
     mesh,
     dimensionedScalar("hrhoU",  dimVelocity*dimMass/(dimVolume),1e-09)
);

// --- Small step for rhoE jacobian
volScalarField hrhoE
(
    IOobject
    (
         "hrhoE",
         runTime.timeName(),
         mesh
    ),
    mesh,
    dimensionedScalar("hrhoE",  dimEnergy/dimVolume, 1e-09)
);

/*--------------------------Definition of variable fields---------------------------*/
Info << "Creating Navier Stokes fields \n" << endl;
// --- Block Coupled Solution Vector
volVector6Field rhoUE
(
    IOobject
    (
        "rhoUE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector6("zero", dimless, vector6::zero)
);

// --- Velocity
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// --- UBoundaryTypes
const volVectorField::GeometricBoundaryField& Ubf = U.boundaryField();
wordList UBoundaryTypes = Ubf.types();

// --- Velocity gradient
surfaceVectorField GradU
(
    IOobject
    (
        "GradU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0,0,-1,0,0,0,0), vector::zero)
);

// --- Temperature
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// --- TBoundaryTypes
const volScalarField::GeometricBoundaryField& Tbf = T.boundaryField();
wordList TBoundaryTypes = Tbf.types();

// --- Temperature gradient
surfaceScalarField GradT
(
    IOobject
    (
        "GradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,-1,0,1,0,0,0), 0)
);

// --- Small step of velocity
volVectorField hU
(
    IOobject
    (
        "hU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    1e-09*U,
    UBoundaryTypes
);

// --- Pressure
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// --- Define a special boundary type for rho
#include "rhoBoundaryTypes.H"

// --- Psi = rho/p
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    M/(R*T),
    TBoundaryTypes
);

// --- rPsi = 1.0/Psi
volScalarField rPsi
(
    IOobject
    (
        "rPsi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    R*T/M,
    TBoundaryTypes
 );

// --- Density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p*psi,
    rhoBoundaryTypes
);

// --- This volScalarField is created to use the vector6 in foam extend only to solve 5 eq
volScalarField zeroVec6
(
    IOobject
    (
        "zeroVec6",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("zeroVec6",  dimless, 0)
);

// --- Momentum
volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U,
    UBoundaryTypes
);
rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

// --- Internal Energy
volScalarField e
(
    IOobject
    (
        "e",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Cv*T,
    TBoundaryTypes
);

// --- Total energy density
volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U)),
    TBoundaryTypes
);
rhoE.boundaryField() == rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));

// --- Total energy
volScalarField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (e + 0.5*magSqr(U)),
    TBoundaryTypes
);

// --- Speed of sound
volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    sqrt(rPsi* Cp/Cv),
    UBoundaryTypes
);

// --- Mach number
volScalarField Mach
(
    IOobject
    (
        "Mach",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mag(U)/c,
    UBoundaryTypes
);

// --- Thermal conductivity
volScalarField lambda = Cp*mu;
lambda.rename("lambda");

// --- Sutherland, W. (1893)
dimensionedScalar S("S", dimensionSet(0,0,0,1,0,0,0), 110.4);
dimensionedScalar Tref("Tref", dimensionSet(0,0,0,1,0,0,0), 273.15);

/*------------------------------Definition of objects--------------------------------*/
ReadTable table;
SimControl sim;
Turbulence turbulence;
Write fields;
Diamond grad;
