/*==============================================================================================================*/
/*----------------------------------------Creating fields for turbulence----------------------------------------*/
/*==============================================================================================================*/

Info << "Creating turbulence fields \n" << endl;

volScalarField yPlus
(
  IOobject
  (
      "yPlus",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("yPlus", dimless, 0.0)
);

volScalarField Aplus
(
  IOobject
  (
    "Aplus",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("Aplus",dimless,26)
);

volScalarField uTau
(
    IOobject
    (
        "uTau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("uTau", dimVelocity, SMALL)
);


volScalarField muTinner
(
    IOobject
    (
        "muTinner",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("muTinner", dimensionSet(1,-1,-1,0,0,0,0), SMALL)
);

volScalarField muTouter
(
    IOobject
    (
        "muTouter",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("muTouter", dimensionSet(1,-3,0,0,0,0,0), SMALL)  //Wrong dimension also
);

volScalarField F
(
    IOobject
    (
        "F",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("F", dimensionSet(0,1,-1,0,0,0,0), SMALL)
);


volScalarField Ckleb
(
  IOobject
  (
    "Ckleb",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("Ckleb",dimensionSet(0,-1,0,0,0,0,0),0.3)  //Wrong dimension need to be improve (Fkleb computation)
);

volScalarField TauWall
(
  IOobject
  (
    "TauWall",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("TauWall",dimensionSet(1,-1,-1,0,0,0,0),0)  //Wrong dimension need to be improve (Fkleb computation)
);

volScalarField magVorticity
(
  IOobject
  (
      "magVorticity",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("magVorticity",dimensionSet(0,0,-1,0,0,0,0),0)
);


volScalarField rhoW = rho;
volScalarField muW = mu;
