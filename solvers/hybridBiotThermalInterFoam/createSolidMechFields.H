// Creating Fields for Solid Mechanics 

Info<< "Reading Solid Velocity Field\n" << endl;
volVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



#include "createSurfaceFluxesSolid.H"

// Defining Solid Transport Model
singlePhaseTransportModel solid(Us, phis);     

// Reading Solid Density
dimensionedScalar rho_s ("rho_s",transportProperties.lookup("rho_s"));
// Reading solid heat conductivity
dimensionedScalar kappa_s ("kappa_s",transportProperties.lookup("kappa_s"));
// Reading solid heat capacity
dimensionedScalar Cp_s ("Cp_s",transportProperties.lookup("Cp_s"));

// Converting dimensionedScalar into volScalarField
volScalarField rhos
(
    IOobject
    (
        "rhos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    rho_s
);
//rhos.oldTime();

volScalarField Cps
(
    IOobject
    (
        "Cps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    Cp_s
);
//Cps.oldTime();

volScalarField kappas
(
    IOobject
    (
        "kappas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    kappa_s
);




Info<< "Creating field kinetic energy Ke_s\n" << endl;
volScalarField Ke_s("Ke_s", 0.5*magSqr(Us));  //solids
 
// Define bulk viscosity (these are not thermal conductivities)
volScalarField kappa_dyn0
(
    IOobject
    (
        "kappa_dyn0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    kappaDyn0
);


volScalarField kappa_dyn
(
    IOobject
    (
        "kappa_dyn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    kappaDyn0  
);

volScalarField kappa_kin
(
    IOobject
    (
        "kappa_kin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    kappaDyn0/rho_s  
);


// Reading solid Mechanics activator
scalar activateSolidMechanics (transportProperties.lookupOrDefault<scalar>("activateSolidMechanics",0.));
activateSolidMechanics = activateSolidMechanics; //prevents unused variable warning

scalar activateSwellingPressure (transportProperties.lookupOrDefault<scalar>("activateSwellingPressure",0.));
activateSwellingPressure = activateSwellingPressure; //prevents unused variable warning

// Reading Maximum Courant Number for Solid 
scalar maxCoUs = runTime.controlDict().lookupOrDefault<scalar>("maxCoUs", 1.0);
maxCoUs=maxCoUs; //prevents unused variable warning

