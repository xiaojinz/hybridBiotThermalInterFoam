
Info<< "Reading field T\n" <<endl;
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


Info<< "Reading field p\n" << endl;
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



Info<< "Reading wetting fluid velocity field Uwetting\n" << endl;
volVectorField Uwetting
(
    IOobject
    (
        "Uwetting",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading wetting fluid velocity field UnonWetting\n" << endl;
volVectorField UnonWetting
(
    IOobject
    (
        "UnonWetting",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading fluid velocity field 'U'\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    UnonWetting + Uwetting
);

Info<< "Reading solid fraction field\n" << endl;     
volScalarField epss
(                                                      
    IOobject
    (              
         "epss", 
         runTime.timeName(),      
         mesh,            
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
    ),           
    mesh
);   

// Define fluid fraction (porosity)
volScalarField eps("eps",(1-epss));

volScalarField epss0("epss0",epss);



// Define relevant parameters for calculating swelling pressure
// Swelling pressure calculation is based on this paper: Liu 2013, Colloids and Surfaces A: Physicochemical and Engineering Aspects, Prediction of swelling pressures of different types of bentonite in dilute solutions
volScalarField h("h",(9.6*1e-10*3/2/0.7667*((((1-epss)/(max(epss,.001)))-0.0333)*(((1-epss)/(max(epss,.001)))-0.0333))+(9.6*1e-10*((1-epss)/(max(epss,.001))))));
volScalarField yhinf("yhinf",4*atanh(-0.8903*exp(-h*3.3143*1e8)));// [10 ppt]
volScalarField yminf("yminf",4*atanh(-0.8903*exp(-h*3.3143*1e8/2)));// [10 ppt]
volScalarField ym("ym",asinh(2*sinh(yminf)+((4/(h*3.3143*1e8))*sin(yhinf/2))));// [10 ppt]
// Swelling pressure contribution from the diffuse double layer, a function of concentration [10 ppt] and temperature
volScalarField PDDL("PDDL",2*298*8.31446*10*(cosh(ym)-1)); 
// Swelling pressure contribution from the crystalline swelling
volScalarField PCSS("PCSS",580000000*exp(-1*h/.0000000002));
// Swelling pressure contribution from the van der Waals
volScalarField PVDW("PVDM",1.1671*1e-21*((1/(h*h*h))-(2/((h+5.04*1e-9)*(h+5.04*1e-9)*(h+5.04*1e-9)))+(1/((h+1.008*1e-8)*(h+1.008*1e-8)*(h+1.008*1e-8)))));
// pswell0 is used for assigning the units to swelling pressure
dimensionedScalar pswell0("pswell0",dimensionSet(1, -1, -2, 0, 0, 0, 0),scalar(1));



// Parameters for compuating pressue and temperature dependent density and viscosity
dimensionedScalar pRef("pRef", dimPressure, 101325.0);   // reference atomospheric pressure
dimensionedScalar rhoMin("rhoMin", dimDensity, SMALL);
dimensionedScalar muMin ("muMin" , dimDynamicViscosity, SMALL);
dimensionedScalar rhoMax("rhoMax", dimDensity, 1e5);            
dimensionedScalar muMax ("muMax" , dimDynamicViscosity, 1e3);


//Info<< "Reading field pswell\n" << endl;
volScalarField pswell
(
    IOobject
    (
        "pswell",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);



// Definition of Solid Indicator
volScalarField Solid
(
    IOobject
    (
        "Solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("Solid", dimensionSet(0,0,0,0,0,0,0), 0) 
);

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());


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


dictionary wetting (transportProperties.subDict(mixture.phase1Name())); 
dictionary nonWetting (transportProperties.subDict(mixture.phase2Name()));

// Define parameters for solid rheology model
const dictionary& rheoDict = transportProperties.subDict("HerschelBulkleyQuemadaCoeffs");
dimensionedScalar tau0("tau0",rheoDict.lookup("tau0"));
dimensionedScalar Dstar("Dstar",rheoDict.lookup("Dstar"));
dimensionedScalar epssMax("epssMax",rheoDict.lookup("epssMax"));

// Define parameters for volume viscosity and swelling correction terms
const dimensionedScalar kappaDyn0("kappaDyn0", dimDynamicViscosity, transportProperties);
const dimensionedScalar mCompaction("mCompaction", dimless, transportProperties);
const dimensionedScalar Ks("Ks", dimPressure, transportProperties);    // bulk modulus of solid

volVectorField SwellTerm
(
    IOobject
    (
        "SwellTerm",                 
        runTime.timeName(),          
        mesh,
        IOobject::NO_READ,          
        IOobject::AUTO_WRITE        
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0,1,-2,0,0,0,0), vector::zero)
);
SwellTerm.write();

volScalarField eps_sw
(
    IOobject
    (
        "eps_sw",                    
        runTime.timeName(),          
        mesh,
        IOobject::NO_READ,           
        IOobject::AUTO_WRITE         
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);



// Define the density and kinematic viscosity for wetting and non-wetting phases
// Coefficients for the density of wetting phase 
dimensionedScalar coeff1("coeff1",transportProperties.lookup("coeff1"));
dimensionedScalar coeff2("coeff2",transportProperties.lookup("coeff2"));
dimensionedScalar coeff3("coeff3",transportProperties.lookup("coeff3"));

// Coefficients for the density of non-wetting phase 
dimensionedScalar coeff4("coeff4",transportProperties.lookup("coeff4"));
dimensionedScalar coeff5("coeff5",transportProperties.lookup("coeff5"));
dimensionedScalar coeff6("coeff6",transportProperties.lookup("coeff6"));

// Coefficients for the dynamic viscosity of wetting phase 
dimensionedScalar coeff7("coeff7",transportProperties.lookup("coeff7"));
dimensionedScalar coeff8("coeff8",transportProperties.lookup("coeff8"));
dimensionedScalar coeff9("coeff9",transportProperties.lookup("coeff9"));

// Coefficients for the dynamic viscosity of non-wetting phase 
dimensionedScalar coeff10("coeff10",transportProperties.lookup("coeff10"));
dimensionedScalar coeff11("coeff11",transportProperties.lookup("coeff11"));
dimensionedScalar coeff12("coeff12",transportProperties.lookup("coeff12"));

// Fluid density and dynamic viscosity at ambient temperatures
dimensionedScalar T0("T0",transportProperties.lookup("T0")); //ambient temperature
dimensionedScalar rho10("rho10",coeff1*T0*T0+coeff2*T0+coeff3);
dimensionedScalar rho20("rho20",coeff4*T0*T0+coeff5*T0+coeff6);
dimensionedScalar mu10("mu10",coeff7*T0*T0+coeff8*T0+coeff9);
dimensionedScalar mu20("mu20",coeff10*T0*T0+coeff11*T0+coeff12);
dimensionedScalar nu10("nu10",mu10/rho10);
dimensionedScalar nu20("nu20",mu20/rho20);
volScalarField mu10v
(
    IOobject
    (
        "mu10v",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    mu10
);
surfaceScalarField mu10f ("mu10f",fvc::interpolate(mu10v,"mu10v"));

volScalarField mu20v
(
    IOobject
    (
        "mu20v",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    mu20
);
surfaceScalarField mu20f ("mu20f",fvc::interpolate(mu20v,"mu20v"));


// Define the temperature dependence coefficient of surface tension 
dimensionedScalar coeffSigma("coeffSigma",transportProperties.lookup("coeffSigma")); 

// Density of wetting phase (fitting results)
volScalarField rho1
(
    IOobject
    (
        "rho1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    coeff1*T*T+coeff2*T+coeff3
);
rho1.oldTime();

// Density of nonwetting phase (fitting results)
volScalarField rho2
(
    IOobject
    (
        "rho2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    coeff4*T*T+coeff5*T+coeff6
);
rho2.oldTime();

// Dynamic viscosity of wetting phase
volScalarField mu1
(
    IOobject
    (
        "mu1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    coeff7*T*T + coeff8*T + coeff9
);
mu1.oldTime();

// Dynamic viscosity of nonwetting phase
volScalarField mu2
(
    IOobject
    (
        "mu2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    coeff10*T*T + coeff11*T + coeff12
);
mu2.oldTime();

// Kinematic viscosity of wetting phase
volScalarField nu1
(
    IOobject
    (
        "nu1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mu1/rho1
);
nu1.oldTime();

// Kinematic viscosity of nonwetting phase
volScalarField nu2
(
    IOobject
    (
        "nu2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mu2/rho2
);
nu2.oldTime();



volScalarField muTotal
(
    IOobject
    (
        "muTotal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1*nu1 + alpha2*rho2*nu2
);
muTotal.oldTime();


volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

volScalarField nuTotal
(
    IOobject
    (
        "nuTotal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    muTotal/(alpha1*rho1 + alpha2*rho2)
);
nuTotal.oldTime();


// Defining Mass Flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


Info<< "Calculating superficial mass flux\n" << endl;
surfaceScalarField rhoPhiByEps
(
    IOobject
    (
        "rhoPhiByEps",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rhoPhi*linearInterpolate(1/(eps+SMALL))
);



#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p.name());
mesh.setFluxRequired(alpha1.name());

// MULES flux from previous time-step
surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    phi*fvc::interpolate(alpha1)
);

// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;


#include "createPorousMediaFields.H"
#include "createSolidMechFields.H"

// Output capillary pressure
volScalarField pcWrite
(
    IOobject
    (
        "pcWrite",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    pc
);
pcWrite.oldTime();

// Output effective permeability
volScalarField k1Write
(
    IOobject
    (
        "k1Write",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,2,0,0,0,0,0), 0.0)
);
k1Write.oldTime();


volScalarField k2Write
(
    IOobject
    (
        "k2Write",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,2,0,0,0,0,0), 0.0)
);
k2Write.oldTime();



Info<< "Creating field kinetic energy Ke\n" << endl;
volScalarField Ke("Ke", 0.5*magSqr(U));  //fluids

#include "createMRF.H"
#include "createFvOptions.H"


