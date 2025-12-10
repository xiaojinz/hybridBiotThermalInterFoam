// Defining Permeability Models

Info<< "Reading intrinsic permeability field field K\n" << endl;
volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);



scalar useKozenyCarman(transportProperties.lookupOrDefault<scalar>("useKozenyCarman",1));
scalar activateTemperatureDependenceOfProperties(transportProperties.lookupOrDefault<scalar>("activateTemperatureDependenceOfProperties",1));
scalar activateTemperatureDependenceOfPermeability(transportProperties.lookupOrDefault<scalar>("activateTemperatureDependenceOfPermeability",1));
scalar activateTemperatureDependenceOfCapillarity(transportProperties.lookupOrDefault<scalar>("activateTemperatureDependenceOfCapillarity",1));
scalar activatePressureDependenceOfDensityForWettingPhase(transportProperties.lookupOrDefault<scalar>("activatePressureDependenceOfDensityForWettingPhase",1));
scalar activatePressureDependenceOfDensityForNonWettingPhase(transportProperties.lookupOrDefault<scalar>("activatePressureDependenceOfDensityForNonWettingPhase",1));

volScalarField k0("k0", K);

if (useKozenyCarman == 1)
{k0 = K*(eps*eps*eps)/(max((1-eps)*(1-eps),SMALL)) ;}

surfaceScalarField k0f("k0f",fvc::interpolate(k0,"k0"));

autoPtr<relativePermeabilityModel> krModel = relativePermeabilityModel::New("krModel",transportProperties,alpha1);
krModel->correct(); 
volScalarField kr1 = krModel->krb();
volScalarField kr2 = krModel->kra();

surfaceScalarField kr1f ("kr1f",fvc::interpolate(kr1,"kr1"));
surfaceScalarField kr2f ("kr2f",fvc::interpolate(kr2,"kr2")); //the change I made
surfaceScalarField mu1f ("mu1f",fvc::interpolate(mu1,"mu1"));
surfaceScalarField mu2f ("mu2f",fvc::interpolate(mu2,"mu2"));
surfaceScalarField rho1f ("rho1f",fvc::interpolate(rho1,"rho1"));
surfaceScalarField rho2f ("rho2f",fvc::interpolate(rho2,"rho2"));

// Reading heat conductivity of wetting phase
dimensionedScalar kappa_1 ("kappa_1",transportProperties.lookup("kappa_1"));
// Reading heat capacity of wetting phase
dimensionedScalar Cp_1 ("Cp_1",transportProperties.lookup("Cp_1"));
// Reading heat conductivity of nonwetting phase
dimensionedScalar kappa_2 ("kappa_2",transportProperties.lookup("kappa_2"));
// Reading heat capacity of nonwetting phase
dimensionedScalar Cp_2 ("Cp_2",transportProperties.lookup("Cp_2"));

volScalarField Cp1
(
    IOobject
    (
        "Cp1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    Cp_1
);

volScalarField Cp2
(
    IOobject
    (
        "Cp2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    Cp_2
);

volScalarField kappa1
(
    IOobject
    (
        "kappa1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    kappa_1
);

volScalarField kappa2
(
    IOobject
    (
        "kappa2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    kappa_2
);

volScalarField Cpf
(
    IOobject
    (
        "Cpf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha1*Cp1 + alpha2*Cp2
);
Cpf.oldTime();

surfaceScalarField Cp1f ("Cp1f",fvc::interpolate(Cp1,"Cp1"));
surfaceScalarField Cp2f ("Cp2f",fvc::interpolate(Cp2,"Cp2"));
surfaceScalarField kappa1f ("kappa1f",fvc::interpolate(kappa1,"kappa1"));
surfaceScalarField kappa2f ("kappa2f",fvc::interpolate(kappa2,"kappa2"));
surfaceScalarField Cpff ("Cpff",fvc::interpolate(Cpf,"Cpf"));




// Mobilities
volScalarField M1 ("M1",k0*kr1/mu1);
volScalarField L1 ("L1",rho1*k0*kr1/mu1);	
volScalarField M2 ("M2",k0*kr2/mu2);
volScalarField L2 ("L2",rho2*k0*kr2/mu2);	
volScalarField M ("M",M1+M2);
volScalarField L ("L",L1+L2);

surfaceScalarField M1f ("M1f",k0f*kr1f/mu1f);
surfaceScalarField L1f ("L1f",rho1f*k0f*kr1f/mu1f);	
surfaceScalarField M2f ("M2f",k0f*kr2f/mu2f);
surfaceScalarField L2f ("L2f",rho2f*k0f*kr2f/mu2f);	
surfaceScalarField Mf ("Mf",M1f+M2f);
surfaceScalarField Lf ("Lf",L1f+L2f);

// Drag Coefficient Calculation
surfaceScalarField Dragf ("Dragf",1/Mf);
volScalarField Drag ("Drag", 1/M);

// Ferrari 2015 2.5 D drag coefficient
dimensionedScalar invFk0("Fk0", transportProperties.lookup("invFk0",0.));
volScalarField DragF("DragF", (1-Solid)*invFk0*rho*nuTotal);

// Defining relative permeability and capillary pressure models
autoPtr<capillarityModel> myCapModel = capillarityModel::New("pc",transportProperties,alpha1);
scalar activateCapillarity(transportProperties.lookupOrDefault<scalar>("activateCapillarity",0.));
myCapModel().correct();
volScalarField pc("pc",myCapModel().pc()*activateCapillarity);
volScalarField dPcdS("dPcdS",myCapModel().dpcdS()*activateCapillarity);
volScalarField PcCoeff("PcCoeff", (1./M)*(M1*(1.-alpha1) - M2*alpha1));

Info<< "Initializing Relative Velocity Field Ur\n" << endl;
volVectorField Ur
(
    IOobject
    (
        "Ur",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U*0
);



// Labeling the places that have porous media
forAll(epss,celli) 
{
	if(epss[celli]>0.01) 
	{	
	Solid[celli]=1; 
	}
	else
	{
	Solid[celli]=0; 
	}
}

forAll(mesh.boundary(), patchi) 
{
	forAll(Solid.boundaryField()[patchi],facei)
	{
		if(epss.boundaryField()[patchi][facei]>0.01) 
		{Solid.boundaryFieldRef()[patchi][facei]=1;}
		else
		{Solid.boundaryFieldRef()[patchi][facei]=0;}
	}
}

surfaceScalarField Solidf(fvc::interpolate(Solid));

// Eliminating intermediate values of Solid Indicator Function
forAll(Solidf,facei)
{
	if(Solidf[facei]<1) // if <1 ==0.  
	{	
	Solidf[facei]=0; 
	}
}




