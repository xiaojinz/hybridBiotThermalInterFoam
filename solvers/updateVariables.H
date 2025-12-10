Info<< "\nUpdating variables\n" << endl; 

// Update swelling pressure 
// Swelling pressure calculation is based on this paper: Liu 2013, Colloids and Surfaces A: Physicochemical and Engineering Aspects, Prediction of swelling pressures of different types of bentonite in dilute solutions
if (activateSwellingPressure == 1)
{
    h=(9.6*1e-10*3/2/0.7667*((((1-epss)/(max(epss,.001)))-0.0333)*(((1-epss)/(max(epss,.001)))-0.0333))+(9.6*1e-10*((1-epss)/(max(epss,.001)))));
    yhinf=4*atanh(-0.8903*exp(-h*3.3143*1e8));// [10 ppt]
    yminf=4*atanh(-0.8903*exp(-h*3.3143*1e8/2));// [10 ppt]
    ym=asinh(2*sinh(yminf)+((4/(h*3.3143*1e8))*sin(yhinf/2)));// [10 ppt]

    PDDL=(2*298*8.31446*10*(cosh(ym)-1)); // [10 ppt]
    PCSS=580000000*exp(-1*h/.0000000002);
    PVDW=1.1671*1e-21*((1/(h*h*h))-(2/((h+5.04*1e-9)*(h+5.04*1e-9)*(h+5.04*1e-9)))+(1/((h+1.008*1e-8)*(h+1.008*1e-8)*(h+1.008*1e-8))));
}
pswell=pswell0*(PCSS+PDDL-PVDW)*activateSwellingPressure;


// Update fluid density and viscosity as a function of temperature 
if (activateTemperatureDependenceOfProperties == 1)
{
    // Temperature-dependent polynomial evaluations
    volScalarField rho1_T = coeff1*T*T + coeff2*T + coeff3;  
    volScalarField mu1_T  = coeff7*T*T + coeff8*T + coeff9;

    volScalarField rho2_T = coeff4*T*T + coeff5*T + coeff6;  
    volScalarField mu2_T  = coeff10*T*T + coeff11*T + coeff12;

    // Clamp to avoid negative densities or viscosities
     rho1 = max(min(rho1_T, rhoMax), rhoMin);
    rho2 = max(min(rho2_T, rhoMax), rhoMin);
    mu1  = max(min(mu1_T , muMax),  muMin);
    mu2  = max(min(mu2_T , muMax),  muMin);

/*
    // Pressure dependence for compressible gas
    if (activatePressureDependenceOfDensityForWettingPhase == 1)
    {
        rho1 = max(min(rho1_T * (p/pRef), rhoMax), rhoMin);  
    }

    if (activatePressureDependenceOfDensityForNonWettingPhase == 1)
    {
        rho2 = max(min(rho2_T * (p/pRef), rhoMax), rhoMin);
    }
*/
}


rho1f = fvc::interpolate(rho1,"rho1");
mu1f = fvc::interpolate(mu1,"mu1");
nu1 = mu1 / max(rho1, rhoMin); 

rho2f = fvc::interpolate(rho2,"rho2");
mu2f = fvc::interpolate(mu2,"mu2");
nu2 = mu2 / max(rho2, rhoMin); 
  
muTotal = alpha1*rho1*nu1 + alpha2*rho2*nu2;
nuTotal = muTotal/(alpha1*rho1 + alpha2*rho2);


// Update relative permeability model
if (useKozenyCarman == 1)
{k0 = K*(1-epss)*(1-epss)*(1-epss)/max(epss*epss,SMALL);}

krModel->correct();
kr1=krModel->krb();
kr2=krModel->kra();
kr1f = fvc::interpolate(kr1,"kr1");
kr2f = fvc::interpolate(kr2,"kr2");

// Update permeability as a function of temperature
if (activateTemperatureDependenceOfPermeability == 1)
{
kr1 = (krModel->krb())*( mu1/mu10);
kr2 = (krModel->kra())*( mu2/mu20);
kr1f = (fvc::interpolate(kr1,"kr1"))*( mu1f/mu10f);
kr2f = (fvc::interpolate(kr2,"kr2"))*( mu2f/mu20f);
}

// Output permeability
k1Write = k0*kr1;
k2Write = k0*kr2;

// Update mobilities
M1 = k0*kr1/mu1;
L1 = rho1*k0*kr1/mu1;
M1f = k0f*kr1f/mu1f;
L1f = rho1f*k0f*kr1f/mu1f;

M2 = k0*kr2/mu2;
L2 = rho2*k0*kr2/mu2;
M2f = k0f*kr2f/mu2f;
L2f = rho2f*k0f*kr2f/mu2f;

M = M1+M2;
L = L1+L2;
Mf = M1f+M2f;
Lf = L1f+L2f;


// Calculate Drag Coefficient
Dragf =1./Mf;
Drag= fvc::average(Dragf);

forAll(Drag,celli)
{
	if (Solid[celli] == 0)
		{Drag[celli]= 0;}
}

// Calculate drag coefficient (Ferrari 2015, 2.5 D term)
rho = (alpha1*rho1 + alpha2*rho2)*(1.-Solid) + Solid*L/M;
DragF = (1-Solid)*invFk0*rho*nuTotal; //I made changes

// Update capillary pressure model
myCapModel().correct();
pc = myCapModel().pc()*activateCapillarity; 
dPcdS = myCapModel().dpcdS()*activateCapillarity; 

if (activateCapillarity == 1)
{
// Update the temperature dependence of capillary pressure
if (activateTemperatureDependenceOfCapillarity == 1)
{
   pc = pc*exp(coeffSigma*(T-T0)) ;
   dPcdS = dPcdS*exp(coeffSigma*(T-T0)) ; 
}
}

PcCoeff = (1./M)*(M1*(1.-alpha1) - M2*(alpha1));

// Output capillary pressure
pcWrite = pc;


