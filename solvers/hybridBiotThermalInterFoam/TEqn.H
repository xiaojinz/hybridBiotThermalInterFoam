{

    Info<< "\nSolving temperature equation\n" << endl;

    //-----------------------------------------------------------------------
    // Define thermal properties with explicit dimensions
    //-----------------------------------------------------------------------

    volScalarField alpha1Clipped("alpha1Clipped", min(max(alpha1, scalar(0.01)), scalar(0.99)));
    volScalarField alpha2Clipped("alpha2Clipped", min(max(alpha2, scalar(0.01)), scalar(0.99)));

    // Fluid thermal capacity (rho*Cp) with porosity
    volScalarField rhoCp_f(
        "rhoCp_f",
        (alpha1Clipped*rho1*Cp1 + alpha2Clipped*rho2*Cp2)
    );

    // Solid thermal capacity with porosity
    volScalarField rhoCp_s(
        "rhoCp_s",
        rhos*Cps
    );

    // Total thermal inertia (time derivative coefficient)
    volScalarField rhoCpTotal(
        "rhoCpTotal",
        eps*rhoCp_f + epss*rhoCp_s
    );

    // Thermal conductivity
    volScalarField kappaTotal(
        "kappaTotal",
        eps*(alpha1Clipped*kappa1 + alpha2Clipped*kappa2) + epss*kappas
	);

    



   // ========terms associated with mechanical energy============
    
    // --- Squared magnitudes (scalars) ---
    volScalarField magSqUf("magSqUf", magSqr(U));
    volScalarField magSqUr("magSqUr", magSqr(Ur));
    volScalarField magSqUs("magSqUs", magSqr(Us));
    volScalarField eps2("eps2", sqr(eps));
    volScalarField epss2("epss2", sqr(epss));
    volScalarField eps3("eps3", eps2*eps); 

    // --- Cell-centered energy-like scalars (time derivative terms) ---
    volScalarField meCoeff1(
        "meCoeff1",
            eps * epss2 * rho * magSqUf
          + 2.0 * eps2 * epss2 * alpha1 * alpha2 * (rho1 - rho2) * (Ur & U)
          + eps3 * epss2 * alpha1 * alpha2 * (rho1*alpha2 + rho2*alpha1) * magSqUr
    );

    volScalarField meCoeff2(
        "meCoeff2",
        eps2 * epss * rhos * magSqUs
    );

    // --- Cell-centered cubic flux coefficients ---
    volVectorField meCoeff3flux(
        "meCoeff3flux",
            epss2 * rho * magSqUf * U
          + 3.0 * eps * epss2 * alpha1 * alpha2 * (rho1 - rho2) * magSqUf * Ur
          + 3.0 * eps2 * epss2 * alpha1 * alpha2 * rho * magSqUr * U
          + eps3 * epss2 * alpha1 * alpha2 * (rho1 * sqr(alpha2) - rho2 * sqr(alpha1)) * magSqUr * Ur
    );

    // --- Face-centered cubic fluxes for divergence operator ---
    surfaceScalarField meCoeff3(
        "meCoeff3",
        fvc::interpolate(meCoeff3flux) & mesh.Sf()
    );

    volVectorField meCoeff4flux(
        "meCoeff4flux",
        eps2 * rhos * magSqUs * Us
    );

    surfaceScalarField meCoeff4(
        "meCoeff4",
        fvc::interpolate(meCoeff4flux) & mesh.Sf()
    );

    volVectorField meCoeff5flux(
        "meCoeff5flux",
        2.0 * eps2 * epss2 * p * U
    );

    surfaceScalarField meCoeff5(
        "meCoeff5",
        fvc::interpolate(meCoeff5flux) & mesh.Sf()
    );

    volVectorField meCoeff6flux(
        "meCoeff6flux",
        2.0 * eps2 * epss2 * p * Us
    );

    surfaceScalarField meCoeff6(
        "meCoeff6",
        fvc::interpolate(meCoeff6flux) & mesh.Sf()
    );

    // ======== End of mechanical-energy terms ==========




    // --- Print Ur statistics ---
    volScalarField magUr("magUr", mag(Ur));
    
    // === Mask to restrict to air-water interface ===
    // Define sigmoid parameters
    scalar steepness = 40;              // steepness
    scalar centerThreshold = 0.1;       // center threshold 
    scalar alphaMixCutoff = 0.01;       //cutoff (0.01*0.99=0.0099)

    volScalarField interfaceMask(
        IOobject("interfaceMask", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        pos(alpha1Clipped * alpha2Clipped - scalar(0.01))  // interface zones
       // Solid * pos(alpha1Clipped * alpha2Clipped - scalar(0.01))  // interface zones //turn this off in solid-free regions
       	//  1.0 / (1.0 + exp(-steepness * (alpha1Clipped * alpha2Clipped - centerThreshold)))
    );


    // --- Coefficient for phiUr (cell-centered) ---
    volScalarField phiUrCoeff
    (
        IOobject("phiUrCoeff", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
       alpha1Clipped * alpha2Clipped * (rho1*Cp1 - rho2*Cp2) * interfaceMask 
    );


    // === Safe Ur computation using clamped alpha1 and alpha2 ===
    volVectorField UrSafe
    (
        IOobject("UrSafe", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
      (Uwetting / alpha1Clipped - UnonWetting / alpha2Clipped) * interfaceMask  //here eps is cancelled out by the eps in phiUrCoeff
    );

    volVectorField gradT = fvc::grad(T);
    volScalarField phiUrSource
    (
        IOobject("phiUrSource", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        phiUrCoeff * (gradT & UrSafe)
    );
    phiUrSource *= interfaceMask;

    
    //-----------------------------------------------------------------------
    // Define explicit heat source (0 W/m³)
    //-----------------------------------------------------------------------
    volScalarField heatSource(
        IOobject(
            "heatSource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(
            "heatSource",
            dimPower/dimVolume,  // [W/m³]
            0  // Adjust this value as needed
        )
    );




    //-----------------------------------------------------------------------
    // Assemble equation
    //-----------------------------------------------------------------------
    fvScalarMatrix TEqn(
        // Time derivative 
       fvm::ddt(rhoCpTotal, T)

        // Fluid convection 
      //  + fvm::div(phiRhoCpFluid, T)
      + rhoCp_f*fvm::div(phi, T)

        // Solid convection 
      //+ fvm::div(phiRhoCpSolid, T)
      //+ rhoCp_s*fvm::div(phis, T)

        // Diffusion 
      - fvm::laplacian(kappaTotal, T)

        // terms associated with mechanical energy 
      + fvc::ddt(meCoeff1)
      + fvc::ddt(meCoeff2)
      + fvc::div(meCoeff3)
      + fvc::div(meCoeff4)
      + fvc::div(meCoeff5)
      + fvc::div(meCoeff6)

        // Explicit source term 
      ==
        heatSource 

        // Treat Ur-induced transfer as a volumetric source
      - phiUrSource

      	);

    
    //-----------------------------------------------------------------------
    // Solve equation
    //-----------------------------------------------------------------------
    TEqn.relax();
    TEqn.solve();
}

