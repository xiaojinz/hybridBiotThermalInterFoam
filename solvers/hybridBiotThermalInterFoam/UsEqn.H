Info<< "\nCalculating Us" << endl;

        // Effective pressure
        volScalarField p_eff = p; //not using pswell for now
        //volScalarField p_eff = p + pswell;

        // Bulk viscosity                           
        volScalarField epssStar("epssStar", epss);
        epssStar = Foam::min(epssStar, epssMax - 0.01);
        epssStar = Foam::max(epssStar, 0.01);
        kappa_dyn = kappa_dyn0 * pow(max(SMALL, (1-epssStar/epssMax) ), -mCompaction);
        kappa_dyn.correctBoundaryConditions();
	      kappa_kin = (kappa_dyn / rho_s);

        // Swelling strain rate 
        //eps_sw = (1.0 / 3.0 /Ks) * fvc::ddt(pswell);
        eps_sw = (1.0 / 3.0 / Ks) * (alpha1 * fvc::ddt(pswell) + pswell * fvc::ddt(alpha1));
       
       	// Swelling-induced volumetric strain term (swelling corrections)
	      SwellTerm = (3.0/rho_s) * fvc::grad(epss * kappa_dyn * eps_sw);

        // Terzaghi term 
        volVectorField Terzaghi = -(epss/rho_s)*( fvc::grad(p_eff) );

        //Calculating Explicit Capillary Interaction force
        volVectorField CapTerm ((Solid/rho_s)*
                         (
                           (eps*PcCoeff*dPcdS - eps*pc)*fvc::grad(alpha1)
                         +
                            eps*PcCoeff*(rho1-rho2)*g
                         )
                       );

fvVectorMatrix UsEqn 
        (

        //Inertial Terms
        fvm::ddt(epss,Us)
      + fvm::div(phisEpss,Us)

       //Shear viscous terms 
      - fvm::laplacian(solid.nu(),Us)
      - fvc::div(solid.nu()*dev2(Foam::T(fvc::grad(Us))))

       // Volumetric viscous term with swelling corrections
      - fvc::grad(kappa_kin * fvc::div(Us))  // Volumetric viscous term (here epss is neglected because of using Darcy like solid velocity)
      + SwellTerm                      // Swelling-induced volumetric strain term (swelling corrections)

        //Terzaghi Term
      - Terzaghi

        //Gravity Term
      - epss*g

        //Capillary presure source
      + CapTerm

        //Momentum Source from fluid to solid
      - fvc::Sp(eps*Drag/rho_s, U)

        //Momentum Transfer from solid to fluid 
      + fvm::Sp(epss*eps*Drag/rho_s, Us)
);

UsEqn.solve();
Us.correctBoundaryConditions();

Ke_s = 0.5*magSqr(Us);

phis = fvc::interpolate(Us) & mesh.Sf();
