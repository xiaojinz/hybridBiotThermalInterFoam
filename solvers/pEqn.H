    Info<< "\nSolving pEqn.H\n" << endl;

    volScalarField rAU("rAU", 1.0/UEqn.A()); 
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU)); 
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );

    adjustPhi(phiHbyA, U, p);

    // Gravity Contribution
    surfaceScalarField phiG 
    (
       rAUf*fvc::interpolate(rho)*(g & mesh.Sf())
    );

    // Interfacial Force Contributions
    surfaceScalarField phiPc
    (
        (
            mixture.surfaceTensionForce()*(1.-Solidf) //activate free fluid interface model
     +  fvc::interpolate(PcCoeff*dPcdS-pc)*fvc::snGrad(alpha1)*Solidf //activate porous media model
        )*rAUf*mesh.magSf()
    );

    phiHbyA += (phiPc + phiG);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAUf);




    // Add the following three terms to consider compressible fluid conditions
    // Dimensional floor to avoid divide-by-zero in 1/rho
    const dimensionedScalar smallRho( "smallRho", rho1.dimensions(), VSMALL );
    
    // Term 1: time-derivative pieces
    tmp<volScalarField> dtrhoW = fvc::ddt(rho1);
    tmp<volScalarField> dtrhoN = fvc::ddt(rho2);
    tmp<volScalarField> term_dtrho = eps*alpha1/(rho1 + smallRho)*dtrhoW() + eps*alpha2/(rho2 + smallRho)*dtrhoN();

    // Gradient pieces for log-density: (∇ρ)/ρ
    tmp<volVectorField> gradRhoW = fvc::grad(rho1);
    tmp<volVectorField> gradRhoN = fvc::grad(rho2);
    tmp<volVectorField> logGradW = gradRhoW()/(rho1 + smallRho); 
    tmp<volVectorField> logGradN = gradRhoN()/(rho2 + smallRho); 

    // Term 2: advection-like scalar term
    tmp<volScalarField> term_adv = ( U & ( alpha1*logGradW() + alpha2*logGradN() ) );

    // Term 3: relative-velocity scalar term
    tmp<volScalarField> term_rel = ( Ur & ( eps*alpha1*alpha2*( logGradW() - logGradN() ) ) );



    while (pimple.correctNonOrthogonal())
    {
       
        fvScalarMatrix pEqn
        (   
            fvm::laplacian(rAUf, p)
            == fvc::div(phiHbyA)
             + fvc::ddt(eps)
             + term_dtrho()   // added term 1
             + term_adv()     // added term 2
             + term_rel()     // added term 3
        );

        pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();

            p.relax();

            U = HbyA + rAU*fvc::reconstruct(((phiPc+phiG) - pEqn.flux())/rAUf);
            U.correctBoundaryConditions();
        #include "correctUBc.H"
        }
    }

    #include "continuityErrs.H"

    if (p.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }


    Ke = 0.5*magSqr(U);



