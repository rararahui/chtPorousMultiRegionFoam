{
    while (pimple.correctNonOrthogonal())
    {        
        p2=p*p;
        
	p2.correctBoundaryConditions();
	
	const volScalarField p_temp=p;
        
        fvScalarMatrix pEqn
        (
            fvm::ddt(p2) + fvm::laplacian(-(k0*p_temp)/(nu*Ef),p2)
        );  

        pEqn.solve();
        
        p=sqrt(p2);

        p.correctBoundaryConditions();

        U = -k0/nu*fvc::grad(p);

        U.correctBoundaryConditions();

        #include "TEqn.H"
        
        rho=p/(R*T);
        
        rho.correctBoundaryConditions();    
    }
}

//thermo.correct();

//Info<< "Min/max T:" << min(thermo.T()).value() << ' '
//    << max(thermo.T()).value() << endl;
