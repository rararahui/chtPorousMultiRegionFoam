// Initialise solid field pointer lists
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size());
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size());
PtrList<volVectorField> USolid(solidRegions.size());
PtrList<volScalarField> TSolid(solidRegions.size());
PtrList<volScalarField> pSolid(solidRegions.size());
PtrList<volScalarField> p2Solid(solidRegions.size());
PtrList<volScalarField> rhoSolid(solidRegions.size());

// Populate solid field pointer lists

Info<< "Reading physics properties\n" << endl;

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        solidRegions[0],
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar k0
(
    "k0",
    dimensionSet(0, 2, 0, 0, 0, 0, 0),
    physicalProperties.lookup("k0")
);

dimensionedScalar Ef
(
    "Ef",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    physicalProperties.lookup("Ef")
);

dimensionedScalar p0
(
    "p0",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    physicalProperties.lookup("p0")
);

dimensionedScalar nu
(
    "nu",
    dimensionSet(1, -1, -1, 0, 0, 0, 0),
    physicalProperties.lookup("nu")
);

dimensionedScalar kappa0
(
    "kappa0",
    dimensionSet(1, 1, -3, -1, 0, 0, 0),
    physicalProperties.lookup("kappa0")
);

dimensionedScalar rho0
(
    "rho0",
    dimensionSet(1, -3, 0, 0, 0, 0, 0),
    physicalProperties.lookup("rho0")
);

dimensionedScalar Cp0
(
    "Cp0",
    dimensionSet(0, 2, -2, -1, 0, 0, 0),
    physicalProperties.lookup("Cp0")
);

dimensionedScalar kappa1
(
    "kappa1",
    dimensionSet(1, 1, -3, -1, 0, 0, 0),
    physicalProperties.lookup("kappa1")
);

dimensionedScalar rho1
(
    "rho1",
    dimensionSet(1, -3, 0, 0, 0, 0, 0),
    physicalProperties.lookup("rho1")
);

dimensionedScalar Cp1
(
    "Cp1",
    dimensionSet(0, 2, -2, -1, 0, 0, 0),
    physicalProperties.lookup("Cp1")
);

dimensionedScalar molWeight
(
    "molWeight",
    dimensionSet(1, 0, 0, 0, -1, 0, 0),
    physicalProperties.lookup("molWeight")
);

dimensionedScalar DT = Ef*kappa0/(rho0*Cp0) + (1-Ef)*kappa1/(rho1*Cp1);

dimensionedScalar RR(dimensionSet(1,2,-2,-1,-1,0,0),1e3*8.314);
        
dimensionedScalar R=RR/molWeight;

forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding to USolid\n" << endl;
    USolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );
    
    Info<< "    Adding to pSolid\n" << endl;
    pSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );
    
    Info<< "    Adding to TSolid\n" << endl;
    TSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );
    
    Info<< "    Adding to p2Solid\n" << endl;
    p2Solid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p2",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );
    
    Info<< "    Adding to rhoSolid\n" << endl;
    rhoSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pSolid[i]/(R*TSolid[i])
        )
    );
        
    Info<< "    Adding to fvModelsSolid\n" << endl;
    fvModelsSolid.set
    (
        i,
        new Foam::fvModels(solidRegions[i])
    );

    Info<< "    Adding fvConstraintsSolid\n" << endl;
    fvConstraintsSolid.set
    (
        i,
        new fvConstraints(solidRegions[i])
    );
}
