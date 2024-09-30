scalar DiNum = -great;

forAll(solidRegions, i)
{
    // Note: do not use setRegionSolidFields.H to avoid double registering Cp
    //#include "setRegionSolidFields.H"
    //const solidThermo& thermo = thermoSolid[i];
    //const fluidThermo& thermo = thermoSolid[i];
    
    volScalarField magKappa
    (
        IOobject
        (
            "magKappa",
            runTime.timeName(),
            solidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidRegions[i],
        dimensionedScalar(dimensionSet(1, 1, -3, -1, 0, 0, 0),kappa0.value())
    );

    volScalarField cp
    (
        IOobject
        (
            "cp",
            runTime.timeName(),
            solidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidRegions[i],
        dimensionedScalar(dimensionSet(0, 2, -2, -1, 0, 0, 0),Cp0.value())
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            solidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidRegions[i],
        dimensionedScalar(dimensionSet(1, -3, 0, 0, 0, 0, 0),rho0.value())
    );

    DiNum = max
    (
        solidRegionDiffNo
        (
            solidRegions[i],
            runTime,
            rho*cp,
            magKappa
        ),
        DiNum
    );

}
