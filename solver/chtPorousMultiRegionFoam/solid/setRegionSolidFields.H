const fvMesh& mesh = solidRegions[i];

volVectorField& U = USolid[i];

volScalarField& T = TSolid[i];

volScalarField& p = pSolid[i];

volScalarField& p2 = p2Solid[i];

volScalarField& rho = rhoSolid[i];

#include "createPhi.H"

const Foam::fvModels& fvModels = fvModelsSolid[i];
Foam::fvConstraints& fvConstraints = fvConstraintsSolid[i];

#include "checkRadiationModel.H"

solidNoLoopControl& pimple = pimples.solid(i);
