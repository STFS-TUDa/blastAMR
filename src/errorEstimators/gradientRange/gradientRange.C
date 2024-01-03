/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gradientRange.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradientEstimator, 0);

addRemovableToRunTimeSelectionTable
(
    errorEstimator,
    gradientEstimator,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gradientEstimator::gradientEstimator
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    errorEstimator(mesh, dict, name),
    gradField_(dict.lookup("gradField")),
    refineThreshold_(dict.lookupOrDefault("refineThreshold", 0.5)),
    unrefineThreshold_(dict.lookupOrDefault("unrefineThreshold", 0.4))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gradientEstimator::~gradientEstimator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gradientEstimator::read(const dictionary& dict)
{
    errorEstimator::read(dict);
    gradField_= word(dict.lookup("gradField"));
    refineThreshold_ = dict.lookupOrDefault("refineThreshold", 0.5);
    unrefineThreshold_ = dict.lookupOrDefault("unrefineThreshold", 0.4);
}


void gradientEstimator::update(const bool scale)
{
    if (updateCurTimeIndex(!scale))
    {
        return;
    }

    error_ = 0.0;
    const auto& vf = mesh_.lookupObject<volScalarField>(gradField_);
    const scalar refineThreshold(dict().lookupOrDefault("refineThreshold", 0.5));
    const scalar unrefineThreshold(dict().lookupOrDefault("unrefineThreshold", 0.4));
    Info << "refineThreshold = " << refineThreshold << endl;
    Info << "unrefineThreshold = " << unrefineThreshold << endl;
    error_ == mag(fvc::grad(vf)) / dimensionedScalar(vf.dimensions()/dimLength,1.0);

    scalar maxGrad = gMax(error_);
    scalar minGrad = gMin(error_);

    lowerRefine_ = minGrad + refineThreshold*(maxGrad-minGrad);
    upperRefine_ =  GREAT;
    lowerUnrefine_ = minGrad + unrefineThreshold*(maxGrad-minGrad);
    upperUnrefine_ =  GREAT;
    if (scale) normalize(error_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// ************************************************************************* //

