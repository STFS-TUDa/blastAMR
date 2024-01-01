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

Description
    Refine and unrefine cells based on the magnitude of the gradient of
    a target scalar field and threshold values for gradient range percentages.
        |------------=======xxxxxxxxxxxxxxxxxxxxx|
    minGrad       0.4    0.5                  maxGrad
     -: unrefine, =: do nothing, x: refine

SourceFiles
    gradientRange.C

\*---------------------------------------------------------------------------*/

#ifndef errorEstimatorTemplate_H
#define errorEstimatorTemplate_H

#include "errorEstimator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                         A templated functionObject
\*---------------------------------------------------------------------------*/

namespace errorEstimators
{

class gradientEstimator
:
    public errorEstimator
{

    //- Target field
    word gradField_;

    //- Gradient range to refine > threshold
    scalar refineThreshold_;

    //- Gradient range to unrefine > threshold
    scalar unrefineThreshold_;

public:

    //- Runtime type information
    TypeName("gradientRange");


    // Constructors

        //- Construct from Time and dictionary
        gradientEstimator
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& name
        );


    //- Destructor
    virtual ~gradientEstimator();


    // Member Functions

        //- Read the dictionary
        virtual void read(const dictionary& dict);

        //- Update
        virtual void update(const bool scale);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

