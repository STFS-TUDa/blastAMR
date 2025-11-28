/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "cellSize.H"
#include "meshSizeObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(cellSize, 0);
    addToRunTimeSelectionTable(errorEstimator, cellSize, dictionary);
}
}

template<>
const char* Foam::NamedEnum<Foam::errorEstimators::cellSize::SizeType, 4>::names[] =
{
    "volume",
    "cmpt",
    "characteristic",
    "mag"
};

const Foam::NamedEnum<Foam::errorEstimators::cellSize::SizeType, 4>
Foam::errorEstimators::cellSize::sizeTypeNames;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::cellSize::cellSize
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    errorEstimator(mesh, dict, name),
    sizeType_(CHARACTERISTIC),
    cmpts_()
{
    this->read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::cellSize::~cellSize()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::errorEstimators::cellSize::readCmpts(Istream& is) const
{
    labelHashSet cmpts;
    while (is.good())
    {
        token t(is);

        if (t.isNumber())
        {
            cmpts.insert(t.labelToken());
        }
        else if (t.isWord())
        {
            word c(t.wordToken());
            if (c == "x" || c == "X")
            {
                cmpts.insert(0);
            }
            else if (c == "y" || c == "Y")
            {
                cmpts.insert(1);
            }
            else if (c == "z" || c == "Z")
            {
                cmpts.insert(2);
            }
            else
            {
                FatalIOErrorInFunction(is)
                    << "Unknown component " << c << endl
                    << "Valid components are: " << nl
                    << "x/X" << nl
                    << "y/Y" << nl
                    << "z/Z" << endl
                    << abort(FatalIOError);
            }
        }
    }
    const Vector<label> geoD(mesh_.geometricD());
    labelList cmptLables;
    forAllConstIter(labelHashSet, cmpts, iter)
    {
        if (geoD[iter.key()] == 1)
        {
            cmptLables.append(iter.key());
        }
    }
    return cmptLables;
}


Foam::labelList Foam::errorEstimators::cellSize::maxRefinement() const
{
    return mesh_.lookupObject<labelIOList>("cellLevel") + 1;
}

void Foam::errorEstimators::cellSize::computeDx() {
    if (sizeType_ == VOLUME)
    {
        this->error_.internalFieldRef().field() = mesh_.V().field();
    }
    else if (sizeType_ == CHARACTERISTIC)
    {
        this->error_.internalFieldRef().field() = meshSizeObject::New(mesh_).dx();
    }
    else if (sizeType_ == MAG)
    {
        const vectorField& dX = meshSizeObject::New(mesh_).dX();
        const Vector<label> geoD(mesh_.geometricD());
        this->error_.internalFieldRef() = 0.0;
        forAll(geoD, cmpti)
        {
            if (geoD[cmpti] == 1)
            {
                this->error_.internalFieldRef().field() += sqr(dX.component(cmpti));
            }
        }
        this->error_.internalFieldRef() = sqrt(this->error_.internalFieldRef());
    }
    else if (sizeType_ == CMPT)
    {
        const vectorField& dX = meshSizeObject::New(mesh_).dX();
        forAll(cmpts_, i)
        {
            label cmpti = cmpts_[i];
            forAll(this->error_.internalFieldRef(), celli)
            {
                if (dX[celli][cmpti] > maxDX_[cmpti])
                {
                    this->error_.internalFieldRef()[celli] = max(1.0, this->error_.internalFieldRef()[celli]);
                }
                else if (dX[celli][cmpti] < minDX_[cmpti])
                {
                    this->error_.internalFieldRef()[celli] = max(-1.0, this->error_.internalFieldRef()[celli]);
                }
                else
                {
                    this->error_.internalFieldRef()[celli] = max(0, this->error_.internalFieldRef()[celli]);
                }
            }
        }
    }
}


void Foam::errorEstimators::cellSize::update(const bool scale)
{
    if (updateCurTimeIndex(!scale))
    {
        return;
    }

    computeDx();

    if (sizeType_ != CMPT)
    {
        forAll(this->error_.internalFieldRef(), celli)
        {
            if (this->error_.internalFieldRef()[celli] < lowerUnrefine_)
            {
                this->error_.internalFieldRef()[celli] = -1.0;
            }
            else if (this->error_.internalFieldRef()[celli] > lowerRefine_)
            {
                this->error_.internalFieldRef()[celli] = 1.0;
            }
            else
            {
                this->error_.internalFieldRef()[celli] = 0.0;
            }
        }
    }
}


void Foam::errorEstimators::cellSize::read(const dictionary& dict)
{
    sizeType_ = sizeTypeNames.read(dict.lookup("type"));
    if (sizeType_ == CMPT)
    {
        cmpts_ = readCmpts(dict.lookup("cmpts"));
        minDX_ = vector(dict.lookup("minDX"));
        maxDX_ = vector(dict.lookup("maxDX"));
    }
    else if (sizeType_ == VOLUME)
    {
        lowerUnrefine_ = readScalar(dict.lookup("minVolume"));
        lowerRefine_ = readScalar(dict.lookup("maxVolume"));
    }
    else
    {
        lowerUnrefine_ = readScalar(dict.lookup("minDx"));
        lowerRefine_ = readScalar(dict.lookup("maxDx"));
    }

    maxLevel_ = dict.lookupOrDefault<label>("maxRefinement", 10);
}

// ************************************************************************* //
