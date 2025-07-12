/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "hexRef.H"

#include "polyMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::hexRef> Foam::hexRef::New
(
    const polyMesh& mesh,
    const bool readHistory
)
{
    word hexRefTypeName;

    // Infer the type of hexRef we need to use from the number of dimensions in
    // the polymesh
    // nSolutionD == 3 && nGeometricD == 3: 3D mesh => hexRef8
    // nSolutionD == 3 && nGeometricD == 2: axisymmetric mesh => hexRef4Axi
    // nSolutionD == 2 && nGeometricD == 2: 2D mesh => hexRef4

    label nSoluD(mesh.nSolutionD());
    label nGeomD(mesh.nGeometricD());

    if (nSoluD == 3 && nGeomD == 3)
    {
        hexRefTypeName = "hexRef3D";
    }
    else if (nSoluD == 3 && nGeomD == 2)
    {
        hexRefTypeName = "hexRef2DAxi";
    }
    else if (nSoluD == 2 && nGeomD == 2)
    {
        hexRefTypeName = "hexRef2D";
    }

    auto* ctorPtr = meshConstructorTable(hexRefTypeName);

    if (ctorPtr == nullptr)
    {
        FatalErrorInFunction
            << "Unsupported mesh number of dimensions for hex refinement" << nl
            << "nSolutionD: " << nSoluD << ", nGeometricD: " << nGeomD << nl
            << "Only 3D, 2D and 2D axisymmetric mesh refinements are supported"
            << exit(FatalError);
    }

    return autoPtr<hexRef>
    (
        ctorPtr(mesh, readHistory)
    );
}

Foam::autoPtr<Foam::hexRef> Foam::hexRef::New
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const hexRefRefinementHistory& history,
    const scalar level0Edge
)
{
    word hexRefTypeName;

    // Infer the type of hexRef we need to use from the number of dimensions in
    // the polymesh
    // nSolutionD == 3 && nGeometricD == 3: 3D mesh => hexRef8
    // nSolutionD == 3 && nGeometricD == 2: axisymmetric mesh => hexRef4Axi
    // nSolutionD == 2 && nGeometricD == 2: 2D mesh => hexRef4

    label nSoluD(mesh.nSolutionD());
    label nGeomD(mesh.nGeometricD());

    if (nSoluD == 3 && nGeomD == 3)
    {
        hexRefTypeName = "hexRef8";
    }
    else if (nSoluD == 3 && nGeomD == 2)
    {
        hexRefTypeName = "hexRef4Axi";
    }
    else if (nSoluD == 2 && nGeomD == 2)
    {
        hexRefTypeName = "hexRef4";
    }

    auto* ctorPtr = levelsHistConstructorTable(hexRefTypeName);

    if (ctorPtr == nullptr)
    {
        FatalErrorInFunction
            << "Unsupported mesh number of dimensions for hex refinement" << nl
            << "nSolutionD: " << nSoluD << ", nGeometricD: " << nGeomD << nl
            << "Only 3D, 2D and 2D axisymmetric mesh refinements are supported"
            << exit(FatalError);
    }

    return autoPtr<hexRef>
    (
        ctorPtr(mesh, cellLevel, pointLevel, history, level0Edge)
    );
}

Foam::autoPtr<Foam::hexRef> Foam::hexRef::New
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar level0Edge
)
{
    word hexRefTypeName;

    // Infer the type of hexRef we need to use from the number of dimensions in
    // the polymesh
    // nSolutionD == 3 && nGeometricD == 3: 3D mesh => hexRef8
    // nSolutionD == 3 && nGeometricD == 2: axisymmetric mesh => hexRef4Axi
    // nSolutionD == 2 && nGeometricD == 2: 2D mesh => hexRef4

    label nSoluD(mesh.nSolutionD());
    label nGeomD(mesh.nGeometricD());

    if (nSoluD == 3 && nGeomD == 3)
    {
        hexRefTypeName = "hexRef8";
    }
    else if (nSoluD == 3 && nGeomD == 2)
    {
        hexRefTypeName = "hexRef4Axi";
    }
    else if (nSoluD == 2 && nGeomD == 2)
    {
        hexRefTypeName = "hexRef4";
    }

    auto* ctorPtr = levelsConstructorTable(hexRefTypeName);

    if (ctorPtr == nullptr)
    {
        FatalErrorInFunction
            << "Unsupported mesh number of dimensions for hex refinement" << nl
            << "nSolutionD: " << nSoluD << ", nGeometricD: " << nGeomD << nl
            << "Only 3D, 2D and 2D axisymmetric mesh refinements are supported"
            << exit(FatalError);
    }

    return autoPtr<hexRef>
    (
        ctorPtr(mesh, cellLevel, pointLevel, level0Edge)
    );
}

// ************************************************************************* //
