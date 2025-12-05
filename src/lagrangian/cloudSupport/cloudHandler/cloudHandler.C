/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "cloudHandler.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cloudHandler, 0);
    defineRunTimeSelectionTable(cloudHandler, cloud);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cloudHandler> Foam::cloudHandler::New(const word& handlerName)
{
    const auto* tablePtr = cloudConstructorTablePtr_;

    if (!tablePtr)
    {
        FatalErrorInFunction
            << "Cloud handler RTS table not initialized. "
            << "No handlers are registered." << nl
            << exit(FatalError);
        return nullptr;
    }

    auto cstrIter = tablePtr->cfind(handlerName);

    if (!cstrIter.good())
    {
        FatalErrorInFunction
            << "Unknown cloudHandler type '" << handlerName << "'" << nl << nl
            << "Valid cloudHandler types are:" << nl
            << tablePtr->sortedToc() << nl
            << exit(FatalError);
        return nullptr;
    }

    return autoPtr<cloudHandler>(cstrIter.val()());
}


Foam::autoPtr<Foam::cloudHandler> Foam::cloudHandler::tryNew
(
    cloud& c,
    const fvMesh& mesh,
    const word& handlerName
)
{
    // If explicit handler name provided, use it directly
    if (!handlerName.empty())
    {
        DebugInfo
            << "cloudHandler::tryNew: Using explicit handler '"
            << handlerName << "' for cloud '" << c.name() << "'" << endl;

        return New(handlerName);
    }

    // Otherwise, search for a handler that can handle this cloud
    const auto* tablePtr = cloudConstructorTablePtr_;

    if (!tablePtr)
    {
        WarningInFunction
            << "Cloud handler RTS table not initialized. "
            << "No handlers are registered." << endl;
        return nullptr;
    }

    if (tablePtr->empty())
    {
        WarningInFunction
            << "Cloud handler RTS table is empty. "
            << "No handlers are registered." << endl;
        return nullptr;
    }

    // Collect all handlers that can handle this cloud, sorted by priority
    DynamicList<std::pair<label, word>> candidates;

    DebugInfo
        << "cloudHandler::tryNew: Testing handlers for cloud type '"
        << c.type() << "' (name: " << c.name() << ")" << nl
        << "  Available handlers: " << tablePtr->size() << endl;

    // Iterate over the hash table using forAllConstIters
    forAllConstIters(*tablePtr, iter)
    {
        // The key (handler type name)
        const word& name = iter.key();

        // The value is a function pointer - call it to construct a handler
        cloudConstructorPtr constructorPtr = iter.val();

        // Construct the handler to test if it can handle this cloud
        autoPtr<cloudHandler> handler(constructorPtr());

        if (handler && handler->canHandle(c))
        {
            DebugInfo
                << "  Handler '" << name << "' can handle cloud" << endl;
            candidates.append(std::make_pair(handler->priority(), name));
        }
    }

    if (candidates.empty())
    {
        return nullptr;
    }

    // Sort by priority (descending) - highest priority first
    std::sort
    (
        candidates.begin(),
        candidates.end(),
        [](const auto& a, const auto& b)
        {
            return a.first > b.first;
        }
    );

    // Return the highest priority handler
    const word& selectedHandler = candidates[0].second;
    auto cstrIter = tablePtr->cfind(selectedHandler);

    if (cstrIter.good())
    {
        return autoPtr<cloudHandler>(cstrIter.val()());
    }

    return nullptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cloudHandler::autoMap(cloud& c, const mapPolyMesh& map)
{
    // Default implementation: just call the cloud's virtual autoMap
    c.autoMap(map);
}


// ************************************************************************* //
