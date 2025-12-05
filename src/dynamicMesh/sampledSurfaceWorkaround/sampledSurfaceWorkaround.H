/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2024
     \\/     M anipulation  | Synthetik Applied Technologies
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

Description
    Workaround for OpenFOAM's surfaceFieldValue not properly expiring its
    sampledSurface cache after mesh topology changes. This uses the "member
    smuggler" technique to access the private sampledPtr_ member.

\*---------------------------------------------------------------------------*/

#ifndef sampledSurfaceWorkaround_H
#define sampledSurfaceWorkaround_H

#include "surfaceFieldValue.H"
#include "sampledSurface.H"
#include "timeControlFunctionObject.H"
#include "Pstream.H"
#include "messageStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Template machinery to access private members
template<class T, class Secret>
struct memberSmuggler
{
    static inline T v{};
};

template<class Secret, auto T>
struct memberSmugglerTag
{
    // The standard says access specifier constraints are lifted for names
    // in **declarations** of explicit instantiation/specialization;
    // so we smuggle the member with a chained assignment
    static const inline auto v = memberSmuggler<decltype(T), Secret>::v = T;
};

// Instantiate to extract surfaceFieldValue::sampledPtr_
using SFV = functionObjects::fieldValues::surfaceFieldValue;
using SampledPtrType = autoPtr<sampledSurface> SFV::*;

// Use a class in the Foam namespace for the secret type
struct SampledPtrSecret {};
template struct memberSmugglerTag<SampledPtrSecret, &SFV::sampledPtr_>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Expire sampled surfaces in all surfaceFieldValue function objects
//  Call this after mesh topology changes (refinement, redistribution)
inline void expireSampledSurfaces(const Time& time, const bool enabled)
{
    auto& funcObjList = const_cast<Time&>(time).functionObjects();

    // Track if we found any surfaceFieldValue with sampledPtr
    bool foundTargets = false;

    for (label i = 0; i < funcObjList.size(); ++i)
    {
        functionObject& fo = funcObjList[i];

        // Check if it's a timeControl wrapper (common case)
        auto* tcPtr =
            dynamic_cast<functionObjects::timeControl*>(&fo);

        functionObject* targetFO =
            tcPtr
          ? const_cast<functionObject*>(&(tcPtr->filter()))
          : &fo;

        // Try to cast to surfaceFieldValue
        auto* sfvPtr = dynamic_cast<SFV*>(targetFO);
        if (sfvPtr)
        {
            // Access the private sampledPtr_ member
            auto& sampledPtr =
                sfvPtr->*memberSmuggler<SampledPtrType, SampledPtrSecret>::v;

            if (sampledPtr)
            {
                foundTargets = true;

                if (enabled)
                {
                    sampledPtr->expire();
                    Info<< "    Expired sampledSurface in surfaceFieldValue: "
                        << fo.name() << endl;
                }
            }
        }
    }

    // Warn on master if we found targets but expiration is disabled
    if (foundTargets && !enabled && Pstream::master())
    {
        WarningInFunction
            << "surfaceFieldValue function objects with sampledSurface "
            << "detected, but expireSampledSurfacesOnLB is disabled." << nl
            << "This may cause garbage values in sampled data after mesh "
            << "topology changes (refinement/load balancing)." << nl
            << "Set 'expireSampledSurfacesOnLB true;' in dynamicMeshDict "
            << "to enable automatic cache invalidation." << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
