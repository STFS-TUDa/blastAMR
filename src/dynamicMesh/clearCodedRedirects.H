/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#ifndef Foam_blastAMR_clearCodedRedirects_H
#define Foam_blastAMR_clearCodedRedirects_H

#include "fvMesh.H"
#include "codedFixedValueFvPatchField.H"
#include "codedMixedFvPatchField.H"

namespace Foam
{

// Reset cached redirect patch fields on coded BCs (codedFixedValue,
// codedMixed) after any topology change (refine, unrefine, redistribute).
//
// OpenFOAM does not override autoMap/rmap on these BCs to clear the cached
// redirectPatchFieldPtr_, so it keeps a pointer at the pre-change patch field
// memory. After autoMap reallocates the parent buffer, the redirect dangles;
// any user "code" block that writes per-face data buffer-overruns the new
// patch field, and destruction of the redirect later double-frees / hits an
// invalid pointer.
//
// Workaround: replace each affected patch field with its clone. The
// codedFixedValueFvPatchField/codedMixedFvPatchField copy constructors
// initialise redirectPatchFieldPtr_ to null, so the clone has a fresh redirect
// that is rebuilt lazily on the next updateCoeffs at the correct post-change
// size.
template<class GeoField>
void clearCodedRedirects(fvMesh& mesh)
{
    using ValueType = typename GeoField::value_type;
    using FixedT    = codedFixedValueFvPatchField<ValueType>;
    using MixedT    = codedMixedFvPatchField<ValueType>;

    HashTable<GeoField*> flds(mesh.template lookupClass<GeoField>());
    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();
        auto& bfld = fld.boundaryFieldRef();
        forAll(bfld, patchi)
        {
            const auto& pf = bfld[patchi];
            if (isA<FixedT>(pf) || isA<MixedT>(pf))
            {
                bfld.set(patchi, pf.clone(fld.internalField()).ptr());
            }
        }
    }
}


// Convenience: clear coded redirects on all common volField types.
inline void clearCodedRedirectsAllVol(fvMesh& mesh)
{
    clearCodedRedirects<volScalarField>(mesh);
    clearCodedRedirects<volVectorField>(mesh);
    clearCodedRedirects<volSphericalTensorField>(mesh);
    clearCodedRedirects<volSymmTensorField>(mesh);
    clearCodedRedirects<volTensorField>(mesh);
}

} // End namespace Foam

#endif
