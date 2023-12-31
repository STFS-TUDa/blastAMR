/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::hexRef2DAxi

Description
    Refinement of (split) hexes using polyTopoChange.

SourceFiles
    hexRef2DAxi.C

Authors
    Original code:
    Ahmad Baniabedalruhman 2015
    Adapted by:
    Yun Ouedraogo 2017

\*---------------------------------------------------------------------------*/

#ifndef hexRef2DAxi_H
#define hexRef2DAxi_H

#include "hexRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class hexRef2DAxi Declaration
\*---------------------------------------------------------------------------*/

class hexRef2DAxi : public hexRef
{

protected:

        //- Get cell added to point of celli (if any)
        label getAnchorCell
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const label celli,
            const label facei,
            const label pointi
        ) const;

        //- Store in maps correspondence from midpoint to anchors and faces.
        label storeMidPointInfo
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& faceMidPoint,
            const labelList& edgeMidPoint,
            const label celli,
            const label facei,
            const bool faceOrder,
            const label midPointI,
            const label anchorPointI,
            const label faceMidPointI,

            Map<edge>& midPointToAnchors,
            Map<edge>& midPointToFaceMids,
            polyTopoChange& meshMod
        ) const;

        //- Create all internal faces to split celli into 4.
        void createInternalFaces
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& faceMidPoint,
            const labelList& faceAnchorLevel,
            const labelList& edgeMidPoint,
            const label celli,
            polyTopoChange& meshMod
        ) const;

        //- Disallow default bitwise copy construct
        hexRef2DAxi(const hexRef2DAxi&);

        //- Disallow default bitwise assignment
        void operator=(const hexRef2DAxi&);


public:

    //- Runtime type information
    TypeName("hexRef2DAxi");


    // Constructors

        //- Construct from mesh, read_if_present refinement data
        //  (from write below). If readHistory is true does read_if_present
        //  of refinement history. If false clears all history
        hexRef2DAxi(const polyMesh& mesh, const bool readHistory = true);

        //- Construct from mesh and un/refinement data and optional size of
        //  starting cells
        hexRef2DAxi
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const hexRefRefinementHistory& history,
            const scalar level0Edge = -GREAT
        );

        //- Construct from mesh and refinement data and optional size of
        //  starting cells
        hexRef2DAxi
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const scalar level0Edge = -GREAT
        );


    // Member Functions

        // Refinement

            //- Insert refinement. All selected cells will be split into 8.
            //  Returns per element in cells the 8 cells they were split into.
            //  Guarantees that the 0th element is the original cell label.
            //  Mapping:
            //  -split cells: 3 new ones get added from original
            //  -split faces: original gets modified; 3 new ones get added
            //               from original
            //  -added internal faces: added from original cell face(if
            //   that was internal) or created out-of-nothing (so will not
            //   get mapped!). Note: could make this inflate from point but
            //   that will allocate interpolation.
            //  -points added to split edge: added from edge start()
            //  -midpoints added: added from cellPoints[0].
            labelListList setRefinement
            (
                const labelList& cells,
                polyTopoChange&
            );

        // Unrefinement
            //- Return the mid elems in top-level split cells
            //  that can be unsplit.
            labelList getSplitElems() const;

            //- Given proposed
            //  splitElems to unrefine according to calculate any clashes
            //  (due to 2:1) and return ok list of mid elems to unrefine.
            labelList consistentUnrefinement
            (
                const labelList& elemsToUnrefine,
                const bool maxSet
            ) const;

            //- Remove some refinement. Needs to be supplied output of
            //  consistentUnrefinement. Only call if undoable set.
            //  All n edgeCells of a split edge will be combined into
            //  the lowest numbered cell of those n.
            void setUnrefinement
            (
                const labelList& splitElemLabels,
                polyTopoChange&
            );

            labelList selectUnrefineElems
            (
                const scalar unrefineLevel,
                const boolList& markedCell,
                const scalarField& vFld
            ) const;

            void calcFaceToSplitPoint
            (
                const labelList& splitElems,
                Map<label>& faceToSplitPoint
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
