/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::errorEstimators::codedErrorEstimator

Description
    An error Estimator for AMR with the code supplied by an on-the-fly compiled C++
    expression.

    The code entries:
    \plaintable
       codeInclude | include files
       codeOptions | compiler line: added to EXE_INC (Make/options)
       codeLibs    | linker line: added to LIB_LIBS (Make/options)
       localCode   | c++; local static functions
       code        | c++; code to select cells to be refined
       codeRead    | c++; read dictionary into class members
    \endplaintable

Usage
    Example:
    \verbatim
    errorEstimator coded;      // Model used to estimate error
    name percentOfMax;
    
    code
    #{
        Info<< "---->! custom error estimator !<----" << endl;
        const auto& T = mesh_.lookupObject<volScalarField>("T");
        error_.interalField() = T.interalField();
    
        const scalar min = gMin(T);
        const scalar max = gMax(T);
        const scalar range = max - min;
    
        lowerRefine_ = min + 0.8*range;
        upperRefine_ = GREAT;
        lowerUnrefine_ = min + 0.8*range;
        upperUnrefine_ =  GREAT;

        if (scale) normalize(error_);
        Info<< "---->! end    error estimator !<----" << endl;
    #};
    \endverbatim

Note

See also
    Foam::dynamicCode
    Foam::errorEstimator

SourceFiles
    codedErrorEstimator.C

\*---------------------------------------------------------------------------*/

#ifndef codedErrorEstimator_H
#define codedErrorEstimator_H

#include "codedBase.H"
#include "errorEstimator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{

/*---------------------------------------------------------------------------*\
                         Class codedErrorEstimator Declaration
\*---------------------------------------------------------------------------*/

class codedErrorEstimator
:
    public errorEstimator,
    protected codedBase
{
    // Private Data

        //- Dictionary contents for the function
        const dictionary dict_;

        word redirectName_;

        // custom code fields
        string codeData_;
        string codeDataConstruct_;
        string codeRead_;

        //- Pointer to the real refinement error estimator object
        mutable autoPtr<errorEstimator> redirectErrorEstimatorPtr_;

    // Private Member Functions

        //- Get reference to the underlying error estimator
        errorEstimator& redirectErrorEstimator();


protected:

    // Protected Member Functions

        //- Mutable access to the loaded dynamic libraries
        virtual dlLibraryTable& libs() const;

        //- Description (type + name) for the output
        virtual string description() const;

        //- Clear redirected object(s)
        virtual void clearRedirect() const;

        //- Additional 'codeContext' dictionary to pass through
        virtual const dictionary& codeContext() const;

        // Get the code (sub)dictionary
        virtual const dictionary& codeDict(const dictionary& dict) const;

        // Get the code dictionary
        virtual const dictionary& codeDict() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;


    // Generated Methods

        //- No copy assignment
        void operator=(const codedErrorEstimator&) = delete;

public:

    // Static Data Members

        //- Name of the C code template to be used
        static constexpr const char* const codeTemplateC
            = "codedErrorEstimatorTemplate.C";

        //- Name of the H code template to be used
        static constexpr const char* const codeTemplateH
            = "codedErrorEstimatorTemplate.H";


    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from mesh, entry name and dictionary

        //- Copy construct
        codedErrorEstimator
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& name
        );



    //- Destructor
    virtual ~codedErrorEstimator() = default;


    // Member Functions

        //- Update error
        virtual void update(const bool scale = true);

        // Integrate etc are not implemented!

        //- Write in dictionary format
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace errorEstimators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
