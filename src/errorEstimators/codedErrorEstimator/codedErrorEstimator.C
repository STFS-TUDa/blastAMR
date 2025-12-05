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

\*---------------------------------------------------------------------------*/

#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "codedErrorEstimator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(codedErrorEstimator, 0);
    addToRunTimeSelectionTable(errorEstimator, codedErrorEstimator, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::dlLibraryTable&
Foam::errorEstimators::codedErrorEstimator::libs() const
{
    return const_cast<Time&>(this->time()).libs();
}


Foam::string
Foam::errorEstimators::codedErrorEstimator::description() const
{
    return "codedErrorEstimator " + redirectName_;
}


void Foam::errorEstimators::codedErrorEstimator::clearRedirect() const
{
    redirectErrorEstimatorPtr_.clear();
}


const Foam::dictionary&
Foam::errorEstimators::codedErrorEstimator::codeContext() const
{
    // What else would make sense?
    return dict_;
}


const Foam::dictionary&
Foam::errorEstimators::codedErrorEstimator::codeDict
(
    const dictionary& dict
) const
{
    // Use named subdictionary if present to provide the code.
    // This allows running with multiple Function1s

    return
    (
        dict.found("code")
      ? dict
      : dict.subDict(redirectName_)
    );
}


const Foam::dictionary&
Foam::errorEstimators::codedErrorEstimator::codeDict() const
{
    return codeDict(dict_);
}


void Foam::errorEstimators::codedErrorEstimator::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    if (context.code().empty())
    {
        FatalIOErrorInFunction(dict_)
            << "No code section in input dictionary for errorEstimator "
            << " name " << redirectName_
            << exit(FatalIOError);
    }

    // Take no chances - typeName must be identical to redirectName_
    dynCode.setFilterVariable("typeName", redirectName_);
    dynCode.setFilterVariable("codeData", codeData_);
    dynCode.setFilterVariable("codeDataConstruct", codeDataConstruct_);
    dynCode.setFilterVariable("codeRead", codeRead_);

    // Set TemplateType and FieldType filter variables
    //dynCode.setFieldTemplates<Type>();

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH);

    // Define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = \\\n"
        "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
        "-I$(AMRLB_PROJECT)/src/errorEstimators/lnInclude "
      + context.options()
      + "\n\nLIB_LIBS = \\\n"
        "    -lOpenFOAM \\\n"
        "    -lfiniteVolume \\\n"
        "    -lmeshTools "
      + context.libs()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::codedErrorEstimator::codedErrorEstimator
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    errorEstimator(mesh, dict, name),
    codedBase(),
    dict_(dict),
    redirectName_(dict.getOrDefault<word>("name", name)),
    codeData_(dict.getOrDefault<string>("codeData", "")),
    codeDataConstruct_(dict.getOrDefault<string>("codeDataConstruct", "")),
    codeRead_(dict.getOrDefault<string>("codeRead", ""))
{
    //this->codedBase::setCodeContext(dict_);
    read(dict_);

    // No additional code chunks...

    updateLibrary(redirectName_);
    redirectErrorEstimator();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::errorEstimator&
Foam::errorEstimators::codedErrorEstimator::redirectErrorEstimator()
{
    if (!redirectErrorEstimatorPtr_)
    {
        dictionary constructDict(dict_);
        constructDict.set("errorEstimator", redirectName_);

        redirectErrorEstimatorPtr_ = errorEstimator::New
        (
            mesh_,
            constructDict,
            redirectName_
        );
    }
    return *redirectErrorEstimatorPtr_;
}


void Foam::errorEstimators::codedErrorEstimator::update(const bool scale)
{
    // Ensure library containing user-defined code is up-to-date
    updateLibrary(redirectName_);
    redirectErrorEstimator().update(scale);
    error_ = redirectErrorEstimator().error();
}


bool Foam::errorEstimators::codedErrorEstimator::writeData
(
    Ostream& os
) const
{
    dict_.writeEntry(this->name(), os);
    return true;
}


// ************************************************************************* //
