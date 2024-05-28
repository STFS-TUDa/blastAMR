/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 AUTHOR,AFFILIATION
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

Application
    updateMesh

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "ReadFields.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption("overwrite", "Overwrite the current polyMesh in <time>");
    argList::addOption("iterations", "n iterations ","How many iterations of refinement should be performed"); 
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    bool overwrite = args.found("overwrite");
    label iterations = args.getOrDefault<label>("iterations",1);
    Info << "Number of refinements: "<<iterations<<endl;
    #include "createTime.H"
    auto Times = runTime.times();
    #include "checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);

    #include "createDynamicFvMesh.H"
    #include "createFields.H"

    for (int i=1; i <= iterations ; i++)
    {
        if (!overwrite)
        {
            runTime++;
        }
	Info << "updating Mesh"<<endl;
	mesh.update();
    }
    runTime.writeNow();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
