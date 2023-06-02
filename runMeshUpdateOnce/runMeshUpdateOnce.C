/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-201X OpenFOAM Foundation
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

Application
  surfRoughGen

Description
  This tool moves the dynamic mesh once by calling mesh.update() function once.
  Then it writes data in current time directory.

Usage
  - runMeshUpdateOnce

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "coupledPatchInterpolation.H"
//#include "pointPatchField.H"
//#include "normalMotionSlipPointPatchVectorField.H"
//#include "fixedValuePointPatchField.H"
//#include "velocityMotionSolver.H"
//#include "motionSolver.H"

/*
 #######################################################################################
 *    Main program body
 #######################################################################################
*/

int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
  double cpuTime = runTime.elapsedCpuTime();
  
  mesh.update();
  
  cpuTime = runTime.elapsedCpuTime() - cpuTime;
  
  Info << nl << "Time statistics:" << nl;
  
  Info << "Running time:                             " << cpuTime << nl << endl;

  Info << "Overwriting points in current time directory." << nl;
  runTime.writeNow();
  
  Info << "End" << nl;
  return 0;
}

// **************************** End of the solver ******************************** //
