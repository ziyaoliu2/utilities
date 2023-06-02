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
  Preprocessing utility to modify the fracture surface represented by
  two parallel plates. It overwrites 0 directory with modified surface.
  
  Limitations now are the next:
    - a geometry should include one or two flat plates

Usage
  - surfRoughGen

Needs dictionary
  system/surfRoughGenDict
    
\*---------------------------------------------------------------------------*/

#include <complex>
#include <vector>
#include <fftw3.h>
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "coupledPatchInterpolation.H"
#include "pointPatchField.H"

#include "normalMotionSlipPointPatchVectorField.H"
#include "fixedValuePointPatchField.H"

#include "velocityMotionSolver.H"

#include "motionSolver.H"

/*
 #######################################################################################
 *    Main program body
 #######################################################################################
*/

class RoughnessGenerator
{
protected:
  // Protected data
  int seed;
  int M;
  int N;
  double majLen;
  double minLen;
  double rgh;
  double dHurst;        // Fractal dimension D = 3 - dHurst
  double cutLen;
  double maxDisp;
  word wayToApply;
  
public:
  // Constructor
  RoughnessGenerator
  (
    int seed_,
    int M_,
    int N_,
    double majLen_,
    double minLen_,
    double rgh_,
    double dHurst,
    double cutLen_,
    double maxDisp_,
    word wayToApply_
  )
  :
    seed(seed_),
    M(M_),
    N(N_),
    majLen(majLen_),
    minLen(minLen_),
    rgh(rgh_),
    dHurst(dHurst),
    cutLen(cutLen_),
    maxDisp(maxDisp_),
    wayToApply(wayToApply_)
  {
  }

  void getSurfaceDisplacement
  (
    dynamicFvMesh& mesh,
    scalarField& wd,
    label& patchID,
    int majDir,
    int minDir
  )
  {
    scalarField sFn(M*N, 0.0);
    fftDisp(sFn);
    scalarField sFp(M*N, 0.0);
    if(wayToApply=="asymmetric")
    {
      // shift the seed to get two different numbers
      seed += 125522;
      fftDisp(sFp);
    }
    else
    {
      sFp = sFn;
    }
    
    pointField pointFace = mesh.boundaryMesh()[patchID].faceCentres();
    pointField normlFace = mesh.boundaryMesh()[patchID].faceNormals();
    scalar maxMaj = max( pointFace.component(majDir) );
    scalar maxLat = max( pointFace.component(minDir) );
    scalar minMaj = min( pointFace.component(majDir) );
    scalar minLat = min( pointFace.component(minDir) );
    double Llat = maxLat - minLat;
    double Lmaj = maxMaj - minMaj;
    
    int mveDir = 3 - (minDir + majDir);
    
    if(wayToApply=="internalCylinder")
    {
      // TODO this is just a quick solution, make it general
      scalar cL = constant::mathematical::pi * 198.0;
      forAll(pointFace, i)
      {
        double x,y;
        x = pointFace[i].x();
        y = pointFace[i].y();
        
        scalar phi = Foam::atan2(y,x);
        phi /= constant::mathematical::twoPi;
        if(phi<0) phi += 1.0;
        
        scalar curMaj = pointFace[i].z();
        
        if(curMaj<cL)
        {
          int curm = std::floor(curMaj / cL * (M-1));
          int curn = std::floor(phi * (N-1));
          int ind = curn + N * curm;
          wd[i] = sFn[ind] * ( 1 - curMaj / cL );
        }
      }
    }
    else
    {
      forAll(pointFace, i)
      {
        scalar curMaj = pointFace[i].component(majDir) - minMaj;
        scalar curLat = pointFace[i].component(minDir) - minLat;

        scalar sign = normlFace[i].component(mveDir) 
                / 
                mag(normlFace[i].component(mveDir));

        int curm = std::floor(curMaj / Lmaj * (M-1));
        int curn = std::floor(curLat / Llat * (N-1));
        int ind = curn + N * curm;

        if(wayToApply=="symmetric" || wayToApply=="oneSurface")
          wd[i] = sFn[ind];

        if(wayToApply=="oneSurfaceDecay")
        {
          double factor = ( 1 - curMaj / Llat );
          if( factor<0.0 ) factor = 0.0;
          wd[i] = sFn[ind] * factor;
        }
        
        if(wayToApply=="synchronous")
          wd[i] = sign * sFn[ind];

        if(wayToApply=="asymmetric")
        {
          if(sign<0)
            wd[i] = sFn[ind];
          else
            wd[i] = sFp[ind];
        }
      }
    }
    Info << "Displacement calculated" << nl;
  }

private:
  // converts indexes
  label index(label m, label n){ return n + N * m; }
  
  scalar power(double ksq)
  {
    if (ksq == 0)  return 0;        // <rad^2> ~ 1/ksq^(1+H)
    if (ksq >  1)  return 0;        // cutoff wavelength = cutLen
    scalar p = Foam::pow(ksq, -(dHurst+1) );
//    p *= Foam::exp(-ksq);
    return std::sqrt(p);
  }

  void fftDisp(scalarField& disp)
  {
    unsigned int MN = M*N;
    std::vector<std::complex<double> > f, F;
    f.resize(MN);
    F.resize(MN);

    Random rnd( seed );
    scalar TwoPi = constant::mathematical::twoPi;

    Info <<  "Displacement calc starts...." << nl;
    /*
     *   --- ---
     *  | 1 | 2 |
     *   --- ---
     *  | 3 | 4 |
     *   --- ---
     */
    // calculating 1 and 4
    for(int m=0; m<M/2+1; ++m)
    {
      for(int n=0; n<N/2+1; ++n)
      {
        scalar p = TwoPi * rnd.sample01<scalar>();
        scalar rad;
        double majk, mink, ksq;
        majk = m*cutLen/majLen;
        mink = n*cutLen/minLen;
        ksq   = majk*majk + mink*mink;
        rad = power(ksq) * rnd.GaussNormal<scalar>();

        f[ index(m,n) ] = 
                rad * std::complex<double>(Foam::cos(p),  Foam::sin(p));
        f[ index(((M-m)%M),(N-n)%N) ] = 
                rad * std::complex<double>(Foam::cos(p), -Foam::sin(p));
      }
    }
    f[ index(M/2,0)   ].imag(0.0);
    f[ index(0,  N/2) ].imag(0.0);
    f[ index(M/2,N/2) ].imag(0.0);

    for(int m=1; m<M/2; ++m)
    {
      for(int n=1; n<N/2; ++n)
      {
        scalar p = TwoPi * rnd.sample01<scalar>();
        scalar rad;
        double majk, mink, ksq;
        majk = TwoPi*m/majLen;
        mink = TwoPi*n/minLen;
        ksq  = majk*majk + mink*mink;
        rad  = power(ksq) * rnd.GaussNormal<scalar>();

        f[ index(  m, N-n) ] = 
                rad * std::complex<double>(Foam::cos(p),  Foam::sin(p));
        f[ index(M-m,   n) ] = 
                rad * std::complex<double>(Foam::cos(p), -Foam::sin(p));
      }
    }

    fftw_plan plan;
    plan = fftw_plan_dft_2d(M, N,
                             reinterpret_cast<fftw_complex*>(&f[0]),
                             reinterpret_cast<fftw_complex*>(&F[0]),
                             FFTW_BACKWARD,
                             FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);


    scalarField sF(MN);
    forAll(sF, ii)
    {
      sF[ii] = F[ii].real();
    }

    scalarField sF2 = sqr(sF);
    scalar avSF     = average(sF);
    scalar avSF2    = average(sF2);

    scalar factor = rgh / Foam::sqrt( mag(avSF2 - sqr(avSF)) );

    sF *= factor;

    disp = sF;
  }
};


int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  // reading dictionary surfRoughGenDict
  IOdictionary surfRoughGenDict
  (
    IOobject
    (
      "surfRoughGenDict",
      runTime.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );
  
  word wayToApply;
  if( !surfRoughGenDict.readIfPresent<word>("apply", wayToApply) ){
    SeriousErrorIn("main")
         << "There is no `synchronous` parameter in dictionary"
         << exit(FatalError);
  }
  word patchName;
  if( !surfRoughGenDict.readIfPresent<word>("patchName", patchName) )
  {
    SeriousErrorIn("main")
         <<  "There is no `patchName` parameter in dictionary"
         <<  exit(FatalError);
  }
  
  int majDir;
  if( !surfRoughGenDict.readIfPresent<int>("majDir", majDir) ){
    SeriousErrorIn("main")
         << "There is no `majDir` parameter in dictionary"
         << exit(FatalError);
  }
  int minDir;
  if( !surfRoughGenDict.readIfPresent<int>("minDir", minDir) ){
    SeriousErrorIn("main")
         << "There is no `minDir` parameter in dictionary"
         << exit(FatalError);
  }
  int majLen;
  if( !surfRoughGenDict.readIfPresent<int>("majLen", majLen) ){
    SeriousErrorIn("main")
         << "There is no `majLen` parameter in dictionary"
         << exit(FatalError);
  }
  int minLen;
  if( !surfRoughGenDict.readIfPresent<int>("minLen", minLen) ){
    SeriousErrorIn("main")
         << "There is no `minLen` parameter in dictionary"
         << exit(FatalError);
  }
  int M;
  if( !surfRoughGenDict.readIfPresent<int>("majNum", M) ){
    SeriousErrorIn("main")
         << "There is no `majNum` parameter in dictionary"
         << exit(FatalError);
  }
  int N;
  if( !surfRoughGenDict.readIfPresent<int>("minNum", N) ){
    SeriousErrorIn("main")
         << "There is no `minNum` parameter in dictionary"
         << exit(FatalError);
  }
  
  
  int seed;
  if( !surfRoughGenDict.readIfPresent<int>("seed", seed) ){
    SeriousErrorIn("main")
         << "There is no `seed` parameter in dictionary"
         << exit(FatalError);
  }
  scalar rgh;
  if( !surfRoughGenDict.readIfPresent<scalar>("roughness", rgh) ){
    SeriousErrorIn("main")
         << "There is no `roughness` parameter in dictionary"
         << exit(FatalError);
  }
  double dHurst;
  if( !surfRoughGenDict.readIfPresent<double>("dHurst", dHurst) ){
    SeriousErrorIn("main")
         << "There is no `dHurst` parameter in dictionary"
         << exit(FatalError);
  }
  double cutLen;
  if( !surfRoughGenDict.readIfPresent<double>("cutLen", cutLen) ){
    SeriousErrorIn("main")
         << "There is no `cutLen` parameter in dictionary"
         << exit(FatalError);
  }
  double maxDisp;
  if( !surfRoughGenDict.readIfPresent<double>("maxDisp", maxDisp) ){
    SeriousErrorIn("main")
         << "There is no `maxDisp` parameter in dictionary"
         << exit(FatalError);
  }
  
  Info <<  "patch:         "  <<  patchName   <<  endl;
  Info <<  "apply (method):"  <<  wayToApply  <<  endl;
  Info <<  "majDir:        "  <<  majDir      <<  endl;
  Info <<  "minDir:        "  <<  minDir      <<  endl;
  Info <<  "majLen:        "  <<  majLen      <<  endl;
  Info <<  "minLen:        "  <<  minLen      <<  endl;
  Info <<  "majNum:        "  <<  M           <<  endl;
  Info <<  "minNum:        "  <<  N           <<  endl;
  Info <<  "seed:          "  <<  seed        <<  endl;
  Info <<  "roughness:     "  <<  rgh         <<  endl;
  Info <<  "dHurst:        "  <<  dHurst      <<  endl;
  Info <<  "cutLen:        "  <<  cutLen      <<  endl;
  Info <<  "maxDisp:       "  <<  maxDisp     <<  endl;
  Info <<  "Setup RoughnessGenerator class"   <<  endl;

  RoughnessGenerator rg( seed, M, N, majLen, minLen, 
                         rgh, dHurst, cutLen, maxDisp, wayToApply );
  
  double cpuTime = runTime.elapsedCpuTime();
  
  // Get patch ID for moving boundaries
  label patchID  = mesh.boundaryMesh().findPatchID(patchName);
  if( patchID==-1 ){
    SeriousErrorIn("main")
         << "patch "  << patchName  << " is missing"  << exit(FatalError);
  }
  
  coupledPatchInterpolation patchInterpolator
    ( mesh.boundaryMesh()[patchID], mesh );
  
  const pointField& boundaryPoints 
    = mesh.boundaryMesh()[patchID].localPoints();
  vectorField pointDispWall(boundaryPoints.size(), vector::zero);
  vectorField pointNface = mesh.boundaryMesh()[patchID].faceNormals();
  vectorField motionN = patchInterpolator.faceToPointInterpolate(pointNface);
  forAll(motionN, ii) motionN[ii]/=mag(motionN[ii]);
  
  scalarField faceDisp(pointNface.size(), 0.0);
  rg.getSurfaceDisplacement(mesh, faceDisp, patchID, majDir, minDir);
  scalarField pointDisp = patchInterpolator.faceToPointInterpolate(faceDisp);
  
  forAll( pointDisp, i )
  {
    pointDisp[i] = std::min(pointDisp[i],  maxDisp);
    pointDisp[i] = std::max(pointDisp[i], -maxDisp);
  }
  
  Info <<  "Maximum and minimum face displacements  " 
       <<  max(faceDisp)  <<  "  " << min(faceDisp) << endl;
  Info <<  "Maximum and minimum point displacements  "
       <<  max(pointDisp) <<  "  " << min(pointDisp) << endl;
  
  forAll( pointDispWall, i )
    pointDispWall[i] = pointDisp[i] * motionN[i];
  
  pointVectorField& pointVelocity = const_cast<pointVectorField&>
  (
    mesh.objectRegistry::lookupObject<pointVectorField>( "pointMotionU" )
  );
  
  pointVelocity.boundaryFieldRef()[patchID] == pointDispWall;
  mesh.update();
  
  cpuTime = runTime.elapsedCpuTime() - cpuTime;
  
  Info << nl << "Time statistics:" << nl;
  
  int wlNP = mesh.boundaryMesh()[patchID].nPoints();
  Info << "Total number of points:                   " << mesh.nPoints() << nl;
  Info << "Number of points on the walls:            " << wlNP << nl;
  Info << "Running time:                             " << cpuTime << nl << endl;

  Info << "Overwriting points in current time directory." << nl;
  runTime.writeNow();
  
  Info << "End" << nl;
  return 0;
}

// **************************** End of the solver ******************************** //
