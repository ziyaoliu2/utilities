# Pre- and post-processing utilities

1) surfRoughGen: OpenFOAM utility to generate a rough surface by Fourier synthesis.

usage: surfRoughGen

It has a dictionary system/surfRoughGenDict - the case dissolFrac uses surfRoughGen

2) runMeshUpdateOnce: OpenFOAM utility for a single mesh update cycle (after surfRoughGen)

usage: runMeshUpdateOnce

3) dissolCalc: OpenFOAM utility to postprocess fracture fields by integrating over the aperture

usage: dissolCalc fieldMap2D all 1000 100

Calculates all the 2D fields using 1000 cells in flow direction and 100 cells in transverse direction. It uses a dictionary fieldMap2Ddict

orderBoundaries: User library adds support for rearranging boundaries after meshing
