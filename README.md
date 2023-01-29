# The AMR library from blastFoam

This repository extracts library parts from [blastFoam](https://github.com/synthetik-technologies/blastfoam)
which are relevant to adaptive mesh refinement while retaining the original commit history.

## Objectives
- Have a reasonable Load-balanced AMR for 2D/3D OpenFOAM meshes in the **ESI fork**.
- Support load balancing of Lagrangian particle clouds.
- Make it easy to include this lib in other projects through submodules/subtrees.

`blastFoam` is GPL licensed; see the included [COPYING](COPYING).

## Notes

- Everything compiles to your `$FOAM_USER_LIBBIN` with a recent ESI ([openfoam.com](https://openfoam.com)) version.
