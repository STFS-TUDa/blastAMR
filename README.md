# A Load-balanced adaptive mesh refinement for OpenFOAM

> [!NOTE]
> If you're using OpenFOAM v2112 or newer, please switch to the `v2212` branch. The `master` branch is for OpenFOAM v2012
> and earlier versions  (v2106 is not supported). See #6. The `Allwmake` does this automatically for you now.

[![Build Status](https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2FSTFS-TUDa%2FblastAMR%2Fbadge%3Fref%3Dmaster%26&style=for-the-badge&label=OF2006%28master%29)](https://actions-badge.atrox.dev/STFS-TUDa/blastAMR/goto?ref=master) [![Build Status](https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2FSTFS-TUDa%2FblastAMR%2Fbadge%3Fref%3Dv2212%26&style=for-the-badge&label=OF2212%28v2212%29)](https://actions-badge.atrox.dev/STFS-TUDa/blastAMR/goto?ref=v2212) [![Build Status](https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2FSTFS-TUDa%2FblastAMR%2Fbadge%3Fref%3Dv2212%26&style=for-the-badge&label=OF2206%28v2212%2CGCC13%29)](https://actions-badge.atrox.dev/STFS-TUDa/blastAMR/goto?ref=v2212) ![Testing framework](https://img.shields.io/badge/tested_with_foamUT-00000?style=for-the-badge) ![GitHub contributors](https://img.shields.io/github/contributors/STFS-TUDa/blastAMR?style=for-the-badge) ![GitHub issues](https://img.shields.io/github/issues/STFS-TUDa/blastAMR?style=for-the-badge) ![GitHub release (with filter)](https://img.shields.io/github/v/release/STFS-TUDa/blastAMR?style=for-the-badge) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.8427734-lightgreen?style=for-the-badge)](https://zenodo.org/record/8427734)

<details open><summary>H2 injection case</summary>

https://github.com/STFS-TUDa/blastAMR/assets/34474472/44cdc485-0c90-41b0-90b7-956e013abf9c

![2023-10-10_14-30](https://github.com/STFS-TUDa/blastAMR/assets/34474472/086dec0e-392f-4cf5-a769-d91f174a7b33)

</details>

<details open><summary>Free-propagating flame case</summary>

https://github.com/STFS-TUDa/blastAMR/assets/34474472/29d1cf64-ade9-4c5c-8044-8e63797d31b7

![2023-10-10_14-29](https://github.com/STFS-TUDa/blastAMR/assets/34474472/cab32728-b497-438b-8150-940a3449a188)

</details>

This repository extracts library parts from [blastFoam](https://github.com/synthetik-technologies/blastfoam) which are relevant to load-balanced adaptive mesh refinement for polyhedral meshes while retaining the original commit history; this is not intended as a full port! Only **tested features** are to be trusted. **!!WIP!!**

To get started, you can visit the one and only [wiki page](https://github.com/STFS-TUDa/blastAMR/wiki).

## Objectives

- Have a reasonable Load-balanced AMR for 2D/3D OpenFOAM meshes in the **ESI fork**.
- Support load balancing of Lagrangian particle clouds (not yet implemented).
- Make it easy to include this lib in other projects through submodules/subtrees.

`blastFoam` is GPL licensed; see the included (original) [COPYING](COPYING).

## Quick Notes

- Everything compiles to your `$FOAM_USER_LIBBIN` with a recent ESI ([openfoam.com](https://openfoam.com)) version.
- Our favorite refinement cell selector is the [`coded`](https://github.com/STFS-TUDa/blastAMR/blob/0e70469d718ee53a5ce892e0e24a4f940bfba369/tutorials/poly_freelyPropFlame2D_H2/constant/dynamicMeshDict#L25). Although a wide range of selectors is provided, they are not well tested. So, use the coded refinement whenever possible.
