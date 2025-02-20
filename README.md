## Code description ##

This repository contains the source code of the latest version of the `darcyInterTransportFoam` solver, as well as that of its complementary pre- and post-processing utilities and boundary conditions developed alongside it. 
`darcyInterTransportFoam` is a fully-coupled, three-dimensional surface water–groundwater (SW-GW) solver implemented in the powerful, open-source CFD toolbox [**OpenFOAM**](https://www.openfoam.com/). 
The model is designed to simulate flow, conservative solute transport and heat transfer processes across both surface and subsurface domains, along with their corresponding interactions.

The hosted source code is organized in the following directories:

* <ins>**applications**</ins>:
  - **solver**: includes the fully-coupled model.
  - **utilities**: contains the complementary utilities of `darcyInterTransportFoam`.

* <ins>**src**</ins>: comprises the newly implemented boundary conditions as well as modified libraries to be used with `darcyInterTransportFoam`.

The installation procedure is the same as for any [**OpenFOAM**](https://www.openfoam.com/) code. 
Therefore, both the solver and its supplements can be installed either as local applications in the user's account (`$FOAM_USER_APPBIN` or `$FOAM_USER_LIBBIN` in the *Make/files* file) 
or as standard applications in the [**OpenFOAM**](https://www.openfoam.com/) installation directory (`$FOAM_APPBIN` or `$FOAM_LIBBIN` in the *Make/files* file). 
In both cases, the code must be compiled using [**OpenFOAM**](https://www.openfoam.com/)'s compilation tool, `wmake`.

## Author ##

**Pardo-Álvarez, Álvaro**  
Ph.D. Student  
[Centre for Hydrogeology and Geothermics](https://www.unine.ch/chyn)  
University of Neuchâtel
