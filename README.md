# miniwasp
This repository contains codes for simulating the optics of the 
eye of a mini wasp. We assume that the electric permitivity and the
magnetic permeability is piecewise constant in each component 
of the eye and assume that the incoming electric and magnetic fields
are time harmonic. 

The library currently depends on FMM3D 
(https://github.com/flatironinstitute/FMM3D), and 
fmm3dbie (https://gitlab.com/fastalgorithms/fmm3dbie).

In order to build the python package, make sure the python packages and
the fortran libraries corresponding to FMM3D and fmm3dbie are installed. 
This can be done by running `make install` and `make python` in the
respective repositories. See
[FMM3D](https://fmm3d.readthedocs.io/en/latest/install.html) and 
[fmm3dbie](https://fmm3dbie.readthedocs.io/en/latest/install.html) 
for detailed installation instructions and optimized compilation of the
respective libraries.

The easiest way to run this library is through it's python interface.
To build the python interface, run `make python`. For a simple demo of
the user callable routines of this library see `python\example1.py`.

The documentation of the python routines can be found
[here](https://miniwaspbie.readthedocs.io/en/latest)

