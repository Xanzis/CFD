A set of tools for CFD steady flow analysis, transient modelling, and eventual combustion modeling.

Current progress: 
Working on developing 1d situation, including the associated tools for manipulating arrays of scalars and velocites.

nrutil.h is taken from numerical recipes in C. The source is currently on git due to lack of licenses.

Status:

Due to misalignment between cfd_model_1d and the new fieldutil and matrixutil header files, development of the former has been halted. Hopefully with this next program I'll be able to keep train of thought intact.

For this purpose, here's a roadmap

------ROADMAP------

The next two todo items:
[PASSABLE DONE]- Write a simple matrixutil.c with a matrix inversion function
[DONE]- Tackle the simplest possible case: one-segment (but arbitrary number of segments) 1d flow-diffusion of a single property phi in a known velocity field, with static densities, pressures, etc.

Future items:
[IN PROGRESS]- add hybrid differencing schemes (QUICK looks unnecessary)
[IN PROGRESS - delayed for now]- basic vector or matrix grapher
- Arbitrary boundary-defined compartments along a 1d line of cells
[IN PROGRESS]- Full 1d SIMPLE - velocity field and pressure solutions 
- 2d: flow-diffusion solutions with a given velocity field in a square grid of cells
- Full SIMPLE in 2d (steady state)
- Optimize matrix solvers for band-diagonal systems
- *configuring boundaries to allow an internal part for fluid to flow around*

Aspirational:
- Temperature and density changes
- Visulaization of flow (image processing in c woo)
- Non-steady state solutions
- Python wrapper integration?
- 3d