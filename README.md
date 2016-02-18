# Mantaflow #

This is a clone of the [Mantaflow framework](http://mantaflow.com) written by Tobias Pfaff and Nils Thuerey.
My goal was to remove all Python bindings to have easy access to the functionality from C++ alone. For testing purposes, the original Python scripts are rewritten as small console applications.

###(incomplete) list of changes:
* removed python wrapper and code generation
* changed directory organization
* replaced include guards with #pragma once directive
* all kernels are expanded to the OpenMP version for now
* Python script files rewritten as console applications
* functions moved to general.h/general.cpp:
	- setDebugLevel
	- buildInfoString
	- printInfoString
* created headers to expose plugin functionality:
	- advection.h
	- extforces.h
	- flip.h
	- initplugins.h
	- pressure.h
* functions exposed in fastmarch.h:
	- extrapolateMACSimple
	- extrapolateMACFromWeight
	- extrapolateLsSimple
	- extrapolateVec3Simple
* class Shape:
	- templated member applyToGrid
* class ParticleDataImpl:
	- additional constructor taking a particle system as an argument
* a gazillion completely unnecessary, OCD-related style changes...
