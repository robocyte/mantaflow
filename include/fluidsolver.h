/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Main class for the fluid solver
 *
 ******************************************************************************/

#pragma once

#include "manta.h"
#include "vectorbase.h"
#include <vector>
#include <map>

namespace Manta
{ 
	
//! Encodes grid size, timstep etc.
class FluidSolver : public PbClass
{
public:
	FluidSolver(Vec3i gridSize, int dim=3);
	virtual ~FluidSolver();
	
	// accessors
	Vec3i getGridSize()         { return mGridSize; }
    inline Real  getDt()        { return mDt; }
	inline Real  getDx()        { return 1.0 / mGridSize.max(); }
	inline Real  getTime()      { return mTimeTotal; }

	//! Check dimensionality
	inline bool is2D() const    { return mDim==2; }
	inline bool is3D() const    { return mDim==3; }
	
	void printMemInfo();
	
	//! Advance the solver one timestep, update GUI if present
	void step();
	
	//! Update the timestep size based on given maximal velocity magnitude 
	void adaptTimestep(Real maxVel);
	
    // temp grid and plugin stuff: you shouldn't call this manually
	template<class T> T* getGridPointer();
	template<class T> void freeGridPointer(T* ptr);    

	//! expose animation time
    Real getTimeStep() const            { return mDt; }
    Real getTimeTotal() const           { return mTimeTotal; }
    int  getFrame() const               { return mFrame; }
    void setTimeStep(Real timestep)     { mDt = timestep; }
    void setTimeTotal(Real timetotal)   { mTimeTotal = timetotal; }
    void setFrame(int frame)            { mFrame = frame; }

    //! parameters for adaptive time stepping
	Real mCflCond;
	Real mDtMin;
	Real mDtMax;
	Real mFrameLength;

protected:
	Vec3i     mGridSize;
	const int mDim;
	Real      mTimePerFrame;
	bool      mLockDt;
	bool      mAdaptDt;

	Real      mDt;
	Real      mTimeTotal;
	int       mFrame;

	// subclass for managing grid memory, stored as a stack to allow fast allocation
	template<class T>
    struct GridStorage
    {
		GridStorage()
            : used(0)
        { }

		T*   get(Vec3i size);
		void free();
		void release(T* ptr);
		
		std::vector<T*> grids;
		int             used;
	};
	
	GridStorage<int>  mGridsInt;
	GridStorage<Real> mGridsReal;
    GridStorage<Vec3> mGridsVec;
};

} //namespace
