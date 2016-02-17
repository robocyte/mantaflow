/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Functions for calculating wavelet turbulence,
 * plus helpers to compute vorticity, and strain rate magnitude
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "shapes.h"
#include "commonkernels.h"
#include "noisefield.h"

using namespace std;

namespace Manta
{

//*****************************************************************************

// first some fairly generic interpolation functions for grids with multiple sizes


void interpolateGrid( Grid<Real>& target, Grid<Real>& source , Vec3 scale=Vec3(1.), Vec3 offset=Vec3(0.), int orderSpace=1 ) {
	Vec3 sourceFactor = calcGridSizeFactor( source.getSize(), target.getSize() );

	// a brief note on a mantaflow specialty: the target grid has to be the first argument here!
	// the parent fluidsolver object is taken from the first grid, and it determines the size of the
	// loop for the kernel call. as we're writing into target, it's important to loop exactly over
	// all cells of the target grid... (note, when calling the plugin in python, it doesnt matter anymore).

	// sourceFactor offset necessary to shift eval points by half a small cell width
	knInterpolateGridTempl<Real>(target, source, sourceFactor*scale, sourceFactor*0.5 + offset, orderSpace);
}

void interpolateGridVec3( Grid<Vec3>& target, Grid<Vec3>& source , Vec3 scale=Vec3(1.), Vec3 offset=Vec3(0.), int orderSpace=1 ) {
	Vec3 sourceFactor = calcGridSizeFactor( source.getSize(), target.getSize() );
	knInterpolateGridTempl<Vec3>(target, source, sourceFactor*scale, sourceFactor*0.5 + offset, orderSpace); 
}

//!interpolate a mac velocity grid from one size to another size
 struct KnInterpolateMACGrid : public KernelBase { KnInterpolateMACGrid(MACGrid& target, MACGrid& source, const Vec3& sourceFactor, const Vec3& off, int orderSpace) :  KernelBase(&target,0) ,target(target),source(source),sourceFactor(sourceFactor),off(off),orderSpace(orderSpace)   { run(); }  inline void op(int i, int j, int k, MACGrid& target, MACGrid& source, const Vec3& sourceFactor, const Vec3& off, int orderSpace )  {
	Vec3 pos = Vec3(i,j,k) * sourceFactor + off;

	Real vx = source.getInterpolatedHi(pos - Vec3(0.5,0,0), orderSpace)[0];
	Real vy = source.getInterpolatedHi(pos - Vec3(0,0.5,0), orderSpace)[1];
	Real vz = 0.f;
	if(source.is3D()) vz = source.getInterpolatedHi(pos - Vec3(0,0,0.5), orderSpace)[2];

	target(i,j,k) = Vec3(vx,vy,vz);
}   inline MACGrid& getArg0() { return target; } typedef MACGrid type0;inline MACGrid& getArg1() { return source; } typedef MACGrid type1;inline const Vec3& getArg2() { return sourceFactor; } typedef Vec3 type2;inline const Vec3& getArg3() { return off; } typedef Vec3 type3;inline int& getArg4() { return orderSpace; } typedef int type4; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,target,source,sourceFactor,off,orderSpace);  } } else { const int k=0; 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,target,source,sourceFactor,off,orderSpace);  } }  } MACGrid& target; MACGrid& source; const Vec3& sourceFactor; const Vec3& off; int orderSpace;   };

void interpolateMACGrid(MACGrid& target, MACGrid& source, Vec3 scale=Vec3(1.), Vec3 offset=Vec3(0.), int orderSpace=1) {
	Vec3 sourceFactor = calcGridSizeFactor( source.getSize(), target.getSize() );

	// see interpolateGrid for why the target grid needs to come first in the parameters!  
	KnInterpolateMACGrid(target, source, sourceFactor*scale, sourceFactor*0.5 + offset, orderSpace);
}

//*****************************************************************************

//! Apply vector noise to grid, this is a simplified version - no position scaling or UVs



 struct knApplySimpleNoiseVec : public KernelBase { knApplySimpleNoiseVec(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, Real scale, Grid<Real>* weight ) :  KernelBase(&flags,0) ,flags(flags),target(target),noise(noise),scale(scale),weight(weight)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, Real scale, Grid<Real>* weight  )  {
	if ( !flags.isFluid(i,j,k) ) return; 
	Real factor = 1;
	if(weight) factor = (*weight)(i,j,k);
	target(i,j,k) += noise.evaluateCurl( Vec3(i,j,k) ) * scale * factor;
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Vec3>& getArg1() { return target; } typedef Grid<Vec3> type1;inline WaveletNoiseField& getArg2() { return noise; } typedef WaveletNoiseField type2;inline Real& getArg3() { return scale; } typedef Real type3;inline Grid<Real>* getArg4() { return weight; } typedef Grid<Real> type4; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,target,noise,scale,weight);  } } else { const int k=0; 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,target,noise,scale,weight);  } }  } FlagGrid& flags; Grid<Vec3>& target; WaveletNoiseField& noise; Real scale; Grid<Real>* weight;   };

void applySimpleNoiseVec3(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, Real scale=1.0 , Grid<Real>* weight=NULL ) {
	// note - passing a MAC grid here is slightly inaccurate, we should evaluate each component separately
	knApplySimpleNoiseVec(flags, target, noise, scale , weight );
}

//! Simple noise for a real grid , follows applySimpleNoiseVec3
 struct knApplySimpleNoiseReal : public KernelBase { knApplySimpleNoiseReal(FlagGrid& flags, Grid<Real>& target, WaveletNoiseField& noise, Real scale, Grid<Real>* weight ) :  KernelBase(&flags,0) ,flags(flags),target(target),noise(noise),scale(scale),weight(weight)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& target, WaveletNoiseField& noise, Real scale, Grid<Real>* weight  )  {
	if ( !flags.isFluid(i,j,k) ) return; 
	Real factor = 1;
	if(weight) factor = (*weight)(i,j,k);
	target(i,j,k) += noise.evaluate( Vec3(i,j,k) ) * scale * factor;
}   inline FlagGrid& getArg0()
{
    return flags;
} typedef FlagGrid type0;inline Grid<Real>& getArg1() { return target; } typedef Grid<Real> type1;inline WaveletNoiseField& getArg2() { return noise; } typedef WaveletNoiseField type2;inline Real& getArg3() { return scale; } typedef Real type3;inline Grid<Real>* getArg4() { return weight; } typedef Grid<Real> type4; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,target,noise,scale,weight);  } } else { const int k=0; 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,target,noise,scale,weight);  } }  } FlagGrid& flags; Grid<Real>& target; WaveletNoiseField& noise; Real scale; Grid<Real>* weight;   };

void applySimpleNoiseReal(FlagGrid& flags, Grid<Real>& target, WaveletNoiseField& noise, Real scale=1.0 , Grid<Real>* weight=NULL ) {
	knApplySimpleNoiseReal(flags, target, noise, scale , weight );
}

//! Apply vector-based wavelet noise to target grid
//! This is the version with more functionality - supports uv grids, and on-the-fly interpolation
//! of input grids.
 struct knApplyNoiseVec : public KernelBase { knApplyNoiseVec(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, Real scale, Real scaleSpatial, Grid<Real>* weight, Grid<Vec3>* uv, bool uvInterpol, const Vec3& sourceFactor ) :  KernelBase(&flags,0) ,flags(flags),target(target),noise(noise),scale(scale),scaleSpatial(scaleSpatial),weight(weight),uv(uv),uvInterpol(uvInterpol),sourceFactor(sourceFactor)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, Real scale, Real scaleSpatial, Grid<Real>* weight, Grid<Vec3>* uv, bool uvInterpol, const Vec3& sourceFactor  )  {
	if ( !flags.isFluid(i,j,k) ) return;

	// get weighting, interpolate if necessary
	Real w = 1;
	if(weight) {
		if(!uvInterpol) {
			w = (*weight)(i,j,k);
		} else {
			w = weight->getInterpolated( Vec3(i,j,k) * sourceFactor );
		}
	}

	// compute position where to evaluate the noise
	Vec3 pos = Vec3(i,j,k);
	if(uv) {
		if(!uvInterpol) {
			pos = (*uv)(i,j,k);
		} else {
			pos = uv->getInterpolated( Vec3(i,j,k) * sourceFactor );
			// uv coordinates are in local space - so we need to adjust the values of the positions
			pos /= sourceFactor;
		}
	}
	pos *= scaleSpatial;

	Vec3 noiseVec3 = noise.evaluateCurl( pos ) * scale * w; 
	//noiseVec3=pos; // debug , show interpolated positions
	target(i,j,k) += noiseVec3;
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Vec3>& getArg1() { return target; } typedef Grid<Vec3> type1;inline WaveletNoiseField& getArg2() { return noise; } typedef WaveletNoiseField type2;inline Real& getArg3() { return scale; } typedef Real type3;inline Real& getArg4() { return scaleSpatial; } typedef Real type4;inline Grid<Real>* getArg5() { return weight; } typedef Grid<Real> type5;inline Grid<Vec3>* getArg6() { return uv; } typedef Grid<Vec3> type6;inline bool& getArg7() { return uvInterpol; } typedef bool type7;inline const Vec3& getArg8() { return sourceFactor; } typedef Vec3 type8; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,target,noise,scale,scaleSpatial,weight,uv,uvInterpol,sourceFactor);  } } else { const int k=0; 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,target,noise,scale,scaleSpatial,weight,uv,uvInterpol,sourceFactor);  } }  } FlagGrid& flags; Grid<Vec3>& target; WaveletNoiseField& noise; Real scale; Real scaleSpatial; Grid<Real>* weight; Grid<Vec3>* uv; bool uvInterpol; const Vec3& sourceFactor;   };

void applyNoiseVec3(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, Real scale=1.0 , Real scaleSpatial=1.0 , Grid<Real>* weight=NULL , Grid<Vec3>* uv=NULL ) {
	// check whether the uv grid has a different resolution
	bool uvInterpol = false; 
	// and pre-compute conversion (only used if uvInterpol==true)
	// used for both uv and weight grid...
	Vec3 sourceFactor = Vec3(1.);
	if(uv) {
		uvInterpol = (target.getSize() != uv->getSize());
		sourceFactor = calcGridSizeFactor( uv->getSize(), target.getSize() );
	} else if(weight) {
		uvInterpol = (target.getSize() != weight->getSize());
		sourceFactor = calcGridSizeFactor( weight->getSize(), target.getSize() );
	}
	if(uv && weight) assertMsg( uv->getSize() == weight->getSize(), "UV and weight grid have to match!");

	// note - passing a MAC grid here is slightly inaccurate, we should evaluate each component separately
	knApplyNoiseVec(flags, target, noise, scale, scaleSpatial, weight , uv,uvInterpol,sourceFactor );
}

//! Compute energy of a staggered velocity field (at cell center)
 struct KnApplyComputeEnergy : public KernelBase { KnApplyComputeEnergy( FlagGrid& flags, MACGrid& vel, Grid<Real>& energy ) :  KernelBase(&flags,0) ,flags(flags),vel(vel),energy(energy)   { run(); }  inline void op(int i, int j, int k,  FlagGrid& flags, MACGrid& vel, Grid<Real>& energy  )  {
	Real e = 0.f;
	if ( flags.isFluid(i,j,k) ) {
		Vec3 v = vel.getCentered(i,j,k);
		e = 0.5 * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	}
	energy(i,j,k) = e;
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline MACGrid& getArg1() { return vel; } typedef MACGrid type1;inline Grid<Real>& getArg2() { return energy; } typedef Grid<Real> type2; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,vel,energy);  } } else { const int k=0; 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,vel,energy);  } }  } FlagGrid& flags; MACGrid& vel; Grid<Real>& energy;   };

void computeEnergy( FlagGrid& flags, MACGrid& vel, Grid<Real>& energy ) {
	KnApplyComputeEnergy( flags, vel, energy );
}

void computeWaveletCoeffs(Grid<Real>& input)
{
	Grid<Real> temp1(input.getParent()), temp2(input.getParent());
	WaveletNoiseField::computeCoefficients(input, temp1, temp2);
}

// note - alomst the same as for vorticity confinement
void computeVorticity(MACGrid& vel, Grid<Vec3>& vorticity, Grid<Real>* norm)
{
	Grid<Vec3> velCenter(vel.getParent());
	GetCentered(velCenter, vel);
	CurlOp(velCenter, vorticity);
	if(norm) GridNorm( *norm, vorticity);
}

// note - very similar to KnComputeProductionStrain, but for use as wavelet turb weighting
 struct KnComputeStrainRateMag : public KernelBase { KnComputeStrainRateMag(const MACGrid& vel, const Grid<Vec3>& velCenter, Grid<Real>& prod ) :  KernelBase(&vel,1) ,vel(vel),velCenter(velCenter),prod(prod)   { run(); }  inline void op(int i, int j, int k, const MACGrid& vel, const Grid<Vec3>& velCenter, Grid<Real>& prod  )  {
	// compute Sij = 1/2 * (dU_i/dx_j + dU_j/dx_i)
	Vec3 diag = Vec3(vel(i+1,j,k).x, vel(i,j+1,k).y, 0. ) - vel(i,j,k);
	if(vel.is3D()) diag[2] += vel(i,j,k+1).z;
	else           diag[2]  = 0.;

	Vec3 ux =         0.5*(velCenter(i+1,j,k)-velCenter(i-1,j,k));
	Vec3 uy =         0.5*(velCenter(i,j+1,k)-velCenter(i,j-1,k));
	Vec3 uz;
	if(vel.is3D()) uz=0.5*(velCenter(i,j,k+1)-velCenter(i,j,k-1));

	Real S12 = 0.5*(ux.y+uy.x);
	Real S13 = 0.5*(ux.z+uz.x);
	Real S23 = 0.5*(uy.z+uz.y);
	Real S2 = square(diag.x) + square(diag.y) + square(diag.z) +
		2.0*square(S12) + 2.0*square(S13) + 2.0*square(S23);
	prod(i,j,k) = S2;
}   inline const MACGrid& getArg0() { return vel; } typedef MACGrid type0;inline const Grid<Vec3>& getArg1() { return velCenter; } typedef Grid<Vec3> type1;inline Grid<Real>& getArg2() { return prod; } typedef Grid<Real> type2; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,vel,velCenter,prod);  } } else { const int k=0; 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,vel,velCenter,prod);  } }  } const MACGrid& vel; const Grid<Vec3>& velCenter; Grid<Real>& prod;   };

void computeStrainRateMag(MACGrid& vel, Grid<Real>& mag) {
	Grid<Vec3> velCenter(vel.getParent());
	GetCentered(velCenter, vel);
	KnComputeStrainRateMag(vel, velCenter, mag);
}

// extrapolate a real grid into a flagged region (based on initial flags)
// by default extrapolates from fluid to obstacle cells
template<class T> 
void extrapolSimpleFlagsHelper (FlagGrid& flags, Grid<T>& val, int distance = 4, 
									int flagFrom=FlagGrid::TypeFluid, int flagTo=FlagGrid::TypeObstacle ) 
{
	Grid<int> tmp( flags.getParent() );
	int dim = (flags.is3D() ? 3:2);
	const Vec3i nb[6] = { 
		Vec3i(1 ,0,0), Vec3i(-1,0,0),
		Vec3i(0,1 ,0), Vec3i(0,-1,0),
		Vec3i(0,0,1 ), Vec3i(0,0,-1) };

	// remove all fluid cells (set to 1)
	tmp.clear();
	bool foundTarget = false;
	FOR_IJK_BND(flags,0) {
		if (flags(i,j,k) & flagFrom) 
			tmp( Vec3i(i,j,k) ) = 1;
		if (!foundTarget && (flags(i,j,k) & flagTo)) foundTarget=true;
	}
	// optimization, skip extrapolation if we dont have any cells to extrapolate to
	if(!foundTarget) {
		debMsg("No target cells found, skipping extrapolation", 1);
		return;
	}

	// extrapolate for given distance
	for(int d=1; d<1+distance; ++d) {

		// TODO, parallelize
		FOR_IJK_BND(flags,1) {
			if (tmp(i,j,k) != 0)          continue;
			if (!(flags(i,j,k) & flagTo)) continue;

			// copy from initialized neighbors
			Vec3i p(i,j,k);
			int nbs = 0;
			T avgVal = 0.;
			for (int n=0; n<2*dim; ++n) {
				if (tmp(p+nb[n]) == d) {
					avgVal += val(p+nb[n]);
					nbs++;
				}
			}

			if(nbs>0) {
				tmp(p) = d+1;
				val(p) = avgVal / nbs;
			}
		}

	} // distance 
}

void extrapolateSimpleFlags(FlagGrid& flags, GridBase* val, int distance = 4, int flagFrom=FlagGrid::TypeFluid, int flagTo=FlagGrid::TypeObstacle ) {
	if (val->getType() & GridBase::TypeReal) {
		extrapolSimpleFlagsHelper<Real>(flags,*((Grid<Real>*) val),distance,flagFrom,flagTo);
	}
	else if (val->getType() & GridBase::TypeInt) {    
		extrapolSimpleFlagsHelper<int >(flags,*((Grid<int >*) val),distance,flagFrom,flagTo);
	}
	else if (val->getType() & GridBase::TypeVec3) {    
		extrapolSimpleFlagsHelper<Vec3>(flags,*((Grid<Vec3>*) val),distance,flagFrom,flagTo);
	}
	else
		errMsg("extrapolateSimpleFlags: Grid Type is not supported (only int, Real, Vec3)");    
}

} // namespace


