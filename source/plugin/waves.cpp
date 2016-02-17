/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Wave equation
 *
 ******************************************************************************/

#include "levelset.h"
#include "commonkernels.h"
#include "particle.h"
#include "conjugategrad.h"
#include <cmath>

using namespace std;

namespace Manta
{


/******************************************************************************
 *
 * explicit integration
 *
 ******************************************************************************/


 struct knCalcSecDeriv2d : public KernelBase { knCalcSecDeriv2d(const Grid<Real>& v, Grid<Real>& ret) :  KernelBase(&v,1) ,v(v),ret(ret)   { run(); }  inline void op(int i, int j, int k, const Grid<Real>& v, Grid<Real>& ret )  {

    ret(i,j,k) = 
		( -4. * v(i,j,k) + v(i-1,j,k) + v(i+1,j,k) + v(i,j-1,k) + v(i,j+1,k) );

}   inline const Grid<Real>& getArg0() { return v; } typedef Grid<Real> type0;inline Grid<Real>& getArg1() { return ret; } typedef Grid<Real> type1; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,v,ret);  } } else { const int k=0; 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,v,ret);  } }  } const Grid<Real>& v; Grid<Real>& ret;   };

void calcSecDeriv2d(const Grid<Real>& v, Grid<Real>& curv) {
	knCalcSecDeriv2d(v,curv);
}

// mass conservation 
 struct knTotalSum : public KernelBase { knTotalSum(Grid<Real>& h) :  KernelBase(&h,1) ,h(h) ,sum(0)  { run(); }  inline void op(int i, int j, int k, Grid<Real>& h ,double& sum)  { sum += h(i,j,k); }   inline Grid<Real>& getArg0() { return h; } typedef Grid<Real> type0; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  double sum = 0; 
#pragma omp for nowait 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,h,sum); 
#pragma omp critical
{this->sum += sum; } } } else { const int k=0; 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  double sum = 0; 
#pragma omp for nowait 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,h,sum); 
#pragma omp critical
{this->sum += sum; } } }  } Grid<Real>& h;  double sum;  };


Real totalSum(Grid<Real>& height) {
	knTotalSum ts(height);
	return ts.sum;
}

void normalizeSumTo(Grid<Real>& height, Real target) {
	knTotalSum ts(height);
	Real factor = target / ts.sum;
	height.multConst(factor);
}

/******************************************************************************
 *
 * implicit time integration
 *
 ******************************************************************************/

//! Kernel: Construct the right-hand side of the poisson equation
 struct MakeRhsWE : public KernelBase { MakeRhsWE(FlagGrid& flags, Grid<Real>& rhs, Grid<Real>& ut, Grid<Real>& utm1, Real s, bool crankNic=false) :  KernelBase(&flags,1) ,flags(flags),rhs(rhs),ut(ut),utm1(utm1),s(s),crankNic(crankNic)   { run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& rhs, Grid<Real>& ut, Grid<Real>& utm1, Real s, bool crankNic=false )  {
	rhs(i,j,k) = ( 2.*ut(i,j,k) - utm1(i,j,k) );
	if(crankNic) {
		rhs(i,j,k) += s * ( -4.*ut(i,j,k) + 1.*ut(i-1,j,k) + 1.*ut(i+1,j,k) + 1.*ut(i,j-1,k) + 1.*ut(i,j+1,k) );
	} 
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return rhs; } typedef Grid<Real> type1;inline Grid<Real>& getArg2() { return ut; } typedef Grid<Real> type2;inline Grid<Real>& getArg3() { return utm1; } typedef Grid<Real> type3;inline Real& getArg4() { return s; } typedef Real type4;inline bool& getArg5() { return crankNic; } typedef bool type5; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,rhs,ut,utm1,s,crankNic);  } } else { const int k=0; 
#pragma omp parallel 
 { this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for 
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,rhs,ut,utm1,s,crankNic);  } }  } FlagGrid& flags; Grid<Real>& rhs; Grid<Real>& ut; Grid<Real>& utm1; Real s; bool crankNic;   };

//! do a CG solve (note, out grid only there for debugging... could be removed)
void cgSolveWE(FlagGrid& flags, Grid<Real>& ut, Grid<Real>& utm1, Grid<Real>& out, bool crankNic = false, Real cSqr = 0.25, Real cgMaxIterFac = 1.5, Real cgAccuracy = 1e-5 ) {
	// reserve temp grids
	FluidSolver* parent = flags.getParent();
	Grid<Real> rhs(parent);
	Grid<Real> residual(parent);
	Grid<Real> search(parent);
	Grid<Real> A0(parent);
	Grid<Real> Ai(parent);
	Grid<Real> Aj(parent);
	Grid<Real> Ak(parent);
	Grid<Real> tmp(parent);
	//Grid<Real> pca0(parent);
	//Grid<Real> pca1(parent);
	//Grid<Real> pca2(parent);
	//Grid<Real> pca3(parent);
	// solution...
	//Grid<Real> pressure(parent);
	out.clear();
		
	// setup matrix and boundaries
	MakeLaplaceMatrix (flags, A0, Ai, Aj, Ak);
	Real dt   = parent->getDt();
	Real s    = dt*dt*cSqr * 0.5;
	FOR_IJK(flags) {
		Ai(i,j,k) *= s;
		Aj(i,j,k) *= s;
		Ak(i,j,k) *= s;
		A0(i,j,k) *= s;
		A0(i,j,k) += 1.;
	}
	
	// compute divergence and init right hand side
	rhs.clear();
	// h=dt
	// rhs:   = 2 ut - ut-1 
   	// A:    (h2 c2/ dx)=s   ,  (1+4s)uij + s ui-1j + ...
	// Cr.Nic.
	// rhs:  cr nic = 2 ut - ut-1 + h^2c^2/2 b 
   	// A:    (h2 c2/2 dx)=s   ,  (1+4s)uij + s ui-1j + ...
	MakeRhsWE kernMakeRhs(flags, rhs, ut,utm1, s, crankNic);
	
	const int maxIter = (int)(cgMaxIterFac * flags.getSize().max()) * (flags.is3D() ? 1 : 4);
	GridCgInterface *gcg;
	if (flags.is3D())
		gcg = new GridCg<ApplyMatrix  >(out, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
	else
		gcg = new GridCg<ApplyMatrix2D>(out, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
	
	gcg->setAccuracy( cgAccuracy ); 

	// optional preconditioning
	//gcg->setPreconditioner( precondition ? GridCgInterface::PC_mICP : GridCgInterface::PC_None, &pca0, &pca1, &pca2, &pca3);

	for (int iter=0; iter<maxIter; iter++) {
		if (!gcg->iterate()) iter=maxIter;
	} 
	debMsg("FluidSolver::solvePressure iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma(), 1);

	utm1.swap( ut );
	ut.copyFrom( out );

	delete gcg;
}

} //namespace
