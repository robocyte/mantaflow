/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugins for pressure correction: solve_pressure, and ghost fluid helpers
 *
 ******************************************************************************/

#include "pressure.h"
#include "kernel.h"

using namespace std;

namespace Manta
{

//! Kernel: Construct the right-hand side of the poisson equation
struct MakeRhs : public KernelBase
{
    MakeRhs(FlagGrid& flags, Grid<Real>& rhs, MACGrid& vel, Grid<Real>* perCellCorr, MACGrid* fractions)
        : KernelBase(&flags,1)
        , flags(flags)
        , rhs(rhs)
        , vel(vel)
        , perCellCorr(perCellCorr)
        , fractions(fractions)
        , cnt(0),sum(0)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& rhs, MACGrid& vel, Grid<Real>* perCellCorr, MACGrid* fractions ,int& cnt,double& sum)
    {
	    if (!flags.isFluid(i,j,k))
        {
		    rhs(i,j,k) = 0;
		    return;
	    }

	    // compute divergence 
	    // no flag checks: assumes vel at obstacle interfaces is set to zero
	    Real set(0);
	    if(!fractions)
        {
		    set =               vel(i,j,k).x - vel(i+1,j,k).x + 
				 			    vel(i,j,k).y - vel(i,j+1,k).y; 
		    if(vel.is3D()) set+=vel(i,j,k).z - vel(i,j,k+1).z;
	    } else
        {
		    set =               (*fractions)(i,j,k).x * vel(i,j,k).x - (*fractions)(i+1,j,k).x * vel(i+1,j,k).x + 
							    (*fractions)(i,j,k).y * vel(i,j,k).y - (*fractions)(i,j+1,k).y * vel(i,j+1,k).y; 
		    if(vel.is3D()) set+=(*fractions)(i,j,k).z * vel(i,j,k).z - (*fractions)(i,j,k+1).z * vel(i,j,k+1).z;
	    }
	
	    // per cell divergence correction (optional)
	    if(perCellCorr) set += perCellCorr->get(i,j,k);
	
	    // obtain sum, cell count
	    sum += set;
	    cnt++;
	
	    rhs(i,j,k) = set;
    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline Grid<Real>& getArg1() { return rhs; }
    typedef Grid<Real> type1;
    inline MACGrid& getArg2() { return vel; }
    typedef MACGrid type2;
    inline Grid<Real>* getArg3() { return perCellCorr; }
    typedef Grid<Real> type3;
    inline MACGrid* getArg4() { return fractions; }
    typedef MACGrid type4;
    
    void run()
    {
        const int _maxX = maxX;
        const int _maxY = maxY;
        if (maxZ > 1)
        {
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
                int cnt = 0;
                double sum = 0;
#pragma omp for nowait
                for (int k=minZ; k < maxZ; k++)
                    for (int j=1; j < _maxY; j++)
                        for (int i=1; i < _maxX; i++)
                            op(i,j,k,flags,rhs,vel,perCellCorr,fractions,cnt,sum);
#pragma omp critical
                {
                    this->cnt += cnt; this->sum += sum;
                }
            }
        } else
        {
            const int k=0;
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
                int cnt = 0;
                double sum = 0;
#pragma omp for nowait
                for (int j=1; j < _maxY; j++)
                    for (int i=1; i < _maxX; i++)
                        op(i,j,k,flags,rhs,vel,perCellCorr,fractions,cnt,sum);
#pragma omp critical
                {
                    this->cnt += cnt; this->sum += sum;
                }
            }
        }
    }
    
    FlagGrid& flags;
    Grid<Real>& rhs;
    MACGrid& vel;
    Grid<Real>* perCellCorr;
    MACGrid* fractions;
    int cnt;
    double sum;
};

//! Kernel: Apply velocity update from poisson equation
struct CorrectVelocity : public KernelBase
{
    CorrectVelocity(FlagGrid& flags, MACGrid& vel, Grid<Real>& pressure)
        : KernelBase(&flags,1)
        , flags(flags)
        , vel(vel)
        , pressure(pressure)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, MACGrid& vel, Grid<Real>& pressure)
    {
	    int idx = flags.index(i,j,k);
	    if (flags.isFluid(idx))
	    {
		    if (flags.isFluid(i-1,j,k)) vel[idx].x -= (pressure[idx] - pressure(i-1,j,k));
		    if (flags.isFluid(i,j-1,k)) vel[idx].y -= (pressure[idx] - pressure(i,j-1,k));
		    if (flags.is3D() && flags.isFluid(i,j,k-1)) vel[idx].z -= (pressure[idx] - pressure(i,j,k-1));
 
		    if (flags.isEmpty(i-1,j,k)) vel[idx].x -= pressure[idx];
		    if (flags.isEmpty(i,j-1,k)) vel[idx].y -= pressure[idx];
		    if (flags.is3D() && flags.isEmpty(i,j,k-1)) vel[idx].z -= pressure[idx];
	    } else if (flags.isEmpty(idx)&&!flags.isOutflow(idx)) // don't change velocities in outflow cells
	    {
		    if (flags.isFluid(i-1,j,k)) vel[idx].x += pressure(i-1,j,k);
		    else                        vel[idx].x  = 0.f;
		    if (flags.isFluid(i,j-1,k)) vel[idx].y += pressure(i,j-1,k);
		    else                        vel[idx].y  = 0.f;
		    if (flags.is3D() ) {
		    if (flags.isFluid(i,j,k-1)) vel[idx].z += pressure(i,j,k-1);
		    else                        vel[idx].z  = 0.f;
		    }
	    }
    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline MACGrid& getArg1() { return vel; }
    typedef MACGrid type1;
    inline Grid<Real>& getArg2() { return pressure; }
    typedef Grid<Real> type2;
    
    void run()
    {
        const int _maxX = maxX;
        const int _maxY = maxY;
        if (maxZ > 1)
        { 
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int k=minZ; k < maxZ; k++)
                    for (int j=1; j < _maxY; j++)
                        for (int i=1; i < _maxX; i++)
                            op(i,j,k,flags,vel,pressure);
            }
        } else
        {
            const int k=0;
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int j=1; j < _maxY; j++)
                    for (int i=1; i < _maxX; i++)
                        op(i,j,k,flags,vel,pressure);
            }
        }
    }
    
    FlagGrid& flags;
    MACGrid& vel;
    Grid<Real>& pressure;
};

// *****************************************************************************
// Ghost fluid helpers

// calculate fraction filled with liquid (note, assumes inside value is < outside!)
inline static Real thetaHelper(Real inside, Real outside)
{
	Real denom = inside-outside;
	if (denom > -1e-04) return 0.5; // should always be neg, and large enough...
	return std::max(Real(0), std::min(Real(1), inside/denom));
}

// calculate ghost fluid factor, cell at idx should be a fluid cell
inline static Real ghostFluidHelper(int idx, int offset, const Grid<Real> &phi, Real gfClamp)
{
	Real alpha = thetaHelper(phi[idx], phi[idx+offset]);
	if (alpha < gfClamp) return alpha = gfClamp;
	return (1-(1/alpha)); 
}

//! Kernel: Adapt A0 for ghost fluid
struct ApplyGhostFluidDiagonal : public KernelBase
{
    ApplyGhostFluidDiagonal(Grid<Real> &A0, const FlagGrid &flags, const Grid<Real> &phi, Real gfClamp)
        : KernelBase(&A0,1)
        , A0(A0)
        , flags(flags)
        , phi(phi)
        , gfClamp(gfClamp)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<Real> &A0, const FlagGrid &flags, const Grid<Real> &phi, Real gfClamp)
    {
	    const int X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	    int idx = flags.index(i,j,k);
	    if (!flags.isFluid(idx)) return;

	    if (flags.isEmpty(i-1,j,k)) A0[idx] -= ghostFluidHelper(idx, -X, phi, gfClamp);
	    if (flags.isEmpty(i+1,j,k)) A0[idx] -= ghostFluidHelper(idx, +X, phi, gfClamp);
	    if (flags.isEmpty(i,j-1,k)) A0[idx] -= ghostFluidHelper(idx, -Y, phi, gfClamp);
	    if (flags.isEmpty(i,j+1,k)) A0[idx] -= ghostFluidHelper(idx, +Y, phi, gfClamp);
	    if (flags.is3D())
        {
		    if (flags.isEmpty(i,j,k-1)) A0[idx] -= ghostFluidHelper(idx, -Z, phi, gfClamp);
		    if (flags.isEmpty(i,j,k+1)) A0[idx] -= ghostFluidHelper(idx, +Z, phi, gfClamp);
	    }
    }
    
    inline Grid<Real> & getArg0() { return A0; }
    typedef Grid<Real> type0;
    inline const FlagGrid& getArg1() { return flags; }
    typedef FlagGrid type1;
    inline const Grid<Real> & getArg2() { return phi; }
    typedef Grid<Real> type2;
    inline Real& getArg3() { return gfClamp; }
    typedef Real type3;
    
    void run()
    {
        const int _maxX = maxX;
        const int _maxY = maxY;
        if (maxZ > 1)
        {
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int k=minZ; k < maxZ; k++)
                    for (int j=1; j < _maxY; j++)
                        for (int i=1; i < _maxX; i++)
                            op(i,j,k,A0,flags,phi,gfClamp);
            }
        } else
        {
            const int k=0;
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int j=1; j < _maxY; j++)
                    for (int i=1; i < _maxX; i++)
                        op(i,j,k,A0,flags,phi,gfClamp);
            }
        }
    }
    
    Grid<Real> & A0;
    const FlagGrid& flags;
    const Grid<Real> & phi;
    Real gfClamp;
};

//! Kernel: Apply velocity update: ghost fluid contribution
struct CorrectVelocityGhostFluid : public KernelBase
{
    CorrectVelocityGhostFluid(MACGrid &vel, const FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp)
        : KernelBase(&vel,1)
        , vel(vel)
        , flags(flags)
        , pressure(pressure)
        , phi(phi)
        , gfClamp(gfClamp)
    {
        run();
    }
    
    inline void op(int i, int j, int k, MACGrid &vel, const FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp)
    {
	    const int X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	    const int idx = flags.index(i,j,k);

	    if (flags.isFluid(idx))
	    {
		    if (flags.isEmpty(i-1,j,k)) vel[idx][0] += pressure[idx] * ghostFluidHelper(idx, -X, phi, gfClamp);
		    if (flags.isEmpty(i,j-1,k)) vel[idx][1] += pressure[idx] * ghostFluidHelper(idx, -Y, phi, gfClamp);
		    if (flags.is3D() && flags.isEmpty(i,j,k-1)) vel[idx][2] += pressure[idx] * ghostFluidHelper(idx, -Z, phi, gfClamp);
	    } else if (flags.isEmpty(idx)&&!flags.isOutflow(idx)) // do not change velocities in outflow cells
	    {
		    if (flags.isFluid(i-1,j,k)) vel[idx][0] -= pressure(i-1,j,k) * ghostFluidHelper(idx-X, +X, phi, gfClamp);
		    else                        vel[idx].x  = 0.f;
		    if (flags.isFluid(i,j-1,k)) vel[idx][1] -= pressure(i,j-1,k) * ghostFluidHelper(idx-Y, +Y, phi, gfClamp);
		    else                        vel[idx].y  = 0.f;
		    if (flags.is3D() )
            {
		        if (flags.isFluid(i,j,k-1)) vel[idx][2] -= pressure(i,j,k-1) * ghostFluidHelper(idx-Z, +Z, phi, gfClamp);
		        else                        vel[idx].z  = 0.f;
		    }
	    }
    }
    
    inline MACGrid& getArg0() { return vel; }
    typedef MACGrid type0;
    inline const FlagGrid& getArg1() { return flags; }
    typedef FlagGrid type1;
    inline const Grid<Real> & getArg2() { return pressure; }
    typedef Grid<Real> type2;
    inline const Grid<Real> & getArg3() { return phi; }
    typedef Grid<Real> type3;
    inline Real& getArg4() { return gfClamp; }
    typedef Real type4;
    
    void run()
    {
        const int _maxX = maxX;
        const int _maxY = maxY; if (maxZ > 1)
        {
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int k=minZ; k < maxZ; k++)
                    for (int j=1; j < _maxY; j++)
                        for (int i=1; i < _maxX; i++)
                            op(i,j,k,vel,flags,pressure,phi,gfClamp);
            }
        } else
        {
            const int k=0;
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int j=1; j < _maxY; j++)
                    for (int i=1; i < _maxX; i++)
                        op(i,j,k,vel,flags,pressure,phi,gfClamp);
            }
        }
    }
    
    MACGrid& vel;
    const FlagGrid& flags;
    const Grid<Real> & pressure;
    const Grid<Real> & phi;
    Real gfClamp;
};

// improve behavior of clamping for large time steps:
inline static Real ghostFluidWasClamped(int idx, int offset, const Grid<Real> &phi, Real gfClamp)
{
	Real alpha = thetaHelper(phi[idx], phi[idx+offset]);
	if (alpha < gfClamp) return true;
	return false;
}

struct ReplaceClampedGhostFluidVels : public KernelBase
{
    ReplaceClampedGhostFluidVels(MACGrid &vel, FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp)
        : KernelBase(&vel,1)
        ,vel(vel)
        ,flags(flags)
        ,pressure(pressure)
        ,phi(phi)
        ,gfClamp(gfClamp)
    {
        run();
    }
    
    inline void op(int i, int j, int k, MACGrid &vel, FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp)
    {
	    const int X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	    const int idx = flags.index(i,j,k);
	    if (!flags.isEmpty(idx)) return;

	    if( (flags.isFluid(i-1,j,k)) && ( ghostFluidWasClamped(idx-X, +X, phi, gfClamp) ) )
		    vel[idx][0] = vel[idx-X][0];
	    if( (flags.isFluid(i,j-1,k)) && ( ghostFluidWasClamped(idx-Y, +Y, phi, gfClamp) ) )
		    vel[idx][1] = vel[idx-Y][1];
	    if( flags.is3D() &&
	      ( (flags.isFluid(i,j,k-1)) && ( ghostFluidWasClamped(idx-Z, +Z, phi, gfClamp) ) ))
		    vel[idx][2] = vel[idx-Z][2];

	    if( (flags.isFluid(i+1,j,k)) && ( ghostFluidWasClamped(idx+X, -X, phi, gfClamp)) )
		    vel[idx][0] = vel[idx+X][0];
	    if( (flags.isFluid(i,j+1,k)) && ( ghostFluidWasClamped(idx+Y, -Y, phi, gfClamp)) )
		    vel[idx][1] = vel[idx+Y][1];
	    if( flags.is3D() && 
	       (flags.isFluid(i,j,k+1))  && ( ghostFluidWasClamped(idx+Z, -Z, phi, gfClamp)) )
		    vel[idx][2] = vel[idx+Z][2];
    }
    
    inline MACGrid& getArg0() { return vel; }
    typedef MACGrid type0;
    inline FlagGrid& getArg1() { return flags; }
    typedef FlagGrid type1;
    inline const Grid<Real> & getArg2() { return pressure; }
    typedef Grid<Real>  type2;
    inline const Grid<Real> & getArg3() { return phi; }
    typedef Grid<Real>  type3;
    inline Real& getArg4() { return gfClamp; }
    typedef Real type4;
    
    void run()
    {
        const int _maxX = maxX;
        const int _maxY = maxY;
        if (maxZ > 1)
        {
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int k=minZ; k < maxZ; k++)
                    for (int j=1; j < _maxY; j++)
                        for (int i=1; i < _maxX; i++)
                            op(i,j,k,vel,flags,pressure,phi,gfClamp);
            }
        } else
        {
            const int k=0;
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int j=1; j < _maxY; j++)
                    for (int i=1; i < _maxX; i++)
                        op(i,j,k,vel,flags,pressure,phi,gfClamp);
            }
        }
    }
    
    MACGrid& vel;
    FlagGrid& flags;
    const Grid<Real> & pressure;
    const Grid<Real> & phi;
    Real gfClamp;
};

//! Kernel: Compute min value of Real grid
struct CountEmptyCells : public KernelBase
{
    CountEmptyCells(FlagGrid& flags)
        : KernelBase(&flags,0)
        , flags(flags)
        , numEmpty(0)
    {
        run();
    }
    
    inline void op(int idx, FlagGrid& flags ,int& numEmpty)
    {
	    if (flags.isEmpty(idx)) numEmpty++;
    }
    
    inline operator int () { return numEmpty; }
    inline int  & getRet() { return numEmpty; }
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
            int numEmpty = 0;
#pragma omp for nowait
            for (int i=0; i < _sz; i++)
                op(i,flags,numEmpty);
#pragma omp critical
            {
                this->numEmpty += numEmpty;
            }
        }
    }
    
    FlagGrid& flags;
    int numEmpty;
};

// *****************************************************************************
// Main pressure solve

//! Perform pressure projection of the velocity grid
void solvePressure(MACGrid& vel, Grid<Real>& pressure, FlagGrid& flags, Grid<Real>* phi, Grid<Real>* perCellCorr, MACGrid* fractions, Real gfClamp, Real cgMaxIterFac, Real cgAccuracy, bool precondition, bool enforceCompatibility, bool useL2Norm, bool zeroPressureFixing, Grid<Real>* retRhs)
{
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
	Grid<Real> pca0(parent);
	Grid<Real> pca1(parent);
	Grid<Real> pca2(parent);
	Grid<Real> pca3(parent);
		
	// setup matrix and boundaries 
	MakeLaplaceMatrix (flags, A0, Ai, Aj, Ak, fractions);

	if (phi) ApplyGhostFluidDiagonal(A0, flags, *phi, gfClamp);
	
	// compute divergence and init right hand side
	MakeRhs kernMakeRhs(flags, rhs, vel, perCellCorr, fractions);
	
	if (enforceCompatibility) rhs += (Real)(-kernMakeRhs.sum / (Real)kernMakeRhs.cnt);
	
	// check whether we need to fix some pressure value...
	// (manually enable, or automatically for high accuracy, can cause asymmetries otherwise)
	int fixPidx = -1;
	if(zeroPressureFixing || cgAccuracy<1e-07) 
	{
		if (FLOATINGPOINT_PRECISION==1) debMsg("Warning - high CG accuracy with single-precision floating point accuracy might not converge...", 2);

		int numEmpty = CountEmptyCells(flags);
		if (numEmpty==0)
        {
			FOR_IJK_BND(flags,1)
            {
				if(flags.isFluid(i,j,k))
                {
					fixPidx = flags.index(i,j,k);
					break;
				}
			}
			//debMsg("No empty cells! Fixing pressure of cell "<<fixPidx<<" to zero",1);
		}

		if (fixPidx >= 0)
        {
			flags[fixPidx] |= FlagGrid::TypeZeroPressure;
			rhs[fixPidx] = 0.; 
			debMsg("Pinning pressure of cell "<<fixPidx<<" to zero", 2);
		}
	}

	// CG setup
	// note: the last factor increases the max iterations for 2d, which right now can't use a preconditioner 
	const int maxIter = (int)(cgMaxIterFac * flags.getSize().max()) * (flags.is3D() ? 1 : 4);
	GridCgInterface *gcg;
	if (vel.is3D())
		gcg = new GridCg<ApplyMatrix>  (pressure, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak);
	else
		gcg = new GridCg<ApplyMatrix2D>(pressure, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak);
	
	gcg->setAccuracy(cgAccuracy); 
	gcg->setUseL2Norm(useL2Norm);

	// optional preconditioning
	gcg->setPreconditioner(precondition ? GridCgInterface::PC_mICP : GridCgInterface::PC_None, &pca0, &pca1, &pca2, &pca3);

	for (int iter=0; iter<maxIter; iter++)
    {
		if (!gcg->iterate()) iter=maxIter;
	} 

	debMsg("FluidSolver::solvePressure iterations:" << gcg->getIterations() << ", res:" << gcg->getSigma(), 1);
	delete gcg;
	
	CorrectVelocity(flags, vel, pressure);

	if (phi)
    {
		CorrectVelocityGhostFluid(vel, flags, pressure, *phi, gfClamp);
		// improve behavior of clamping for large time steps:
		ReplaceClampedGhostFluidVels(vel, flags, pressure, *phi, gfClamp);
	}

	if (fixPidx>=0) flags[fixPidx] &= ~FlagGrid::TypeZeroPressure;

	// optionally , return RHS
	if (retRhs) retRhs->copyFrom( rhs );
}

} // end namespace
