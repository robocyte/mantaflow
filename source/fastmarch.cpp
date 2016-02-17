/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Fast marching and extrapolation
 *
 ******************************************************************************/

#include "fastmarch.h"
#include "levelset.h"
#include "kernel.h"
#include <algorithm>

using namespace std;

namespace Manta
{

/*****************************************************************************/
//! FastMarch members
/*****************************************************************************/
template<class COMP, int TDIR>
FastMarch<COMP,TDIR>::FastMarch(FlagGrid& flags, Grid<int>& fmFlags, LevelsetGrid& levelset, Real maxTime, MACGrid* velTransport)
	: mLevelset(levelset)
    , mFlags(flags)
    , mFmFlags(fmFlags)
{
	if (velTransport)
		mVelTransport.initMarching(velTransport, &flags);
	
	mMaxTime = maxTime * TDIR;
}

// helper for individual components to calculateDistance
template<class COMP, int TDIR> template<int C>
Real FastMarch<COMP,TDIR>::calcWeights(int& okcnt, int& invcnt, Real* v, const Vec3i& idx)
{
	Real val = 0.;
	Vec3i idxPlus(idx), idxMinus(idx);
	idxPlus[C]++;
	idxMinus[C]--;
	
	mWeights[C*2] = mWeights[C*2+1] = 0.;
	if (mFmFlags(idxPlus)==FlagInited)
    {
		// somewhat arbitrary - choose +1 value over -1 ...
		val = mLevelset(idxPlus);
		v[okcnt] = val; okcnt++;
		mWeights[C*2] = 1.;
	} else if (mFmFlags(idxMinus)==FlagInited)
    {
		val = mLevelset(idxMinus);
		v[okcnt] = val; okcnt++;
		mWeights[C*2+1] = 1.;
	} else
    {
		invcnt++;
	}

	return val;
}

template<class COMP, int TDIR>
inline Real FastMarch<COMP,TDIR>::calculateDistance(const Vec3i& idx)
{
	//int invflag = 0;
	int invcnt = 0;
	Real v[3];
	int okcnt = 0;
	
	Real aVal = calcWeights<0>(okcnt, invcnt, v, idx);
	Real bVal = calcWeights<1>(okcnt, invcnt, v, idx);
	Real cVal = 0.;
	if (mLevelset.is3D())   cVal = calcWeights<2>(okcnt, invcnt, v, idx);
	else					{ invcnt++; mWeights[4] = mWeights[5] = 0.; }

	Real ret = InvalidTime();
	switch(invcnt)
    {
	case 0:
    {
		// take all values
		const Real ca=v[0], cb=v[1], cc=v[2];
		const Real csqrt = max(0. , 
				-2.*(ca*ca+cb*cb- cb*cc + cc*cc - ca*(cb+cc)) + 3 );
		// clamp to make sure the sqrt is valid
		ret = 0.333333*( ca+cb+cc+ TDIR*sqrt(csqrt) );

		// weights needed for transport (transpTouch)
		mWeights[0] *= fabs(ret-ca);
		mWeights[1] *= fabs(ret-ca);
		mWeights[2] *= fabs(ret-cb);
		mWeights[3] *= fabs(ret-cb);
		mWeights[4] *= fabs(ret-cc);
		mWeights[5] *= fabs(ret-cc);

		Real norm = 0.0; // try to force normalization
		for(int i=0;i<6;i++)
        { 
			norm += mWeights[i]; 
		}

		norm = 1.0/norm;
		for(int i=0;i<6;i++)
        {
            mWeights[i] *= norm;
        } 

	} break; 
	case 1:
    {
		// take just the 2 ok values
		// t=0.5*( a+b+ (2*g*g-(b-a)*(b-a))^0.5) 
		const Real csqrt = max(0. , 2.-(v[1]-v[0])*(v[1]-v[0]) );
		// clamp to make sure the sqrt is valid
		ret = 0.5*( v[0]+v[1]+ TDIR*sqrt(csqrt) );

		// weights needed for transport (transpTouch)
		mWeights[0] *= fabs(ret-aVal);
		mWeights[1] *= fabs(ret-aVal);
		mWeights[2] *= fabs(ret-bVal);
		mWeights[3] *= fabs(ret-bVal);
		mWeights[4] *= fabs(ret-cVal);
		mWeights[5] *= fabs(ret-cVal);

		Real norm = 0.0; // try to force normalization
		for(int i=0;i<6;i++)
        { 
			norm += mWeights[i]; 
		}

		norm = 1.0/norm;
		for(int i=0;i<6;i++)
        {
            mWeights[i] *= norm;
        } 

	} break; 
	case 2:
    {
		// just use the one remaining value
		ret = v[0]+ (Real)(TDIR) ; // direction = +- 1
	} break; 
	default:
		throw Error("FastMarch :: Invalid invcnt");
		break;
	}

	return ret;
}

template<class COMP, int TDIR>
void FastMarch<COMP,TDIR>::addToList(const Vec3i& p, const Vec3i& src)
{
	if (!mLevelset.isInBounds(p,1)) return;
	const int idx = mLevelset.index(p);
	
	// already known value, value alreay set to valid value? skip cell...
	if(mFmFlags[idx] == FlagInited) return;

	// discard by source time now , TODO do instead before calling all addtolists?
	Real srct = mLevelset(src);
	if(COMP::compare(srct, mMaxTime)) return;

	Real ttime = calculateDistance(p);
	
	// remove old entry if larger
	bool found=false;

	Real oldt = mLevelset[idx];
	if (mFmFlags[idx] == FlagIsOnHeap)
    {
		found = true;
		// is old time better?
		if(COMP::compare(ttime,oldt)) return;        
	}

	// update field
	mFmFlags[idx] = FlagIsOnHeap;
	mLevelset[idx] = ttime;
	// debug info std::cout<<"set "<< idx <<","<< ttime <<"\n";
	
	if (mVelTransport.isInitialized())
		mVelTransport.transpTouch(p.x, p.y, p.z, mWeights, ttime);

	// the following adds entries to the heap of active cells
	// current: (!found) , previous: always add, might lead to duplicate
	//     entries, but the earlier will be handled earlier, the second one will skip to the FlagInited check above
	if(!found) 
	{
		// add list entry with source value
		COMP entry;
		entry.p    = p;
		entry.time = mLevelset[idx];

		mHeap.push( entry );
		// debug info std::cout<<"push "<< entry.p <<","<< entry.time <<"\n";
	}	

}

//! Enforce delta_phi = 0 on boundaries
struct SetLevelsetBoundaries : public KernelBase
{
    SetLevelsetBoundaries(LevelsetGrid& phi)
        : KernelBase(&phi,0)
        , phi(phi)
    {
        run();
    }
    
    inline void op(int i, int j, int k, LevelsetGrid& phi)
    {
	    if (i==0)      phi(i,j,k) = phi(1,j,k);
	    if (i==maxX-1) phi(i,j,k) = phi(i-1,j,k);

	    if (j==0)      phi(i,j,k) = phi(i,1,k);
	    if (j==maxY-1) phi(i,j,k) = phi(i,j-1,k);

	    if(phi.is3D())
        {
		    if (k==0)      phi(i,j,k) = phi(i,j,1);
		    if (k==maxZ-1) phi(i,j,k) = phi(i,j,k-1);
	    }
    }
    
    inline LevelsetGrid& getArg0() { return phi; }
    typedef LevelsetGrid type0;
    
    void run()
    {
        const int _maxX = maxX;
        const int _maxY = maxY;
        for (int k=minZ; k< maxZ; k++)
            for (int j=0; j< _maxY; j++)
                for (int i=0; i< _maxX; i++)
                    op(i,j,k, phi);
    }
    
    LevelsetGrid& phi;
};

template<class COMP, int TDIR>
void FastMarch<COMP,TDIR>::performMarching()
{
	mReheapVal = 0.0;
	while(mHeap.size() > 0)
    {
		
		const COMP& ce = mHeap.top();
		Vec3i p = ce.p; 
		mFmFlags(p) = FlagInited;
		mHeap.pop();
		// debug info std::cout<<"pop "<< ce.p <<","<< ce.time <<"\n";

		addToList(Vec3i(p.x-1,p.y,p.z), p);
		addToList(Vec3i(p.x+1,p.y,p.z), p);
		addToList(Vec3i(p.x,p.y-1,p.z), p);
		addToList(Vec3i(p.x,p.y+1,p.z), p);
		if (mLevelset.is3D())
        {
			addToList(Vec3i(p.x,p.y,p.z-1), p);
			addToList(Vec3i(p.x,p.y,p.z+1), p);        
		}
	}
	
	// set boundary for plain array
	SetLevelsetBoundaries setls(mLevelset); 
	setls.getArg0(); // get rid of compiler warning...
}

// explicit instantiation
template class FastMarch<FmHeapEntryIn, -1>;
template class FastMarch<FmHeapEntryOut, +1>;



/*****************************************************************************/
// simpler extrapolation functions (primarily for FLIP)
/*****************************************************************************/
struct knExtrapolateMACSimple : public KernelBase
{
    knExtrapolateMACSimple(MACGrid& vel, int distance, Grid<int>& tmp, const int d, const int c)
        : KernelBase(&vel,1)
        , vel(vel)
        , distance(distance)
        , tmp(tmp)
        , d(d)
        , c(c)
    {
        run();
    }
    
    inline void op(int i, int j, int k, MACGrid& vel, int distance, Grid<int>& tmp, const int d, const int c)
    {
	    static const Vec3i nb[6] =
        { 
		    Vec3i(1 ,0,0), Vec3i(-1,0,0),
		    Vec3i(0,1 ,0), Vec3i(0,-1,0),
		    Vec3i(0,0,1 ), Vec3i(0,0,-1)
        };

	    const int dim = (vel.is3D() ? 3:2);

	    if (tmp(i,j,k) != 0) return;

	    // copy from initialized neighbors
	    Vec3i p(i,j,k);
	    int nbs = 0;
	    Real avgVel = 0.;
	    for (int n=0; n<2*dim; ++n)
        {
		    if (tmp(p+nb[n]) == d)
            {
			    //vel(p)[c] = (c+1.)*0.1;
			    avgVel += vel(p+nb[n])[c];
			    nbs++;
		    }
	    }

	    if (nbs>0)
        {
		    tmp(p)    = d+1;
		    vel(p)[c] = avgVel / nbs;
	    }
    }
    
    inline MACGrid& getArg0() { return vel; }
    typedef MACGrid type0;
    inline int& getArg1() { return distance; }
    typedef int type1;
    inline Grid<int>& getArg2() { return tmp; }
    typedef Grid<int> type2;
    inline const int& getArg3() { return d; }
    typedef int type3;
    inline const int& getArg4() { return c; }
    typedef int type4;
    
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
                            op(i,j,k,vel,distance,tmp,d,c);
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
                        op(i,j,k,vel,distance,tmp,d,c);
            }
        }
    }
    
    MACGrid& vel;
    int distance;
    Grid<int>& tmp;
    const int d;
    const int c;
};

// TODO, make parallel by only modifying ijk (ie, turn push into pull)
struct knExtrapolateIntoBnd : public KernelBase
{
    knExtrapolateIntoBnd(FlagGrid& flags, MACGrid& vel)
        : KernelBase(&flags,0)
        , flags(flags)
        , vel(vel)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, MACGrid& vel)
    {
	    int c=0;
	    Vec3 v(0,0,0);

	    if (i == 0)
        { 
		    v = vel(i+1,j,k);
		    if(v[0] < 0.) v[0] = 0.;
		    c++;
	    } else if (i == (flags.getSizeX() - 1))
        { 
		    v = vel(i-1,j,k);
		    if(v[0] > 0.) v[0] = 0.;
		    c++;
	    }

	    if (j == 0)
        { 
		    v = vel(i,j+1,k);
		    if(v[1] < 0.) v[1] = 0.;
		    c++;
	    } else if (j == (flags.getSizeY() - 1))
        { 
		    v = vel(i,j-1,k);
		    if(v[1] > 0.) v[1] = 0.;
		    c++;
	    }

	    if (flags.is3D())
        {
	        if (k == 0)
            { 
		        v = vel(i,j,k+1);
		        if (v[2] < 0.) v[2] = 0.;
		        c++;
	        } else if (k == (flags.getSizeY() - 1))
            { 
		        v = vel(i,j,k-1);
		        if(v[2] > 0.) v[2] = 0.;
		        c++;
	        }
        }

	    if (c > 0)
        {
		    vel(i,j,k) = v/(Real)c;
	    }
    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline MACGrid& getArg1() { return vel; }
    typedef MACGrid type1;
    
    void run()
    {
        const int _maxX = maxX;
        const int _maxY = maxY;
        for (int k=minZ; k< maxZ; k++)
            for (int j=0; j< _maxY; j++)
                for (int i=0; i< _maxX; i++)
                    op(i,j,k, flags,vel);
    }
    
    FlagGrid& flags;
    MACGrid& vel;
};

inline Vec3 getNormal(const Grid<Real>& data, int i, int j, int k)
{
	if (i > data.getSizeX()-2) i= data.getSizeX()-2;
	if (i < 1) i = 1;
	if (j > data.getSizeY()-2) j= data.getSizeY()-2;
	if (j < 1) j = 1;

	int kd = 1;
	if(data.is3D())
    {
	    if (k > data.getSizeZ()-2) k= data.getSizeZ()-2;
	    if (k < 1) k = 1;
	} else
    {
        kd=0;
    }

	return Vec3( data(i+1,j  ,k   ) - data(i-1,j  ,k   ) ,
				 data(i  ,j+1,k   ) - data(i  ,j-1,k   ) ,
				 data(i  ,j  ,k+kd) - data(i  ,j  ,k-kd) );
}

struct knUnprojectNormalComp : public KernelBase
{
    knUnprojectNormalComp(FlagGrid& flags, MACGrid& vel, LevelsetGrid& phi, Real maxDist)
        : KernelBase(&flags,1)
        , flags(flags)
        , vel(vel)
        , phi(phi)
        , maxDist(maxDist)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, MACGrid& vel, LevelsetGrid& phi, Real maxDist)
    {
	    // apply inside, within range near obstacle surface
	    if(phi(i,j,k)>0. || phi(i,j,k)<-maxDist) return;

	    Vec3 n = getNormal(phi, i,j,k);
	    Vec3 v = vel(i,j,k);
	    if(dot(n,v) < 0.)
        { 
		    normalize(n);
		    Real l = dot(n,v);
		    vel(i,j,k) -= n*l;
	    }
    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline MACGrid& getArg1() { return vel; }
    typedef MACGrid type1;
    inline LevelsetGrid& getArg2() { return phi; }
    typedef LevelsetGrid type2;
    inline Real& getArg3() { return maxDist; }
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
                            op(i,j,k,flags,vel,phi,maxDist);
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
                        op(i,j,k,flags,vel,phi,maxDist);
            }
        }
    }
    
    FlagGrid& flags;
    MACGrid& vel;
    LevelsetGrid& phi;
    Real maxDist;
};

// a simple extrapolation step, used for cases where there's no levelset
// (note, less accurate than fast marching extrapolation!)
// into obstacle is a special mode for second order obstable boundaries (extrapolating
// only fluid velocities, not those at obstacles)
void extrapolateMACSimple(FlagGrid& flags, MACGrid& vel, int distance, LevelsetGrid* phiObs, bool intoObs)
{
	Grid<int>tmp( flags.getParent());
	int dim = (flags.is3D() ? 3:2);

	for (int c=0; c<dim; ++c)
    {
		Vec3i dir = 0;
		dir[c] = 1;
		tmp.clear();

		// remove all fluid cells (not touching obstacles)
		FOR_IJK_BND(flags,1)
        {
			Vec3i p(i,j,k);
			bool mark = false;
			if(!intoObs)
            {
				if( flags.isFluid(p) || flags.isFluid(p-dir) ) mark = true;
			} else
            {
				if( (flags.isFluid(p) || flags.isFluid(p-dir) ) && 
					(!flags.isObstacle(p)) && (!flags.isObstacle(p-dir)) ) mark = true;
			}

			if(mark) tmp(p) = 1;
		}
		
		// extrapolate for distance
		for(int d=1; d<1+distance; ++d)
        {
			knExtrapolateMACSimple(vel, distance, tmp, d, c);
		} // d
	}

	if(phiObs) knUnprojectNormalComp( flags, vel, *phiObs, distance );

	// copy tangential values into sides
	knExtrapolateIntoBnd(flags, vel);
}

struct knExtrapolateMACFromWeight : public KernelBase
{
    knExtrapolateMACFromWeight(MACGrid& vel, Grid<Vec3>& weight, int distance, const int d, const int c)
        : KernelBase(&vel,1)
        , vel(vel)
        , weight(weight)
        , distance(distance)
        , d(d)
        , c(c)
    {
        run();
    }
    
    inline void op(int i, int j, int k, MACGrid& vel, Grid<Vec3>& weight, int distance, const int d, const int c)
    {
	    static const Vec3i nb[6] =
        { 
		    Vec3i(1 ,0,0), Vec3i(-1,0,0),
		    Vec3i(0,1 ,0), Vec3i(0,-1,0),
		    Vec3i(0,0,1 ), Vec3i(0,0,-1)
        };

	    const int dim = (vel.is3D() ? 3:2);

	    if (weight(i,j,k)[c] != 0) return;

	    // copy from initialized neighbors
	    Vec3i p(i,j,k);
	    int nbs = 0;
	    Real avgVel = 0.;
	    for (int n=0; n<2*dim; ++n)
        {
		    if (weight(p+nb[n])[c] == d)
            {
			    avgVel += vel(p+nb[n])[c];
			    nbs++;
		    }
	    }

	    if(nbs>0)
        {
		    weight(p)[c]    = d+1;
		    vel(p)[c] = avgVel / nbs;
	    }
    }
    
    inline MACGrid& getArg0() { return vel; }
    typedef MACGrid type0;
    inline Grid<Vec3>& getArg1() { return weight; }
    typedef Grid<Vec3> type1;
    inline int& getArg2() { return distance; }
    typedef int type2;
    inline const int& getArg3() { return d; }
    typedef int type3;
    inline const int& getArg4() { return c; }
    typedef int type4;
    
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
                            op(i,j,k,vel,weight,distance,d,c);
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
                        op(i,j,k,vel,weight,distance,d,c);
            }
        }
    }
    
    MACGrid& vel;
    Grid<Vec3>& weight;
    int distance;
    const int d;
    const int c;
};

// same as extrapolateMACSimple, but uses weight vec3 grid instead of flags to check
// for valid values (to be used in combination with mapParticlesToMAC)
// note - the weight grid values are destroyed! the function is necessary due to discrepancies
// between velocity mapping on surface-levelset / fluid-flag creation. With this
// extrapolation we make sure the fluid region is covered by initial velocities
void extrapolateMACFromWeight(MACGrid& vel, Grid<Vec3>& weight, int distance)
{
	const int dim = (vel.is3D() ? 3:2);

	for(int c=0; c<dim; ++c)
    {
		Vec3i dir = 0;
		dir[c] = 1;

		// reset weight values to 0 (uninitialized), and 1 (initialized inner values)
		FOR_IJK_BND(vel,1)
        {
			Vec3i p(i,j,k);
			if(weight(p)[c]>0.) weight(p)[c] = 1.0;
		}
		
		// extrapolate for distance
		for(int d=1; d<1+distance; ++d)
        {
			knExtrapolateMACFromWeight(vel, weight, distance, d, c);
		} // d

	}
}

// simple extrapolation functions for levelsets
static const Vec3i nb[6] =
{ 
	Vec3i(1 ,0,0), Vec3i(-1,0,0),
	Vec3i(0,1 ,0), Vec3i(0,-1,0),
	Vec3i(0,0,1 ), Vec3i(0,0,-1)
};

// larger neighborhood (>1 dx away)
static const Vec3i nbExt2d[4] =
{ 
	Vec3i(1, 1,0), Vec3i(-1, 1,0),
	Vec3i(1,-1,0), Vec3i(-1,-1,0)
};

template <class S>
struct knExtrapolateLsSimple : public KernelBase
{
    knExtrapolateLsSimple(Grid<S>& val, int distance, Grid<int>& tmp, const int d, S direction)
        : KernelBase(&val,1)
        , val(val)
        , distance(distance)
        , tmp(tmp)
        , d(d)
        , direction(direction)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<S>& val, int distance, Grid<int>& tmp, const int d, S direction)
    {
	    const int dim = (val.is3D() ? 3:2); 
	    if (tmp(i,j,k) != 0) return;

	    // copy from initialized neighbors
	    Vec3i p(i,j,k);
	    int   nbs = 0;
	    S     avg(0.);
	    for (int n=0; n<2*dim; ++n)
        {
		    if (tmp(p+nb[n]) == d)
            {
			    avg += val(p+nb[n]);
			    nbs++;
		    }
	    }

	    if(nbs>0)
        {
		    tmp(p) = d+1;
		    val(p) = avg / nbs + direction;
	    } 
    }
    
    inline Grid<S>& getArg0() { return val; }
    typedef Grid<S> type0;
    inline int& getArg1() { return distance; }
    typedef int type1;
    inline Grid<int>& getArg2() { return tmp; }
    typedef Grid<int> type2;
    inline const int& getArg3() { return d; }
    typedef int type3;
    inline S& getArg4() { return direction; }
    typedef S type4;
    
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
                            op(i,j,k,val,distance,tmp,d,direction);
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
                        op(i,j,k,val,distance,tmp,d,direction);
            }
        }
    }
    
    Grid<S>& val;
    int distance;
    Grid<int>& tmp;
    const int d;
    S direction;
};

template <class S>
struct knSetRemaining : public KernelBase
{
    knSetRemaining(Grid<S>& phi, Grid<int>& tmp, S distance)
        : KernelBase(&phi,1)
        , phi(phi)
        , tmp(tmp)
        , distance(distance)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<S>& phi, Grid<int>& tmp, S distance)
    {
	    if (tmp(i,j,k) != 0) return;
	    phi(i,j,k) = distance;
    }
    
    inline Grid<S>& getArg0() { return phi; }
    typedef Grid<S> type0;
    inline Grid<int>& getArg1() { return tmp; }
    typedef Grid<int> type1;
    inline S& getArg2() { return distance; }
    typedef S type2;
    
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
                            op(i,j,k,phi,tmp,distance);
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
                        op(i,j,k,phi,tmp,distance);
            }
        }
    }
    
    Grid<S>& phi;
    Grid<int>& tmp;
    S distance;
};

void extrapolateLsSimple(Grid<Real>& phi, int distance, bool inside)
{
	Grid<int> tmp( phi.getParent() );
	tmp.clear();
	const int dim = (phi.is3D() ? 3:2);

	// by default, march outside
	Real direction = 1.;
	if (!inside)
    { 
		// mark all inside
		FOR_IJK_BND(phi,1)
        {
			if (phi(i,j,k) < 0.)
            {
                tmp(i,j,k) = 1;
            }
		} 
	} else
    {
		direction = -1.;
		FOR_IJK_BND(phi,1)
        {
			if ( phi(i,j,k) > 0. ) { tmp(i,j,k) = 1; }
		} 
	}

	// + first layer around
	FOR_IJK_BND(phi,1)
    {
		Vec3i p(i,j,k);
		if (tmp(p)) continue;
		for (int n = 0; n < 2 * dim; ++n)
        {
			if (tmp(p+nb[n]) == 1)
            {
				tmp(i,j,k) = 2; n=2*dim;
			}
		}
	} 

	// extrapolate for distance
	for (int d=2; d<1+distance; ++d)
    {
		knExtrapolateLsSimple<Real>(phi, distance, tmp, d, direction);
	} 

	// set all remaining cells to max
	knSetRemaining<Real>(phi, tmp, Real(direction * (distance+2)));
}

// extrapolate centered vec3 values from marked fluid cells
void extrapolateVec3Simple(Grid<Vec3>& vel, Grid<Real>& phi, int distance)
{
	Grid<int>tmp(vel.getParent());
	tmp.clear();
	const int dim = (vel.is3D() ? 3:2);

	// mark all inside
	FOR_IJK_BND(vel,1)
    {
		if (phi(i,j,k) < 0.)
        {
            tmp(i,j,k) = 1;
        }
	} 
	// + first layer outside
	FOR_IJK_BND(vel,1)
    {
		Vec3i p(i,j,k);
		if (tmp(p)) continue;

		for (int n = 0; n < 2 * dim; ++n)
        {
			if (tmp(p + nb[n]) == 1)
            {
				tmp(i,j,k) = 2;
                n = 2 * dim;
			}
		}
	} 

	for(int d=2; d<1+distance; ++d)
    {
		knExtrapolateLsSimple<Vec3>(vel, distance, tmp, d, Vec3(0.));
	}

	knSetRemaining<Vec3>(vel, tmp, Vec3(0.));
}

} // namespace
