/******************************************************************************
 *
 * MantaFlow fluid solver framework 
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * FLIP (fluid implicit particles)
 * for use with particle data fields
 *
 ******************************************************************************/

#include "flip.h"

using namespace std;

namespace Manta
{

//! note - this is a simplified version , sampleLevelsetWithParticles has more functionality
void sampleFlagsWithParticles(FlagGrid& flags, BasicParticleSystem& parts, int discretization, Real randomness)
{
	bool is3D = flags.is3D();
	Real jlen = randomness / discretization;
	Vec3 disp (1.0 / discretization, 1.0 / discretization, 1.0/discretization);
	RandomStream mRand(9832);
 
	FOR_IJK_BND(flags, 0)
    {
		if (flags.isObstacle(i,j,k)) continue;
		if (flags.isFluid(i,j,k))
        {
			Vec3 pos (i,j,k);
			for (int dk=0; dk<(is3D ? discretization : 1); dk++)
			for (int dj=0; dj<discretization; dj++)
			for (int di=0; di<discretization; di++)
            {
				Vec3 subpos = pos + disp * Vec3(0.5+di, 0.5+dj, 0.5+dk);
				subpos += jlen * (Vec3(1,1,1) - 2.0 * mRand.getVec3());
				if(!is3D) subpos[2] = 0.5; 
				parts.addBuffered( subpos);
			}
		}
	}

	parts.insertBufferedParticles();
}

//! sample a level set with particles, use reset to clear the particle buffer,
//! and skipEmpty for a continuous inflow (in the latter case, only empty cells will
//! be re-filled once they empty when calling sampleLevelsetWithParticles during 
//! the main loop).
void sampleLevelsetWithParticles(LevelsetGrid& phi, FlagGrid& flags, BasicParticleSystem& parts, int discretization, Real randomness, bool reset, bool refillEmpty)
{
	bool is3D = phi.is3D();
	Real jlen = randomness / discretization;
	Vec3 disp (1.0 / discretization, 1.0 / discretization, 1.0/discretization);
	RandomStream mRand(9832);
 
	if(reset)
    {
		parts.clear(); 
		parts.doCompress();
	}

	FOR_IJK_BND(phi, 0)
    {
		if (flags.isObstacle(i,j,k)) continue;
		if (refillEmpty && flags.isFluid(i,j,k)) continue;
		if (phi(i,j,k) < 1.733)
        {
			Vec3 pos (i,j,k);
			for (int dk=0; dk<(is3D ? discretization : 1); dk++)
			for (int dj=0; dj<discretization; dj++)
			for (int di=0; di<discretization; di++)
            {
				Vec3 subpos = pos + disp * Vec3(0.5+di, 0.5+dj, 0.5+dk);
				subpos += jlen * (Vec3(1,1,1) - 2.0 * mRand.getVec3());
				if(!is3D) subpos[2] = 0.5; 
				if(phi.getInterpolated(subpos) > 0.) continue; 
				parts.addBuffered( subpos);
			}
		}
	}

	parts.insertBufferedParticles();
}

//! mark fluid cells and helpers
struct knClearFluidFLags : public KernelBase
{
    knClearFluidFLags(FlagGrid& flags, int dummy=0)
        : KernelBase(&flags,0)
        , flags(flags)
        , dummy(dummy)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, int dummy=0 )
    {
	    if (flags.isFluid(i,j,k))
        {
		    flags(i,j,k) = (flags(i,j,k) | FlagGrid::TypeEmpty) & ~FlagGrid::TypeFluid;
	    }
    }
    
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
                    for (int j=0; j < _maxY; j++)
                        for (int i=0; i < _maxX; i++)
                            op(i,j,k,flags,dummy);
            }
        } else
        {
            const int k=0;
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int j=0; j < _maxY; j++)
                    for (int i=0; i < _maxX; i++)
                        op(i,j,k,flags,dummy);
            }
        }
    }
    
    FlagGrid& flags;
    int dummy;
};

struct knSetNbObstacle : public KernelBase
{
    knSetNbObstacle(FlagGrid& flags, Grid<Real>* phiObs)
        : KernelBase(&flags,1)
        , flags(flags)
        , phiObs(phiObs)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>* phiObs)
    {
	    if ((*phiObs)(i,j,k)>0.) return;
	    if (flags.isEmpty(i,j,k))
        {
		    bool set=false;
		    if( (flags.isFluid(i-1,j,k)) && (flags.isObstacle(i+1,j,k)) ) set=true;
		    if( (flags.isFluid(i+1,j,k)) && (flags.isObstacle(i-1,j,k)) ) set=true;
		    if( (flags.isFluid(i,j-1,k)) && (flags.isObstacle(i,j+1,k)) ) set=true;
		    if( (flags.isFluid(i,j+1,k)) && (flags.isObstacle(i,j-1,k)) ) set=true;
		    if(flags.is3D()) {
		    if( (flags.isFluid(i,j,k-1)) && (flags.isObstacle(i,j,k+1)) ) set=true;
		    if( (flags.isFluid(i,j,k+1)) && (flags.isObstacle(i,j,k-1)) ) set=true;
		    }
		    if(set) flags(i,j,k) = (flags(i,j,k) | FlagGrid::TypeFluid) & ~FlagGrid::TypeEmpty;
	    }
    }
    
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
                            op(i,j,k,flags,phiObs);
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
                        op(i,j,k,flags,phiObs);
            }
        }
    }
    
    FlagGrid& flags;
    Grid<Real>* phiObs;
};

void markFluidCells(BasicParticleSystem& parts, FlagGrid& flags, Grid<Real>* phiObs)
{
	// remove all fluid cells
	knClearFluidFLags(flags, 0);
	
	// mark all particles in flaggrid as fluid
	for(int idx=0;idx<parts.size();idx++)
    {
		if (!parts.isActive(idx)) continue;
		Vec3i p = toVec3i(parts.getPos(idx));
		if (flags.isInBounds(p) && flags.isEmpty(p))
			flags(p) = (flags(p) | FlagGrid::TypeFluid) & ~FlagGrid::TypeEmpty;
	}

	// special for second order obstacle BCs, check empty cells in boundary region
	if(phiObs) knSetNbObstacle( flags, phiObs );
}

void testInitGridWithPos(Grid<Real>& grid)
{
	FOR_IJK(grid)
    {
        grid(i,j,k) = norm(Vec3(i,j,k));
    }
}

//! helper to calculate particle radius factor to cover the diagonal of a cell in 2d/3d
inline Real calculateRadiusFactor(Grid<Real>& grid, Real factor)
{
	return (grid.is3D() ? sqrt(3.) : sqrt(2.) ) * (factor+.01); // note, a 1% safety factor is added here
} 

//! re-sample particles based on an input levelset 
// optionally skip seeding new particles in "exclude" SDF
void adjustNumber(BasicParticleSystem& parts, MACGrid& vel, FlagGrid& flags, int minParticles, int maxParticles, LevelsetGrid& phi, Real radiusFactor, Real narrowBand, Grid<Real>* exclude)
{
	// which levelset to use as threshold
	const Real SURFACE_LS = -1.0 * calculateRadiusFactor(phi, radiusFactor);
	Grid<int> tmp(vel.getParent());
	std::ostringstream out;

	// count particles in cells, and delete excess particles
	for (int idx=0; idx<(int)parts.size(); idx++)
    {
		if (parts.isActive(idx))
        {
			Vec3i p = toVec3i(parts.getPos(idx));
			if (!tmp.isInBounds(p))
            {
				parts.kill(idx); // out of domain, remove
				continue;
			}

			Real phiv = phi.getInterpolated(parts.getPos(idx));
			if (narrowBand>0. && phiv < -narrowBand)
            {
                parts.kill(idx);
                continue;
            }

			bool atSurface = false;
			if (phiv > SURFACE_LS) atSurface = true;
			int num = tmp(p);
			
			// dont delete particles in non fluid cells here, the particles are "always right"
			if (num > maxParticles && (!atSurface))
            {
				parts.kill(idx);
			} else
            {
				tmp(p) = num+1;
			}
		}
	}

	// seed new particles
	RandomStream mRand(9832);
	FOR_IJK(tmp)
    {
		int cnt = tmp(i,j,k);
		
		// skip cells near surface
		if (phi(i,j,k) > SURFACE_LS) continue;
		if (narrowBand>0. && phi(i,j,k) < -narrowBand)
        {
            continue;
        }
		if (exclude && ( (*exclude)(i,j,k) < 0.))
        {
            continue;
        }

		if (flags.isFluid(i,j,k) && cnt < minParticles)
        {
			for (int m=cnt; m < minParticles; m++)
            { 
				Vec3 pos = Vec3(i,j,k) + mRand.getVec3();
				//Vec3 pos (i + 0.5, j + 0.5, k + 0.5); // cell center
				parts.addBuffered( pos ); 
			}
		}
	}

	parts.doCompress();
	parts.insertBufferedParticles();
}

// simple and slow helper conversion to show contents of int grids like a real grid in the ui
// (use eg to quickly display contents of the particle-index grid)
void debugIntToReal( Grid<int>& source, Grid<Real>& dest, Real factor=1.)
{
	FOR_IJK(source)
    {
        dest(i,j,k) = (Real)source(i,j,k) * factor;
    }
}

void gridParticleIndex(BasicParticleSystem& parts, ParticleIndexSystem& indexSys, FlagGrid& flags, Grid<int>& index, Grid<int>* counter)
{
	bool delCounter = false;
	if(!counter) { counter = new Grid<int>(  flags.getParent() ); delCounter=true; }
	else         { counter->clear(); }
	
	// count particles in cells, and delete excess particles
	index.clear();
	int inactive = 0;
	for (int idx=0; idx<(int)parts.size(); idx++)
    {
		if (parts.isActive(idx))
        {
			// check index for validity...
			Vec3i p = toVec3i(parts.getPos(idx));
			if (!index.isInBounds(p)) { inactive++; continue; }

			index(p)++;
		} else
        {
			inactive++;
		}
	}

	// note - this one might be smaller...
	indexSys.resize( parts.size()-inactive );

	// convert per cell number to continuous index
	int idx=0;
	FOR_IJK(index)
    {
		int num = index(i,j,k);
		index(i,j,k) = idx;
		idx += num;
	}

	// add particles to indexed array, we still need a per cell particle counter
	for (int idx=0; idx<(int)parts.size(); idx++)
    {
		if (!parts.isActive(idx)) continue;
		Vec3i p = toVec3i(parts.getPos(idx));
		if (!index.isInBounds(p)) continue;

		// initialize position and index into original array
		//indexSys[ index(p)+(*counter)(p) ].pos        = parts[idx].pos;
		indexSys[ index(p)+(*counter)(p) ].sourceIndex = idx;
		(*counter)(p)++;
	}

	if(delCounter) delete counter;
}

struct ComputeUnionLevelsetPindex : public KernelBase
{
    ComputeUnionLevelsetPindex(Grid<int>& index, BasicParticleSystem& parts, ParticleIndexSystem& indexSys, LevelsetGrid& phi, Real radius=1.)
        : KernelBase(&index,0)
        , index(index)
        , parts(parts)
        , indexSys(indexSys)
        , phi(phi)
        , radius(radius)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<int>& index, BasicParticleSystem& parts, ParticleIndexSystem& indexSys, LevelsetGrid& phi, Real radius=1.)
    {
	    const Vec3 gridPos = Vec3(i,j,k) + Vec3(0.5); // shifted by half cell
	    Real phiv = radius * 1.0;  // outside

	    int r  = int(radius) + 1;
	    int rZ = phi.is3D() ? r : 0;
	    for(int zj=k-rZ; zj<=k+rZ; zj++) 
	    for(int yj=j-r ; yj<=j+r ; yj++) 
	    for(int xj=i-r ; xj<=i+r ; xj++)
        {
		    if (!phi.isInBounds(Vec3i(xj,yj,zj))) continue;

		    // note, for the particle indices in indexSys the access is periodic (ie, dont skip for eg inBounds(sx,10,10)
		    int isysIdxS = index.index(xj,yj,zj);
		    int pStart = index(isysIdxS), pEnd=0;
		    if(phi.isInBounds(isysIdxS+1)) pEnd = index(isysIdxS+1);
		    else                           pEnd = indexSys.size();

		    // now loop over particles in cell
		    for(int p=pStart; p<pEnd; ++p)
            {
			    const int psrc = indexSys[p].sourceIndex;
			    const Vec3 pos = parts[psrc].pos; 
			    phiv = std::min(phiv, fabs(norm(gridPos-pos)) - radius);
		    }
	    }

	    phi(i,j,k) = phiv;
    }
    
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
                    for (int j=0; j < _maxY; j++)
                        for (int i=0; i < _maxX; i++)
                            op(i,j,k,index,parts,indexSys,phi,radius);
            }
        } else
        {
            const int k=0;
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int j=0; j < _maxY; j++)
                    for (int i=0; i < _maxX; i++)
                        op(i,j,k,index,parts,indexSys,phi,radius);
            }
        }
    }
    
    Grid<int>& index;
    BasicParticleSystem& parts;
    ParticleIndexSystem& indexSys;
    LevelsetGrid& phi;
    Real radius;
};

void unionParticleLevelset(BasicParticleSystem& parts, ParticleIndexSystem& indexSys, FlagGrid& flags, Grid<int>& index, LevelsetGrid& phi, Real radiusFactor)
{
	// use half a cell diagonal as base radius
	const Real radius = 0.5 * calculateRadiusFactor(phi, radiusFactor);

	// no reset of phi necessary here 
	ComputeUnionLevelsetPindex(index, parts, indexSys, phi, radius);

	phi.setBound(0.5, 0);
}

struct ComputeAveragedLevelsetWeight : public KernelBase
{
    ComputeAveragedLevelsetWeight(BasicParticleSystem& parts, Grid<int>& index, ParticleIndexSystem& indexSys, LevelsetGrid& phi, Real radius=1.)
        : KernelBase(&index,0)
        , parts(parts)
        , index(index)
        , indexSys(indexSys)
        , phi(phi)
        , radius(radius)
    {
        run();
    }
    
    inline void op(int i, int j, int k, BasicParticleSystem& parts, Grid<int>& index, ParticleIndexSystem& indexSys, LevelsetGrid& phi, Real radius=1.)
    {
	    const Vec3 gridPos = Vec3(i,j,k) + Vec3(0.5); // shifted by half cell
	    Real phiv = radius * 1.0; // outside 

	    // loop over neighborhood, similar to ComputeUnionLevelsetPindex
	    const Real sradiusInv = 1. / (4. * radius * radius);
	    int   r = int(1. * radius) + 1;
	    int   rZ = phi.is3D() ? r : 0;

	    // accumulators
	    Real  wacc = 0.;
	    Vec3  pacc = Vec3(0.);
	    Real  racc = 0.;

	    for(int zj=k-rZ; zj<=k+rZ; zj++) 
	    for(int yj=j-r ; yj<=j+r ; yj++) 
	    for(int xj=i-r ; xj<=i+r ; xj++)
        {
		    if (!phi.isInBounds(Vec3i(xj,yj,zj))) continue;

		    int isysIdxS = index.index(xj,yj,zj);
		    int pStart = index(isysIdxS), pEnd=0;
		    if(phi.isInBounds(isysIdxS+1)) pEnd = index(isysIdxS+1);
		    else                           pEnd = indexSys.size();
		    for(int p=pStart; p<pEnd; ++p)
            {
			    int   psrc = indexSys[p].sourceIndex;
			    Vec3  pos  = parts[psrc].pos; 
			    Real  s    = normSquare(gridPos-pos) * sradiusInv;
			    //Real  w = std::max(0., cubed(1.-s) );
			    Real  w = std::max(0., (1.-s) ); // a bit smoother
			    wacc += w;
			    racc += radius * w;
			    pacc += pos    * w;
		    } 
	    }

	    if(wacc > VECTOR_EPSILON)
        {
		    racc /= wacc;
		    pacc /= wacc;
		    phiv = fabs( norm(gridPos-pacc) )-racc;
	    }

	    phi(i,j,k) = phiv;
    }
    
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
                    for (int j=0; j < _maxY; j++)
                        for (int i=0; i < _maxX; i++)
                            op(i,j,k,parts,index,indexSys,phi,radius);
            }
        } else
        {
            const int k=0;
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num();
                this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int j=0; j < _maxY; j++)
                    for (int i=0; i < _maxX; i++)
                        op(i,j,k,parts,index,indexSys,phi,radius);
            }
        }
    }
    
    BasicParticleSystem& parts;
    Grid<int>& index;
    ParticleIndexSystem& indexSys;
    LevelsetGrid& phi;
    Real radius;
};

template<class T>
T smoothingValue(Grid<T> val, int i, int j, int k, T center)
{
	return val(i,j,k);
}

// smoothing, and  
template <class T>
struct knSmoothGrid : public KernelBase
{
    knSmoothGrid(Grid<T>& me, Grid<T>& tmp, Real factor)
        : KernelBase(&me,1)
        , me(me)
        , tmp(tmp)
        , factor(factor)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<T>& me, Grid<T>& tmp, Real factor)
    {
	    T val = me(i,j,k) + 
			    me(i+1,j,k) + me(i-1,j,k) + 
			    me(i,j+1,k) + me(i,j-1,k);

	    if(me.is3D())
        {
		    val += me(i,j,k+1) + me(i,j,k-1);
	    }

	    tmp(i,j,k) = val * factor;
    }
    
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
                            op(i,j,k,me,tmp,factor);
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
                        op(i,j,k,me,tmp,factor);
            }
        }
    }
    
    Grid<T>& me;
    Grid<T>& tmp;
    Real factor;
};

template <class T>
struct knSmoothGridNeg : public KernelBase
{
    knSmoothGridNeg(Grid<T>& me, Grid<T>& tmp, Real factor)
        :  KernelBase(&me,1)
        ,me(me)
        ,tmp(tmp)
        ,factor(factor)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<T>& me, Grid<T>& tmp, Real factor)
    {
	    T val = me(i,j,k) + 
			    me(i+1,j,k) + me(i-1,j,k) + 
			    me(i,j+1,k) + me(i,j-1,k);

	    if(me.is3D())
        {
		    val += me(i,j,k+1) + me(i,j,k-1);
	    }

	    val *= factor;
	    if(val<tmp(i,j,k)) tmp(i,j,k) = val;
	    else               tmp(i,j,k) = me(i,j,k);
    }
    
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
                            op(i,j,k,me,tmp,factor);
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
                        op(i,j,k,me,tmp,factor);
            }
        }
    }
    
    Grid<T>& me;
    Grid<T>& tmp;
    Real factor;
};

void averagedParticleLevelset(BasicParticleSystem& parts, ParticleIndexSystem& indexSys, FlagGrid& flags, Grid<int>& index, LevelsetGrid& phi, Real radiusFactor, int smoothen, int smoothenNeg)
{
	// use half a cell diagonal as base radius
	const Real radius = 0.5 * calculateRadiusFactor(phi, radiusFactor); 
	ComputeAveragedLevelsetWeight(parts,  index, indexSys, phi, radius);

	// post-process level-set
	for(int i=0; i<std::max(smoothen,smoothenNeg); ++i)
    {
		LevelsetGrid tmp(flags.getParent());
		if(i<smoothen)
        { 
			knSmoothGrid    <Real> (phi,tmp, 1./(phi.is3D() ? 7. : 5.) );
			phi.swap(tmp);
		}

		if(i<smoothenNeg)
        { 
			knSmoothGridNeg <Real> (phi,tmp, 1./(phi.is3D() ? 7. : 5.) );
			phi.swap(tmp);
		}
	}

	phi.setBound(0.5, 0);
}

struct knPushOutofObs : public KernelBase
{
    knPushOutofObs(BasicParticleSystem& parts, FlagGrid& flags, Grid<Real>& phiObs, Real shift=0.05, Real thresh=0.)
        : KernelBase(parts.size())
        , parts(parts)
        , flags(flags)
        , phiObs(phiObs)
        , shift(shift)
        , thresh(thresh)
    {
        run();
    }
    
    inline void op(int idx, BasicParticleSystem& parts, FlagGrid& flags, Grid<Real>& phiObs, Real shift=0.05, Real thresh=0.)
    {
	    if (!parts.isActive(idx)) return;
	    Vec3i p = toVec3i(parts.getPos(idx));

	    if (!flags.isInBounds(p)) return;
	    Real v = phiObs.getInterpolated(parts.getPos(idx));
	    if(v < thresh)
        {
		    Vec3 grad = getGradient( phiObs, p.x,p.y,p.z );
		    if( normalize(grad) < VECTOR_EPSILON ) return;
		    parts.setPos(idx, parts.getPos(idx) + shift * grad );
	    }
    }
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,parts,flags,phiObs,shift,thresh);
        }
    }
    
    BasicParticleSystem& parts;
    FlagGrid& flags;
    Grid<Real>& phiObs;
    Real shift;
    Real thresh;
};

//! slightly push particles out of obstacle levelset
void pushOutofObs(BasicParticleSystem& parts, FlagGrid& flags, Grid<Real>& phiObs, Real shift, Real thresh)
{
	knPushOutofObs(parts, flags, phiObs, shift, thresh);
}

//******************************************************************************
// grid interpolation functions
template <class T>
struct knSafeDivReal : public KernelBase
{
    knSafeDivReal(Grid<T>& me, const Grid<Real>& other, Real cutoff=VECTOR_EPSILON)
        : KernelBase(&me,0)
        , me(me)
        , other(other)
        , cutoff(cutoff)
    {
        run();
    }
    
    inline void op(int idx, Grid<T>& me, const Grid<Real>& other, Real cutoff=VECTOR_EPSILON)
    { 
	    if(other[idx]<cutoff)
        {
		    me[idx] = 0.;
	    } else
        {
		    T div( other[idx] );
		    me[idx] = safeDivide(me[idx], div ); 
	    }
    }
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,other,cutoff);
        }
    }
    
    Grid<T>& me;
    const Grid<Real>& other;
    Real cutoff;
};

// Set velocities on the grid from the particle system
struct knStompVec3PerComponent : public KernelBase
{
    knStompVec3PerComponent(Grid<Vec3>& grid, Real threshold)
        : KernelBase(&grid,0)
        , grid(grid)
        , threshold(threshold)
    {
        run();
    }
    
    inline void op(int idx, Grid<Vec3>& grid, Real threshold)
    {
	    if(grid[idx][0] < threshold) grid[idx][0] = 0.;
	    if(grid[idx][1] < threshold) grid[idx][1] = 0.;
	    if(grid[idx][2] < threshold) grid[idx][2] = 0.;
    }
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for
                (int i=0; i < _sz; i++)
                op(i,grid,threshold);
        }
    }
    
    Grid<Vec3>& grid;
    Real threshold;
};

struct knMapLinearVec3ToMACGrid : public KernelBase
{
    knMapLinearVec3ToMACGrid(BasicParticleSystem& p, FlagGrid& flags, MACGrid& vel, Grid<Vec3>& tmp, ParticleDataImpl<Vec3>& pvel)
        : KernelBase(p.size())
        , p(p)
        , flags(flags)
        , vel(vel)
        , tmp(tmp)
        , pvel(pvel)
    {
        run();
    }
    
    inline void op(int idx,  BasicParticleSystem& p, FlagGrid& flags, MACGrid& vel, Grid<Vec3>& tmp, ParticleDataImpl<Vec3>& pvel)
    {
	    unusedParameter(flags);
	    if (!p.isActive(idx)) return;
	    vel.setInterpolated(p[idx].pos, pvel[idx], &tmp[0]);
    }
    
    void run()
    {
        const int _sz = size;
        for (int i=0; i < _sz; i++)
            op(i, p,flags,vel,tmp,pvel);
    }
    
    BasicParticleSystem& p;
    FlagGrid& flags;
    MACGrid& vel;
    Grid<Vec3>& tmp;
    ParticleDataImpl<Vec3>& pvel;
};

// optionally, this function can use an existing vec3 grid to store the weights
// this is useful in combination with the simple extrapolation function
void mapParticlesToMAC(FlagGrid& flags, MACGrid& vel, MACGrid& velOld, BasicParticleSystem& parts, ParticleDataImpl<Vec3>& partVel, Grid<Vec3>* weight)
{
	// interpol -> grid. tmpgrid for particle contribution weights
	bool freeTmp = false;
	if(!weight)
    {
		weight = new Grid<Vec3>(flags.getParent());
		freeTmp = true;
	} else
    {
		weight->clear(); // make sure we start with a zero grid!
	}
	vel.clear();
	knMapLinearVec3ToMACGrid(parts, flags, vel, *weight, partVel);

	// stomp small values in weight to zero to prevent roundoff errors
	knStompVec3PerComponent(*weight, VECTOR_EPSILON);
	vel.safeDivide(*weight);
	
	// store original state
	velOld.copyFrom(vel);
	if(freeTmp) delete weight;
}

struct knCombineVels : public KernelBase
{
    knCombineVels(MACGrid& vel, Grid<Vec3>& w, MACGrid& combineVel)
        : KernelBase(&vel,0)
        , vel(vel)
        , w(w)
        , combineVel(combineVel)
    {
        run();
    }
    
    inline void op(int idx, MACGrid& vel, Grid<Vec3>& w, MACGrid& combineVel)
    {
	    for(int c=0; c<3; ++c) if(w[idx][c]>0.1) combineVel[idx][c] = vel[idx][c];
    }
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++) op(i,vel,w,combineVel);
        }
    }
    
    MACGrid& vel;
    Grid<Vec3>& w;
    MACGrid& combineVel;
};

void combineGridVel(MACGrid& vel, Grid<Vec3>& weight, MACGrid& combineVel)
{
	knCombineVels(vel, weight, combineVel);
}

template <class T>
struct knMapLinear : public KernelBase
{
    knMapLinear(BasicParticleSystem& p, FlagGrid& flags, Grid<T>& target, Grid<Real>& gtmp, ParticleDataImpl<T>& psource)
        : KernelBase(p.size())
        , p(p)
        , flags(flags)
        , target(target)
        , gtmp(gtmp)
        , psource(psource)
    {
        run();
    }
    
    inline void op(int idx,  BasicParticleSystem& p, FlagGrid& flags, Grid<T>& target, Grid<Real>& gtmp, ParticleDataImpl<T>& psource)
    {
	    unusedParameter(flags);
	    if (!p.isActive(idx)) return;
	    target.setInterpolated(p[idx].pos, psource[idx], gtmp);
    }
    
    void run()
    {
        const int _sz = size;
        for (int i=0; i < _sz; i++)
            op(i, p,flags,target,gtmp,psource);
    }
    
    BasicParticleSystem& p;
    FlagGrid& flags;
    Grid<T>& target;
    Grid<Real>& gtmp;
    ParticleDataImpl<T>& psource;
}; 

template<class T>
void mapLinearRealHelper(FlagGrid& flags, Grid<T>& target, BasicParticleSystem& parts, ParticleDataImpl<T>& source) 
{
	Grid<Real> tmp(flags.getParent());
	target.clear();
	knMapLinear<T>(parts, flags, target, tmp, source); 
	knSafeDivReal<T>(target, tmp);
}

void mapPartsToGrid(FlagGrid& flags, Grid<Real>& target , BasicParticleSystem& parts , ParticleDataImpl<Real>& source)
{
	mapLinearRealHelper<Real>(flags,target,parts,source);
}

void mapPartsToGridVec3(FlagGrid& flags, Grid<Vec3>& target , BasicParticleSystem& parts , ParticleDataImpl<Vec3>& source)
{
	mapLinearRealHelper<Vec3>(flags,target,parts,source);
}

template <class T>
struct knMapFromGrid : public KernelBase
{
    knMapFromGrid(BasicParticleSystem& p, Grid<T>& gsrc, ParticleDataImpl<T>& target)
        : KernelBase(p.size())
        , p(p)
        , gsrc(gsrc)
        , target(target)
    {
        run();
    }
    
    inline void op(int idx,  BasicParticleSystem& p, Grid<T>& gsrc, ParticleDataImpl<T>& target)
    {
	    if (!p.isActive(idx)) return;
	    target[idx] = gsrc.getInterpolated( p[idx].pos );
    }
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,p,gsrc,target);
        }
    }
    
    BasicParticleSystem& p;
    Grid<T>& gsrc;
    ParticleDataImpl<T>& target;
};
 
void mapGridToParts(Grid<Real>& source , BasicParticleSystem& parts , ParticleDataImpl<Real>& target)
{
	knMapFromGrid<Real>(parts, source, target);
}

void mapGridToPartsVec3(Grid<Vec3>& source , BasicParticleSystem& parts , ParticleDataImpl<Vec3>& target)
{
	knMapFromGrid<Vec3>(parts, source, target);
}

// Get velocities from grid
struct knMapLinearMACGridToVec3_PIC : public KernelBase
{
    knMapLinearMACGridToVec3_PIC(BasicParticleSystem& p, FlagGrid& flags, MACGrid& vel, ParticleDataImpl<Vec3>& pvel)
        : KernelBase(p.size())
        , p(p)
        , flags(flags)
        , vel(vel),pvel(pvel)
    {
        run();
    }
    
    inline void op(int idx,  BasicParticleSystem& p, FlagGrid& flags, MACGrid& vel, ParticleDataImpl<Vec3>& pvel)
    {
	    if (!p.isActive(idx)) return;
	    // pure PIC
	    pvel[idx] = vel.getInterpolated( p[idx].pos );
    }
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,p,flags,vel,pvel);
        }
    }
    
    BasicParticleSystem& p;
    FlagGrid& flags;
    MACGrid& vel;
    ParticleDataImpl<Vec3>& pvel;
};

void mapMACToParts(FlagGrid& flags, MACGrid& vel , BasicParticleSystem& parts , ParticleDataImpl<Vec3>& partVel)
{
	knMapLinearMACGridToVec3_PIC( parts, flags, vel, partVel );
}

// with flip delta interpolation 
struct knMapLinearMACGridToVec3_FLIP : public KernelBase
{
    knMapLinearMACGridToVec3_FLIP(BasicParticleSystem& p, FlagGrid& flags, MACGrid& vel, MACGrid& oldVel, ParticleDataImpl<Vec3>& pvel, Real flipRatio)
        : KernelBase(p.size())
        , p(p)
        , flags(flags)
        , vel(vel)
        , oldVel(oldVel)
        , pvel(pvel)
        , flipRatio(flipRatio)
    {
        run();
    }
    
    inline void op(int idx,  BasicParticleSystem& p, FlagGrid& flags, MACGrid& vel, MACGrid& oldVel, ParticleDataImpl<Vec3>& pvel, Real flipRatio)
    {
	    if (!p.isActive(idx)) return; 
	    Vec3 v     =        vel.getInterpolated(p[idx].pos);
	    Vec3 delta = v - oldVel.getInterpolated(p[idx].pos); 
	    pvel[idx] = flipRatio * (pvel[idx] + delta) + (1.0 - flipRatio) * v;    
    }
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,p,flags,vel,oldVel,pvel,flipRatio);
        }
    }
    
    BasicParticleSystem& p;
    FlagGrid& flags;
    MACGrid& vel;
    MACGrid& oldVel;
    ParticleDataImpl<Vec3>& pvel;
    Real flipRatio;
};

void flipVelocityUpdate(FlagGrid& flags, MACGrid& vel, MACGrid& velOld, BasicParticleSystem& parts, ParticleDataImpl<Vec3>& partVel, Real flipRatio)
{
	knMapLinearMACGridToVec3_FLIP( parts, flags, vel, velOld, partVel, flipRatio );
}

// NT_DEBUG
void floatify(BasicParticleSystem& parts)
{
	for (int idx=0; idx<(int)parts.size(); idx++)
    {
		float f = parts.getPos(idx).x;
		float g = parts.getPos(idx).y;
		float h = parts.getPos(idx).z;
		parts.setPos(idx, Vec3(f,g,h) );
	}
}

} // namespace
