/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Levelset
 *
 ******************************************************************************/

#include "levelset.h"
#include "fastmarch.h"
#include "kernel.h"
#include "mcubes.h"
#include "mesh.h"

using namespace std;

namespace Manta
{

//************************************************************************
// Helper functions and kernels for marching
//************************************************************************
static const int FlagInited = FastMarch<FmHeapEntryOut, +1>::FlagInited;

// neighbor lookup vectors
static const Vec3i neighbors[6] =
{
    Vec3i(-1,0,0), Vec3i(1,0,0), Vec3i(0,-1,0), Vec3i(0,1,0), Vec3i(0,0,-1), Vec3i(0,0,1)
};

struct InitFmIn : public KernelBase
{
    InitFmIn(FlagGrid& flags, Grid<int>& fmFlags, LevelsetGrid& phi, bool ignoreWalls, int obstacleType)
        : KernelBase(&flags,1)
        , flags(flags)
        , fmFlags(fmFlags)
        , phi(phi)
        , ignoreWalls(ignoreWalls)
        , obstacleType(obstacleType)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, Grid<int>& fmFlags, LevelsetGrid& phi, bool ignoreWalls, int obstacleType)
    {
	    const int idx = flags.index(i,j,k);
	    const Real v = phi[idx];
	    if (ignoreWalls)
        {
		    if (v>=0. && ((flags[idx] & obstacleType) == 0))
			    fmFlags[idx] = FlagInited;
		    else
			    fmFlags[idx] = 0;
	    } else
        {
		    if (v>=0) fmFlags[idx] = FlagInited;
		    else      fmFlags[idx] = 0;
	    }
    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline Grid<int>& getArg1() { return fmFlags; }
    typedef Grid<int> type1;
    inline LevelsetGrid& getArg2() { return phi; }
    typedef LevelsetGrid type2;
    inline bool& getArg3() { return ignoreWalls; }
    typedef bool type3;
    inline int& getArg4() { return obstacleType; }
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
                            op(i,j,k,flags,fmFlags,phi,ignoreWalls,obstacleType);
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
                        op(i,j,k,flags,fmFlags,phi,ignoreWalls,obstacleType);
            }
        }
    }
    
    FlagGrid& flags;
    Grid<int>& fmFlags;
    LevelsetGrid& phi;
    bool ignoreWalls;
    int obstacleType;
};

struct InitFmOut : public KernelBase
{
    InitFmOut(FlagGrid& flags, Grid<int>& fmFlags, LevelsetGrid& phi, bool ignoreWalls, int obstacleType)
        : KernelBase(&flags,1)
        , flags(flags)
        , fmFlags(fmFlags)
        , phi(phi)
        , ignoreWalls(ignoreWalls)
        , obstacleType(obstacleType)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, Grid<int>& fmFlags, LevelsetGrid& phi, bool ignoreWalls, int obstacleType)
    {
	    const int idx = flags.index(i,j,k);
	    const Real v = phi[idx];

	    if (ignoreWalls)
        {
		    fmFlags[idx] = (v<0) ? FlagInited : 0;
		    if ((flags[idx] & obstacleType) != 0)
            {
			    fmFlags[idx] = 0;
			    phi[idx] = 0;
		    }
	    } else
        {
		    fmFlags[idx] = (v<0) ? FlagInited : 0;
	    }
    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline Grid<int>& getArg1() { return fmFlags; }
    typedef Grid<int> type1;
    inline LevelsetGrid& getArg2() { return phi; }
    typedef LevelsetGrid type2;
    inline bool& getArg3() { return ignoreWalls; }
    typedef bool type3;
    inline int& getArg4() { return obstacleType; }
    typedef int type4;
    
    void run()
    {  const int _maxX = maxX;
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
                        op(i,j,k,flags,fmFlags,phi,ignoreWalls,obstacleType);
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
                    op(i,j,k,flags,fmFlags,phi,ignoreWalls,obstacleType);
        }
    }
    }
    
    FlagGrid& flags;
    Grid<int>& fmFlags;
    LevelsetGrid& phi;
    bool ignoreWalls;
    int obstacleType;
};

struct SetUninitialized : public KernelBase
{
    SetUninitialized(Grid<int>& flags, Grid<int>& fmFlags, LevelsetGrid& phi, const Real val, int ignoreWalls, int obstacleType)
        : KernelBase(&flags,1)
        , flags(flags)
        , fmFlags(fmFlags)
        , phi(phi)
        , val(val)
        , ignoreWalls(ignoreWalls)
        , obstacleType(obstacleType)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<int>& flags, Grid<int>& fmFlags, LevelsetGrid& phi, const Real val, int ignoreWalls, int obstacleType)
    {
	    if (ignoreWalls)
        {
		    if ((fmFlags(i,j,k) != FlagInited) && ((flags(i,j,k) & obstacleType) == 0))
            {
			    phi(i,j,k) = val;
            }
	    } else
        {
		    if ((fmFlags(i,j,k) != FlagInited)) phi(i,j,k) = val;
	    }
    }
    
    inline Grid<int>& getArg0() { return flags; }
    typedef Grid<int> type0;
    inline Grid<int>& getArg1() { return fmFlags; }
    typedef Grid<int> type1;
    inline LevelsetGrid& getArg2() { return phi; }
    typedef LevelsetGrid type2;
    inline const Real& getArg3() { return val; }
    typedef Real type3;
    inline int& getArg4() { return ignoreWalls; }
    typedef int type4;
    inline int& getArg5() { return obstacleType; }
    typedef int type5;
    
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
                            op(i,j,k,flags,fmFlags,phi,val,ignoreWalls,obstacleType);
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
                        op(i,j,k,flags,fmFlags,phi,val,ignoreWalls,obstacleType);
            }
        }
    }
    
    Grid<int>& flags;
    Grid<int>& fmFlags;
    LevelsetGrid& phi;
    const Real val;
    int ignoreWalls;
    int obstacleType;
};

template<bool inward>
inline bool isAtInterface(Grid<int>& fmFlags, LevelsetGrid& phi, const Vec3i& p)
{
	// check for interface
	int max = phi.is3D() ? 6 : 4;

	for (int nb=0; nb<max; nb++)
    {
		const Vec3i pn(p + neighbors[nb]);
		if (!fmFlags.isInBounds(pn)) continue;
		
		if (fmFlags(pn) != FlagInited) continue;
		if (( inward && phi(pn) >= 0.) || 
			(!inward && phi(pn) <  0.)) return true;
	}

	return false;
}



//************************************************************************
// Levelset class members
//************************************************************************
LevelsetGrid::LevelsetGrid(FluidSolver* parent, bool show) 
	: Grid<Real>(parent, show) 
{ 
	mType = (GridType)(TypeLevelset | TypeReal);    
}    

Real LevelsetGrid::invalidTimeValue()
{
	return FastMarch<FmHeapEntryOut, 1>::InvalidTime();
}

//! Kernel: perform levelset union
struct KnJoin : public KernelBase
{
    KnJoin(Grid<Real>& a, const Grid<Real>& b)
        : KernelBase(&a,0)
        , a(a)
        , b(b)
    {
        run();
    }
    
    inline void op(int idx, Grid<Real>& a, const Grid<Real>& b)
    {
	    a[idx] = min(a[idx], b[idx]);
    }
    
    inline Grid<Real>& getArg0() { return a; }
    typedef Grid<Real> type0;
    inline const Grid<Real>& getArg1() { return b; }
    typedef Grid<Real> type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,a,b);
        }
    }
    
    Grid<Real>& a;
    const Grid<Real>& b;
};
 
void LevelsetGrid::join(const LevelsetGrid& o)
{
    KnJoin(*this, o);
}

//! subtract b, note does not preserve SDF!
struct KnSubtract : public KernelBase
{
    KnSubtract(Grid<Real>& a, const Grid<Real>& b)
        : KernelBase(&a,0)
        , a(a)
        , b(b)
    {
        run();
    }
    
    inline void op(int idx, Grid<Real>& a, const Grid<Real>& b)
    {
	    if(b[idx]<0.) a[idx] = b[idx] * -1.;
    }
    
    inline Grid<Real>& getArg0() { return a; }
    typedef Grid<Real> type0;
    inline const Grid<Real>& getArg1() { return b; }
    typedef Grid<Real> type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,a,b);
        }
    }
    
    Grid<Real>& a;
    const Grid<Real>& b;
};
 
void LevelsetGrid::subtract(const LevelsetGrid& o)
{
    KnSubtract(*this, o);
}

//! re-init levelset and extrapolate velocities (in & out)
//  note - uses flags to identify border (could also be done based on ls values)
void LevelsetGrid::reinitMarching(FlagGrid& flags, Real maxTime, MACGrid* velTransport,bool ignoreWalls, bool correctOuterLayer, int obstacleType)
{
	const int dim = (is3D() ? 3 : 2);
	
	Grid<int> fmFlags(mParent);
	LevelsetGrid& phi = *this;
	FastMarch<FmHeapEntryIn,  -1> marchIn (flags, fmFlags, phi, maxTime, NULL);
	
	// march inside
	InitFmIn (flags, fmFlags, phi, ignoreWalls, obstacleType);
	
	FOR_IJK_BND(flags, 1)
    {
		if (fmFlags(i,j,k) == FlagInited) continue;
		if (ignoreWalls && ((flags(i,j,k) & obstacleType) != 0)) continue;
		const Vec3i p(i,j,k);
				
		if(isAtInterface<true>(fmFlags, phi, p))
        {
			// set value
			fmFlags(p) = FlagInited;
			
			// add neighbors that are not at the interface
			for (int nb = 0; nb < 2 * dim; nb++)
            {
				const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
				if (ignoreWalls && ((flags.get(pn) & obstacleType) != 0)) continue;
				
				// check neighbors of neighbor
				if (phi(pn) < 0. && !isAtInterface<true>(fmFlags, phi, pn))
                {
					marchIn.addToList(pn, p); 
				}
			}            
		}
	}

	marchIn.performMarching();     
	// done with inwards marching
   
	// now march out...    
	
	// set un initialized regions
	SetUninitialized(flags, fmFlags, phi, -maxTime - 1., ignoreWalls, obstacleType); 

	InitFmOut(flags, fmFlags, phi, ignoreWalls, obstacleType);
	
	FastMarch<FmHeapEntryOut, +1>marchOut(flags, fmFlags, phi, maxTime, velTransport);

	// by default, correctOuterLayer is on
	if (correctOuterLayer)
    {
		// normal version, inwards march is done, now add all outside values (0..2] to list
		// note, this might move the interface a bit! but keeps a nice signed distance field...        
		FOR_IJK_BND(flags, 1)
        {
			if (ignoreWalls && ((flags(i,j,k) & obstacleType) != 0)) continue;
			const Vec3i p(i,j,k);
			
			// check nbs
			for (int nb = 0; nb < 2 * dim; nb++)
            {
				const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
				
				if (fmFlags(pn) != FlagInited) continue;
				if (ignoreWalls && ((flags.get(pn) & obstacleType)) != 0) continue;
				
				const Real nbPhi = phi(pn);
				
				// only add nodes near interface, not e.g. outer boundary vs. invalid region                
				if (nbPhi < 0 && nbPhi >= -2)
					marchOut.addToList(p, pn); 
			}
		}         
	} else
    {
		// alternative version, keep interface, do not distort outer cells
		// add all ouside values, but not those at the IF layer
		FOR_IJK_BND(flags, 1)
        {
			if (ignoreWalls && ((flags(i,j,k) & obstacleType) != 0)) continue;
			
			// only look at ouside values
			const Vec3i p(i,j,k);
			if (phi(p) < 0) continue;
			
			if (isAtInterface<false>(fmFlags, phi, p))
            {
				// now add all non, interface neighbors
				fmFlags(p) = FlagInited;
				
				// add neighbors that are not at the interface
				for (int nb = 0; nb < 2 * dim; nb++)
                {
					const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
					if (ignoreWalls && ((flags.get(pn) & obstacleType) != 0)) continue;
				
					// check neighbors of neighbor
					if (phi(pn) > 0. && !isAtInterface<false>(fmFlags, phi, pn))
                    {
						marchOut.addToList(pn, p);
					}
				}            
			}
		}
	}    

	marchOut.performMarching();

	// set un initialized regions
	SetUninitialized(flags, fmFlags, phi, +maxTime + 1., ignoreWalls, obstacleType);    
}

void LevelsetGrid::initFromFlags(FlagGrid& flags, bool ignoreWalls)
{
	FOR_IDX(*this)
    {
		if (flags.isFluid(idx) || (ignoreWalls && flags.isObstacle(idx)))
			mData[idx] = -0.5;
		else
			mData[idx] = 0.5;
	}
}

//! run marching cubes to create a mesh for the 0-levelset
void LevelsetGrid::createMesh(Mesh& mesh)
{
	assertMsg(is3D(), "Only 3D grids supported so far");
	
	mesh.clear();
		
	const Real invalidTime = invalidTimeValue();
	const Real isoValue = 1e-4;
	
	// create some temp grids
	Grid<int> edgeVX(mParent);
	Grid<int> edgeVY(mParent);
	Grid<int> edgeVZ(mParent);
	
	for(int i=0; i<mSize.x-1; i++)
	for(int j=0; j<mSize.y-1; j++)
	for(int k=0; k<mSize.z-1; k++)
    {
		 Real value[8] =
         {
             get(i,j,k),   get(i+1,j,k),   get(i+1,j+1,k),   get(i,j+1,k),
			 get(i,j,k+1), get(i+1,j,k+1), get(i+1,j+1,k+1), get(i,j+1,k+1)
         };
		
		// build lookup index, check for invalid times
		bool skip = false;
		int cubeIdx = 0;
		for (int l = 0; l < 8 ;l++)
        {
			value[l] *= -1;
			if (-value[l] <= invalidTime)
				skip = true;
			if (value[l] < isoValue) 
				cubeIdx |= 1<<l;
		}

		if (skip || (mcEdgeTable[cubeIdx] == 0)) continue;
		
		// where to look up if this point already exists
		int triIndices[12];
		int *eVert[12] = { &edgeVX(i,j,k),   &edgeVY(i+1,j,k),   &edgeVX(i,j+1,k),   &edgeVY(i,j,k), 
						   &edgeVX(i,j,k+1), &edgeVY(i+1,j,k+1), &edgeVX(i,j+1,k+1), &edgeVY(i,j,k+1), 
						   &edgeVZ(i,j,k),   &edgeVZ(i+1,j,k),   &edgeVZ(i+1,j+1,k), &edgeVZ(i,j+1,k) };
		
		const Vec3 pos[9] = { Vec3(i,j,k),   Vec3(i+1,j,k),   Vec3(i+1,j+1,k),   Vec3(i,j+1,k),
						      Vec3(i,j,k+1), Vec3(i+1,j,k+1), Vec3(i+1,j+1,k+1), Vec3(i,j+1,k+1) };
		
		for (int e = 0; e < 12; e++)
        {
			if (mcEdgeTable[cubeIdx] & (1<<e))
            {
				// vertex already calculated ?
				if (*eVert[e] == 0)
                {
					// interpolate edge
					const int e1 = mcEdges[e*2  ];
					const int e2 = mcEdges[e*2+1];
					const Vec3 p1 = pos[ e1  ];    // scalar field pos 1
					const Vec3 p2 = pos[ e2  ];    // scalar field pos 2
					const float valp1  = value[ e1  ];  // scalar field val 1
					const float valp2  = value[ e2  ];  // scalar field val 2
					const float mu = (isoValue - valp1) / (valp2 - valp1);

					// init isolevel vertex
					Node vertex;
					vertex.pos = p1 + (p2-p1)*mu;
					vertex.normal = getNormalized( 
										getGradient( *this, i+cubieOffsetX[e1], j+cubieOffsetY[e1], k+cubieOffsetZ[e1]) * (1.0-mu) +
										getGradient( *this, i+cubieOffsetX[e2], j+cubieOffsetY[e2], k+cubieOffsetZ[e2]) * (    mu)) ;
					
					triIndices[e] = mesh.addNode(vertex) + 1;
					
					// store vertex 
					*eVert[e] = triIndices[e];
				} else
                {
					// retrieve  from vert array
					triIndices[e] = *eVert[e];
				}
			}
		}
		
		// Create the triangles... 
		for(int e = 0; mcTriTable[cubeIdx][e] != -1; e+=3)
        {
			mesh.addTri(Triangle(triIndices[ mcTriTable[cubeIdx][e+0]] - 1,
								 triIndices[ mcTriTable[cubeIdx][e+1]] - 1,
								 triIndices[ mcTriTable[cubeIdx][e+2]] - 1));
		}
	}
	
	//mesh.rebuildCorners();
	//mesh.rebuildLookup();
}

} //namespace
