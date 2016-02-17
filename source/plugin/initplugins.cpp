/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Tools to setup fields and inflows
 *
 ******************************************************************************/

#include "initplugins.h"
#include "particle.h"
#include "simpleimage.h"
#include "mesh.h"

using namespace std;

namespace Manta
{
	
//! Apply noise to grid
struct KnApplyNoiseInfl : public KernelBase
{
    KnApplyNoiseInfl(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>& sdf, Real scale, Real sigma)
        : KernelBase(&flags,0)
        ,flags(flags)
        ,density(density)
        ,noise(noise)
        ,sdf(sdf)
        ,scale(scale)
        ,sigma(sigma)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>& sdf, Real scale, Real sigma)
    {
	    if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
	    Real factor = clamp(1.0-0.5/sigma * (sdf(i,j,k)+sigma), 0.0, 1.0);
	
	    Real target = noise.evaluate(Vec3(i,j,k)) * scale * factor;
	    if (density(i,j,k) < target)
		    density(i,j,k) = target;
    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline Grid<Real>& getArg1() { return density; }
    typedef Grid<Real> type1;
    inline WaveletNoiseField& getArg2() { return noise; }
    typedef WaveletNoiseField type2;
    inline Grid<Real>& getArg3() { return sdf; }
    typedef Grid<Real> type3;
    inline Real& getArg4() { return scale; }
    typedef Real type4;
    inline Real& getArg5() { return sigma; }
    typedef Real type5;
    
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
                            op(i,j,k,flags,density,noise,sdf,scale,sigma);
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
                        op(i,j,k,flags,density,noise,sdf,scale,sigma);
            }
        }
    }
    
    FlagGrid& flags;
    Grid<Real>& density;
    WaveletNoiseField& noise;
    Grid<Real>& sdf;
    Real scale;
    Real sigma;
};

//! Apply noise to real grid based on an SDF
struct KnAddNoise : public KernelBase
{
    KnAddNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>* sdf, Real scale)
        : KernelBase(&flags,0)
        , flags(flags)
        , density(density)
        , noise(noise)
        , sdf(sdf)
        , scale(scale)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>* sdf, Real scale)
    {
	    if (!flags.isFluid(i,j,k) || (sdf && (*sdf)(i,j,k) > 0.)) return;
	    density(i,j,k) += noise.evaluate(Vec3(i,j,k)) * scale;
    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline Grid<Real>& getArg1() { return density; }
    typedef Grid<Real> type1;
    inline WaveletNoiseField& getArg2() { return noise; }
    typedef WaveletNoiseField type2;
    inline Grid<Real>* getArg3() { return sdf; }
    typedef Grid<Real> type3;
    inline Real& getArg4() { return scale; }
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
                    for (int j=0; j < _maxY; j++)
                        for (int i=0; i < _maxX; i++)
                            op(i,j,k,flags,density,noise,sdf,scale);
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
                        op(i,j,k,flags,density,noise,sdf,scale);
            }
        }
    }
    
    FlagGrid& flags;
    Grid<Real>& density;
    WaveletNoiseField& noise;
    Grid<Real>* sdf;
    Real scale;
};

void addNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>* sdf=NULL, Real scale=1.0)
{
	KnAddNoise(flags, density, noise, sdf, scale );
}

//! sample noise field and set pdata with its values (for convenience, scale the noise values)
template <class T>
struct knSetPdataNoise : public KernelBase
{
    knSetPdataNoise(BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale)
        : KernelBase(parts.size())
        , parts(parts)
        , pdata(pdata)
        , noise(noise)
        , scale(scale)
    {
        run();
    }
    
    inline void op(int idx, BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale)
    {
	    pdata[idx] = noise.evaluate( parts.getPos(idx) ) * scale;
    }
    
    inline BasicParticleSystem& getArg0() { return parts; }
    typedef BasicParticleSystem type0;
    inline ParticleDataImpl<T>& getArg1() { return pdata; }
    typedef ParticleDataImpl<T> type1;
    inline WaveletNoiseField& getArg2() { return noise; }
    typedef WaveletNoiseField type2;
    inline Real& getArg3() { return scale; }
    typedef Real type3;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,parts,pdata,noise,scale);
        }
    }
    
    BasicParticleSystem& parts;
    ParticleDataImpl<T>& pdata;
    WaveletNoiseField& noise;
    Real scale;
};

template <class T>
struct knSetPdataNoiseVec : public KernelBase
{
    knSetPdataNoiseVec(BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale)
        : KernelBase(parts.size())
        , parts(parts)
        , pdata(pdata)
        , noise(noise)
        , scale(scale)
    {
        run();
    }
    
    inline void op(int idx, BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale)
    {
	    pdata[idx] = noise.evaluateVec( parts.getPos(idx) ) * scale;
    }
    
    inline BasicParticleSystem& getArg0() { return parts; }
    typedef BasicParticleSystem type0;
    inline ParticleDataImpl<T>& getArg1() { return pdata; }
    typedef ParticleDataImpl<T> type1;
    inline WaveletNoiseField& getArg2() { return noise; }
    typedef WaveletNoiseField type2;
    inline Real& getArg3() { return scale; }
    typedef Real type3;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,parts,pdata,noise,scale);
        }
    }
    
    BasicParticleSystem& parts;
    ParticleDataImpl<T>& pdata;
    WaveletNoiseField& noise;
    Real scale;
};

void setNoisePdata(BasicParticleSystem& parts, ParticleDataImpl<Real>& pd, WaveletNoiseField& noise, Real scale=1.)
{
    knSetPdataNoise<Real>(parts, pd,noise,scale);
}

void setNoisePdataVec3(BasicParticleSystem& parts, ParticleDataImpl<Vec3>& pd, WaveletNoiseField& noise, Real scale=1.)
{
    knSetPdataNoiseVec<Vec3>(parts, pd,noise,scale);
}

void setNoisePdataInt(BasicParticleSystem& parts, ParticleDataImpl<int >& pd, WaveletNoiseField& noise, Real scale=1.)
{
    knSetPdataNoise<int>(parts, pd,noise,scale);
}

//! SDF gradient from obstacle flags
Grid<Vec3> obstacleGradient(FlagGrid& flags)
{
	LevelsetGrid levelset(flags.getParent(),false);
	Grid<Vec3> gradient(flags.getParent());
	
	// rebuild obstacle levelset
	FOR_IDX(levelset)
    {
		levelset[idx] = flags.isObstacle(idx) ? -0.5 : 0.5;
	}

	levelset.reinitMarching(flags, 6.0, 0, true, false, FlagGrid::TypeReserved);
	
	// build levelset gradient
	GradientOp(gradient, levelset);
	
	FOR_IDX(levelset)
    {
		Vec3 grad = gradient[idx];
		Real s = normalize(grad);
		if (s <= 0.1 || levelset[idx] >= 0) 
			grad=Vec3(0.);        
		gradient[idx] = grad * levelset[idx];
	}
	
	return gradient;
}

LevelsetGrid obstacleLevelset(FlagGrid& flags)
{
    LevelsetGrid levelset(flags.getParent(),false);
    Grid<Vec3> gradient(flags.getParent());

    // rebuild obstacle levelset
    FOR_IDX(levelset)
    {
	    levelset[idx] = flags.isObstacle(idx) ? -0.5 : 0.5;
    }

    levelset.reinitMarching(flags, 6.0, 0, true, false, FlagGrid::TypeReserved);

    return levelset;
}



//*****************************************************************************
// blender init functions 
//*****************************************************************************
struct KnApplyEmission : public KernelBase
{
    KnApplyEmission(FlagGrid& flags, Grid<Real>& density, Grid<Real>& emission, bool isAbsolute)
        : KernelBase(&flags,0)
        , flags(flags)
        , density(density)
        , emission(emission)
        , isAbsolute(isAbsolute)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& density, Grid<Real>& emission, bool isAbsolute)
    {
	    if (!flags.isFluid(i,j,k) || emission(i,j,k) == 0.) return;
	    if (isAbsolute)
		    density(i,j,k) = emission(i,j,k);
	    else
		    density(i,j,k) += emission(i,j,k);
    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline Grid<Real>& getArg1() { return density; }
    typedef Grid<Real> type1;
    inline Grid<Real>& getArg2() { return emission; }
    typedef Grid<Real> type2;
    inline bool& getArg3() { return isAbsolute; }
    typedef bool type3;
    
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
                            op(i,j,k,flags,density,emission,isAbsolute);
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
                        op(i,j,k,flags,density,emission,isAbsolute);
            }
        }
    }
    
    FlagGrid& flags;
    Grid<Real>& density;
    Grid<Real>& emission;
    bool isAbsolute;
};

//! Add emission values
//isAbsolute: whether to add emission values to existing, or replace
void applyEmission(FlagGrid& flags, Grid<Real>& density, Grid<Real>& emission, bool isAbsolute)
{
	KnApplyEmission(flags, density, emission, isAbsolute);
}

// blender init functions for meshes
struct KnApplyDensity : public KernelBase
{
    KnApplyDensity(FlagGrid& flags, Grid<Real>& density, Grid<Real>& sdf, Real value, Real sigma)
        : KernelBase(&flags,0)
        , flags(flags)
        , density(density)
        , sdf(sdf)
        , value(value)
        , sigma(sigma)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& density, Grid<Real>& sdf, Real value, Real sigma)
    {
	    if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
	    density(i,j,k) = value;
    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline Grid<Real>& getArg1() { return density; }
    typedef Grid<Real> type1;
    inline Grid<Real>& getArg2() { return sdf; }
    typedef Grid<Real> type2;
    inline Real& getArg3() { return value; }
    typedef Real type3;
    inline Real& getArg4() { return sigma; }
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
                    for (int j=0; j < _maxY; j++)
                        for (int i=0; i < _maxX; i++)
                            op(i,j,k,flags,density,sdf,value,sigma);
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
                        op(i,j,k,flags,density,sdf,value,sigma);
            }
        }
    }
    
    FlagGrid& flags;
    Grid<Real>& density;
    Grid<Real>& sdf;
    Real value;
    Real sigma;
};

void densityInflow(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Shape* shape, Real scale, Real sigma)
{
	Grid<Real> sdf = shape->computeLevelset();
	KnApplyNoiseInfl(flags, density, noise, sdf, scale, sigma);
}

void densityInflowMeshNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Mesh* mesh, Real scale, Real sigma)
{
	LevelsetGrid sdf(density.getParent(), false);
	mesh->computeLevelset(sdf, 1.);
	KnApplyNoiseInfl(flags, density, noise, sdf, scale, sigma);
}

void densityInflowMesh(FlagGrid& flags, Grid<Real>& density, Mesh* mesh, Real value, Real cutoff, Real sigma)
{
	LevelsetGrid sdf(density.getParent(), false);
	mesh->computeLevelset(sdf, 2., cutoff);
	KnApplyDensity(flags, density, sdf, value, sigma);
}

//! check for symmetry , optionally enfore by copying
void checkSymmetry(Grid<Real>& a, Grid<Real>* err = nullptr, bool symmetrize = false, int axis = 0, int bound = 0)
{
	const int c  = axis; 
	const int s = a.getSize()[c];
	FOR_IJK(a)
    { 
		Vec3i idx(i,j,k), mdx(i,j,k);
		mdx[c] = s-1-idx[c];
		if (bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound)))) continue;

		if (err) (*err)(idx) = fabs((double)(a(idx) - a(mdx))); 
		if (symmetrize && (idx[c]<s/2))
        {
			a(idx) = a(mdx);
		}
	}
}

//! check for symmetry, mac grid version
void checkSymmetryVec3( Grid<Vec3>& a, Grid<Real>* err = nullptr, bool symmetrize = false, int axis = 0, int bound = 0, int disable = 0)
{
	if(err) err->setConst(0.);

	// each dimension is measured separately for flexibility (could be combined)
	const int c  = axis;
	const int o1 = (c+1)%3;
	const int o2 = (c+2)%3;

	// x
	if (!(disable&1))
    {
		const int s = a.getSize()[c]+1; 
		FOR_IJK(a)
        { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if (mdx[c] >= a.getSize()[c]) continue; 
			if (bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound)))) continue;

			// special case: center "line" of values , should be zero!
			if (mdx[c] == idx[c])
            {
				if (err) (*err)(idx) += fabs((double)( a(idx)[c]));
				if (symmetrize) a(idx)[c] = 0.;
				continue; 
			}

			// note - the a(mdx) component needs to be inverted here!
			if (err) (*err)(idx) += fabs((double)( a(idx)[c]- (a(mdx)[c]*-1.)));
			if (symmetrize && (idx[c]<s/2))
            {
				a(idx)[c] = a(mdx)[c] * -1.;
			}
		}
	}

	// y
	if (!(disable&2))
    {
		const int s = a.getSize()[c];
		FOR_IJK(a)
        { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if (bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound)))) continue;

			if (err) (*err)(idx) += fabs((double)( a(idx)[o1]-a(mdx)[o1])); 
			if (symmetrize && (idx[c]<s/2))
            {
				a(idx)[o1] = a(mdx)[o1];
			}
		}
	} 

	// z
	if (!(disable&4))
    {
		const int s = a.getSize()[c];
		FOR_IJK(a)
        { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if (bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound)))) continue;

			if (err) (*err)(idx) += fabs((double)( a(idx)[o2]-a(mdx)[o2]));
			if (symmetrize && (idx[c]<s/2))
            {
				a(idx)[o2] = a(mdx)[o2];
			}
		}
	} 

}

//! project data onto a plane, and write ppm
void projectPpmOut(Grid<Real>& val, string name, int axis = 2, Real scale = 1.)
{
	const int c  = axis; 
	const int o1 = (c+1)%3;
	const int o2 = (c+2)%3;

	SimpleImage img;
	img.init(val.getSize()[o1],  val.getSize()[o2]);
	Real s = 1. / (Real)val.getSize()[c];

	FOR_IJK(val)
    { 
		Vec3i idx(i,j,k); 
		//if(idx[c]==val.getSize()[c]/2) img( idx[o1], idx[o2] ) = val(idx); 
		img( idx[o1], idx[o2] ) += s * val(idx);
	}

	img.mapRange( 1./scale );
	img.writePpm( name );
}

struct KnPrecompLight : public KernelBase
{
    KnPrecompLight(Grid<Real>& density, Grid<Real>& L, Vec3 light = Vec3(1, 1, 1))
        : KernelBase(&density,0)
        , density(density)
        , L(L)
        , light(light)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<Real>& density, Grid<Real>& L, Vec3 light = Vec3(1, 1, 1))
    {
	    Vec3 n = getGradient( density, i,j,k ) * -1.; 
	    normalize(n);

	    Real d = dot( light, n );
	    L(i,j,k) = d;
    }
    
    inline Grid<Real>& getArg0() { return density; }
    typedef Grid<Real> type0;
    inline Grid<Real>& getArg1() { return L; }
    typedef Grid<Real> type1;
    inline Vec3& getArg2() { return light; }
    typedef Vec3 type2;
    
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
                            op(i,j,k,density,L,light);
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
                        op(i,j,k,density,L,light);
            }
        }
    }
    
    Grid<Real>& density;
    Grid<Real>& L;
    Vec3 light;
};

// simple shading with pre-computed gradient
static inline void shadeCell(Vec3& dst, int shadeMode, Real src, Real light, int depthPos, Real depthInv)
{	
	switch(shadeMode)
    {

	case 1:
    {
		// surfaces
		Vec3 ambient = Vec3(0.1,0.1,0.1);
		Vec3 diffuse = Vec3(0.9,0.9,0.9); 
		Real alpha = src; 

		// different color for depth?
		diffuse[0] *= ((Real)depthPos * depthInv) * 0.7 + 0.3;
		diffuse[1] *= ((Real)depthPos * depthInv) * 0.7 + 0.3;

		Vec3 col = ambient + diffuse * light; 

		//img( 0+i, j ) = (1.-alpha) * img( 0+i, j ) + alpha * col;
		dst = (1.-alpha) * dst + alpha * col;
	} break;

	default:
    {
		// volumetrics / smoke
		dst += depthInv * Vec3(src,src,src);
	} break;

	}
}

//! output shaded (all 3 axes at once for 3D)
//! shading modes: 0 smoke, 1 surfaces
void projectPpmFull( Grid<Real>& val, string name, int shadeMode=0, Real scale=1.)
{
	Vec3i s  = val.getSize();
	Vec3  si = Vec3( 1. / (Real)s[0], 1. / (Real)s[1], 1. / (Real)s[2] );

	SimpleImage img;
	int imgSx = s[0];
	if (val.is3D()) imgSx += s[2]+s[0]; // mult views in 3D
	img.init( imgSx, std::max(s[0], std::max( s[1],s[2])) );

	// precompute lighting
	Grid<Real> L(val);
	KnPrecompLight( val, L , Vec3(1,1,1) );

	FOR_IJK(val)
    { 
		Vec3i idx(i,j,k);
		// img( 0+i, j ) += si[2] * val(idx); // averaging
		shadeCell( img( 0+i, j ) , shadeMode, val(idx), L(idx), k, si[2]);
	}

	if (val.is3D())
    {
	    FOR_IJK(val)
        { 
		    Vec3i idx(i,j,k);
		    //img( s[0]+k, j ) += si[0] * val(idx);
		    shadeCell( img( s[0]+k, j ) , shadeMode, val(idx), L(idx), i, si[0]);
	    }

	    FOR_IJK(val)
        { 
		    Vec3i idx(i,j,k);
		    //img( s[0]+s[2]+i, k ) += si[1] * val(idx);
		    shadeCell( img( s[0]+s[2]+i, k ) , shadeMode, val(idx), L(idx), j, si[1]);
	    }
	} // 3d

	img.mapRange( 1./scale );
	img.writePpm( name );
}

// helper functions for pdata operator tests
//! init some test particles at the origin
void addTestParts(BasicParticleSystem& parts, int num)
{
	for (int i = 0; i < num; ++i)
		parts.addBuffered(Vec3(0,0,0));

	parts.doCompress();
	parts.insertBufferedParticles();
}

// calculate the difference between two pdata fields (note - slow!, not parallelized)
Real pdataMaxDiff(ParticleDataBase* a, ParticleDataBase* b)
{    
	double maxVal = 0.;
	//debMsg(" PD "<< a->getType()<<"  as"<<a->getSizeSlow()<<"  bs"<<b->getSizeSlow() , 1);
	assertMsg(a->getType()     == b->getType()    , "pdataMaxDiff problem - different pdata types!");
	assertMsg(a->getSizeSlow() == b->getSizeSlow(), "pdataMaxDiff problem - different pdata sizes!");
	
	if (a->getType() & ParticleDataBase::TypeReal) 
	{
		ParticleDataImpl<Real>& av = *dynamic_cast<ParticleDataImpl<Real>*>(a);
		ParticleDataImpl<Real>& bv = *dynamic_cast<ParticleDataImpl<Real>*>(b);
		FOR_PARTS(av) 
        {
			maxVal = std::max(maxVal, (double)fabs(av[idx]-bv[idx]));
		}
	} else if (a->getType() & ParticleDataBase::TypeInt) 
	{
		ParticleDataImpl<int>& av = *dynamic_cast<ParticleDataImpl<int>*>(a);
		ParticleDataImpl<int>& bv = *dynamic_cast<ParticleDataImpl<int>*>(b);
		FOR_PARTS(av)
        {
			maxVal = std::max(maxVal, (double)fabs((double)av[idx]-bv[idx]));
		}
	} else if (a->getType() & ParticleDataBase::TypeVec3)
    {
		ParticleDataImpl<Vec3>& av = *dynamic_cast<ParticleDataImpl<Vec3>*>(a);
		ParticleDataImpl<Vec3>& bv = *dynamic_cast<ParticleDataImpl<Vec3>*>(b);
		FOR_PARTS(av)
        {
			double d = 0.;
			for(int c=0; c<3; ++c)
            { 
				d += fabs((double)av[idx][c] - (double)bv[idx][c]);
			}
			maxVal = std::max(maxVal, d);
		}
	} else
    {
		errMsg("pdataMaxDiff: Grid Type is not supported (only Real, Vec3, int)");    
	}

	return maxVal;
}



//*****************************************************************************
// helper functions for volume fractions
//*****************************************************************************
struct kninitVortexVelocity : public KernelBase
{
    kninitVortexVelocity(Grid<Real> &phiObs, MACGrid& vel, const Vec3 &center, const Real &radius)
        : KernelBase(&phiObs,0)
        , phiObs(phiObs)
        , vel(vel)
        , center(center)
        , radius(radius)
    
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<Real> &phiObs, MACGrid& vel, const Vec3 &center, const Real &radius)
    {
	    if(phiObs(i,j,k) >= -1.)
        {
		    Real dx = i - center.x; if(dx>=0) dx -= .5; else dx += .5;
		    Real dy = j - center.y;
		    Real r = std::sqrt(dx*dx+dy*dy);
		    Real alpha = atan2(dy,dx);

		    vel(i,j,k).x = -std::sin(alpha)*(r/radius);

		    dx = i - center.x;
		    dy = j - center.y; if(dy>=0) dy -= .5; else dy += .5;
		    r = std::sqrt(dx*dx+dy*dy);
		    alpha = atan2(dy,dx);

		    vel(i,j,k).y = std::cos(alpha)*(r/radius);
	    }
    }

    inline Grid<Real> & getArg0() { return phiObs; }
    typedef Grid<Real>  type0;
    inline MACGrid& getArg1() { return vel; }
    typedef MACGrid type1;
    inline const Vec3& getArg2() { return center; }
    typedef Vec3 type2;
    inline const Real& getArg3() { return radius; }
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
                    for (int j=0; j < _maxY; j++)
                        for (int i=0; i < _maxX; i++)
                            op(i,j,k,phiObs,vel,center,radius);
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
                        op(i,j,k,phiObs,vel,center,radius);
            }
        }
    }
    
    Grid<Real> & phiObs;
    MACGrid& vel;
    const Vec3& center;
    const Real& radius;
};

void initVortexVelocity(Grid<Real> &phiObs, MACGrid& vel, const Vec3 &center, const Real &radius)
{
	kninitVortexVelocity(phiObs,  vel, center, radius);
}

inline static Real calcFraction(Real phi1, Real phi2)
{
	if (phi1>0. && phi2>0.) return 1.;
	if (phi1<0. && phi2<0.) return 0.;

	// make sure phi1 < phi2
	if (phi2<phi1)
    {
        Real t = phi1;
        phi1= phi2;
        phi2 = t;
    }

	Real denom = phi1-phi2;
	if (denom > -1e-04) return 0.5; 

	Real frac = 1. - phi1/denom;
	if (frac<0.01) frac = 0.; // skip , dont mark as fluid
	return std::max(Real(0), std::min(Real(1), frac));
}

struct KnUpdateFractions : public KernelBase
{
    KnUpdateFractions(FlagGrid& flags, Grid<Real>& phiObs, MACGrid& fractions, const int &boundaryWidth)
        : KernelBase(&flags,1)
        , flags(flags)
        , phiObs(phiObs)
        , fractions(fractions)
        , boundaryWidth(boundaryWidth)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& phiObs, MACGrid& fractions, const int &boundaryWidth)
    {
	    // walls at domain bounds and inner objects
	    fractions(i,j,k).x = calcFraction( phiObs(i,j,k) , phiObs(i-1,j,k));
	    fractions(i,j,k).y = calcFraction( phiObs(i,j,k) , phiObs(i,j-1,k));
        if (phiObs.is3D())
        {
	        fractions(i,j,k).z = calcFraction( phiObs(i,j,k) , phiObs(i,j,k-1));
	    }

	    // remaining BCs at the domain boundaries 
	    const int w = boundaryWidth;
	    // only set if not in obstacle
 	    if (phiObs(i,j,k)<0.) return;

	    // x-direction boundaries
	    if (i <= w+1)
        {                     //min x
		    if( (flags.isInflow(i-1,j,k))  ||
			    (flags.isOutflow(i-1,j,k)) ||
			    (flags.isOpen(i-1,j,k)))
            {
                fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		    }
	    }

	    if (i >= flags.getSizeX()-w-2)
        {    //max x
		    if ( (flags.isInflow(i+1,j,k))  ||
			     (flags.isOutflow(i+1,j,k)) ||
			     (flags.isOpen(i+1,j,k)))
            {
			    fractions(i+1,j,k).x = fractions(i+1,j,k).y = 1.; if(flags.is3D()) fractions(i+1,j,k).z = 1.;
		    }
	    }

	    // y-direction boundaries
 	    if (j <= w+1)
        {                     //min y
		    if(	(flags.isInflow(i,j-1,k))  ||
			    (flags.isOutflow(i,j-1,k)) ||
			    (flags.isOpen(i,j-1,k)))
            {
			    fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		    }
 	    }

 	    if (j >= flags.getSizeY()-w-2)
        {      //max y
		    if(	(flags.isInflow(i,j+1,k))  ||
			    (flags.isOutflow(i,j+1,k)) ||
			    (flags.isOpen(i,j+1,k)))
            {
			    fractions(i,j+1,k).x = fractions(i,j+1,k).y = 1.; if(flags.is3D()) fractions(i,j+1,k).z = 1.;
		    }
 	    }

	    // z-direction boundaries
	    if (flags.is3D()) 
        {
	        if(k <= w+1)
            {                 //min z
		        if(	(flags.isInflow(i,j,k-1))  ||
			        (flags.isOutflow(i,j,k-1)) ||
			        (flags.isOpen(i,j,k-1)) )
                {
			        fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		        }
	        }

	        if (j >= flags.getSizeZ()-w-2)
            { //max z
		        if(	(flags.isInflow(i,j,k+1))  ||
			        (flags.isOutflow(i,j,k+1)) ||
			        (flags.isOpen(i,j,k+1)) )
                {
			        fractions(i,j,k+1).x = fractions(i,j,k+1).y = 1.; if(flags.is3D()) fractions(i,j,k+1).z = 1.;
		        }
	        }
	    }

    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline Grid<Real>& getArg1() { return phiObs; }
    typedef Grid<Real> type1;
    inline MACGrid& getArg2() { return fractions; }
    typedef MACGrid type2;
    inline const int& getArg3() { return boundaryWidth; }
    typedef int type3;
    
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
                            op(i,j,k,flags,phiObs,fractions,boundaryWidth);
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
                        op(i,j,k,flags,phiObs,fractions,boundaryWidth);
            }
        }
    }
    
    FlagGrid& flags;
    Grid<Real>& phiObs;
    MACGrid& fractions;
    const int& boundaryWidth;
};

void updateFractions(FlagGrid& flags, Grid<Real>& phiObs, MACGrid& fractions, const int &boundaryWidth)
{
	fractions.setConst(Vec3(0.));
	KnUpdateFractions(flags, phiObs, fractions, boundaryWidth);
}

struct KnUpdateFlags : public KernelBase
{
    KnUpdateFlags(FlagGrid& flags, MACGrid& fractions, Grid<Real>& phiObs)
        : KernelBase(&flags,1)
        , flags(flags)
        , fractions(fractions)
        , phiObs(phiObs)
    {
        run();
    }
    
    inline void op(int i, int j, int k, FlagGrid& flags, MACGrid& fractions, Grid<Real>& phiObs)
    {
	    Real test = 0.;
	    test += fractions.get(i  ,j,k).x;
	    test += fractions.get(i+1,j,k).x;
	    test += fractions.get(i,j  ,k).y;
	    test += fractions.get(i,j+1,k).y;

	    if (flags.is3D())
        {
	        test += fractions.get(i,j,k  ).z;
	        test += fractions.get(i,j,k+1).z;
        }

	    if (test==0. && phiObs(i,j,k) < 0.) flags(i,j,k) = FlagGrid::TypeObstacle; 
	    else flags(i,j,k) = FlagGrid::TypeEmpty; 
    }
    
    inline FlagGrid& getArg0() { return flags; }
    typedef FlagGrid type0;
    inline MACGrid& getArg1() { return fractions; }
    typedef MACGrid type1;
    inline Grid<Real>& getArg2() { return phiObs; }
    typedef Grid<Real> type2;
    
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
                            op(i,j,k,flags,fractions,phiObs);
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
                        op(i,j,k,flags,fractions,phiObs);
            }
        }
    }
    
    FlagGrid& flags;
    MACGrid& fractions;
    Grid<Real>& phiObs;
};

void setObstacleFlags(FlagGrid& flags, MACGrid& fractions, Grid<Real>& phiObs)
{
	KnUpdateFlags(flags,fractions, phiObs);
}

} // namespace
