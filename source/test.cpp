/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Use this file to test new functionality
 *
 ******************************************************************************/

#include "levelset.h"
#include "commonkernels.h"
#include "particle.h"
#include <cmath>

using namespace std;

namespace Manta
{

struct reductionTest : public KernelBase
{
    reductionTest(const Grid<Real>& v) : KernelBase(&v,0) ,v(v) ,sum(0)
    {
        run();
    }
    
    inline void op(int idx, const Grid<Real>& v ,double& sum)
    {
	    sum += v[idx];
    }
    
    inline operator double () { return sum; }
    inline double  & getRet() { return sum; }
    inline const Grid<Real>& getArg0() { return v; }
    typedef Grid<Real> type0;
    
    void run()
    {
        const int _sz = size; 
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
            double sum = 0;
#pragma omp for nowait
            for (int i=0; i < _sz; i++)
                op(i,v,sum);
#pragma omp critical
            {
                this->sum += sum;
            }
        }
    }
    
    const Grid<Real>& v;
    double sum;
};

struct minReduction : public KernelBase
{
    minReduction(const Grid<Real>& v) 
        : KernelBase(&v,0)
        , v(v)
        , sum(0)
    {
        run();
    }
    
    inline void op(int idx, const Grid<Real>& v ,double& sum)
    {
        if (sum < v[idx])
            sum = v[idx];   
    }

    inline operator double () { return sum; }
    inline double  & getRet() { return sum; }
    inline const Grid<Real>& getArg0() { return v; }
    typedef Grid<Real> type0;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
            double sum = 0;
#pragma omp for nowait
            for (int i=0; i < _sz; i++)
                op(i,v,sum);
#pragma omp critical
            {
                this->sum = min(sum, this->sum);
            }
        }
    }
    
    const Grid<Real>& v;
    double sum;
};

void getCurl(MACGrid& vel, Grid<Real>& vort, int comp)
{
	Grid<Vec3> velCenter(vel.getParent()), curl(vel.getParent());
	
	GetCentered(velCenter, vel);
	CurlOp(velCenter, curl);
	GetComponent(curl, vort, comp);
}

void setinflow(FlagGrid& flags, MACGrid& vel, LevelsetGrid& phi, Real h)
{
	FOR_IJK(vel)
    {
		if (i<=2)
        {
			if (j < h*flags.getSizeY())
            {
				vel(i,j,k).x = 1;            
				if (!flags.isObstacle(i,j,k))
                { 
					flags(i,j,k) = 1;        
					phi(i,j,k) = -1;
				}                
			} else
            {
				vel(i,j,k).x = 0;                            
				if (!flags.isObstacle(i,j,k))
                { 
					flags(i,j,k) = 4;
					phi(i,j,k) = 1;
				}
			}
		} else if (i>=flags.getSizeX()-2)
        {
			vel(i,j,k).x = 1;            
		}
	}
}
	
void testDiscardNth(BasicParticleSystem& parts, int skip=1)
{ 
	//knSetPdataConst<Real>(pd,value); 
	for (int i = 0; i < parts.size(); ++i)
    {
		if (i%(skip+1) == skip)
        { // keep 
		} else
        {
			parts.setPos(i, Vec3(-100000) );
		}
	}
}

} //namespace
