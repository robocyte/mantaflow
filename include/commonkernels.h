/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Common grid kernels
 *
 ******************************************************************************/

#pragma once

#include "general.h"
#include "kernel.h"
#include "grid.h"

namespace Manta
{
   
//! Kernel: Invert real values, if positive and fluid
struct InvertCheckFluid : public KernelBase
{
    InvertCheckFluid(FlagGrid& flags, Grid<Real>& grid) : KernelBase(&flags,0) ,flags(flags),grid(grid)
    {
        run();
    }

    inline void op(int idx, FlagGrid& flags, Grid<Real>& grid )
    {
        if (flags.isFluid(idx) && grid[idx] > 0)
            grid[idx] = 1.0 / grid[idx];
    }

    void run()
    {
        const int _sz = size; 
#pragma omp parallel 
        {
            this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  
#pragma omp for
            for (int i=0; i < _sz; i++) op(i,flags,grid);
        }
    }
     
    FlagGrid& flags;
    Grid<Real>& grid;
 };

//! Kernel: Squared sum over grid
struct GridSumSqr : public KernelBase
{
    GridSumSqr(Grid<Real>& grid) : KernelBase(&grid,0) ,grid(grid) ,sum(0)
    {
        run();
    }
    
    inline void op(int idx, Grid<Real>& grid ,double& sum)
    {
        sum += square((double)grid[idx]);
    }
    
    void run()
    {
        const int _sz = size; 
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();  double sum = 0;
#pragma omp for nowait
            for (int i=0; i < _sz; i++) op(i,grid,sum);
#pragma omp critical
            {
                this->sum += sum;
            }
        }
    }
    
    Grid<Real>& grid;  double sum;
};

//! Kernel: rotation operator \nabla x v for centered vector fields
struct CurlOp : public KernelBase
{
    CurlOp(const Grid<Vec3>& grid, Grid<Vec3>& dst) : KernelBase(&grid,1) ,grid(grid),dst(dst)
    {
        run();
    }
    
    inline void op(int i, int j, int k, const Grid<Vec3>& grid, Grid<Vec3>& dst )
    {
        Vec3 v = Vec3(0. , 0. , 0.5*((grid(i+1,j,k).y - grid(i-1,j,k).y) - (grid(i,j+1,k).x - grid(i,j-1,k).x)) );
        if(dst.is3D())
        {
            v[0] = 0.5*((grid(i,j+1,k).z - grid(i,j-1,k).z) - (grid(i,j,k+1).y - grid(i,j,k-1).y));
            v[1] = 0.5*((grid(i,j,k+1).x - grid(i,j,k-1).x) - (grid(i+1,j,k).z - grid(i-1,j,k).z));
        }
        dst(i,j,k) = v;
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
                        for (int i=1; i < _maxX; i++) op(i,j,k,grid,dst);
            }
        } else
        {
            const int k=0;
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int j=1; j < _maxY; j++)
                    for (int i=1; i < _maxX; i++) op(i,j,k,grid,dst);
            }
        }
    }
    
    const Grid<Vec3>& grid; Grid<Vec3>& dst;
};

//! Kernel: divergence operator (from MAC grid)
struct DivergenceOpMAC : public KernelBase
{
    DivergenceOpMAC(Grid<Real>& div, const MACGrid& grid) : KernelBase(&div,1) ,div(div),grid(grid)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<Real>& div, const MACGrid& grid )
    {
        Vec3 del = Vec3(grid(i+1,j,k).x, grid(i,j+1,k).y, 0.) - grid(i,j,k); 
        if(grid.is3D()) del[2] += grid(i,j,k+1).z;
        else            del[2]  = 0.;
        div(i,j,k) = del.x + del.y + del.z;
    }
    
    void run()
    {
        const int _maxX = maxX;
        const int _maxY = maxY;
        if (maxZ > 1)
        {
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int k=minZ; k < maxZ; k++)
                    for (int j=1; j < _maxY; j++)
                        for (int i=1; i < _maxX; i++)
                            op(i,j,k,div,grid);
            }
        } else
        {
            const int k=0;
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int j=1; j < _maxY; j++) 
                    for (int i=1; i < _maxX; i++) op(i,j,k,div,grid);
            }
        }
    }
    
    Grid<Real>& div;
    const MACGrid& grid;
};

//! Kernel: gradient operator for MAC grid
struct GradientOpMAC : public KernelBase
{
    GradientOpMAC(MACGrid& gradient, const Grid<Real>& grid) : KernelBase(&gradient,1) ,gradient(gradient),grid(grid)
    {
        run();
    }
    
    inline void op(int i, int j, int k, MACGrid& gradient, const Grid<Real>& grid )
    {
        Vec3 grad = (Vec3(grid(i,j,k)) - Vec3(grid(i-1,j,k), grid(i,j-1,k), 0. ));
        if(grid.is3D()) grad[2] -= grid(i,j,k-1);
        else            grad[2]  = 0.;
        gradient(i,j,k) = grad;
    }
    
    void run()
    {
        const int _maxX = maxX;
        const int _maxY = maxY;
        
        if (maxZ > 1)
        { 
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int k=minZ; k < maxZ; k++)
                    for (int j=1; j < _maxY; j++)
                        for (int i=1; i < _maxX; i++)
                            op(i,j,k,gradient,grid);
            }
        } else
        {
            const int k=0;
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int j=1; j < _maxY; j++)
                    for (int i=1; i < _maxX; i++) op(i,j,k,gradient,grid);
            }
        }
    }
    
    MACGrid& gradient;
    const Grid<Real>& grid;
};

//! Kernel: centered gradient operator 
struct GradientOp : public KernelBase
{
    GradientOp(Grid<Vec3>& gradient, const Grid<Real>& grid) : KernelBase(&gradient,1)
        , gradient(gradient)
        , grid(grid)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<Vec3>& gradient, const Grid<Real>& grid )
    {
        Vec3 grad = 0.5 * Vec3( grid(i+1,j,k)-grid(i-1,j,k), 
								grid(i,j+1,k)-grid(i,j-1,k), 0.);
        if(grid.is3D()) grad[2]= 0.5*( grid(i,j,k+1)-grid(i,j,k-1) );
        gradient(i,j,k) = grad;
    }
    
    void run()
    {
        const int _maxX = maxX;
        const int _maxY = maxY;
        if (maxZ > 1)
        {
#pragma omp parallel
            {
                this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();
#pragma omp for
                for (int k=minZ; k < maxZ; k++)
                    for (int j=1; j < _maxY; j++)
                        for (int i=1; i < _maxX; i++)
                            op(i,j,k,gradient,grid);
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
                        op(i,j,k,gradient,grid);
            }
        }
    }
    
    Grid<Vec3>& gradient;
    const Grid<Real>& grid;
};

//! Kernel: Laplace operator
struct LaplaceOp : public KernelBase
{
    LaplaceOp(Grid<Real>& laplace, const Grid<Real>& grid) : KernelBase(&laplace,1)
        , laplace(laplace)
        , grid(grid)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<Real>& laplace, const Grid<Real>& grid)
    {
        laplace(i,j,k) = -(6.0*grid(i,j,k)-grid(i+1,j,k)-grid(i-1,j,k)-grid(i,j+1,k)-grid(i,j-1,k)-grid(i,j,k+1)-grid(i,j,k-1));
    }
    
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
                        for (int i=1; i < _maxX; i++) op(i,j,k,laplace,grid);
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
                        op(i,j,k,laplace,grid);
            }
        }
    }
    
    Grid<Real>& laplace;
    const Grid<Real>& grid;
};

//! Kernel: get component at MAC positions
struct GetShiftedComponent : public KernelBase
{
    GetShiftedComponent(const Grid<Vec3>& grid, Grid<Real>& comp, int dim) : KernelBase(&grid,1)
        , grid(grid)
        , comp(comp)
        , dim(dim)
    {
        run();
    }
    
    inline void op(int i, int j, int k, const Grid<Vec3>& grid, Grid<Real>& comp, int dim )
    {
        Vec3i ishift(i,j,k);
        ishift[dim]--;
        comp(i,j,k) = 0.5*(grid(i,j,k)[dim] + grid(ishift)[dim]);
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
                            op(i,j,k,grid,comp,dim);
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
                        op(i,j,k,grid,comp,dim);
            }
        }
    }
    
    const Grid<Vec3>& grid;
    Grid<Real>& comp;
    int dim;
};

//! Kernel: get component (not shifted)
struct GetComponent : public KernelBase
{
    GetComponent(const Grid<Vec3>& grid, Grid<Real>& comp, int dim) : KernelBase(&grid,0)
        , grid(grid)
        , comp(comp)
        , dim(dim)
    {
        run();
    }
    
    inline void op(int idx, const Grid<Vec3>& grid, Grid<Real>& comp, int dim )
    {
        comp[idx] = grid[idx][dim];    
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
                op(i,grid,comp,dim);
        }
    }
    
    const Grid<Vec3>& grid;
    Grid<Real>& comp;
    int dim;
};

//! Kernel: get norm of centered grid
struct GridNorm : public KernelBase
{
    GridNorm(Grid<Real>& n, const Grid<Vec3>& grid) : KernelBase(&n,0)
        , n(n)
        , grid(grid)
    {
        run();
    }
    
    inline void op(int idx, Grid<Real>& n, const Grid<Vec3>& grid )
    {
        n[idx] = norm(grid[idx]);
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
                op(i,n,grid);
        }
    }
    
    Grid<Real>& n;
    const Grid<Vec3>& grid;
};

//! Kernel: set component (not shifted)
struct SetComponent : public KernelBase
{
    SetComponent(Grid<Vec3>& grid, const Grid<Real>& comp, int dim) : KernelBase(&grid,0)
        , grid(grid)
        , comp(comp)
        , dim(dim)
    {
        run();
    }
    
    inline void op(int idx, Grid<Vec3>& grid, const Grid<Real>& comp, int dim)
    {
        grid[idx][dim] = comp[idx];
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
                op(i,grid,comp,dim);
        }
    }
    
    Grid<Vec3>& grid;
    const Grid<Real>& comp;
    int dim;
};

//! Kernel: compute centered velocity field from MAC
struct GetCentered : public KernelBase
{
    GetCentered(Grid<Vec3>& center, const MACGrid& vel) : KernelBase(&center,1)
        , center(center)
        , vel(vel)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<Vec3>& center, const MACGrid& vel)
    {
        Vec3 v = 0.5 * ( vel(i,j,k) + Vec3(vel(i+1,j,k).x, vel(i,j+1,k).y, 0. ) );
        if(vel.is3D()) v[2] += 0.5 * vel(i,j,k+1).z;
        else           v[2]  = 0.;
        center(i,j,k) = v;
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
                            op(i,j,k,center,vel);
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
                        op(i,j,k,center,vel);
            }
        }
    }
    
    Grid<Vec3>& center;
    const MACGrid& vel;
};

//! Kernel: compute MAC from centered velocity field
struct GetMAC : public KernelBase
{
    GetMAC(MACGrid& vel, const Grid<Vec3>& center) : KernelBase(&vel,1)
        , vel(vel)
        , center(center)
    {
        run();
    }
    
    inline void op(int i, int j, int k, MACGrid& vel, const Grid<Vec3>& center )
    {
        Vec3 v = 0.5*(center(i,j,k) + Vec3(center(i-1,j,k).x, center(i,j-1,k).y, 0. ));
        if(vel.is3D()) v[2] += 0.5 * center(i,j,k-1).z;
        else           v[2]  = 0.;
        vel(i,j,k) = v;
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
                            op(i,j,k,vel,center);
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
                        op(i,j,k,vel,center);
            }
        }
    }
    
    MACGrid& vel;
    const Grid<Vec3>& center;
};

//! Fill in the domain boundary cells (i,j,k=0/size-1) from the neighboring cells
struct FillInBoundary : public KernelBase
{
    FillInBoundary(Grid<Vec3>& grid, int g) : KernelBase(&grid,0)
        , grid(grid)
        , g(g)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<Vec3>& grid, int g)
    {
        if (i==0) grid(i,j,k) = grid(i+1,j,k);
        if (j==0) grid(i,j,k) = grid(i,j+1,k);
        if (k==0) grid(i,j,k) = grid(i,j,k+1);
        if (i==grid.getSizeX()-1) grid(i,j,k) = grid(i-1,j,k);
        if (j==grid.getSizeY()-1) grid(i,j,k) = grid(i,j-1,k);
        if (k==grid.getSizeZ()-1) grid(i,j,k) = grid(i,j,k-1);
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
                            op(i,j,k,grid,g);
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
                        op(i,j,k,grid,g);
            }
        }
    }
    
    Grid<Vec3>& grid;
    int g;
};

} // namespace
