/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * shapes classes
 *
 ******************************************************************************/

#pragma once

#include "manta.h"
#include "vectorbase.h"
#include "levelset.h"

namespace Manta
{

class Mesh;



//! Base class for all shapes
class Shape : public PbClass
{
public:
	enum GridType { TypeNone = 0, TypeBox = 1, TypeSphere = 2, TypeCylinder = 3, TypeSlope = 4 };
	
	Shape(FluidSolver* parent);
	
	//! Get the type of grid
	inline GridType getType() const { return mType; }
	
	//! Apply shape to flag grid, set inside cells to <value>
    template<typename T>
	void applyToGrid(Grid<T>* grid, T value, FlagGrid* respectFlags = nullptr)
    {
        ApplyShapeToGrid<T>(grid, this, value, respectFlags);
    }

    void applyToMACGrid(MACGrid* grid, Vec3 value, FlagGrid* respectFlags);
    void applyToGridSmooth(GridBase* grid, Real sigma = 1.0, Real shift = 0, FlagGrid* respectFlags = 0);

    LevelsetGrid computeLevelset();
	void collideMesh(Mesh& mesh);

	virtual Vec3 getCenter() const              { return Vec3::Zero; }
	virtual void setCenter(const Vec3& center)  { }
	virtual Vec3 getExtent() const              { return Vec3::Zero; }
	
	//! Inside test of the shape
	virtual bool isInside(const Vec3& pos) const;
	inline  bool isInsideGrid(int i, int j, int k) const { return isInside(Vec3(i+0.5,j+0.5,k+0.5)); };
	
	virtual void generateMesh(Mesh* mesh) {} ;    
	virtual void generateLevelset(Grid<Real>& phi) {};    
	
protected:
    GridType mType;
};

//! Dummy shape
class NullShape : public Shape
{   
public:
	NullShape(FluidSolver* parent)
        : Shape(parent)
    {}
	
	virtual bool isInside(const Vec3& pos) const { return false; }
	virtual void generateMesh(Mesh* mesh) {}
	
protected:
    virtual void generateLevelset(Grid<Real>& phi) { gridSetConst<Real>(phi, 1000.0f); }
};

//! Box shape
class Box : public Shape
{   
public:
	Box(FluidSolver* parent, Vec3 center = Vec3::Invalid, Vec3 p0 = Vec3::Invalid, Vec3 p1 = Vec3::Invalid, Vec3 size = Vec3::Invalid);
	
	inline Vec3 getSize() const { return mP1-mP0; }
	inline Vec3 getP0() const { return mP0; }
	inline Vec3 getP1() const { return mP1; }
	virtual void setCenter(const Vec3& center) { Vec3 dh=0.5*(mP1-mP0); mP0 = center-dh; mP1 = center+dh;}
	virtual Vec3 getCenter() const { return 0.5*(mP1+mP0); }
	virtual Vec3 getExtent() const { return getSize(); }
	virtual bool isInside(const Vec3& pos) const;
	virtual void generateMesh(Mesh* mesh);
	virtual void generateLevelset(Grid<Real>& phi);
	
protected:
    Vec3 mP0, mP1;
};

//! Spherical shape
class Sphere : public Shape
{   
public:
	Sphere(FluidSolver* parent, Vec3 center, Real radius, Vec3 scale=Vec3(1,1,1));
	
	virtual void setCenter(const Vec3& center) { mCenter = center; }
	virtual Vec3 getCenter() const { return mCenter; }
	inline Real getRadius() const { return mRadius; }
	virtual Vec3 getExtent() const { return Vec3(2.0*mRadius); }    
	virtual bool isInside(const Vec3& pos) const;
	virtual void generateMesh(Mesh* mesh);
	virtual void generateLevelset(Grid<Real>& phi);
	
protected:
	Vec3 mCenter, mScale;
    Real mRadius;
};

//! Cylindrical shape
class Cylinder : public Shape
{   
public:
	Cylinder(FluidSolver* parent, Vec3 center, Real radius, Vec3 z);
	
	void setRadius(Real r) { mRadius = r; }
	void setZ(Vec3 z) { mZDir=z; mZ=normalize(mZDir); }
	
	virtual void setCenter(const Vec3& center) { mCenter=center; }
	virtual Vec3 getCenter() const { return mCenter; }
	inline Real getRadius() const { return mRadius; }
	inline Vec3 getZ() const { return mZ*mZDir; }
	virtual Vec3 getExtent() const { return Vec3(2.0*sqrt(square(mZ)+square(mRadius))); }    
	virtual bool isInside(const Vec3& pos) const;
	virtual void generateMesh(Mesh* mesh);
	virtual void generateLevelset(Grid<Real>& phi);

protected:
	Vec3 mCenter, mZDir;
    Real mRadius, mZ;
};

//! Slope shape
// generates a levelset based on a plane
// plane is specified by two angles and an offset on the y axis in (offset vector would be ( 0, offset, 0) )
// the two angles are specified in degrees, between: y-axis and x-axis 
//                                                   y-axis and z-axis
class Slope : public Shape
{
public:
	Slope(FluidSolver* parent, Real anglexy, Real angleyz, Real origin, Vec3 gs);

	virtual void setOrigin (const Real& origin)  { mOrigin=origin; }
	virtual void setAnglexy(const Real& anglexy) { mAnglexy=anglexy; }
	virtual void setAngleyz(const Real& angleyz) { mAnglexy=angleyz; }

	inline Real getOrigin()   const { return mOrigin; }
	inline Real getmAnglexy() const { return mAnglexy; }
	inline Real getmAngleyz() const { return mAngleyz; }
	virtual bool isInside(const Vec3& pos) const;
	virtual void generateMesh(Mesh* mesh);
	virtual void generateLevelset(Grid<Real>& phi);

protected:
	Real mAnglexy, mAngleyz;
	Real mOrigin;
    Vec3 mGs;
};



//! Kernel: Apply a shape to a grid, setting value inside
template <class T>
struct ApplyShapeToGrid : public KernelBase
{
    ApplyShapeToGrid(Grid<T>* grid, Shape* shape, T value, FlagGrid* respectFlags) 
        : KernelBase(grid,0)
        , grid(grid)
        , shape(shape)
        , value(value)
        , respectFlags(respectFlags)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<T>* grid, Shape* shape, T value, FlagGrid* respectFlags)
    {
	    if (respectFlags && respectFlags->isObstacle(i,j,k)) return;
	    if (shape->isInsideGrid(i,j,k))                      (*grid)(i,j,k) = value;
    }
    
    inline Grid<T>* getArg0() { return grid; }
    typedef Grid<T> type0;
    inline Shape* getArg1() { return shape; }
    typedef Shape type1;
    inline T& getArg2() { return value; }
    typedef T type2;
    inline FlagGrid* getArg3() { return respectFlags; }
    typedef FlagGrid type3;
    
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
                            op(i,j,k,grid,shape,value,respectFlags);
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
                        op(i,j,k,grid,shape,value,respectFlags);
            }
        }
    }
    
    Grid<T>* grid;
    Shape* shape;
    T value;
    FlagGrid* respectFlags;
};

//! Kernel: Apply a shape to a MAC grid, setting value inside
struct ApplyShapeToMACGrid : public KernelBase
{
    ApplyShapeToMACGrid(MACGrid* grid, Shape* shape, Vec3 value, FlagGrid* respectFlags)
        : KernelBase(grid,0)
        , grid(grid)
        , shape(shape)
        , value(value)
        , respectFlags(respectFlags)
    {
        run();
    }
    
    inline void op(int i, int j, int k, MACGrid* grid, Shape* shape, Vec3 value, FlagGrid* respectFlags )
    {
	    if (respectFlags && respectFlags->isObstacle(i,j,k))
		    return;    
	    if (shape->isInside(Vec3(i,j+0.5,k+0.5))) (*grid)(i,j,k).x = value.x;
	    if (shape->isInside(Vec3(i+0.5,j,k+0.5))) (*grid)(i,j,k).y = value.y;
	    if (shape->isInside(Vec3(i+0.5,j+0.5,k))) (*grid)(i,j,k).z = value.z;
    }
    
    inline MACGrid* getArg0() { return grid; }
    typedef MACGrid type0;
    inline Shape* getArg1() { return shape; }
    typedef Shape type1;
    inline Vec3& getArg2() { return value; }
    typedef Vec3 type2;
    inline FlagGrid* getArg3() { return respectFlags; }
    typedef FlagGrid type3;
    
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
                    for (int j=0; j < _maxY; j++)
                        for (int i=0; i < _maxX; i++)
                            op(i,j,k,grid,shape,value,respectFlags);
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
                        op(i,j,k,grid,shape,value,respectFlags);
            }
        }
    }
    
    MACGrid* grid;
    Shape* shape;
    Vec3 value;
    FlagGrid* respectFlags;
};

//! Kernel: Apply a shape to a grid, setting value inside (scaling by SDF value)
template <class T>
struct ApplyShapeToGridSmooth : public KernelBase
{
    ApplyShapeToGridSmooth(Grid<T>* grid, Grid<Real>& phi, Real sigma, Real shift, T value, FlagGrid* respectFlags)
        : KernelBase(grid,0)
        , grid(grid)
        , phi(phi)
        , sigma(sigma)
        , shift(shift)
        , value(value)
        , respectFlags(respectFlags)
    {
        run();
    }
    
    inline void op(int i, int j, int k, Grid<T>* grid, Grid<Real>& phi, Real sigma, Real shift, T value, FlagGrid* respectFlags)
    {
	    if (respectFlags && respectFlags->isObstacle(i,j,k)) return;

	    const Real p = phi(i,j,k) - shift;
	    if (p < -sigma)
		    (*grid)(i,j,k) = value;
	    else if (p < sigma)
		    (*grid)(i,j,k) = value*(0.5f*(1.0f-p/sigma));
    }
    
    inline Grid<T>* getArg0() { return grid; }
    typedef Grid<T> type0;
    inline Grid<Real>& getArg1() { return phi; }
    typedef Grid<Real> type1;
    inline Real& getArg2() { return sigma; }
    typedef Real type2;
    inline Real& getArg3() { return shift; }
    typedef Real type3;
    inline T& getArg4() { return value; }
    typedef T type4;
    inline FlagGrid* getArg5() { return respectFlags; }
    
    typedef FlagGrid type5;
    
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
                            op(i,j,k,grid,phi,sigma,shift,value,respectFlags);
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
                        op(i,j,k,grid,phi,sigma,shift,value,respectFlags);
            }
        }
    }
    
    Grid<T>* grid;
    Grid<Real>& phi;
    Real sigma;
    Real shift;
    T value;
    FlagGrid* respectFlags;
};

} //namespace
