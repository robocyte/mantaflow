/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * moving obstacles
 *
 ******************************************************************************/

#pragma once

#include "shapes.h"
#include "particle.h"

namespace Manta
{

//! Moving obstacle composed of basic shapes
class MovingObstacle : public PbClass
{
public:
	MovingObstacle(FluidSolver* parent, int emptyType=FlagGrid::TypeEmpty);
	
	void add(Shape* shape);

    //! If t in [t0,t1], apply linear motion path from p0 to p1
	void moveLinear(Real t, Real t0, Real t1, Vec3 p0, Vec3 p1, FlagGrid& flags, MACGrid& vel, bool smooth=true);
	
    //! Compute levelset, and project FLIP particles outside obstacles
	void projectOutside(FlagGrid& flags, BasicParticleSystem& flip);
	
protected:
	std::vector<Shape*> mShapes;
	int mEmptyType;
	int mID;
    static int sIDcnt;
};

} //namespace
