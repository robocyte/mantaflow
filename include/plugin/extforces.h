/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Set boundary conditions, gravity
 *
 ******************************************************************************/

#pragma once

#include "vectorbase.h"
#include "grid.h"
#include "commonkernels.h"
#include "particle.h"

#include <string>

namespace Manta
{ 

void setWallBcs(FlagGrid& flags, MACGrid& vel, MACGrid* fractions = 0, Grid<Real>* phiObs = 0, int boundaryWidth = 0);

void addGravity(FlagGrid& flags, MACGrid& vel, Vec3 gravity);

void addBuoyancy(FlagGrid& flags, Grid<Real>& density, MACGrid& vel, Vec3 gravity);

void setOpenBound(FlagGrid& flags, int bWidth, std::string openBound = "", int type = FlagGrid::TypeOutflow | FlagGrid::TypeEmpty);

void resetOutflow(FlagGrid& flags, Grid<Real>* phi = 0, BasicParticleSystem* parts = 0, Grid<Real>* real = 0, Grid<int>* index = 0, ParticleIndexSystem* indexSys = 0);

} // namespace
