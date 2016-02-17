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

#pragma once

#include "particle.h"
#include "grid.h"
#include "randomstream.h"
#include "levelset.h"

namespace Manta
{

void sampleFlagsWithParticles(FlagGrid& flags, BasicParticleSystem& parts, int discretization, Real randomness);

void sampleLevelsetWithParticles(LevelsetGrid& phi, FlagGrid& flags, BasicParticleSystem& parts, int discretization, Real randomness, bool reset = false, bool refillEmpty = false);

void mapParticlesToMAC(FlagGrid& flags, MACGrid& vel, MACGrid& velOld, BasicParticleSystem& parts, ParticleDataImpl<Vec3>& partVel, Grid<Vec3>* weight = NULL);

void markFluidCells(BasicParticleSystem& parts, FlagGrid& flags, Grid<Real>* phiObs = NULL);

void adjustNumber(BasicParticleSystem& parts, MACGrid& vel, FlagGrid& flags, int minParticles, int maxParticles, LevelsetGrid& phi, Real radiusFactor = 1.0, Real narrowBand = -1.0, Grid<Real>* exclude = nullptr);

void flipVelocityUpdate(FlagGrid& flags, MACGrid& vel, MACGrid& velOld, BasicParticleSystem& parts, ParticleDataImpl<Vec3>& partVel, Real flipRatio);

void mapGridToPartsVec3(Grid<Vec3>& source , BasicParticleSystem& parts , ParticleDataImpl<Vec3>& target);

void pushOutofObs(BasicParticleSystem& parts, FlagGrid& flags, Grid<Real>& phiObs, Real shift = 0.05, Real thresh = 0.0);

// for testing purposes only...
void testInitGridWithPos(Grid<Real>& grid);

// build a grid that contains indices for a particle system
// the particles in a cell i,j,k are particles[index(i,j,k)] to particles[index(i+1,j,k)-1]
// (ie,  particles[index(i+1,j,k)] alreadu belongs to cell i+1,j,k)
void gridParticleIndex( BasicParticleSystem& parts, ParticleIndexSystem& indexSys, FlagGrid& flags, Grid<int>& index, Grid<int>* counter = nullptr);

void unionParticleLevelset(BasicParticleSystem& parts, ParticleIndexSystem& indexSys, FlagGrid& flags, Grid<int>& index, LevelsetGrid& phi, Real radiusFactor = 1.0);

void averagedParticleLevelset(BasicParticleSystem& parts, ParticleIndexSystem& indexSys, FlagGrid& flags, Grid<int>& index, LevelsetGrid& phi, Real radiusFactor=1.0, int smoothen = 1 , int smoothenNeg = 1);

} // namespace
