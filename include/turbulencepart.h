/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Turbulence particles
 *
 ******************************************************************************/

#pragma once

#include "particle.h"
#include "noisefield.h"

namespace Manta
{
class Shape;
	
struct TurbulenceParticleData
{
	TurbulenceParticleData()
        : pos(_0)
        , color(1.)
        , tex0(_0)
        , tex1(_0)
        , flag(0)
    {}

	TurbulenceParticleData(const Vec3& p, const Vec3& color = Vec3(1.))
        : pos(p)
        , color(color)
        , tex0(p)
        , tex1(p)
        , flag(0)
    {}

	Vec3 pos, color;
	Vec3 tex0, tex1;
	int flag;

	static ParticleBase::SystemType getType()
    {
        return ParticleBase::TURBULENCE;
    }
};

//! Turbulence particles
class TurbulenceParticleSystem : public ParticleSystem<TurbulenceParticleData>
{
public:
	TurbulenceParticleSystem(FluidSolver* parent, WaveletNoiseField& noise);
  
	void resetTexCoords(int num, const Vec3& inflow);
	void seed(Shape* source, int num);
	void synthesize(FlagGrid& flags, Grid<Real>& k, int octaves=2, Real switchLength=10.0, Real L0=0.1, Real scale=1.0, Vec3 inflowBias=_0);
	void deleteInObstacle(FlagGrid& flags);
		
	virtual ParticleBase* clone();
	
private:
    WaveletNoiseField& noise;
};

} // namespace
