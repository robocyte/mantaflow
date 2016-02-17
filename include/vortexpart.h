/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Vortex particles
 *
 ******************************************************************************/

#pragma once

#include "particle.h"

namespace Manta
{
class Mesh;
	
struct VortexParticleData
{
	VortexParticleData()
        : pos(_0)
        , vorticity(_0)
        , sigma(0)
        , flag(0)
    {}

	VortexParticleData(const Vec3& p, const Vec3& v, Real sig)
        : pos(p)
        , vorticity(v)
        , sigma(sig)
        , flag(0)
    {}

	Vec3 pos, vorticity;
	Real sigma;
	int flag;    

	static ParticleBase::SystemType getType()
    {
        return ParticleBase::VORTEX;
    }
};

//! Vortex particles
class VortexParticleSystem : public ParticleSystem<VortexParticleData>
{
public:
	VortexParticleSystem(FluidSolver* parent);
  
	void advectSelf(Real scale=1.0, int integrationMode=IntRK4);
	void applyToMesh(Mesh& mesh, Real scale=1.0, int integrationMode=IntRK4);
	virtual ParticleBase* clone();
};

} // namespace
