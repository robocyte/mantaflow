#include "manta.h"
#include "grid.h"
#include "particle.h"
#include "shapes.h"
#include "flip.h"
#include "fastmarch.h"
#include "extforces.h"
#include "pressure.h"

int main()
{
    int  dimension       = 3;
    int  particle_number = 2;
    auto resolution      = Manta::Vec3i(64);

    if (dimension == 2)
    {
        resolution.z = 1;
        particle_number = 3;
    }

    auto main_solver = Manta::FluidSolver(resolution, dimension);
    main_solver.setTimeStep(0.5f);

    // solver grids
    auto flags      = Manta::FlagGrid(&main_solver);
    auto vel        = Manta::MACGrid(&main_solver);
    auto vel_old    = Manta::MACGrid(&main_solver);
    auto pressure   = Manta::Grid<Real>(&main_solver);
    auto tmp_vec3   = Manta::Grid<Manta::Vec3>(&main_solver);

    // particles
    auto pp         = Manta::BasicParticleSystem(&main_solver);
    auto pVel       = Manta::PdataVec3(&pp);

    flags.initDomain();

    auto fluid_box = Manta::Box(&main_solver, Manta::Vec3::Invalid,
                                              Manta::Vec3(resolution.x * 0.4f, resolution.y * 0.72f, resolution.z * 0.4f),
                                              Manta::Vec3(resolution.x * 0.6f, resolution.y * 0.92f, resolution.z * 0.6f));

    auto phi_init = fluid_box.computeLevelset();

    flags.updateFromLevelset(phi_init);

    Manta::sampleFlagsWithParticles(flags, pp, particle_number, 0.2f);

    for (int i = 0; i < 100; i++)
    {
        // FLIP
        pp.advectInGrid(flags, vel, Manta::IntegrationMode::IntRK4, false);
        Manta::mapParticlesToMAC(flags, vel, vel_old, pp, pVel, &tmp_vec3);
        Manta::extrapolateMACFromWeight(vel, tmp_vec3, 2);
        Manta::markFluidCells(pp, flags);

        Manta::addGravity(flags, vel, Manta::Vec3(0.0f, -0.002f, 0.0f));

        // pressure solve
        Manta::setWallBcs(flags, vel);
        Manta::solvePressure(vel, pressure, flags);
        Manta::setWallBcs(flags, vel);

	    // we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
	    Manta::extrapolateMACSimple(flags, vel);
	
	    // FLIP velocity update
	    Manta::flipVelocityUpdate(flags, vel, vel_old, pp, pVel, 0.97f);
	
	    main_solver.step();
    }
}
