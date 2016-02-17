#include <sstream>
#include <iomanip>

#include "manta.h"
#include "grid.h"
#include "levelset.h"
#include "particle.h"
#include "mesh.h"
#include "shapes.h"
#include "flip.h"
#include "fastmarch.h"
#include "extforces.h"
#include "pressure.h"
#include "initplugins.h"

int main()
{
    int  dimension  = 3;
    auto resolution = Manta::Vec3i(64);

    if (dimension == 2) resolution.z = 1;

    auto main_solver = Manta::FluidSolver(resolution, dimension);
    main_solver.setTimeStep(0.8f);
    Real minParticles = pow(2, dimension);

    // size of particles
    Real radiusFactor = 1.0f;

    // solver grids
    auto flags      = Manta::FlagGrid(&main_solver);
    auto phi        = Manta::LevelsetGrid(&main_solver);
    auto phiObs     = Manta::LevelsetGrid(&main_solver);

    auto vel        = Manta::MACGrid(&main_solver);
    auto velOld     = Manta::MACGrid(&main_solver);
    auto pressure   = Manta::Grid<Real>(&main_solver);
    auto fractions  = Manta::MACGrid(&main_solver);
    auto tmpVec3    = Manta::Grid<Manta::Vec3>(&main_solver);
    auto phiWalls   = Manta::LevelsetGrid(&main_solver);

    // particles
    auto pp         = Manta::BasicParticleSystem(&main_solver);
    auto pVel       = Manta::PdataVec3(&pp);
    auto mesh       = Manta::Mesh(&main_solver);

    // acceleration data for particle nbs
    auto pindex     = Manta::ParticleIndexSystem(&main_solver);
    auto gpi        = Manta::Grid<int>(&main_solver);

    // scene setup
    int bWidth = 1;
    flags.initDomain(bWidth, "xXyYzZ", "      ", "      ", "      ", &phiWalls);
    phi.setConst(999.0f);
    phiObs.setConst(999.0f);

    // standing dam
    auto fluidbox   = Manta::Box(&main_solver, Manta::Vec3::Invalid,
                                               Manta::Vec3(0.0f, 0.0f, 0.0f),
                                               Manta::Vec3(resolution.x * 1.0f, resolution.y * 0.3f, resolution.z * 1.0f));
    phi.join(fluidbox.computeLevelset());
    auto fluidbox2  = Manta::Box(&main_solver, Manta::Vec3::Invalid,
                                               Manta::Vec3(resolution.x * 0.1f, 0.0f, 0.0f),
                                               Manta::Vec3(resolution.x * 0.2f, resolution.y * 0.75f, resolution.z * 1.0f));
    phi.join(fluidbox2.computeLevelset());

    phiObs.join(phiWalls);

    auto sphere     = Manta::Sphere(&main_solver, Manta::Vec3(resolution.x * 0.66f, resolution.y * 0.3f, resolution.z * 0.5f), resolution.x * 0.2f);
    phiObs.join(sphere.computeLevelset());

    flags.updateFromLevelset(phi);
    phi.subtract(phiObs);
    Manta::sampleLevelsetWithParticles(phi, flags, pp, 2, 0.05f);

    // also set boundary flags for phiObs
    Manta::updateFractions(flags, phiObs, fractions, bWidth);
    Manta::setObstacleFlags(flags, fractions, phiObs);

    for (int t = 0; t < 250; t++)
    {
        // FLIP
        pp.advectInGrid(flags, vel, Manta::IntegrationMode::IntRK4, false, false);
        Manta::pushOutofObs(pp, flags, phiObs);

        // make sure we have velocities throught liquid region
        Manta::mapParticlesToMAC(flags, vel, velOld, pp, pVel, &tmpVec3);
        Manta::extrapolateMACFromWeight(vel, tmpVec3, 2);
        Manta::markFluidCells(pp, flags, &phiObs);

        // create approximate surface level set, resample particles
        Manta::gridParticleIndex(pp, pindex, flags, gpi);
        Manta::unionParticleLevelset(pp, pindex, flags, gpi, phi, radiusFactor);

	    // extend levelset somewhat, needed by particle resampling in adjustNumber
	    Manta::extrapolateLsSimple(phi, 4, true);

        // forces and pressure solve
        Manta::addGravity(flags, vel, Manta::Vec3(0.0f, -0.001f, 0.0f));
        Manta::extrapolateMACSimple(flags, vel, 2, nullptr, true);
        Manta::setWallBcs(flags, vel, &fractions, &phiObs);
        Manta::solvePressure(vel, pressure, flags, &phi, nullptr, &fractions);
        Manta::extrapolateMACSimple(flags, vel, 4, nullptr, true);
        Manta::setWallBcs(flags, vel, &fractions, &phiObs);

        // set source grids for resampling, used in adjustNumber!
        pVel.setSource(&vel, true);
        Manta::adjustNumber(pp, vel, flags, 1 * minParticles, 2 * minParticles, phi, radiusFactor, -1.0f, &phiObs);

        // make sure we have proper velocities
        Manta::flipVelocityUpdate(flags, vel, velOld, pp, pVel, 0.97f);

        if (dimension == 3) phi.createMesh(mesh);

        main_solver.step();

    	// generate data directly and for flip03_gen.py surface generation scene
        std::stringstream filename1, filename2;
        filename1 << "simulation data\\flip05_secOrderBnd\\fluidsurface_final_" << std::setw(4) << std::setfill('0') << t << ".bobj.gz";
        filename2 << "simulation data\\flip05_secOrderBnd\\flipParts_"          << std::setw(4) << std::setfill('0') << t << ".uni";

        mesh.save(filename1.str());
        //pp.save(filename2.str());
    }

}
