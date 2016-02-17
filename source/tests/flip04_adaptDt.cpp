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

int main()
{
    int  dimension  = 3;
    auto resolution = Manta::Vec3i(80);
    if (dimension == 2) resolution.z = 1;

    auto main_solver = Manta::FluidSolver(resolution, dimension);

    // how many frames to calculate
    int frames = 250;

    // adaptive time stepping
    main_solver.mFrameLength    = 0.6f;     // length of one frame (in "world time")
    main_solver.mDtMin          = 0.1f;     // time step range
    main_solver.mDtMax          = 2.0f;
    main_solver.mCflCond        = 1.5f;     // maximal velocity per cell
    main_solver.setTimeStep(0.5f * (main_solver.mDtMax + main_solver.mDtMin));

    Real minParticles = pow(2, dimension);

    // size of particles
    Real radiusFactor = 1.0f;

    // prepare grids
    auto flags      = Manta::FlagGrid(&main_solver);
    auto phi        = Manta::LevelsetGrid(&main_solver);
    auto vel        = Manta::MACGrid(&main_solver);
    auto velOld     = Manta::MACGrid(&main_solver);
    auto pressure   = Manta::Grid<Real>(&main_solver);
    auto tmpVec3    = Manta::Grid<Manta::Vec3>(&main_solver);
    auto tstGrid    = Manta::Grid<Real>(&main_solver);
    auto phiObs     = Manta::LevelsetGrid(&main_solver);

    // prepare particles
    auto pp         = Manta::BasicParticleSystem(&main_solver);
    auto pVel       = Manta::PdataVec3(&pp);
    auto pTest      = Manta::PdataReal(&pp);
    auto mesh       = Manta::Mesh(&main_solver);

    // acceleration data for particle nbs
    auto pindex     = Manta::ParticleIndexSystem(&main_solver);
    auto gpi        = Manta::Grid<int>(&main_solver);

    int setup = 0;
    flags.initDomain();
    auto fluidVel = Manta::Sphere(&main_solver, Manta::Vec3(0.0f), 1.0f);
    Manta::Vec3 fluidSetVel;

    if (setup == 0)
    {
        // breakin dam
        auto fluidbox   = Manta::Box(&main_solver, Manta::Vec3::Invalid,
                                                   Manta::Vec3(0.0f, 0.0f, 0.0f),
                                                   Manta::Vec3(resolution.x * 0.4f, resolution.y * 0.6f, resolution.z * 1.0f));
        phi.copyFrom(fluidbox.computeLevelset());

    } else if (setup == 1)
    {
        // drop into box
        auto fluidbasin = Manta::Box(&main_solver, Manta::Vec3::Invalid,
                                                   Manta::Vec3(0.0f, 0.0f, 0.0f),
                                                   Manta::Vec3(resolution.x * 1.0f, resolution.y * 0.2f, resolution.z * 1.0f));
        auto dropCenter = Manta::Vec3(0.5f, 0.5f, 0.5f);
        Real dropRadius = 0.15f;
        fluidSetVel     = Manta::Vec3(0.0f, -0.03f, 0.0f);
        auto fluidDrop  = Manta::Sphere(&main_solver, Manta::Vec3(resolution.x * dropCenter.x, resolution.y * dropCenter.y, resolution.z * dropCenter.z),
                                                      resolution.x * dropRadius);
        fluidVel        = Manta::Sphere(&main_solver, Manta::Vec3(resolution.x * dropCenter.x, resolution.y * dropCenter.y, resolution.z * dropCenter.z),
                                                      resolution.x * (dropRadius + 0.05f));
        phi.copyFrom(fluidbasin.computeLevelset());
        phi.join(fluidDrop.computeLevelset());
    }

    flags.updateFromLevelset(phi);

    if (dimension == 3)
    {
        auto obsBox     = Manta::Box(&main_solver, Manta::Vec3::Invalid,
                                                   Manta::Vec3(resolution.x * 0.7f, resolution.y * 0.0f, resolution.z * 0.5f),
                                                   Manta::Vec3(resolution.x * 0.8f, resolution.y * 1.0f, resolution.z * 0.8f));
        obsBox.applyToGrid(&flags, (int)Manta::FlagGrid::CellType::TypeObstacle);
    }

    Manta::sampleLevelsetWithParticles(phi, flags, pp, 2, 0.05f);
    Manta::mapGridToPartsVec3(vel, pp, pVel);

    if (setup == 1)
    {
        fluidVel.applyToGrid(&vel, fluidSetVel);
        Manta::mapGridToPartsVec3(vel, pp, pVel);
    }

    Manta::testInitGridWithPos(tstGrid);
    pTest.setConst(0.1f);

    int t = 0;
    while (main_solver.getFrame() < frames)
    {
        auto maxvel = vel.getMaxValue();
        main_solver.adaptTimestep(maxvel);

        // FLIP
        pp.advectInGrid(flags, vel, Manta::IntegrationMode::IntRK4, false);

	    // make sure we have velocities throught liquid region
	    Manta::mapParticlesToMAC(flags, vel, velOld, pp, pVel, &tmpVec3);
	    Manta::extrapolateMACFromWeight(vel, tmpVec3, 2); // note, tmpVec3 could be free'd now...
        Manta::markFluidCells(pp, flags);

	    // create approximate surface level set, resample particles
        Manta::gridParticleIndex(pp, pindex, flags, gpi);
        Manta::unionParticleLevelset(pp, pindex, flags, gpi, phi, radiusFactor);
	    Manta::extrapolateLsSimple(phi, 4, true);
	    // note - outside levelset doesnt matter...

        // forces and pressure solve
        Manta::addGravity(flags, vel, Manta::Vec3(0.0f, -0.003f, 0.0f));
        Manta::setWallBcs(flags, vel);
        Manta::solvePressure(vel, pressure, flags, &phi);
        Manta::setWallBcs(flags, vel);

        // set source grids for resampling, used in adjustNumber!
        pVel.setSource(&vel, true);
        pTest.setSource(&tstGrid);
        Manta::adjustNumber(pp, vel, flags, 1 * minParticles, 2 * minParticles, phi, radiusFactor);

        // make sure we have proper velocities
        Manta::extrapolateMACSimple(flags, vel, (int)(maxvel * 1.5f) + 2);
        Manta::flipVelocityUpdate(flags, vel, velOld, pp, pVel, 0.97f);

        if (dimension == 3) phi.createMesh(mesh);

        main_solver.step();

        // generate data directly and for flip03_gen.exe surface generation scene
        // note: 1 file corresponds to 1 simulation step, not "world time"
        std::stringstream filename1, filename2;
        filename1 << "simulation data\\flip04_adaptDt\\fluidsurface_final_" << std::setw(4) << std::setfill('0') << t << ".bobj.gz";
        filename2 << "simulation data\\flip04_adaptDt\\flipParts_"          << std::setw(4) << std::setfill('0') << t << ".uni";

        mesh.save(filename1.str());
        pp.save(filename2.str());

        t++;
    }
}
