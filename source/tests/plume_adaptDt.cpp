#include "manta.h"
#include "grid.h"
#include "noisefield.h"
#include "extforces.h"
#include "shapes.h"
#include "initplugins.h"
#include "advection.h"
#include "pressure.h"
//#include "fileio.h"

int main()
{
    int dim = 3;
    int res = 64;
    auto gs = Manta::Vec3i(res, 1.5f * res, res);
    if (dim == 2) gs.z = 1;
    auto main_solver = Manta::FluidSolver(gs, dim);

    int frames = 100;

    // set time step range
    main_solver.mFrameLength    = 1.2f;
    main_solver.mDtMin          = 0.2;
    main_solver.mDtMax          = 2.0f;
    main_solver.mCflCond        = 3.0f;
    main_solver.setTimeStep((main_solver.mDtMax + main_solver.mDtMin) * 0.5f);

    // prepare grids
    auto flags      = Manta::FlagGrid(&main_solver);
    auto vel        = Manta::MACGrid(&main_solver);
    auto density    = Manta::Grid<Real>(&main_solver);
    auto pressure   = Manta::Grid<Real>(&main_solver);
    
    // noise field
    auto noise = Manta::WaveletNoiseField(&main_solver, -1, 1);
    noise.mPosScale = Manta::Vec3(45);
    noise.mClamp = true;
    noise.mClampNeg = 0;
    noise.mClampPos = 1;
    noise.mValScale = 1;
    noise.mValOffset = 0.75f;
    noise.mTimeAnim = 0.2f;

    flags.initDomain();
    flags.fillGrid();

    auto source = Manta::Cylinder(&main_solver, Manta::Vec3(64.0f * 0.5f,
                                                            92.0f * 0.1f,
                                                            64.0f * 0.5f),
                                                            res * 0.14f, Manta::Vec3(0.0f,
                                                                                     92.0f * 0.02f,
                                                                                     0.0f));

    while (main_solver.getFrame() < frames)
    {
        auto maxvel = vel.getMaxValue();
        main_solver.adaptTimestep(maxvel);

        if (main_solver.getTimeTotal() < 50.0f)
        {
            Manta::densityInflow(flags, density, noise, &source, 1.0f, 0.5f);
        }

        Manta::advectSemiLagrange(&flags, &vel, &density, 2);
        Manta::advectSemiLagrange(&flags, &vel, &vel, 2);

        Manta::setWallBcs(flags, vel);
        Manta::addBuoyancy(flags, density, vel, Manta::Vec3(0.0f, -0.006f, 0.0f));

        Manta::solvePressure(vel, pressure, flags);
        Manta::setWallBcs(flags, vel);

        //timings.display();
        main_solver.step();

        //if (i == 150) Manta::writeGridTxt("test150.txt", &density);
    }
}
