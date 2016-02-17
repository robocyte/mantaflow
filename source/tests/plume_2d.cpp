#include "manta.h"
#include "timing.h"
#include "grid.h"
#include "extforces.h"
#include "shapes.h"
#include "advection.h"
#include "pressure.h"
#include "fileio.h"

int main()
{
    int res = 64;
    auto gs = Manta::Vec3i(res, res, 1);
    auto main_solver = Manta::FluidSolver(gs, 2);
    main_solver.setTimeStep(0.5f);
    auto timings = Manta::Timings();
    timings.setParent(&main_solver);

    // prepare grids
    auto flags      = Manta::FlagGrid(&main_solver);
    auto vel        = Manta::MACGrid(&main_solver);
    auto density    = Manta::Grid<Real>(&main_solver);
    auto pressure   = Manta::Grid<Real>(&main_solver);
    
    flags.initDomain();
    flags.fillGrid();

    Manta::setOpenBound(flags, 1, "yY", Manta::FlagGrid::TypeOutflow | Manta::FlagGrid::TypeEmpty);

    auto source = Manta::Cylinder(&main_solver, Manta::Vec3(64.0f * 0.5f,
                                                            64.0f * 0.1f,
                                                            1.0f * 0.5f),
                                                            res * 0.14f, Manta::Vec3(0.0f,
                                                                                     64.0f * 0.02f,
                                                                                     0.0f));

    for (int i = 0; i < 400; i++)
    {

        if (i < 300) source.applyToGrid<Real>(&density, 1.0f);

        Manta::advectSemiLagrange(&flags, &vel, &density, 2);
        Manta::advectSemiLagrange(&flags, &vel, &vel, 2, 1.0f, 1, true, 1);
        Manta::resetOutflow(flags, nullptr, nullptr, &density);

        Manta::setWallBcs(flags, vel);
        Manta::addBuoyancy(flags, density, vel, Manta::Vec3(0.0f, -0.003f, 0.0f));

        Manta::solvePressure(vel, pressure, flags);

        //timings.display();
        main_solver.step();

        //if (i == 150) Manta::writeGridTxt("test150.txt", &density);
    }
}
