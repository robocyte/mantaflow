/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugins for pressure correction: solve_pressure, and ghost fluid helpers
 *
 ******************************************************************************/

#pragma once

#include "conjugategrad.h"

namespace Manta
{

void solvePressure(MACGrid& vel, Grid<Real>& pressure, FlagGrid& flags, Grid<Real>* phi = 0, Grid<Real>* perCellCorr = 0, MACGrid* fractions = 0, Real gfClamp = 1e-04, Real cgMaxIterFac = 1.5, Real cgAccuracy = 1e-3, bool precondition = true, bool enforceCompatibility = false, bool useL2Norm = false, bool zeroPressureFixing = false, Grid<Real>* retRhs = NULL);

} // namespace
