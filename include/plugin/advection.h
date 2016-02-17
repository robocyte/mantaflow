/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2015 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugins for advection
 *
 ******************************************************************************/

#pragma once

#include "vectorbase.h"
#include "grid.h"
#include "kernel.h"
#include <limits>

namespace Manta
{

void advectSemiLagrange(FlagGrid* flags, MACGrid* vel, GridBase* grid, int order = 1, Real strength = 1.0, int orderSpace = 1, bool openBounds = false, int boundaryWidth = 1);

} // namespace
