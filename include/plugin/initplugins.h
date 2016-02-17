/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Tools to setup fields and inflows
 *
 ******************************************************************************/

#pragma once

#include "vectorbase.h"
#include "grid.h"
#include "commonkernels.h"
#include "noisefield.h"
#include "shapes.h"

namespace Manta
{
    
//! Init noise-modulated density inside shape
void densityInflow(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Shape* shape, Real scale = 1.0, Real sigma = 0);

void densityInflowMeshNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Mesh* mesh, Real scale = 1.0, Real sigma = 0);

//! Init constant density inside mesh
void densityInflowMesh(FlagGrid& flags, Grid<Real>& density, Mesh* mesh, Real value = 1., Real cutoff = 7, Real sigma = 0);

void updateFractions(FlagGrid& flags, Grid<Real>& phiObs, MACGrid& fractions, const int &boundaryWidth = 0);

void setObstacleFlags(FlagGrid& flags, MACGrid& fractions, Grid<Real>& phiObs);

} //namesopace
