/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 ******************************************************************************/

#include "kernel.h"
#include "grid.h"

using namespace std;
namespace Manta {

	KERNEL 
	void levelsetToDensityKN(Grid<Real>& phi, Grid<Real>& density){
		// map phi values [-5,5] to density values [0,1]
	}

	PYTHON void levelsetToDensity(Grid<Real>& phi, Grid<Real>& density){
		levelsetToDensityKN(phi, density);
	}
}; 
