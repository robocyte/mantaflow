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
#include <thread>
#include "kernel.h"
#include "grid.h"
#include "pcgsolver.h"

using namespace std;
namespace Manta {
	
// 3.2 
PYTHON void diffuseTemperatureExplicit(...){
	// don't overwrite values in T that will be read again
	// write a KERNEL and make sure that the temperature in boundary cells stays zero
}

// 3.3
KERNEL void setupB(...){ 
	
}

KERNEL void fillT(...){
	// make sure that temperature values in boundary cells are zero!
}

// use the single argument to prevent multithreading (multiple threads might interfere during the matrix setup)
KERNEL (single) void setupA(...){  
	// set with:  A.set_element( index1, index2 , value );
	// if needed, read with: A(index1, index2);

	// avoid zero rows in A -> set the diagonal value for boundary cells to 1.0
}

PYTHON void diffuseTemperatureImplicit(...){
		// solve A T = b
		const int N = ...;
		SparseMatrix<Real> A(N);
		std::vector<Real> b(N);

		setupA(...);
		setupB(...);
		
		// perform solve
		Real pcg_target_residual = 1e-05;
		Real pcg_max_iterations = 1000;
		Real ret_pcg_residual = 1e10;
		int  ret_pcg_iterations = -1;

		SparsePCGSolver<Real> solver;
		solver.set_solver_parameters(pcg_target_residual, pcg_max_iterations, 0.97, 0.25);

		std::vector<Real> x(N);
		for (int j = 0; j<N; ++j) { x[j] = 0.; }

		// preconditioners: 0 off, 1 diagonal, 2 incomplete cholesky
		solver.solve(A, b, x, ret_pcg_residual, ret_pcg_iterations, 0);

		// x contains the new temperature values
		fillT(...);
}
}; 
