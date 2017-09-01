/*
 * reachabilitySequential_GPU_MatrixVector_Multiply.cuh
 *
 *  Created on: 16-Apr-2015
 *      Author: amit
 */

#ifndef REACHABILITYSEQUENTIAL_GPU_MATRIXVECTOR_MULTIPLY_CUH_
#define REACHABILITYSEQUENTIAL_GPU_MATRIXVECTOR_MULTIPLY_CUH_

#include "Utilities/invariantBoundaryCheck.h"
#include "core_system/math/matrix.h"
#include "Utilities/Template_Polyhedra.h"

//Mixing CPU with GPU for Matrix-Vector Multiplication
template_polyhedra::ptr reachabilitySequential_GPU_MatrixVector_Multiply(
		Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen);



#endif /* REACHABILITYSEQUENTIAL_GPU_MATRIXVECTOR_MULTIPLY_CUH_ */
