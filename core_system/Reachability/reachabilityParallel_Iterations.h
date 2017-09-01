/*
 * reachabilityParallel_Iterations.h
 *
 *  Created on: 11-Nov-2014
 *      Author: amit
 */

#ifndef REACHABILITYPARALLEL_ITERATIONS_H_
#define REACHABILITYPARALLEL_ITERATIONS_H_


#include "core_system/math/glpk_lp_solver/glpk_lp_solver.h"
#include "Utilities/invariantBoundaryCheck.h"
#include "Utilities/Template_Polyhedra.h"
#include "core_system/math/matrix.h"
#include <fstream>
#include <omp.h>



const template_polyhedra reachabilityParallelIters(Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial, ReachabilityParameters& ReachParameters,
		polytope::ptr invariant, bool isInvariantExist);




#endif /* REACHABILITYPARALLEL_ITERATIONS_H_ */
