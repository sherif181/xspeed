/*
 * reachParallelExplore.h
 *
 *  Created on: 16-Nov-2014
 *      Author: amit
 */

#ifndef REACHPARALLELEXPLORE_H_
#define REACHPARALLELEXPLORE_H_

#include "core_system/math/glpk_lp_solver/glpk_lp_solver.h"
//#include "matrixOperation.h"
#include "core_system/math/matrix.h"
#include <fstream>
#include <omp.h>
#include "core_system/continuous/ConvexSet/transMinkPoly.h"
#include "Utilities/Template_Polyhedra.h"
#include "application/DataStructureDirections.h"


const template_polyhedra::ptr reachabilityParallel(
			unsigned int NewTotalIteration, Dynamics& SystemDynamics,
			supportFunctionProvider::ptr Initial,
			ReachabilityParameters& ReachParameters, polytope::ptr invariant,
			bool isInvariantExist, int lp_solver_type_choosen);






/**
 * Reachability function which explores the state space in parallel starting at different initial sets.
 *
 * CORES: Its the number of parallel compute cores available in the hardware architecture.
 */
const template_polyhedra::ptr reachParallelExplore(unsigned int NewTotalIteration, Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, const polytope::ptr invariant,
		bool isInvariantExist, int CORES, unsigned int Algorithm_Type,
		int lp_solver_type_choosen);


#endif /* REACHPARALLELEXPLORE_H_ */
