/*
 * reachabilityParallel_Process.h
 *
 *  Created on: 22-May-2014
 *      Author: gurung
 */

#ifndef REACHABILITYPARALLEL_PROCESS_H_
#define REACHABILITYPARALLEL_PROCESS_H_

#include "core_system/math/glpk_lp_solver/glpk_lp_solver.h"
#include "application/DataStructureDirections.h"
#include "core_system/continuous/Polytope/Polytope.h"
#include "Utilities/invariantBoundaryCheck.h"
/*#include <fstream>*/

#include <string.h>
#include <stdio.h>
#include <errno.h>

#include <iostream>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdlib.h>   // Declaration for exit()
#include "core_system/math/matrix.h"
#include "Utilities/Template_Polyhedra.h"
#include "application/sf_utility.h"
#include "application/CopyArray.h"

//#include <sys/types.h>	//required for SharedMemory Operation
#include <sys/ipc.h>		//required for SharedMemory Operation
#include <sys/shm.h>		//required for SharedMemory Operation
#define SHMSZ     27		//required for SharedMemory Operation

typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;

/**
 * coeffMatrix : All facets of the Polytope I
 * bMatrix : All Bound-Value of the Polytope I
 * boundBGreaterThenExpr : All Bound-Sign of the Polytope I  ... 1 for <= and 0 for >=
 *
 * VcoeffMatrix : All facets of the Polytope V
 * VbMatrix : All Bound-Value of the Polytope V
 * VboundBGreaterThenExpr : All Bound-Sign of the Polytope V  ... 1 for <= and 0 for >=
 *
 * AMatrix : is the flow matrix
 * Vector_R : is the matrix (m x n) consisting the list of 'm' directions of 'n' variables each.
 *
 * TotalIteration : Number of iterations of the reachability algorithm.
 * sf_vals : is a row of the resultant optimal reachability matrix computed by each forked process.
 */

void reachFunction(unsigned int eachDirection, Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		lp_solver &s_per_thread_I, lp_solver &s_per_thread_U,
		double* sf_vals);

const template_polyhedra::ptr reachabilityParallel_Process(Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen);

#endif /* REACHABILITYPARALLEL_PROCESS_H_ */
