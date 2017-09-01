/*
 * reachability.h
 *
 *  Created on: 16-Apr-2014
 *      Author: amit
 */

#ifndef REACHABILITY_SEQ_H_
#define REACHABILITY_SEQ_H_

//#include "core_system/math/glpk_lp_solver/glpk_lp_solver.h"
#include "core_system/math/lp_solver/lp_solver.h"
#include <fstream>
#include "Utilities/invariantBoundaryCheck.h"
#include "core_system/math/matrix.h"
#include "Utilities/Template_Polyhedra.h"
#include "application/sf_utility.h"

using namespace std;

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
 * M_Matrix : is the matrix with the optimal reachability result
 * iterationNum : Number of iterations of the reachability algorithm.
 */

// Called By pure Sequential Algorithm with no critical section in this Algorithm
template_polyhedra::ptr reachabilitySequential(unsigned int NewTotalIteration, Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen);

template_polyhedra::ptr reachabilitySequential_For_Parallel_Iterations(unsigned int NewTotalIteration, Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen);

#endif /* REACHABILITY_SEQ_H_ */
