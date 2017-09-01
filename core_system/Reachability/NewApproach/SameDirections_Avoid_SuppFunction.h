/*
 * SameDirections_Avoid_SuppFunction.h
 *
 *  Created on: 12-Mar-2015
 *      Author: amit
 */

#ifndef SAMEDIRECTIONS_AVOID_SUPPFUNCTION_H_
#define SAMEDIRECTIONS_AVOID_SUPPFUNCTION_H_

#include "core_system/math/lp_solver/lp_solver.h"
#include <fstream>
#include "Utilities/invariantBoundaryCheck.h"
#include "core_system/math/matrix.h"
#include "Utilities/Template_Polyhedra.h"
#include "core_system/Reachability/NewApproach/utility_functions.h"

using namespace std;

typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;

/*
 * Called By pure Sequential Algorithm with no critical section in this Algorithm
 */
template_polyhedra::ptr reachabilitySameDirection(Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen);

#endif /* SAMEDIRECTIONS_AVOID_SUPPFUNCTION_H_ */
