/*
 * reachability_Sequential_AllDirections.h
 *
 *  Created on: 23-Mar-2015
 *      Author: amit
 */

#ifndef REACHABILITY_SEQUENTIAL_ALLDIRECTIONS_H_
#define REACHABILITY_SEQUENTIAL_ALLDIRECTIONS_H_

#include "core_system/math/lp_solver/lp_solver.h"
//#include <fstream>
#include "Utilities/invariantBoundaryCheck.h"
#include "core_system/math/matrix.h"
#include "Utilities/Template_Polyhedra.h"
#include "application/sf_utility.h"

using namespace std;

typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;

/*
 * Sequential Algorithm with all the Directions(including transposed directions)
 */

template_polyhedra reachabilitySequentialAllDirections(Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen);

#endif /* REACHABILITY_SEQUENTIAL_ALLDIRECTIONS_H_ */
