/*
 * Partition_BoundingPolytope.h
 *
 *  Created on: 03-Mar-2015
 *      Author: amit
 */

#ifndef PARTITION_BOUNDINGPOLYTOPE_H_
#define PARTITION_BOUNDINGPOLYTOPE_H_

#include "application/sf_directions.h"
#include "application/sf_utility.h"
#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/math/matrix.h"
#include <list>

/*
 * Function takes a polytope as input and returns a list of polytopes after partitioning
 * Inputs parameters:
 * 		polytope::ptr S : initial input polytope S
 * 		Partition number of variables : integer N
 * 		(eg 1 indicate variable x1, 2 indicate variable x1 and x2, 3 indicate variable x1, x2 and x3)
 * Returns : a list of polytopes after partitioning the polytope S.
 */
std::list<polytope::ptr> Partition_BoundingPolytope(polytope::ptr S, unsigned int nVar, unsigned int p) ;



#endif /* PARTITION_BOUNDINGPOLYTOPE_H_ */
