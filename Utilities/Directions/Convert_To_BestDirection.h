/*
 * Convert_To_BestDirection.h
 *
 *  Created on: 07-Dec-2014
 *      Author: amit
 */

#ifndef CONVERT_TO_BESTDIRECTION_H_
#define CONVERT_TO_BESTDIRECTION_H_

#include "core_system/continuous/Polytope/Polytope.h"


/*
 * Returns the direction that is best instead of supplied dir
 *
 */
std::vector<double> getBestDirection(std::vector<double> dir);
std::vector<double> getBestDirection(polytope::ptr initial, std::vector<double> dir, int lp_solver_type_choosen);
bool Issame_direction(std::vector<double> dir1, std::vector<double> dir2);
/*
 * Returns a set of directions may be less than the
 * total number of inputs ( if there are duplicates or nearby direction(limit of near given)
 */
math::matrix<double> getBestTemplateDirections(polytope::ptr initial, math::matrix<double> directions, int lp_solver_type_choosen);


/*
 * Returns a set of directions equal to the total number of inputs ( if there are duplicates
 * or nearby direction(limit of near given) already existing then return the same direction
 */
math::matrix<double> getBestEqualTemplateDirections(polytope::ptr initial, math::matrix<double> directions, int lp_solver_type_choosen);


#endif /* CONVERT_TO_BESTDIRECTION_H_ */
