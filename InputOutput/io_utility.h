/*
 * io_utility.h
 *
 *  Created on: 19-Jan-2016
 *      Author: amit
 */

#ifndef IO_UTILITY_H_
#define IO_UTILITY_H_

#include <vector>
#include <list>
#include <utility>
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"
#include "core_system/continuous/Polytope/Polytope.h"
#include "Utilities/StandardVector.h"

typedef std::vector<std::pair<double, double> > Intervals;

void Interval_Generator(std::list<symbolic_states::ptr>& Symbolic_states_list,
		std::list<std::pair<int, Intervals> > & location_interval_outputs, initial_state::ptr& init_state);

polytope::ptr convertBounding_Box(template_polyhedra::ptr sfm);

#endif /* IO_UTILITY_H_ */

