/*
 * Rotation_Circle_FourLocation.h
 *
 *  Created on: 24-Feb-2015
 *      Author: amit
 */

#ifndef ROTATION_CIRCLE_FOURLOCATION_H_
#define ROTATION_CIRCLE_FOURLOCATION_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"
#include "core_system/math/matrix.h"

void SetRotationCircle4Location_Parameters(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

#endif /* ROTATION_CIRCLE_FOURLOCATION_H_ */
