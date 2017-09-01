/*
 * Rotation_Circle.h
 *
 *	Model with Two Locations First the Upper-Half and Second the Lower-Half
 *
 *  Created on: 17-Sep-2014
 *      Author: amit
 */

#ifndef ROTATION_TIMED_CIRCLE_H_
#define ROTATION_TIMED_CIRCLE_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"
#include "core_system/math/matrix.h"

void SetRotationTimedCircle_Parameters(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

void SetRotationTimedCircle_Parameters(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

#endif /* ROTATION_TIMED_CIRCLE_H_ */
