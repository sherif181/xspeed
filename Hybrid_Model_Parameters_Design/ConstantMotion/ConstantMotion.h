/*
 * ConstantMotion.h
 *
 *  Created on: 06-Jun-2016
 *      Author: amit
 */

#ifndef CONSTANTMOTION_H_
#define CONSTANTMOTION_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"

void SetConstantMotion(hybrid_automata& Hybrid_Automata,
		initial_state::ptr& init_state,
		ReachabilityParameters& reach_parameters);

#endif /* CONSTANTMOTION_H_ */
