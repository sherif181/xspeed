/*
 * TimedBouncingBall.h
 *
 *  Created on: 17-Sep-2014
 *      Author: amit
 */

#ifndef TIMEDBOUNCINGBALL_H_
#define TIMEDBOUNCINGBALL_H_

/*
 * Basic Model Design requirements are ---
 *
 * Initial Polytope
 * Dynamics : Matrix A, Matrix B and Polytope U
 * Gaurd : Polytope
 * Invariant : Polytope *
 * Transition Dynamics  Rx + w : where R is a Assignment_mapping matrix and w a constant vector or a Polytope
 *
 */

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/math/matrix.h"
#include "application/DataStructureDirections.h"

#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"

void SetTimedBouncingBall_ParametersOurOutput(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);
void SetTimedBouncingBall_ParametersHystOutput(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

//Testing for 2 initial set
void SetTimedBouncingBall_2initSet(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

#endif /* TIMEDBOUNCINGBALL_H_ */
