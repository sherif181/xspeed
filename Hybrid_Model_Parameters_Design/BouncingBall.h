/*
 * BouncingBall.h
 *
 *  Created on: 17-Sep-2014
 *      Author: amit
 */

#ifndef BOUNCINGBALL_H_
#define BOUNCINGBALL_H_

/*
 * Basic Model Design requirements are ---
 *
 * Initial Polytope
 * Dynamics : Matrix A, Matrix B and Polytope U
 * Gaurd : Polytope
 * Invariant : Polytope *
 * Invariant Directions : represented in A <= b format as L1 direction and its negative will be L2 direction
 * Transition Dynamics  Rx + w : where R is a matrix and w a constant vector or a Polytope
 *
 */

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"
#include "core_system/math/matrix.h"

/*(ReachabilityParameters& reach_parameters, polytope::ptr initial_polytope_I,
 Dynamics& system_dynamics, polytope& invariant,
 polytope::ptr gaurd_polytope, TransitionDynamics& Rw)*/
void SetBouncingBall_Parameters(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

#endif /* BOUNCINGBALL_H_ */

