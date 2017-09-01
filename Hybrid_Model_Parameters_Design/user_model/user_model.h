/*
 * user_model.h
 *
 *  Created on: 26-Nov-2015
 *      Author: amit
 */

#ifndef USER_MODEL_H_
#define USER_MODEL_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"
#include "core_system/math/matrix.h"
#include "application/All_PP_Definition.h"
#include "core_system/math/uni_sphere.h"	//for obtaining uniformly distributed directions
#include "application/sf_directions.h"
#include "application/sf_utility.h"
#include "application/userOptions.h"


//void user_model(hybrid_automata& Hybrid_Automata,initial_state::ptr& init_state,ReachabilityParameters& reach_parameters);

void user_model(hybrid_automata& ha,
		std::list<initial_state::ptr>& init_state,
		ReachabilityParameters& reach_parameters,
		userOptions& op);

#endif /* USER_MODEL_H_ */
