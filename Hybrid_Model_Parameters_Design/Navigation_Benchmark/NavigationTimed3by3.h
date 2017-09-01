/*
 * NavigationTimed3by3.h
 *
 *  Created on: 22-May-2016
 *      Author: rajarshi
 */

/**
 * Navigation 3 by 3 model with time evolution added
 * to the dynamics.
 */

#ifndef NAVIGATIONTIMED3BY3_H_
#define NAVIGATIONTIMED3BY3_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"

void Set_NavTimed_Parameters(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

#endif /* NAVIGATIONTIMED3BY3_H_ */
