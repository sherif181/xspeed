#ifndef FIVE_DIM_SYS_H_
#define FIVE_DIM_SYS_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"
#include "core_system/math/matrix.h"

void setSysParams(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

#endif		/* FIVE_DIM_SYS_H_	*/
