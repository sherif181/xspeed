/*
 * TTEthernetModel.h
 *
 *  Created on: 3-Oct-2016
 *      Author: amit
 */

#ifndef TTETHERNET_H_
#define TTETHERNET_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"
#include "core_system/math/matrix.h"


/*
 * * Reference for model is http://cps-vo.org/node/26547
 *
 * Model created with Hyst after converting SpaceEx model to flat xml model using Hyst and then to that model to .cpp
 */
void setTTEthernetModel2(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);


#endif /* TTETHERNET_H_ */
