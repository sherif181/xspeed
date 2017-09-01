/*
 * Oscillator.h
 *
 *  Created on: 10-May-2016
 *      Author: amit
 */

#ifndef OSCILLATOR_H_
#define OSCILLATOR_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"
#include "core_system/math/matrix.h"

/*
 * Reference for model is https://ths.rwth-aachen.de/research/projects/hypro/filtered-oscillator/
 *
 */
void SetOscillatorParameters(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

#endif /* OSCILLATOR_H_ */
