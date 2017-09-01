/*
 * reachabilityCaller.h
 *
 *  Created on: 27-Oct-2016
 *      Author: amit
 */

#ifndef REACHABILITYCALLER_H_
#define REACHABILITYCALLER_H_

#include "userOptions.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/Reachability/reachability.h"
#include "counterExample/abstractCE.h"
#include "application/userOptions.h"

#include <list>

#include "core_system/Reachability/SequentialSF.h"
#include "core_system/Reachability/ParallelSF.h"
#include "core_system/Reachability/TimeSliceSF.h"
#include "core_system/Reachability/AGJH.h"
#include "core_system/Reachability/TPBFS.h"
#include "core_system/Reachability/AGJHGPU.cuh"

void reachabilityCaller(hybrid_automata& H, std::list<initial_state::ptr>& I,
		ReachabilityParameters& reach_parameters, userOptions& user_options,
		int lp_solver_type_choosen, int Solver_GLPK_Gurobi_GPU, std::pair<int, polytope::ptr> forbidden_set,
		std::list<symbolic_states::ptr>& Symbolic_states_list,
		std::list<abstractCE::ptr>& ce_candidates);

#endif /* REACHABILITYCALLER_H_ */
