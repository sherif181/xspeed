/*
 * NavigationBenchmark.h
 *
 *  Created on: 24-April-2016
 *      Author: amit
 *      Navigation Benchmark:: deals with an object that moves in the R^2 plane with
 *      velocity V_d is determined by the position of the object in an n x m grid, and the desired velocities may
 *      take values [sin(i * pi/4), cos(i * pi/4)], for i=0, . . . , 7.
 *      The Grid label can as obtained from Author's website
 *
 *   detail can be obtained from the paper title "Benchmarks for Hybrid System Verification" by Ansgar Fehnker and Franjo Ivancic
 */

/*
 * Basic Model Design requirements are ---
 *
 * Initial Polytope
 * Dynamics : Matrix A, Matrix B and Polytope U :: x' = Ax + Bu and u belongs to U. where x and u are convex polytope
 * Gaurd : Polytope
 * Invariant : Polytope *
 * Transition Dynamics  Rx + w : where R is a Assignment_mapping matrix and w a constant vector or a Polytope
 *
 */

#ifndef NAVIGATIONBENCHMARK4Var_H_
#define NAVIGATIONBENCHMARK4Var_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "core_system/symbolic_states/initial_state.h"

void SetNavigationBenchMark4Var(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

//Only the desired-velocities of  Loc 3 and Loc 4 have been changed
void SetNavigationModel2(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

/*
 * Only the desired-velocities of  Loc 1, Loc 2, Loc 3 and Loc 4 have been changed
 * Standard model used by most researchers by the Model name : NAV04
 */
void SetNavigationModel4(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

void SetNavigationModel5by5(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters);

//Model (5 x 5) containing 81-Locations with 280 transitions
/*void SetNavigationModel9by9(hybrid_automata& Hybrid_Automata,
 initial_state::ptr& init_state,
 ReachabilityParameters& reach_parameters); */

#endif /* NAVIGATIONBENCHMARK4Var_H_ */
