/*
 * reachability.h
 *
 *  Created on: 23-Feb-2016
 *      Author: amit
 */

#ifndef REACHABILITY_H_
#define REACHABILITY_H_

#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/initial_state.h"
#include "application/DataStructureDirections.h"
#include "counterExample/abstractCE.h"
#include <utility>	//for std::pair
#include <set>
//***************** From Sequential BFS *****************************
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "counterExample/abstract_symbolic_state.h"
#include "InputOutput/io_utility.h"

//#include "core_system/PWL/pwl.h"
#include "core_system/PWL/pwlist.h"
#include "core_system/symbolic_states/initial_state.h"
#include "core_system/symbolic_states/symbolic_states_utility.h"

#include "Utilities/Template_Polyhedra.h"
#include "Utilities/testPolytopePlotting.h"
#include "Utilities/Post_Assignment.h"
#include <boost/lexical_cast.hpp>
#include "core_system/Reachability/reachabilitySequential.h"
#include "core_system/Reachability/reachParallelExplore.h"
#include "core_system/Reachability/reachabilityParallel_Process.h"
#include "core_system/Reachability/NewApproach/SameDirections_Avoid_SuppFunction.h"
#include "Utilities/computeInitialPolytope.h"
//#include "core_system/math/gurobi_lp_solver/gurobi_lp_solver.h"
#include "application/All_PP_Definition.h"
//#include "core_system/Reachability/reachabilitySequential_GPU_MatrixVector_Multiply.cuh"
#include "core_system/Reachability/GPU_Reach/reach_Sequential_GPU.cuh"
#include "boost/timer/timer.hpp"
#include "counterExample/abstractCE.h"
#include <omp.h>
#include "LockAvoidanceUtility.h"
#include "reachDataStructure.h"
#include <vector>
#include "application/sf_directions.h"
#include "application/sf_utility.h"
//***************** End Sequential BFS *****************************

using namespace std;

class reachability {

public:
/*	reachability(){
	//	cout<< "Calling from reachability Class\n";
	}*/
	void setReachParameter(hybrid_automata& H, std::list<initial_state::ptr>& I,
			ReachabilityParameters& reach_parameters, int bound,
			unsigned int Algorithm_Type, unsigned int Total_Partition,
			int lp_solver_type_choosen, unsigned int number_of_streams,
			int Solver_GLPK_Gurobi_GPU,
			std::pair<int, polytope::ptr> forbidden_set);

	//bound is the maximum number of transitions or jumps permitted.
	//reach_parameters includes the different parameters needed in the computation of reachability.
	//I is the initial symbolic state
	//Although this interface can be pushed in a separate sequential class but it can also be used to call par_SF and time_slice algorithms.
	std::list<symbolic_states::ptr> computeSeqentialBFSReach(std::list<abstractCE::ptr>& ce_candidates);


	//1)  Parallel Breadth First Search for Discrete Jumps with critical section for adding newSymbolicState
	std::list<symbolic_states::ptr> computeParallelBFSReach(
			std::list<abstractCE::ptr>& ce_candidates);

	//2)  Lock Avoidance for adding newSymbolicState:: Parallel Breadth First Search for Discrete Jumps
	//separate Read and Write Queue (pwlist.WaitingList)
		std::list<symbolic_states::ptr> computeParallelBFSReachLockAvoid(std::list<abstractCE::ptr>& ce_candidates);


	/*
	 * List of private variables now converted into public due to class inheritance framework
	 */
		std::list<initial_state::ptr> I; //converted to a list of initial state
		int bound;
		ReachabilityParameters reach_parameters;
		hybrid_automata H; //todo:: have to change it to boost::ptr
		int lp_solver_type_choosen;
		std::pair<int, polytope::ptr> forbidden_set;
		int Solver_GLPK_Gurobi_GPU;

		void sequentialReachSelection(unsigned int NewTotalIteration, location::ptr current_location, polytope::ptr continuous_initial_polytope,
						template_polyhedra::ptr& reach_region);

private:

	//initial_state::ptr I;


	unsigned int Algorithm_Type;
	unsigned int Total_Partition;
	unsigned int number_of_streams;



	void parallelReachSelection(unsigned int NewTotalIteration, location::ptr current_location, polytope::ptr continuous_initial_polytope,
			ReachabilityParameters& reach_parameters, std::vector<symbolic_states::ptr>& S, unsigned int id);


	/*Returns True, if safety has been violated on the current computed Symbolic States and sets/creates the counterExample class
	Returns False, if safety not violated the counterExample list remains empty*/
	bool safetyVerify(symbolic_states::ptr& computedSymStates, std::list<symbolic_states::ptr>& Reachability_Region, std::list<abstractCE::ptr>& ce);

	/*Performs the task of computing support functions on the specified polytopes for the given directions
		and returns the results of support function computation on the same data-structure as be symbolic-state basis
		SEQUENTIAL*/
	void computeBIG_Task(std::vector<LoadBalanceData>& LoadBalanceDS);

	/*
	 * Performs the task of computing support functions in parallel using the approach of task_scheduling
	 */

	void parallelBIG_Task(std::vector<LoadBalanceData>& LoadBalanceDS);

	double boxLPSolver(polytope::ptr poly, std::vector<double> dir);

	double LPSolver(polytope::ptr poly, std::vector<double> dirs);

};

#endif /* REACHABILITY_H_ */
