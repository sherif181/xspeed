#include "application/reachabilityCaller.h"


void reachabilityCaller(hybrid_automata& Hybrid_Automata, std::list<initial_state::ptr>& init_state,
		ReachabilityParameters& reach_parameters, userOptions& user_options,
		int lp_solver_type_choosen, int Solver_GLPK_Gurobi_GPU, std::pair<int, polytope::ptr> forbidden_set,
		std::list<symbolic_states::ptr>& Symbolic_states_list,
		std::list<abstractCE::ptr>& ce_candidates) {

	/*
	 * For algorithm (Seq-SF, Par-SF and Time-Slice) only the Algorithm-Type/user_options.get_algorithm decides to execute
	 */
	if ((user_options.get_algorithm() == 1)
			|| (user_options.get_algorithm() == 2)
			|| (user_options.get_algorithm() == 3)) { //Sequential PostC and PostD
		//sequentialSF reach;
		reachability reach_SEQ_BFS;	//user_options.get_algorithm decides to choose Sequential Algorithm
		reach_SEQ_BFS.setReachParameter(Hybrid_Automata, init_state, reach_parameters,
				user_options.get_bfs_level(), user_options.get_algorithm(),
				user_options.getTotalSliceSize(), lp_solver_type_choosen, user_options.getStreamSize(),
				Solver_GLPK_Gurobi_GPU, forbidden_set);
		if (user_options.get_algorithm() == 1) {
			std::cout << "\nRunning sequential BFS.\n";
		} else if (user_options.get_algorithm() == 2) { //ParallelSF reach;
			//Parallel PostC using Lazy SF algorithm and Sequential PostD
			std::cout
					<< "\nRunning parallel PostC using lazy SF algorithm and sequential PostD.\n";
		} else if (user_options.get_algorithm() == 3) { //Parallel PostC using Time-Slice algorithm and Sequential PostD
			std::cout
					<< "\nRunning parallel PostC using Time-Slice algorithm and sequential PostD.\n";
		}
		Symbolic_states_list = reach_SEQ_BFS.computeSeqentialBFSReach(ce_candidates);

	} else if (user_options.get_algorithm() == 4) { //Adaptation of Gerard J. Holzmann's algorithm (Seq PostC and PBFS)
		agjh reach_AGJH;
		reach_AGJH.setReachParameter(Hybrid_Automata, init_state, reach_parameters,
						user_options.get_bfs_level(), user_options.get_algorithm(),
						user_options.getTotalSliceSize(), lp_solver_type_choosen, user_options.getStreamSize(),
						Solver_GLPK_Gurobi_GPU, forbidden_set);
				std::cout << "\nRunning adaptation of Gerard J. Holzmann's algorithm.\n";
				Symbolic_states_list = reach_AGJH.ParallelBFS_GH(); // Holzmann algorithm adaptation
			} else if (user_options.get_algorithm() == 5) { // TPBFS -- Task parallel BFS (Load Balancing Algorithm)
				tpbfs reach;
				reach.setReachParameter(Hybrid_Automata, init_state, reach_parameters,
								user_options.get_bfs_level(), user_options.get_algorithm(),
								user_options.getTotalSliceSize(), lp_solver_type_choosen, user_options.getStreamSize(),
								Solver_GLPK_Gurobi_GPU, forbidden_set);
						std::cout<< "\nRunning Task parallel (Load Balancing) BFS algorithm.\n";
						Symbolic_states_list = reach.LoadBalanceAll(
								ce_candidates);
			} else if (user_options.get_algorithm() == 6) { //gpu-postc -- PostC in GPU and Sequential BFS
						reachability reach;
						reach.setReachParameter(Hybrid_Automata, init_state, reach_parameters,
										user_options.get_bfs_level(), user_options.get_algorithm(),
										user_options.getTotalSliceSize(), lp_solver_type_choosen, user_options.getStreamSize(),
										Solver_GLPK_Gurobi_GPU, forbidden_set);

							std::cout << "\nRunning PostC in GPU and Sequential BFS.\n";
							//Symbolic_states_list = reach.LoadBalanceAll(ce_candidates);
			}

	if(user_options.get_algorithm() == 7){
		AGJHGPU aghjGpu;
		aghjGpu.setReachParameter(Hybrid_Automata, init_state, reach_parameters, user_options.get_bfs_level(), user_options.get_algorithm(), user_options.getTotalSliceSize(), lp_solver_type_choosen, user_options.getStreamSize(), Solver_GLPK_Gurobi_GPU, forbidden_set);
		Symbolic_states_list = aghjGpu.gpuReachability(ce_candidates);
	}

}

