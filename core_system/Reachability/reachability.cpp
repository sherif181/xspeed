#include "reachability.h"
#include "InputOutput/cpu_utilities/cpu_utilities.h"
#include "Utilities/flow_cost_estimate.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include <ctime>

using namespace std;


void reachability::setReachParameter(hybrid_automata& h, std::list<initial_state::ptr>& i,
		ReachabilityParameters& reach_param, int bound_limit,
		unsigned int algorithm_type, unsigned int total_partition,
		int lp_solver_type, unsigned int streams_size,
		int solver_GLPK_Gurobi_for_GPU,
		std::pair<int, polytope::ptr> forbidden) {
	//std::cout<<"setReachParameter function called"<<std::endl;
	H = h;
	I = i;
	reach_parameters = reach_param;
	bound = bound_limit;
	Algorithm_Type = algorithm_type;	//Important parameter to decide to select an algorithm to execute
	Total_Partition = total_partition;
	lp_solver_type_choosen = lp_solver_type;
	number_of_streams = streams_size;
	Solver_GLPK_Gurobi_GPU = solver_GLPK_Gurobi_for_GPU; //todo:: used for comparing GLPK solver vs GPU. Can be removed
	forbidden_set = forbidden;
}



//bound is the maximum number of transitions or jumps permitted.
//reach_parameters includes the different parameters needed in the computation of reachability.
//I is the initial symbolic state
std::list<symbolic_states::ptr> reachability::computeSeqentialBFSReach(std::list<abstractCE::ptr>& ce_candidates) {

	//std::cout<<"Hello from Sequential Reach!!"<<std::endl;
	std::list < symbolic_states::ptr > Reachability_Region;
	template_polyhedra::ptr reach_region;

	int number_times = 0, BreadthLevel = 0, previous_level = -1;
	std::list<int> queue; // data structure to keep track of the number of transitions
	//discrete_set discrete_state;

	pwlist pw_list; //list of initial_state
	for (std::list<initial_state::ptr>::iterator i=I.begin();i!=I.end();i++){
		pw_list.WaitingList_insert(*(i));
		queue.push_back(BreadthLevel); //insert at REAR, first Location
		//std::cout<<"initial set pushed in pw_list.WaitingList()\n";
	}

	//queue.push_back(BreadthLevel); //insert at REAR, first Location
	bool starting_location = true;
	bool saftey_violated = false;

	polytope::ptr continuous_initial_polytope;

	unsigned int num_flowpipe_computed=0;	//keeping track of number of flowpipe computed
	while (!pw_list.isEmpty_WaitingList()) {
		symbolic_states::ptr S = symbolic_states::ptr(new symbolic_states()); //required to be pushed into the Reachability_Region
		initial_state::ptr U;
		U = pw_list.WaitingList_delete_front();
		int levelDeleted = queue.front(); //get FRONT element
		queue.pop_front(); //delete from FRONT
		if (levelDeleted > bound)
			break; //stopping due to number of transitions exceeds the bound

		int location_id;
	//	cout<<"\nTesting 2 a 2\n";
		//discrete_state = U.getDiscreteSet();
		location_id = U->getLocationId();
		//std::cout<<"location_id from U = "<<location_id<<std::endl;
		discrete_set discrete_state;
		discrete_state.insert_element(location_id); //creating discrete_state

		//continuous_initial_polytope = U.getInitial_ContinousSetptr();
		continuous_initial_polytope = U->getInitialSet();
		reach_parameters.X0 = continuous_initial_polytope;

		S->setDiscreteSet(discrete_state);
		S->setParentPtrSymbolicState(U->getParentPtrSymbolicState()); //keeps track of parent pointer to symbolic_states
		S->setTransitionId(U->getTransitionId()); //keeps track of originating transition_ID


		pw_list.PassedList_insert(U);

		location::ptr current_location;

		current_location = H.getLocation(location_id);
		string name = current_location->getName();
		//cout << "\nLocation ID = "<<location_id<<", Location Name = " << name << "\n";
		if ((name.compare("GOOD") == 0) || (name.compare("BAD") == 0)
				|| (name.compare("UNSAFE") == 0)
				|| (name.compare("FINAL") == 0))
			continue; //do not compute the continuous reachability algorithm
		bool foundSkippingLocation = false;
		for (int StopLoc = 0; StopLoc < reach_parameters.Stop_locID.size();
				StopLoc++) {
			if (location_id == reach_parameters.Stop_locID[StopLoc]) {
				foundSkippingLocation = true;
			}
		}
		if (foundSkippingLocation) //do not compute the continuous reachability algorithm
			continue;

		// ******************* Computing Parameters *******************************
		//In this current_location the parameters will change for eg., alfa, beta and phi_trans  have to be re-computed
		/*
		 * Computing the parameters to avoid multiple computation in the child process
		 * Items Required :: time_step, phi_trans , B_trans, compute_alfa,compute_beta
		 * Computation of compute_alfa depends on initial set. For algorithm PAR_BY_PARTS where the
		 * initial set in divided into parts. Compute_alfa should be computed for each initial sets.
		 * */
	//	cout<<"\nTesting 2 c\n";
		double result_alfa = compute_alfa(reach_parameters.time_step,
				current_location->getSystem_Dynamics(),
				continuous_initial_polytope, lp_solver_type_choosen); //2 glpk object created here

		double result_beta = compute_beta(current_location->getSystem_Dynamics(),
				reach_parameters.time_step, lp_solver_type_choosen); // NO glpk object created here

		reach_parameters.result_alfa = result_alfa;
		reach_parameters.result_beta = result_beta;
		// Intialised the transformation and its transpose matrix
		math::matrix<double> phi_matrix, phi_trans;

		if (!current_location->getSystem_Dynamics().isEmptyMatrixA) { //if A not Empty
			current_location->getSystem_Dynamics().MatrixA.matrix_exponentiation(
					phi_matrix, reach_parameters.time_step);
			phi_matrix.transpose(phi_trans);
			reach_parameters.phi_trans = phi_trans;
		}
		math::matrix<double> B_trans;
		// transpose to be done once and stored in the structure of parameters
		if (!current_location->getSystem_Dynamics().isEmptyMatrixB) { //if B not Empty
			current_location->getSystem_Dynamics().MatrixB.transpose(B_trans);
			reach_parameters.B_trans = B_trans;
		}
	//	cout<<"\nTesting 2 d\n";
		// ******************* Computing Parameters *******************************


		unsigned int NewTotalIteration = reach_parameters.Iterations;

		// ************ Compute flowpipe_cost:: estimation Starts **********************************

		if (current_location->isInvariantExists()){
			if (current_location->getSystem_Dynamics().isEmptyMatrixB==true && current_location->getSystem_Dynamics().isEmptyC==true){
				//Approach of Coarse-time-step and Fine-time-step
				jumpInvariantBoundaryCheck(current_location->getSystem_Dynamics(), continuous_initial_polytope, reach_parameters,
					current_location->getInvariant(), lp_solver_type_choosen, NewTotalIteration);
				cout<<"Running Approach of Coarse-time-step and Fine-time-step\n";
			}else{
				//Approach of Sequential invariant check will work for all case
				InvariantBoundaryCheck(current_location->getSystem_Dynamics(), continuous_initial_polytope,
					reach_parameters, current_location->getInvariant(), lp_solver_type_choosen, NewTotalIteration);
			}
		}
		// ************ Compute flowpipe_cost:: estimation Ends **********************************

		//cout<<"\nTesting 2 e\n";

		sequentialReachSelection(NewTotalIteration, current_location, continuous_initial_polytope, reach_region);
		//cout<<"\nTesting 2 e 1\n";
		num_flowpipe_computed++;//computed one Flowpipe
		//cout<<"\nTesting 2 e 2\n";
		//	*********************************************** Reach or Flowpipe Computed ************************************
		if (previous_level != levelDeleted) {
			previous_level = levelDeleted;
			BreadthLevel++;
		}

		//reach_region->Testing_print();
		//cout<<"\nTesting 2 e 3\n";
		if (reach_region->getTotalIterations() != 0) {
			S->setContinuousSetptr(reach_region);
			Reachability_Region.push_back(S);
		}
		//  ******************************** Safety Verification section ********************************
		std::list < symbolic_states::ptr > list_sym_states;
		std::list < abstract_symbolic_state::ptr > list_abstract_sym_states;
		polytope::ptr abs_flowpipe; //bounding_box Polytope
		polytope::ptr polyI; //initial polytope of the abstract flowpipe
		std::list < transition::ptr > list_transitions;
//		cout<<"\nTesting 2 f\n";
		if (reach_region->getTotalIterations() != 0 && forbidden_set.second != NULL) { //flowpipe exists
				//so perform intersection with forbidden set provided locID matches
			int locID = current_location->getLocId();
			cout<<"Running Safety Check for Loc = "<<locID<<std::endl;

			if (forbidden_set.first==-1 || locID == forbidden_set.first) { //forbidden locID matches
				polytope::ptr forbid_poly = forbidden_set.second;
				std::list < template_polyhedra::ptr > forbid_intersects;
				forbid_intersects = reach_region->polys_intersectionSequential(forbid_poly, lp_solver_type_choosen);
				if (forbid_intersects.size() == 0) {
					std::cout << "\nThe model does NOT violates SAFETY property!!!\n";
				} else {
					std::cout << "\nThe model violates SAFETY property!!!\n";
					symbolic_states::ptr current_forbidden_state;
					current_forbidden_state = S;
					std::cout << "\nReverse Path Trace =>\n";
					int cc = 0;
					do {
						int locationID, locationID2;
						discrete_set ds, ds2;
						ds = current_forbidden_state->getDiscreteSet();
						//insert discrete_set in the abstract_symbolic_state
						// ***********insert bounding_box_polytope as continuousSet in the abstract_symbolic_state***********
						abs_flowpipe = convertBounding_Box(
								current_forbidden_state->getContinuousSetptr());
						// The below gets the first polytope from the SFM. Could be made even better
						// by taking the exact initial set X0.

						polyI = current_forbidden_state->getInitial_ContinousSetptr();

						// Here create a new abstract_symbolic_state
						abstract_symbolic_state::ptr curr_abs_sym_state;
						curr_abs_sym_state = abstract_symbolic_state::ptr(
								new abstract_symbolic_state(ds, abs_flowpipe,
										polyI));

						// ***********insert bounding_box_polytope as continuousSet in the abstract_symbolic_state***********
						for (std::set<int>::iterator it =
								ds.getDiscreteElements().begin();
								it != ds.getDiscreteElements().end(); ++it)
							locationID = (*it); //Assuming only a single element exist in the discrete_set
						int transID =
								current_forbidden_state->getTransitionId(); //a)
						//   **********************************************************
						//create an object of abstractCE[1)list_of_symbolic_states 2)list_of_transition and 3) length]
						//1) ******************** list_of_symbolic_states ********************
						list_sym_states.push_front(current_forbidden_state); //pushing the bad symbolic_state first(at the top)
						list_abstract_sym_states.push_front(curr_abs_sym_state);
						//2) list_of_transition
						//a) current sym_state only have trans_ID but to retrieve this transition I have to
						//b) get the parent to this sym_state using getParentPtrSymbolicState and then in
						//c) that sym_state get the discrete state and then get its loc_ID.
						//d) Now from the class hybrid_automata get an object of the class location
						//e) Now for this location's object get the transition with transID as transition_ID from
						// the data member out_going_transitions.
						//3) length: number of transitions
						//   **********************************************************
						if (cc != 0) {
							std::cout << " --> ";
						}
						std::cout << "(" << locationID << ", " << transID << ")";
						if (current_forbidden_state->getParentPtrSymbolicState() != NULL) { //searching only if not NULL
							//cout<<"check if parentPtr not NULL\n";
							current_forbidden_state = searchSymbolic_state(Reachability_Region,
											current_forbidden_state->getParentPtrSymbolicState()); //b)
							//2) ******************* list_transitions ********************
							ds2 = current_forbidden_state->getDiscreteSet(); //c)
							for (std::set<int>::iterator it =
									ds2.getDiscreteElements().begin();
									it != ds2.getDiscreteElements().end(); ++it)
								locationID2 = (*it); //c)

							location::ptr object_location;
							object_location = H.getLocation(locationID2); //d)

							transition::ptr temp =
									object_location->getTransition(transID); //e)
							list_transitions.push_front(temp); //pushing the transition in the stack
							//2) ******************* list_transitions Ends ********************
						}
						cc++;
					} while (current_forbidden_state->getParentPtrSymbolicState()!= NULL);

					if ((cc >= 1) && (current_forbidden_state->getParentPtrSymbolicState()== NULL)) { //root is missed
						int locationID;
						discrete_set ds;
						ds = current_forbidden_state->getDiscreteSet();
						abs_flowpipe = convertBounding_Box(
								current_forbidden_state->getContinuousSetptr());
						// The below gets the first polytope from the SFM. Could be made even better
						// by taking the exact initial set X0.

						polyI = current_forbidden_state->getInitial_ContinousSetptr();

						// Create a new abstract_symbolic_state
						abstract_symbolic_state::ptr curr_abs_sym_state;

						curr_abs_sym_state = abstract_symbolic_state::ptr(new abstract_symbolic_state(ds, abs_flowpipe, polyI));

						for (std::set<int>::iterator it =
								ds.getDiscreteElements().begin();
								it != ds.getDiscreteElements().end(); ++it)
							locationID = (*it); //Assuming only a single element exist in the discrete_set

						int transID =
								current_forbidden_state->getTransitionId();
						list_sym_states.push_front(current_forbidden_state); //1) pushing the initial/root bad symbolic_state at the top
						list_abstract_sym_states.push_front(curr_abs_sym_state);
						std::cout << " -->  (" << locationID << ", " << transID << ")\n";
					}
					saftey_violated = true;
					abstractCE::ptr ce = abstractCE::ptr(new abstractCE());
					ce->set_length(cc);
					ce->set_sym_states(list_abstract_sym_states);
					ce->set_transitions(list_transitions);
					hybrid_automata::ptr ha = hybrid_automata::ptr(
							new hybrid_automata(H));
					ce->set_automaton(ha);
					ce->set_forbid_poly(forbidden_set.second);
					ce_candidates.push_back(ce); // ce added to the abstract ce candidates list.
				} // end of condition when forbidden state intersects with the flowpipe set
			} //end of condition when forbidden state loc id matches with flowpipe loc id
		} //computed flowpipe is not empty

//		if (saftey_violated) {
//			break; //no need to compute rest of the locations
//		}
		//  ******************************** Safety Verification section Ends********************************
		//  ******* ---POST_D Begins--- ******* Check to see if Computed FlowPipe is Empty  **********
		if (reach_region->getTotalIterations() != 0 && BreadthLevel <= bound) {
			//computed reach_region is empty and optimize transition BreadthLevel-wise
			for (std::list<transition::ptr>::iterator t = current_location->getOut_Going_Transitions().begin();
					t != current_location->getOut_Going_Transitions().end(); t++) {
				// get each destination_location_id and push into the pwl.waiting_list
				int transition_id = (*t)->getTransitionId();
				location::ptr current_destination;
				Assign current_assignment;
				polytope::ptr gaurd_polytope;

				polytope::ptr intersectedRegion; //created two objects here
				discrete_set ds;

				current_destination = H.getLocation((*t)->getDestination_Location_Id());

				string locName = current_destination->getName();
				std::list<polytope::ptr> polys;
				gaurd_polytope = (*t)->getGaurd(); //	GeneratePolytopePlotter(gaurd_polytope);
				if (!gaurd_polytope->getIsUniverse())	//Todo guard and invariants in the model: True is universal and False is unsatisfiable/empty
				{
					//std::cout<<"Not universal\n";
					polys = reach_region->flowpipe_intersectionSequential(gaurd_polytope, lp_solver_type_choosen);
				} else {	//the guard polytope is universal
					//std::cout<<"Guard is universal\n";
					polytope::ptr templated_hull_flowpipe; //bounding_box Polytope
					templated_hull_flowpipe = convertBounding_Box(reach_region);
					polys.push_back(templated_hull_flowpipe);
				}
//std::cout<<"\nTesting 2 g\n";
				//Todo to make is even procedure with Sequential procedure.... so intersection is done first and then decide to skip this loc
				if ((locName.compare("BAD") == 0) || (locName.compare("GOOD") == 0)
						|| (locName.compare("FINAL") == 0) || (locName.compare("UNSAFE") == 0))
					continue; //do not push into the waitingList

				current_assignment = (*t)->getAssignT();
				// *** interesected_polyhedra included with invariant_directions also ******
				int destination_locID = (*t)->getDestination_Location_Id();
				ds.insert_element(destination_locID);
				for (std::list<polytope::ptr>::iterator i = polys.begin(); i != polys.end(); i++) {

					//intersectedRegion = (*i)->getTemplate_approx(lp_solver_type_choosen);
					intersectedRegion = (*i);

					//Returns a single over-approximated polytope from the list of intersected polytopes
					//	GeneratePolytopePlotter(intersectedRegion);
					polytope::ptr newShiftedPolytope, newPolytope; //created an object here
					newPolytope = intersectedRegion->GetPolytope_Intersection(gaurd_polytope);
					//Returns the intersected region as a single newpolytope.
//std::cout<<"\nTesting 2 h\n";
					math::matrix<double> test(current_assignment.Map.size1(),current_assignment.Map.size2());

					if (current_assignment.Map.inverse(test))	//invertible?
					{
						//std::cout<<"Exact Post Assignment\n";
						newShiftedPolytope = post_assign_exact(newPolytope,
								current_assignment.Map, current_assignment.b);
					} else {
						//std::cout<<"Approximate Post Assignment\n";
						newShiftedPolytope = post_assign_approx_deterministic(newPolytope,
								current_assignment.Map, current_assignment.b, reach_parameters.Directions,lp_solver_type_choosen);
					}
					// @Rajarshi: the newShifted satisfy the destination location invariant
					newShiftedPolytope = newShiftedPolytope->GetPolytope_Intersection(H.getLocation(destination_locID)->getInvariant());

					initial_state::ptr newState = initial_state::ptr(new initial_state(destination_locID, newShiftedPolytope));
					newState->setTransitionId(transition_id); // keeps track of the transition_ID
					newState->setParentPtrSymbolicState(S);
					pw_list.WaitingList_insert(newState);
					queue.push_back(BreadthLevel); //insert at REAR first Location
				}
			} //end of multiple transaction
		}
	} //end of while loop checking waiting_list != empty

	cout << "\n ***************************************************************************\n";
	cout << "\nMaximum Iterations Completed = " << num_flowpipe_computed << "\n";
	cout << "\n ***************************************************************************\n";

	return Reachability_Region;
}



void reachability::sequentialReachSelection(unsigned int NewTotalIteration, location::ptr current_location,
		polytope::ptr continuous_initial_polytope,
		template_polyhedra::ptr& reach_region) {

	if (Algorithm_Type == SEQ_SF) { //Continuous Sequential Algorithm
		/*boost::timer::cpu_timer AllReach_time;
		AllReach_time.start();*/
		//std::cout << "\nFlowpipe" << std::endl;
		reach_region = reachabilitySequential(NewTotalIteration, current_location->getSystem_Dynamics(),
				continuous_initial_polytope, reach_parameters, current_location->getInvariant(),
				current_location->isInvariantExists(), lp_solver_type_choosen);
		//std::cout << "\nFlowpipe computed" << std::endl;
		/*AllReach_time.stop();
		double wall_clock1;
		wall_clock1 = AllReach_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
		double return_Time1 = wall_clock1 / (double) 1000;
		std::cout << "\nFlowpipe:Time:Wall(Seconds) = " << return_Time1 << std::endl;*/
	}

	if (Algorithm_Type == PAR_SF) {
		//Parallel implementation using OMP parallel			//	double wall_timer = omp_get_wtime();
		//	cout << "\nRunning Parallel Using OMP Thread\n";
		/*boost::timer::cpu_timer AllReach_time;
		AllReach_time.start();*/
		reach_region = reachabilityParallel(NewTotalIteration,
				current_location->getSystem_Dynamics(),
				continuous_initial_polytope, reach_parameters,
				current_location->getInvariant(),
				current_location->isInvariantExists(), lp_solver_type_choosen);
		/*AllReach_time.stop();
		double wall_clock1;
		wall_clock1 = AllReach_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
		double return_Time1 = wall_clock1 / (double) 1000;
		std::cout << "\nFlowpipe:Time:Wall(Seconds) = " << return_Time1 << std::endl;*/
		//	std::cout << "Parallel Done\n";
		//	std::cout << "Time seen from mop wall timer: "<< omp_get_wtime() - wall_timer << std::endl;
	}


	if (Algorithm_Type == TIME_SLICE) { //Continuous Parallel Algorithm parallelizing the Iterations :: to be debugged (compute initial polytope(s))
		cout << "\nRunning Parallel-over-Iterations(PARTITIONS/Time-Sliced) Using OMP Thread\n";
		bool invertible;
		if (!current_location->getSystem_Dynamics().isEmptyMatrixA) {
			math::matrix<double> A_inv(current_location->getSystem_Dynamics().MatrixA.size1(),
					current_location->getSystem_Dynamics().MatrixA.size2());
			invertible = current_location->getSystem_Dynamics().MatrixA.inverse(A_inv); //size of A_inv must be declared else error
			if (!invertible)
				std::cout << "\nDynamics Matrix A is not invertible\n"; //later can give option to stop or select different algorithm
			else
				reach_parameters.A_inv = A_inv;
			int NCores = Total_Partition; //Number of Partitions (number of threads)
			reach_region = reachParallelExplore(NewTotalIteration, current_location->getSystem_Dynamics(),
					continuous_initial_polytope, reach_parameters, current_location->getInvariant(),
					current_location->isInvariantExists(), NCores, TIME_SLICE, lp_solver_type_choosen);
			std::cout << "Finished computing reachable states\n";
		} else {
			std::cout << "\nFlow Dynamics Matrix A does not exists!!!\n";
			template_polyhedra::ptr poly_emptyp;
			reach_region = poly_emptyp; //returning an empty reach_region for this location

		}
	}

	/*	if (Algorithm_Type == GPU_SF) { //computing all support function in GPU
			cout << "\nRunning GPU Sequential\n";
			boost::timer::cpu_timer AllReachGPU_time;
			AllReachGPU_time.start();
			reachabilitySequential_GPU(NewTotalIteration, current_location->getSystem_Dynamics(), continuous_initial_polytope, reach_parameters,
					current_location->getInvariant(), current_location->isInvariantExists(),
					lp_solver_type_choosen, number_of_streams, Solver_GLPK_Gurobi_GPU, reach_region);
			std::cout << "Out from GPU_reach\n";
			AllReachGPU_time.stop();
			double wall_clock1;
			wall_clock1 = AllReachGPU_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
			double return_Time1 = wall_clock1 / (double) 1000;
			std::cout << "\nAllReach_time: Boost Time:Wall(Seconds) = " << return_Time1 << std::endl;

		}*/
}

//bound is the maximum number of transitions or jumps permitted.
//reach_parameters includes the different parameters needed in the computation of reachability.
//I is the initial symbolic state

std::list<symbolic_states::ptr> reachability::computeParallelBFSReach(
		std::list<abstractCE::ptr>& ce_candidates) {

	std::list < symbolic_states::ptr > Reachability_Region;
	//	template_polyhedra::ptr reach_region;

	pwlist pw_list;
	for (std::list<initial_state::ptr>::iterator i=I.begin();i!=I.end();i++)
		pw_list.WaitingList_insert(*(i));
		//pw_list.WaitingList_insert(I);
	int number_times = 0;
	unsigned int iter_max = 1;

	bool stop_loop = false;

	bool saftey_violated = false;

	while (!pw_list.isEmpty_WaitingList()) {

		//To avoid write-contention for reach_region, Vector/List of reach_region_list created for each threads
		// write in its respective index. So need a unique ID for each thread which can be obtained from
		// the size of the PWList at each iteration

		unsigned int count = pw_list.getWaitingListSize(); //get the size of PWList
		std::vector < symbolic_states::ptr > S(count);

		//Create a sublist of initial_state and work with it inside the parallel region(each thread accesses uniquely)

		vector < initial_state::ptr > list_U(count); //SubList for parallel

		for (int i = 0; i < count; i++) {
			list_U[i] = pw_list.WaitingList_delete_front();
			pw_list.PassedList_insert(list_U[i]);
		}
		//All initial_state have been deleted
		// ********************************* BFS Starts **********************************************************

		omp_set_nested(1); //enable nested parallelism
#pragma omp parallel for // num_threads(2)
		for (unsigned int id = 0; id < count; id++) {
			initial_state::ptr U; //local
			U = list_U[id]; //independent symbolic state to work with
			discrete_set discrete_state; //local
			polytope::ptr continuous_initial_polytope; //local
			ReachabilityParameters reach_parameter_local; //local

			int location_id = U->getLocationId();
			discrete_state.insert_element(location_id);
			continuous_initial_polytope = U->getInitialSet();
			reach_parameter_local = reach_parameters;
			reach_parameter_local.X0 = continuous_initial_polytope; //	cout<<"\nInside for Loop";
			//create an instance of Symbolic_states S
			S[id] = symbolic_states::ptr(new symbolic_states());

			S[id]->setDiscreteSet(discrete_state);
			S[id]->setParentPtrSymbolicState(U->getParentPtrSymbolicState()); //keeps track of parent pointer to symbolic_states
			S[id]->setTransitionId(U->getTransitionId()); //keeps track of originating transition_ID


			location::ptr current_location;
			current_location = H.getLocation(location_id);
			string name = current_location->getName();
			if ((name.compare("GOOD") == 0) || (name.compare("BAD") == 0)
					|| (name.compare("UNSAFE") == 0)
					|| (name.compare("FINAL") == 0))
				continue; //do not compute the continuous reachability algorithm

			// ******************* Computing Parameters *******************************
			double result_alfa = compute_alfa(reach_parameter_local.time_step,
					current_location->getSystem_Dynamics(),
					continuous_initial_polytope, lp_solver_type_choosen); //2 glpk object created here

			double result_beta = compute_beta(
					current_location->getSystem_Dynamics(),
					reach_parameter_local.time_step, lp_solver_type_choosen); // NO glpk object created here

			reach_parameter_local.result_alfa = result_alfa;
			reach_parameter_local.result_beta = result_beta;

			math::matrix<double> phi_matrix, phi_trans;
			if (!current_location->getSystem_Dynamics().isEmptyMatrixA) { //if A not Empty
				current_location->getSystem_Dynamics().MatrixA.matrix_exponentiation(
						phi_matrix, reach_parameter_local.time_step);
				phi_matrix.transpose(phi_trans);
				reach_parameter_local.phi_trans = phi_trans;
			}
			math::matrix<double> B_trans;
			if (!current_location->getSystem_Dynamics().isEmptyMatrixB) { //if B not Empty
				current_location->getSystem_Dynamics().MatrixB.transpose(
						B_trans);
				reach_parameter_local.B_trans = B_trans;
			}
			// ******************* Computing Parameters Done *******************************
			unsigned int NewTotalIteration;
			if (current_location->isInvariantExists()) {
				InvariantBoundaryCheck(current_location->getSystem_Dynamics(), continuous_initial_polytope,
						reach_parameter_local, current_location->getInvariant(), lp_solver_type_choosen, NewTotalIteration);
				std::cout << "NewTotalIteration = " << NewTotalIteration << std::endl;
			}
			//  ********************* FlowPipe or Reach Computation *************************
			parallelReachSelection(NewTotalIteration, current_location, continuous_initial_polytope, reach_parameter_local, S, id);
			// Returns the Flow_Pipe in reach_region_list[id]
			//  ********************* FlowPipe or Reach Computation Done ********************

			//  ************** Check to see if Computed FlowPipe is Empty  **********
			template_polyhedra::ptr t_poly = S[id]->getContinuousSetptr();

			if (t_poly->getTotalIterations() != 0 && number_times < bound) { //computed reach_region is empty && optimize computation

				for (std::list<transition::ptr>::iterator t =
						current_location->getOut_Going_Transitions().begin();
						t != current_location->getOut_Going_Transitions().end();
						t++) { // get each destination_location_id and push into the pwl.waiting_list
					int transition_id = (*t)->getTransitionId();
					if (transition_id == -1) { //Indicates empty transition or no transition exists
						break; //out from transition for-loop as there is no transition for this location
					}
					location::ptr current_destination;
					Assign current_assignment;
					polytope::ptr gaurd_polytope;
					std::list < template_polyhedra::ptr > intersected_polyhedra;
					polytope::ptr intersectedRegion; //created two objects here
					discrete_set ds;
					current_destination = H.getLocation(
							(*t)->getDestination_Location_Id());
					string locName = current_destination->getName();

					gaurd_polytope = (*t)->getGaurd();
					current_assignment = (*t)->getAssignT();

					boost::timer::cpu_timer t100;
					//this intersected_polyhedra will have invariant direction added in it
					string trans_name = (*t)->getLabel();
					//	intersected_polyhedra = t_poly->polys_intersectionParallel(gaurd_polytope, lp_solver_type_choosen); //, intersection_start_point);
					t100.start();
					intersected_polyhedra =
							t_poly->polys_intersectionSequential(gaurd_polytope,
									lp_solver_type_choosen); //, intersection_start_point);
					t100.stop();
					if (intersected_polyhedra.size() > 0) { //there is intersection so new symbolic state will be inserted into the waitingList
						//std::cout << "Intersected = " < intersected_polyhedra.size() << std::endl;
#pragma omp critical
						{
							iter_max++; //
						}
					}
					//todo Handle this later as In SpaceEx model we did not specified BAD or GOOD
					if ((locName.compare("BAD") == 0)
							|| (locName.compare("GOOD") == 0)
							|| (locName.compare("FINAL") == 0)
							|| (locName.compare("UNSAFE") == 0)) {
						continue; //do not push into the waitingList
					}
					//			std::cout << "Before calling getTemplate_approx\n";
					double clock100;
					clock100 = t100.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
					double return100 = clock100 / (double) 1000;
					std::cout << "\nIntersection Time:Wall(Seconds) = "
							<< return100 << std::endl;

					int destination_locID = (*t)->getDestination_Location_Id();
					ds.insert_element(destination_locID);
					for (std::list<template_polyhedra::ptr>::iterator i =
							intersected_polyhedra.begin();
							i != intersected_polyhedra.end(); i++) {

						intersectedRegion = (*i)->getTemplate_approx(
								lp_solver_type_choosen);
						//Returns a single over-approximated polytope from the list of intersected polytopes
						polytope::ptr newShiftedPolytope, newPolytope; //created an object here
						newPolytope =
								intersectedRegion->GetPolytope_Intersection(
										gaurd_polytope); //Retuns only the intersected region as a single newpolytope. ****** with added directions

						newShiftedPolytope = post_assign_exact(newPolytope,
								current_assignment.Map, current_assignment.b); //initial_polytope_I = post_assign_exact(newPolytope, R, w);

						initial_state::ptr newState = initial_state::ptr(
								new initial_state(destination_locID,
										newShiftedPolytope));
						newState->setTransitionId(transition_id); // keeps track of the transition_ID
						newState->setParentPtrSymbolicState(S[id]);

#pragma omp critical
						{
							/*
							 * Todo:: we insert all the symbolic_states into the pwlist for computing FlowPipe iff 1) and 2) holds
							 * Step 1) the new current_destination.getLocID() is NOT in the Passed List and
							 * Step 2) the "New Initial Polytope" or "the newShiftedPolytope" is NOT contained in the FlowPipe of LocationID of step 1
							 *  as the FlowPipe in the same location will be same if the newShiftedPolytope is IN FlowPipe
							 *  But if Step 1 holds and Step 2 does not then it will be inserted in pwlist even if the LocationID is in Passed List
							 */
							pw_list.WaitingList_insert(newState); //RACE CONDITION HERE

						}
					}
				} //end of multiple transaction
				  //Here we have again populated the pwlist for next-round's parallel process
			} // End-if  ******* Check for Empty FlowPipe Done *********
		} //END of parallel FOR-LOOP

		//:: Can be optimized if we can count number_times inside the parallel loop per breadth then we can avoid transaction and intersection
		//:: computation for next transition if number_times exceeds bound ....
		number_times++; //One Level or one Breadth Search over

		// ************************* BFS Ends *************************************

		//Creating a list of objects of "Reachability Set"/Symbolic_states
		bool foundUnSafe = false;
		for (unsigned int index = 0; index < count; index++) {

			if (S[index]->getContinuousSetptr()->getTotalIterations() != 0) //computed reach_region is NOT empty
				Reachability_Region.push_back(S[index]);
			//  ******************************** Safety Verification section ********************************
			saftey_violated = safetyVerify(S[index], Reachability_Region, ce_candidates);
			if (saftey_violated) {
				foundUnSafe = true; //have to insert all flowpipe of same breadth even if foundUnSafe=true for some flowpipe, So not breaking
			}
			//  ******************************** Safety Verification section ********************************
		} //end-for pushing all computed flowpipe
//		if (foundUnSafe) {
//			break; //no need to compute rest of the locations
//		}
		if (number_times > bound) //check to see how many jumps have been made(i.e., number of discrete transitions made)
			break;
	} //end of while loop checking waiting_list != empty
	cout << "\n ***************************************************************************\n";
	cout << "\nMaximum Iterations Completed = " << iter_max << "\n";
	cout << "\n ***************************************************************************\n";

	return Reachability_Region;

}

//Lock Avoidance:: Parallel Breadth First Search for Discrete Jumps
//separate Read and Write Queue (pwlist.WaitingList)
std::list<symbolic_states::ptr> reachability::computeParallelBFSReachLockAvoid(std::list<abstractCE::ptr>& ce_candidates) {

	std::list < symbolic_states::ptr > Reachability_Region;
	//	template_polyhedra::ptr reach_region;

	//Shared Read and Write Queue(pwList)
	int t = 0; //0 for Read and 1 for Write
	std::vector < std::vector<pwlist::ptr> > Qpw_list(2); // QpwList[0] for read and QpwList[1] for write
	//pwlist pw_list;	//Shared Queue(pwList)
	//cout << "Test 1\n";
	unsigned int nos_initial_states=1;
	nos_initial_states=I.size();	//retrieve the size of the number of initial states supplied by the user/model

	Qpw_list[t].resize(nos_initial_states); //resize for the first symbolic_state
	//cout << "Test 2\n";
	int initialState_index=0;// first initial state
	for (std::list<initial_state::ptr>::iterator i=I.begin();i!=I.end();i++){

		Qpw_list[t][initialState_index] = pwlist::ptr(new pwlist()); //have to instantiate it
		Qpw_list[t][initialState_index]->WaitingList_insert(*(i));

		initialState_index++; //next initial state
	}
	//cout << "Test 3\n";
//Since Qpw_list is a temporary working list is it switches between Read and Write so we can not maintain the passedList in it.
	pwlist::ptr allPassedList; //so we create a permanent pwlist for storing only the passedList;
	allPassedList = pwlist::ptr(new pwlist()); //have to instantiate it
	int number_times = 0;
	unsigned int iter_max = 1;

	bool stop_loop = false;

	bool saftey_violated = false;
	//cout << "Test 4\n";
	while (!isEmpty_Qpw_list(Qpw_list[t]) && (number_times <= bound)) {
		cout<<"Breadth - Level === "<<number_times<<"\n";
		//To avoid write-contention for reach_region, Vector/List of reach_region_list created for each threads
		// write in its respective index. So need a unique ID for each thread which can be obtained from
		// the size of the PWList at each iteration
		//	cout << "Test 5\n";
		unsigned int count = getSize_Qpw_list(Qpw_list[t]); //get the size of PWList
		std::vector < symbolic_states::ptr > S(count);

		std::vector <location::ptr> list_currLocation(count);	//Data structure required to separate PostD from PostC
		std::vector <int> list_invBounaryValue(count);	//Data structure required to separate invCheck from PostC
		std::vector <LoadBalanceData> SymDataStruct(count);

		//Create a sublist of initial_state and work with it inside the parallel region(each thread accesses uniquely)
		//vector<symbolic_states> list_U(count); //SubList for parallel
		vector < initial_state::ptr > list_U(count); //SubList for parallel

		list_U = getAllpw_list(Qpw_list, t, count, allPassedList); //All initial_state have been deleted

		Qpw_list[t].resize(0);
		cout<<"Qpw_list[t].size() = "<<Qpw_list[t].size()<<std::endl;

		Qpw_list[1 - t].resize(count); //resize to accommodate
		for (int i = 0; i < count; i++) {
			Qpw_list[1 - t][i] = pwlist::ptr(new pwlist()); //have to instantiate it
		}
		// ********************************* BFS Starts **********************************************************
		//	cout << "Test 7\n";
		omp_set_dynamic(0);	//handles dynamic adjustment of the number of threads within a team
		omp_set_nested(1); //enable nested parallelism
		//omp_set_num_threads(10);
		//omp_set_max_active_levels(2);

		boost::timer::cpu_timer t705;
		t705.start();
#pragma omp parallel for  num_threads(count)
		for (unsigned int id = 0; id < count; id++) {
			boost::timer::cpu_timer t702;
			t702.start();
			if (id==0){
				std::cout<<"\nMax Thread in Outer Level = "<< omp_get_num_threads();

			}
			//there will be different current_location, continuous_initial_polytope, reach_parameters
			initial_state::ptr U; //local
			//	symbolic_states::ptr S = symbolic_states::ptr(new symbolic_states());
			U = list_U[id]; //independent symbolic state to work with
			discrete_set discrete_state; //local
			polytope::ptr continuous_initial_polytope; //local
			ReachabilityParameters reach_parameter_local; //local

			int location_id = U->getLocationId();
			discrete_state.insert_element(location_id);
			continuous_initial_polytope = U->getInitialSet();
			reach_parameter_local = reach_parameters;
			reach_parameter_local.X0 = continuous_initial_polytope; //	cout<<"\nInside for Loop";
			//create an instance of Symbolic_states S
			S[id] = symbolic_states::ptr(new symbolic_states());

			S[id]->setDiscreteSet(discrete_state);
			S[id]->setParentPtrSymbolicState(U->getParentPtrSymbolicState()); //keeps track of parent pointer to symbolic_states
			S[id]->setTransitionId(U->getTransitionId()); //keeps track of originating transition_ID

			location::ptr current_location;
			current_location = H.getLocation(location_id);
			string name = current_location->getName();
			if ((name.compare("GOOD") == 0) || (name.compare("BAD") == 0)
					|| (name.compare("UNSAFE") == 0) || (name.compare("FINAL") == 0))
				continue; //do not compute the continuous reachability algorithm

			// ******************* Computing Parameters *******************************
			//current_location:: parameters alfa, beta and phi_trans have to be re-computed
			//cout<<"\nBefore Compute Alfa";
			double result_alfa = compute_alfa(reach_parameter_local.time_step,
					current_location->getSystem_Dynamics(), continuous_initial_polytope, lp_solver_type_choosen); //2 glpk object created here
			//cout<<"\nCompute Alfa Done";
			double result_beta = compute_beta(current_location->getSystem_Dynamics(), reach_parameter_local.time_step, lp_solver_type_choosen); // NO glpk object created here
			//			cout<<"\nCompute Beta Done";
			reach_parameter_local.result_alfa = result_alfa;
			reach_parameter_local.result_beta = result_beta;

			math::matrix<double> phi_matrix, phi_trans;
			if (!current_location->getSystem_Dynamics().isEmptyMatrixA) { //if A not Empty
				current_location->getSystem_Dynamics().MatrixA.matrix_exponentiation(phi_matrix, reach_parameter_local.time_step);
				phi_matrix.transpose(phi_trans);
				reach_parameter_local.phi_trans = phi_trans;
			}
			math::matrix<double> B_trans;
			if (!current_location->getSystem_Dynamics().isEmptyMatrixB) { //if B not Empty
				current_location->getSystem_Dynamics().MatrixB.transpose(B_trans);
				reach_parameter_local.B_trans = B_trans;
			}

			t702.stop();
			double clock702;
			clock702 = t702.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
			double return702 = clock702 / (double) 1000;
			std::cout << "\nReachParameters computation Time:Wall(Seconds) = "<< return702 << std::endl;

			SymDataStruct[id].X0 = continuous_initial_polytope;
			SymDataStruct[id].current_location = current_location;
			SymDataStruct[id].reach_param = reach_parameter_local;
			// ******************* Computing Parameters Done *******************************
			boost::timer::cpu_timer boundCheck;
			boundCheck.start();
			unsigned int NewTotalIteration;
			if (current_location->isInvariantExists()){
				InvariantBoundaryCheck(SymDataStruct[id].current_location->getSystem_Dynamics(), continuous_initial_polytope,
						SymDataStruct[id].reach_param, SymDataStruct[id].current_location->getInvariant(), lp_solver_type_choosen, NewTotalIteration);
				std::cout << "NewTotalIteration = " << NewTotalIteration << std::endl;
				//list_invBounaryValue[id] = NewTotalIteration;
			}
			boundCheck.stop();
			double wall_clockboundcheck;
			wall_clockboundcheck = boundCheck.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
			double return_TimeBC = wall_clockboundcheck / (double) 1000;
			std::cout<< "\nBoundary Check for Iterations Number Time taken:Wall  (in Seconds) = " << return_TimeBC << std::endl;

			SymDataStruct[id].newIteration = NewTotalIteration;
		}


		boost::timer::cpu_timer t3;
		double cpu_usage;
		init_cpu_usage();
		t3.start();
#pragma omp parallel for  num_threads(count)
		for (unsigned int id = 0; id < count; id++) {
		/*if (shm_NewTotalIteration <= 1) {
				template_polyhedra::ptr poly_emptyp;
				return poly_emptyp;
			}*/
			//  ********************* FlowPipe or Reach Computation *************************

			parallelReachSelection(SymDataStruct[id].newIteration, SymDataStruct[id].current_location,
					SymDataStruct[id].X0, SymDataStruct[id].reach_param, S, id);

			//parallelReachSelection(NewTotalIteration, current_location, continuous_initial_polytope, reach_parameter_local, S, id);
			// Returns the Flow_Pipe in reach_region_list[id]
		//	list_currLocation[id] = current_location;

			/*double res= return702 + return3;
			std::cout << "****************************************************************************\n";
			std::cout << "\nEach Computation (ReachParameter + Flow-pipe) Time:Wall(Seconds) = " << res << std::endl;
			std::cout << "****************************************************************************\n";*/
		}
		t3.stop();
		cpu_usage = getCurrent_ProcessCPU_usage();
		std::cout << "\nCPU Usage:(%) = " << cpu_usage<< std::endl;
		double clock3;
		clock3 = t3.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
		double return3 = clock3 / (double) 1000;
		std::cout << "\nFlow-pipe Computed Time:Wall(Seconds) = " << return3 << std::endl;

		//  ********************* POST_C computation Done ********************
		t705.stop();
		double clock705;
		clock705 = t705.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
		double return705 = clock705 / (double) 1000;
		std::cout << "\nFlow-pipe Computation Done Time:Wall(Seconds) = " << return705 << std::endl;
		std::cout << "****************************************************************************\n";


		/*number_times++; //One Level or one Breadth Search over
		if (number_times > bound) {
			break; //check to see how many jumps have been made(i.e., number of discrete transitions made)
		}*/
//  ************************************** POST_D computation Begins **********************************************************

		omp_set_dynamic(0);	//handles dynamic adjustment of the number of threads within a team
		omp_set_nested(1); //enable nested parallelism
		boost::timer::cpu_timer t72;
		t72.start();
#pragma omp parallel for // num_threads(2)
		for (unsigned int id = 0; id < count; id++) {

			location::ptr current_location;
			//current_location = list_currLocation[id];
			current_location = SymDataStruct[id].current_location;

			/*Removed from here as kept in the WHILE-loop condition
			 * if (number_times > bound) //check to see how many jumps have been made(i.e., number of discrete transitions made)
				break;	//try to remove from here as this may not work in #pragma omp parallel for
			*/

			//  ************** Check to see if Computed FlowPipe is Empty  **********
			template_polyhedra::ptr t_poly = S[id]->getContinuousSetptr();

			if (t_poly->getTotalIterations() != 0 && number_times < bound) { //computed reach_region is empty && optimize computation
			//cout << "\nLoc ID = " << current_location.getLocId() << " Location Name = " << name << "\n";
				for (std::list<transition::ptr>::iterator trans = current_location->getOut_Going_Transitions().begin();
						trans != current_location->getOut_Going_Transitions().end(); trans++) { // get each destination_location_id and push into the pwl.waiting_list
					int transition_id = (*trans)->getTransitionId();
					if (transition_id == -1) { //Indicates empty transition or no transition exists
						break; //out from transition for-loop as there is no transition for this location
					}
					location::ptr current_destination;
					Assign current_assignment;
					polytope::ptr gaurd_polytope;
					std::list < template_polyhedra::ptr > intersected_polyhedra;
					polytope::ptr intersectedRegion; //created two objects here
					discrete_set ds;
					current_destination = H.getLocation((*trans)->getDestination_Location_Id());
					string locName = current_destination->getName();
					//	cout << "\nNext Loc ID = " << current_destination.getLocId() << " Location Name = " << locName << "\n";

					gaurd_polytope = (*trans)->getGaurd();
					current_assignment = (*trans)->getAssignT();

					boost::timer::cpu_timer t100;
					//this intersected_polyhedra will have invariant direction added in it
					string trans_name = (*trans)->getLabel();
					t100.start();
					intersected_polyhedra = t_poly->polys_intersectionParallel(gaurd_polytope, lp_solver_type_choosen); //, intersection_start_point);
					//intersected_polyhedra = t_poly->polys_intersectionSequential(gaurd_polytope, lp_solver_type_choosen); //, intersection_start_point);
					t100.stop();
					if (intersected_polyhedra.size() > 0) { //there is intersection so new symbolic state will be inserted into the waitingList
						//std::cout << "Intersected = " < intersected_polyhedra.size() << std::endl;
#pragma omp critical
						{
							iter_max++; //
						}
					}
					// Handle this later as In SpaceEx model we did not specified BAD or GOOD
					if ((locName.compare("BAD") == 0) || (locName.compare("GOOD") == 0)
							|| (locName.compare("FINAL") == 0) || (locName.compare("UNSAFE") == 0)) {
						continue; //do not push into the waitingList
					}
					//			std::cout << "Before calling getTemplate_approx\n";
					double clock100;
					clock100 = t100.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
					double return100 = clock100 / (double) 1000;
					std::cout << "\nIntersection Time:Wall(Seconds) = " << return100 << std::endl;
					int destination_locID = (*trans)->getDestination_Location_Id();
					ds.insert_element(destination_locID);
					for (std::list<template_polyhedra::ptr>::iterator it = intersected_polyhedra.begin(); it != intersected_polyhedra.end(); it++) {
						//cout << "\nNumber of Intersections #1\n";
						intersectedRegion = (*it)->getTemplate_approx(lp_solver_type_choosen);
						//Returns a single over-approximated polytope from the list of intersected polytopes
						polytope::ptr newShiftedPolytope, newPolytope; //created an object here
						newPolytope = intersectedRegion->GetPolytope_Intersection(
										gaurd_polytope); //Retuns only the intersected region as a single newpolytope. ****** with added directions
						//std::cout << "Before calling post_assign_exact\n";
						newShiftedPolytope = post_assign_exact(newPolytope, current_assignment.Map, current_assignment.b); //initial_polytope_I = post_assign_exact(newPolytope, R, w);
						initial_state::ptr newState = initial_state::ptr(new initial_state(destination_locID,newShiftedPolytope));
						newState->setTransitionId(transition_id); // keeps track of the transition_ID
						newState->setParentPtrSymbolicState(S[id]);
//#pragma omp critical
//					{
						/*
						 * Todo:: we insert all the symbolic_states into the pwlist for computing FlowPipe iff 1) and 2) holds
						 * Step 1) the new current_destination.getLocID() is NOT in the Passed List and
						 * Step 2) the "New Initial Polytope" or "the newShiftedPolytope" is NOT contained in the FlowPipe of LocationID of step 1
						 *  as the FlowPipe in the same location will be same if the newShiftedPolytope is IN FlowPipe
						 *  But if Step 1 holds and Step 2 does not then it will be inserted in pwlist even if the LocationID is in Passed List
						 */
						Qpw_list[1 - t][id]->WaitingList_insert(newState); //RACE CONDITION HERE

	//				}
					} //end of multiple intersected region with guard
					  //cout<<"Size = "<< pwlist.getWaitingList().size()<<endl;
				} //end of multiple transaction
				  //Here we have again populated the pwlist for next-round's parallel process
			} // End-if  ******* Check for Empty FlowPipe Done *********
		} //END of parallel FOR-LOOP
		//  ************************************** POST_D computation Ends **********************************************************
		t72.stop();
		double clock72;
		clock72 = t72.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
		double return72 = clock72 / (double) 1000;
		std::cout << "\nDiscrete Post_D computation Time:Wall(Seconds) = " << return72 << std::endl;

//  ********** Safety Verification Should be after Flowpipe computation but due to omp-parallel-for NOT POSSIBLE to use break*************
		//Creating a list of objects of "Reachability Set"/Symbolic_states
		bool foundUnSafe = false;
		//Having is sequential to avoid critical section moreover this time is very very small
		for (unsigned int index = 0; index < count; index++) {
			if (S[index]->getContinuousSetptr()->getTotalIterations() != 0) //computed reach_region is NOT empty
				Reachability_Region.push_back(S[index]);
		}
#pragma omp parallel for
		for (unsigned int index = 0; index < count; index++) {
			//  ******************************** Safety Verification section ********************************
			saftey_violated = safetyVerify(S[index], Reachability_Region, ce_candidates);
#pragma omp critical
		{
			if (saftey_violated) {
				foundUnSafe = true; //have to insert all flowpipe of same breadth even if foundUnSafe=true for some flowpipe, So not breaking
			}
		}
			//  ******************************** Safety Verification section ********************************
		} //end-for pushing all computed flowpipe
//		if (foundUnSafe) {
//			break; // removed from inner-loop does not work in #pragma omp parallel for
//		}

		t = 1 - t; //Switching Read/Write options for Qpw_list[1-t]
//:: Can be optimize if we can count number_times inside the parallel loop per breadth then we can avoid transaction and intersection
//:: computation for next transition if number_times exceeds bound ....
		number_times++; //One Level or one Breadth Search over
		//cout << "\nnumber_times = " << number_times << "  Bound = " << bound << "\n";
		// ************************* BFS Ends *************************************
	} //end of while loop checking waiting_list != empty
	cout
			<< "\n ***************************************************************************\n";
	cout << "\nMaximum Iterations Completed = " << iter_max << "\n";
	cout
			<< "\n ***************************************************************************\n";

	return Reachability_Region;
}

/*** TODO: Have to optimize invariant_boundary_check() for support function computation ***/

void reachability::parallelReachSelection(unsigned int NewTotalIteration, location::ptr current_location,
		polytope::ptr continuous_initial_polytope,
		ReachabilityParameters& reach_parameters,
		std::vector<symbolic_states::ptr>& S, unsigned int id) {

	template_polyhedra::ptr reach_region;

	if (Algorithm_Type == SEQ_SF) { //Continuous Sequential Algorithm
		boost::timer::cpu_timer AllReach_time;
		AllReach_time.start();

		reach_region = reachabilitySequential(NewTotalIteration,
				current_location->getSystem_Dynamics(),
				continuous_initial_polytope, reach_parameters,
				current_location->getInvariant(),
				current_location->isInvariantExists(), lp_solver_type_choosen);
		S[id]->setContinuousSetptr(reach_region);

		AllReach_time.stop();
		double wall_clock1;
		wall_clock1 = AllReach_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
		double return_Time1 = wall_clock1 / (double) 1000;
		std::cout << "\nFlowPipe Time:Wall(Seconds) = " << return_Time1
				<< std::endl;
	}

	if (Algorithm_Type == PAR_SF) {
		boost::timer::cpu_timer AllReach_time;
		AllReach_time.start();

		reach_region = reachabilityParallel(NewTotalIteration,
				current_location->getSystem_Dynamics(),
				continuous_initial_polytope, reach_parameters,
				current_location->getInvariant(),
				current_location->isInvariantExists(), lp_solver_type_choosen);
		AllReach_time.stop();
		S[id]->setContinuousSetptr(reach_region);
		double wall_clock1;
		wall_clock1 = AllReach_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
		double return_Time1 = wall_clock1 / (double) 1000;
		std::cout << "\nFlowpipe Time:Wall(Seconds) = " << return_Time1
				<< std::endl;
		//	std::cout << "Time seen from mop wall timer: "<< omp_get_wtime() - wall_timer << std::endl;
	}

/*	if (Algorithm_Type == GPU_SF) { //computing all support function in GPU
		cout << "\nRunning GPU Sequential\n";
		boost::timer::cpu_timer AllReachGPU_time;
		AllReachGPU_time.start();
		reachabilitySequential_GPU(NewTotalIteration, current_location->getSystem_Dynamics(),
				continuous_initial_polytope, reach_parameters,
				current_location->getInvariant(),
				current_location->isInvariantExists(), lp_solver_type_choosen,
				number_of_streams, Solver_GLPK_Gurobi_GPU, reach_region);
		S[id]->setContinuousSetptr(reach_region);

		AllReachGPU_time.stop();
		double wall_clock1;
		wall_clock1 = AllReachGPU_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
		double return_Time1 = wall_clock1 / (double) 1000;
		std::cout << "\nAllReach_time: Boost Time:Wall(Seconds) = "
				<< return_Time1 << std::endl;
	}*/

	if (Algorithm_Type == TIME_SLICE) { //Continuous Parallel Algorithm parallelizing the Iterations :: to be debugged (compute initial polytope(s))
		cout
				<< "\nRunning Parallel-over-Iterations(PARTITIONS/Time-Sliced) Using OMP Thread\n";
		bool invertible;
		if (!current_location->getSystem_Dynamics().isEmptyMatrixA) {
			math::matrix<double> A_inv(current_location->getSystem_Dynamics().MatrixA.size1(),
					current_location->getSystem_Dynamics().MatrixA.size2());
			invertible = current_location->getSystem_Dynamics().MatrixA.inverse(A_inv); //size of A_inv must be declared else error
			if (!invertible)
				std::cout << "\nDynamics Matrix A is not invertible\n"; //later can give option to stop or select different algorithm
			else
				reach_parameters.A_inv = A_inv;
			int NCores = Total_Partition; //Number of Partitions (number of threads)
			reach_region = reachParallelExplore(NewTotalIteration, current_location->getSystem_Dynamics(),
					continuous_initial_polytope, reach_parameters, current_location->getInvariant(),
					current_location->isInvariantExists(), NCores, TIME_SLICE, lp_solver_type_choosen);
			S[id]->setContinuousSetptr(reach_region);
			//		std::cout << "Finished computing reachable states\n";
		} else {
			std::cout << "\nFlow Dynamics Matrix A does not exists!!!\n";
			template_polyhedra::ptr poly_emptyp;
			reach_region = poly_emptyp; //returning an empty reach_region for this location
			S[id]->setContinuousSetptr(reach_region);
		}
	}
}



bool reachability::safetyVerify(symbolic_states::ptr& computedSymStates,
		std::list<symbolic_states::ptr>& Reachability_Region,
		std::list<abstractCE::ptr>& ce_candidates) {

	std::list < symbolic_states::ptr > list_sym_states;
	std::list < abstract_symbolic_state::ptr > list_abstract_sym_states;
	polytope::ptr abs_flowpipe; //bounding_box Polytope
	polytope::ptr polyI; // initial polytope of the abstract flowpipe.
	bool saftey_violated = false;
	std::list < transition::ptr > list_transitions;

	if (computedSymStates->getContinuousSetptr()->getTotalIterations() != 0) { //flowpipe exists
		//so perform intersection with forbidden set provided locID matches
		int locID;
		discrete_set ds;
		ds = computedSymStates->getDiscreteSet();
		for (std::set<int>::iterator it = ds.getDiscreteElements().begin();
				it != ds.getDiscreteElements().end(); ++it)
			locID = (*it); //Assuming only a single element exist in the discrete_set
		/*location current_location;
		 current_location = H.getLocation(locID);*/

		//int forbid_locID = current_location.getLocId();
		int forbid_locID = locID;
		if (forbid_locID == forbidden_set.first) { //forbidden locID matches
			polytope::ptr forbid_poly = forbidden_set.second;
			//check intersection with flowpipe/reach_region
			//GeneratePolytopePlotter(forbid_poly);
			std::list < template_polyhedra::ptr > forbid_intersects;
			//forbid_intersects = computedSymStates->getContinuousSetptr()->polys_intersectionParallel(forbid_poly, lp_solver_type_choosen);

			forbid_intersects = computedSymStates->getContinuousSetptr()->polys_intersectionSequential(
							forbid_poly, lp_solver_type_choosen); //TODO:: CHECK RUNNING BOTH PARALLEL AND SEQUENTIAL

			if (forbid_intersects.size() == 0) {
				std::cout << "\nThe model does NOT violates SAFETY property!!!\n";
			} else {
				std::cout << "\nThe model violates SAFETY property!!!\n";
				symbolic_states::ptr current_forbidden_state;
				current_forbidden_state = computedSymStates;
			//				abstract_symbolic_state::ptr curr_abs_sym_state;
			//				curr_abs_sym_state = abstract_symbolic_state::ptr(new abstract_symbolic_state());
				std::cout << "\nReverse Path Trace =>\n";
				int cc = 0;
				do {
					int locationID, locationID2;
					discrete_set ds, ds2;
					ds = current_forbidden_state->getDiscreteSet();

					abs_flowpipe = convertBounding_Box(
							current_forbidden_state->getContinuousSetptr());
					polyI =
							current_forbidden_state->getInitial_ContinousSetptr();
					// Here create a new abstract_symbolic_state
					abstract_symbolic_state::ptr curr_abs_sym_state =
							abstract_symbolic_state::ptr(
									new abstract_symbolic_state(ds,
											abs_flowpipe, polyI));

					// ***********insert bounding_box_polytope as continuousSet in the abstract_symbolic_state***********

					for (std::set<int>::iterator it =
							ds.getDiscreteElements().begin();
							it != ds.getDiscreteElements().end(); ++it)
						locationID = (*it); //Assuming only a single element exist in the discrete_set

					int transID = current_forbidden_state->getTransitionId(); //a)
					//   **********************************************************
					//todo::create an object of abstractCE[1)list_of_symbolic_states 2)list_of_transition and 3) length]
					//1) ******************** list_of_symbolic_states ********************
					list_sym_states.push_front(current_forbidden_state); //pushing the bad symbolic_state first(at the top)
					list_abstract_sym_states.push_front(curr_abs_sym_state);
					//2) list_of_transition
					//a) current sym_state only have trans_ID but to retrieve this transition I have to
					//b) get the parent to this sym_state using getParentPtrSymbolicState and then in
					//c) that sym_state get the discrete state and then get its loc_ID.
					//d) Now from the class hybrid_automata get an object of the class location
					//e) Now for this location's object get the transition with transID as transition_ID from
					// the data member out_going_transitions.
					//3) length: number of transitions
					//   **********************************************************

					if (cc != 0) {
						std::cout << " --> ";
					}
					std::cout << "(" << locationID << ", " << transID << ")";
					current_forbidden_state =
							searchSymbolic_state(Reachability_Region,
									current_forbidden_state->getParentPtrSymbolicState()); //b)

					//2) ******************* list_transitions ********************
					ds2 = current_forbidden_state->getDiscreteSet(); //c)
					for (std::set<int>::iterator it =
							ds2.getDiscreteElements().begin();
							it != ds2.getDiscreteElements().end(); ++it)
						locationID2 = (*it); //c)
					location::ptr object_location;
					object_location = H.getLocation(locationID2); //d)
					transition::ptr temp = object_location->getTransition(
							transID); //e)
					list_transitions.push_front(temp); //pushing the transition in the stack
					//2) ******************* list_transitions Ends ********************
					cc++;
				} while (current_forbidden_state->getParentPtrSymbolicState()
						!= NULL);

				if ((cc >= 1)
						&& (current_forbidden_state->getParentPtrSymbolicState()
								== NULL)) { //root is missed
					ds = current_forbidden_state->getDiscreteSet();

					int locationID, locationID2;
					discrete_set ds, ds2;
					ds = current_forbidden_state->getDiscreteSet();

					abs_flowpipe = convertBounding_Box(
							current_forbidden_state->getContinuousSetptr());
					polyI =
							current_forbidden_state->getInitial_ContinousSetptr();
					// Here create a new abstract_symbolic_state
					abstract_symbolic_state::ptr curr_abs_sym_state =
							abstract_symbolic_state::ptr(
									new abstract_symbolic_state(ds,
											abs_flowpipe, polyI));

					// ***********insert bounding_box_polytope as continuousSet in the abstract_symbolic_state***********

					for (std::set<int>::iterator it =
							ds.getDiscreteElements().begin();
							it != ds.getDiscreteElements().end(); ++it)
						locationID = (*it); //Assuming only a single element exist in the discrete_set

					int transID = current_forbidden_state->getTransitionId();
					list_sym_states.push_front(current_forbidden_state); //1) pushing the initial/root bad symbolic_state at the top
					list_abstract_sym_states.push_front(curr_abs_sym_state);

					std::cout << " -->  (" << locationID << ", " << transID
							<< ")\n";
				}
				saftey_violated = true;
#pragma omp critical
				{
					abstractCE::ptr ce = abstractCE::ptr(new abstractCE());
					ce->set_length(cc);
					ce->set_sym_states(list_abstract_sym_states);
					ce->set_transitions(list_transitions);
					hybrid_automata::ptr h = hybrid_automata::ptr(
							new hybrid_automata(H));
					ce->set_automaton(h);
					ce->set_forbid_poly(forbid_poly);
					ce_candidates.push_back(ce); // ce added to the candidates list
				}
			} // end of condition when forbidden state intersects with the flowpipe set
		} //end of condition when forbidden state loc id matches with flowpipe loc id
	} //computed flowpipe is not empty
	return saftey_violated;
}

void reachability::computeBIG_Task(std::vector<LoadBalanceData>& LoadBalanceDS) {
	for (int i = 0; i < LoadBalanceDS.size(); i++) { //for each symbolic-states
		assert(i == LoadBalanceDS[i].symState_ID);
		int dimension = LoadBalanceDS[i].List_dir_X0.size2();
		lp_solver lp(this->lp_solver_type_choosen), lp_U(
				this->lp_solver_type_choosen);
		lp.setMin_Or_Max(2); //2 for Maximization
		//cout << "Testing A1\n";
		lp.setConstraints(LoadBalanceDS[i].reach_param.X0->getCoeffMatrix(),
				LoadBalanceDS[i].reach_param.X0->getColumnVector(),
				LoadBalanceDS[i].reach_param.X0->getInEqualitySign());
		//cout << "Testing A2\n";
		lp_U.setMin_Or_Max(2);
		bool U_empty;
		U_empty =
				LoadBalanceDS[i].current_location->getSystem_Dynamics().U->getIsEmpty();
		if (!U_empty) {
			LoadBalanceDS[i].sf_U.resize(LoadBalanceDS[i].List_dir_U.size1());
			lp_U.setConstraints(LoadBalanceDS[i].U->getCoeffMatrix(),
					LoadBalanceDS[i].U->getColumnVector(),
					LoadBalanceDS[i].U->getInEqualitySign());
			for (unsigned int j = 0; j < LoadBalanceDS[i].List_dir_U.size1();
					j++) { //for all directions ** Can be done in loop **** A ******
				std::vector<double> dirs(dimension);
				for (int index = 0; index < dimension; index++) {
					dirs[index] = LoadBalanceDS[i].List_dir_U(j, index);
				}
				LoadBalanceDS[i].sf_U[j] = lp_U.Compute_LLP(dirs);
			}
		}
		//cout << "Testing A3\n";
		//unsigned int tot_dirs = LoadBalanceDS[i].List_dir_X0.size1();
		LoadBalanceDS[i].sf_X0.resize(LoadBalanceDS[i].List_dir_X0.size1());
		LoadBalanceDS[i].sf_dotProduct.resize(
				LoadBalanceDS[i].List_dir_X0.size1());
		LoadBalanceDS[i].sf_UnitBall.resize(
				LoadBalanceDS[i].List_dir_X0.size1());
//#pragma omp parallel for
		for (unsigned int j = 0; j < LoadBalanceDS[i].List_dir_X0.size1();
				j++) { //for all directions  ****** A ******
			std::vector<double> dirs(dimension);
			for (int index = 0; index < dimension; index++) {
				dirs[index] = LoadBalanceDS[i].List_dir_X0(j, index);
			}
			//cout << "Testing A4\n";
			LoadBalanceDS[i].sf_X0[j] = lp.Compute_LLP(dirs);
			//cout << "Testing A5\n";
			// ******DotProduction and Support Function of UnitBall  *******
			LoadBalanceDS[i].sf_dotProduct[j] = dot_product(LoadBalanceDS[i].current_location->getSystem_Dynamics().C, dirs);
			LoadBalanceDS[i].sf_UnitBall[j] = support_unitball_infnorm(dirs);
		}
	} //end-for each symbolic-states
}

//can be parallelized to task-based approach
void reachability::parallelBIG_Task(std::vector<LoadBalanceData>& LoadBalanceDS) {

	int dimension = LoadBalanceDS[0].List_dir_X0.size2();
	unsigned int countTotal_X = 0, countTotal_U = 0;
	//cout <<"erro 1\n";
	for (int i = 0; i < LoadBalanceDS.size(); i++) { //for each symbolic-states
		countTotal_X = countTotal_X + LoadBalanceDS[i].List_dir_X0.size1();
	//	cout<<"   = "<<LoadBalanceDS[i].List_dir_X0.size1()<<std::endl;
		countTotal_U = countTotal_U + LoadBalanceDS[i].List_dir_U.size1();
		// *********** resize all result vector  *********************
		LoadBalanceDS[i].sf_X0.resize(LoadBalanceDS[i].List_dir_X0.size1()); // resize
		LoadBalanceDS[i].sf_U.resize(LoadBalanceDS[i].List_dir_U.size1()); // resize
		LoadBalanceDS[i].sf_UnitBall.resize(LoadBalanceDS[i].List_dir_X0.size1()); // resize
		LoadBalanceDS[i].sf_dotProduct.resize(LoadBalanceDS[i].List_dir_X0.size1()); // resize
	} //getCountTotal(LoadBalanceDS, countTotal_X, countTotal_U);
	//cout << "countTotal_X = " << countTotal_X << std::endl;
//cout <<"erro 2\n";
/*	std::ofstream myfile,myfile2;
	myfile.open("./sf_X.txt");
	myfile2.open("./dir_X.txt");*/
//#pragma omp parallel for
	for (int i = 0; i < countTotal_X; i++) {
		int index;
		unsigned int j;
		//cout <<"erro 3\n";
		search_SymState_dirsX0Index(i, LoadBalanceDS, index, j);
		//	cout <<"erro index = "<<index <<"  j = "<<j<<"\n";
		std::vector<double> dirs(dimension);
		for (int ind = 0; ind < dimension; ind++) {
			dirs[ind] = LoadBalanceDS[index].List_dir_X0(j, ind);
		//	myfile2<<dirs[ind]<< "\t";
		}
		//myfile2<< "\n";
		LoadBalanceDS[index].sf_X0[j] = LPSolver(LoadBalanceDS[index].X0, dirs);
	//	myfile <<"i = "<<i<<"   ";
	//	myfile <<LoadBalanceDS[index].sf_X0[j]<<"\n";
		//LoadBalanceDS[index].sf_X0[j] = boxLPSolver(LoadBalanceDS[index].X0, dirs);
		// ******DotProduction and Support Function of UnitBall  *******
		if (!LoadBalanceDS[index].current_location->getSystem_Dynamics().isEmptyC) {
			LoadBalanceDS[index].sf_dotProduct[j] = dot_product(LoadBalanceDS[index].current_location->getSystem_Dynamics().C,dirs);
		}
		LoadBalanceDS[index].sf_UnitBall[j] = support_unitball_infnorm(dirs);
		// ******DotProduction and Support Function of UnitBall  *******
	}
	//myfile.close();
	//myfile2.close();
	//cout <<"erro 3\n";
	bool U_empty;
	U_empty = LoadBalanceDS[0].current_location->getSystem_Dynamics().U->getIsEmpty();
	//todo:: assuming all symbolic states has same setup for polytope U
	if (!U_empty) {
		cout<<"polytope U is NOT empty!!!!\n";
#pragma omp parallel for
		for (int i = 0; i < countTotal_U; i++) {
			int index;
			unsigned int j;
			search_SymState_dirsUIndex(i, LoadBalanceDS, index, j);
			std::vector<double> dirs(dimension);
			for (int ind = 0; ind < dimension; ind++) {
				dirs[ind] = LoadBalanceDS[index].List_dir_U(j, ind);
			}
			LoadBalanceDS[index].sf_U[j] = LPSolver(LoadBalanceDS[index].U, dirs);
		}
	}
}


double reachability::LPSolver(polytope::ptr poly, std::vector<double> dirs) {
	double res;
	if (poly->getIsEmpty())
		res = 0.0;
	else {
		lp_solver lp(this->lp_solver_type_choosen);
		lp.setMin_Or_Max(2); //for maximizing
		lp.setConstraints(poly->getCoeffMatrix(), poly->getColumnVector(), poly->getInEqualitySign());
		res = lp.Compute_LLP(dirs);
	}
	return res;
}

/*
 * A hyperbox LP solver
 */

double reachability::boxLPSolver(polytope::ptr poly, std::vector<double> dir) {
	double res = 0.0;
	math::matrix<double> M = poly->getCoeffMatrix();
	math::vector<double> b = poly->getColumnVector();
	double lb, ub;
	for (unsigned int i = 0; i < dir.size(); i++) {
		if (dir[i] < 0) {
			//iterate over the ith component of each row of M
			for (unsigned int row = 0; row < M.size1(); row++) {
				if (M(row, i) == -1) { // lower bound constraint
					lb = -1 * b[row];
					break;
				}
			}
			res += dir[i] * lb;
		} else {
			//iterate over the ith component of each row of M
			for (unsigned int row = 0; row < M.size1(); row++) {
				if (M(row, i) == 1) { // lower bound constraint
					ub = b[row];
					break;
				}
			}
			res += dir[i] * ub;
		}
	}
	return res;
}

