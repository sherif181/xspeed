/*
 * TPBFS.cpp
 *
 *  Created on: 28-Oct-2016
 *      Author: amit
 */

#include "TPBFS.h"

/*
 * Aggregrate All Flowpipe computation work into one BIG task and will run that in parallel either by multi-core CPU or GPU
 * Also the Post_D is done in the same manner
 */
std::list<symbolic_states::ptr> tpbfs::LoadBalanceAll(std::list<abstractCE::ptr>& ce_candidates) {

	std::list < symbolic_states::ptr > Reachability_Region; //	template_polyhedra::ptr reach_region;
	int t = 0; //0 for Read and 1 for Write

	std::vector < std::vector<pwlist::ptr> > Qpw_list(2); // QpwList[0] for read and QpwList[1] for write 	//cout << "Test 1\n";

	/*Qpw_list[t].resize(1); //resize for the first symbolic_state 	//cout << "Test 2\n";
	Qpw_list[t][0] = pwlist::ptr(new pwlist()); //have to instantiate it
	Qpw_list[t][0]->WaitingList_insert(I); //cout << "Test 3\n";*/

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

	pwlist::ptr allPassedList; //so we create a permanent pwlist for storing only the passedList;
	allPassedList = pwlist::ptr(new pwlist()); //have to instantiate it
	int number_times = 0;
	bool levelcompleted = false;
	unsigned int iter_max = 1;
	int numCores = omp_get_num_procs(); //get the number of cores
	/*boost::timer::cpu_timer t70;
	t70.start();*/

	while (!isEmpty_Qpw_list(Qpw_list[t]) && (number_times <= bound)) { //	cout << "Test 5\n";
	//	cout<<"Breadth - Level === "<<number_times<<"\n";
//		boost::timer::cpu_timer t73;
//		t73.start();
		unsigned int count = getSize_Qpw_list(Qpw_list[t]); //get the size of PWList
		cout << "\nNumber of Flowpipes to be Computed (per Breadths) = " << count << "\n";
		std::vector < symbolic_states::ptr > S(count);
		vector < initial_state::ptr > list_U(count); //SubList for parallel		//	cout << "Test 6\n";
		list_U = getAllpw_list(Qpw_list, t, count, allPassedList); //All initial_state have been deleted //cout<<"Identifed pwList Done\n";
		//cout<<"list_U = "<<list_U.size();

		Qpw_list[t].resize(0);
		//cout<<"Qpw_list[t].size() = "<<Qpw_list[t].size()<<std::endl;

		Qpw_list[1 - t].resize(count); //resize to accommodate
		for (int i = 0; i < count; i++) {
			Qpw_list[1 - t][i] = pwlist::ptr(new pwlist()); //have to instantiate it
		}
//		t73.stop();
//		double clock73;
//		clock73 = t73.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
//		double return73 = clock73 / (double) 1000;
//		std::cout << "\nHow Much queue delete  Time:Wall(Seconds) = "<< return73 << std::endl;

		// ********************************* BFS Starts **********************************************************
		//	cout << "Test 7\n";
// ***************** DIRECTION and BIG_Task CREATION:: per breadth-level *****************
		/*boost::timer::cpu_timer t705;
		t705.start();*/
		std::vector < LoadBalanceData > LoadBalanceDS(count);
// ***************** DIRECTION and BIG_Task CREATION *****************

		//omp_set_nested(1);	//enables nested parallelization
		/*boost::timer::cpu_timer t702;
		t702.start();*/
#pragma omp parallel for // num_threads(count)
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
			S[id] = symbolic_states::ptr(new symbolic_states());
			S[id]->setDiscreteSet(discrete_state);
			S[id]->setParentPtrSymbolicState(U->getParentPtrSymbolicState()); //keeps track of parent pointer to symbolic_states
			S[id]->setTransitionId(U->getTransitionId()); //keeps track of originating transition_ID
			location::ptr current_location;
			current_location = H.getLocation(location_id);
			string name = current_location->getName();
			//	cout<<"Location Name = "<< name<<"\n";
			if ((name.compare("GOOD") == 0) || (name.compare("BAD") == 0)
					|| (name.compare("UNSAFE") == 0) || (name.compare("FINAL") == 0))
				continue; //do not compute the continuous reachability algorithm
			// ******************* Computing Parameters ************************ //current_location:: parameters alfa, beta and phi_trans
			//cout<<"\nBefore Compute Alfa";
			double result_alfa = compute_alfa(reach_parameter_local.time_step,
					current_location->getSystem_Dynamics(), continuous_initial_polytope, lp_solver_type_choosen); //cout<<"\nCompute Alfa Done";
			//cout<<"\nCompute alfa  = " << result_alfa;
			double result_beta = compute_beta(current_location->getSystem_Dynamics(),
					reach_parameter_local.time_step, lp_solver_type_choosen); // NO glpk object created here
			//	cout<<"\nCompute Beta  = " << result_beta;
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
			// ******************* Computing Parameters Done *******************************
// *************** POST_C computation ********** in 3 steps ***************
			bool U_empty = false;
			if (current_location->getSystem_Dynamics().U->getIsEmpty()) { //polytope U can be empty set
				U_empty = true;
			}
			LoadBalanceDS[id].X0 = continuous_initial_polytope;
			LoadBalanceDS[id].U = current_location->getSystem_Dynamics().U;
			LoadBalanceDS[id].current_location = current_location;
			LoadBalanceDS[id].symState_ID = id;
			LoadBalanceDS[id].reach_param = reach_parameter_local;
		}	//END of count FOR-LOOP  -- reach parameters computed

#pragma omp parallel for // num_threads(count)
		for (unsigned int id = 0; id < count; id++) { //separate parameters assignment with Invariant check and preLoadBalanceReachCompute task
			unsigned int NewTotalIteration;

			if (LoadBalanceDS[id].current_location->isInvariantExists()) {

				if (LoadBalanceDS[id].current_location->getSystem_Dynamics().isEmptyMatrixB == true
						&& LoadBalanceDS[id].current_location->getSystem_Dynamics().isEmptyC == true) {
					//Approach of Coarse-time-step and Fine-time-step
					jumpInvariantBoundaryCheck(
							LoadBalanceDS[id].current_location->getSystem_Dynamics(),
							LoadBalanceDS[id].X0, LoadBalanceDS[id].reach_param,
							LoadBalanceDS[id].current_location->getInvariant(),
							lp_solver_type_choosen, NewTotalIteration);
				} else {
					//Approach of Sequential invariant check will work for all case
					InvariantBoundaryCheck(
							LoadBalanceDS[id].current_location->getSystem_Dynamics(),
							LoadBalanceDS[id].X0, LoadBalanceDS[id].reach_param,
							LoadBalanceDS[id].current_location->getInvariant(),
							lp_solver_type_choosen, NewTotalIteration);
				}
				LoadBalanceDS[id].newIteration = NewTotalIteration; //Important to take care
				//cout<<"Invariant setting Done\n";
			}else
				LoadBalanceDS[id].newIteration = LoadBalanceDS[id].reach_param.Iterations; //Important to take care
			//std::cout << "NewTotalIteration = " << NewTotalIteration << std::endl;
		}	//END of count FOR-LOOP  -- invariant boundary check done

		int numCoreAvail = 1;
		if (numCores >= count)
			numCoreAvail = numCores - count;

		#pragma omp parallel for //num_threads(count)
		for (unsigned int id = 0; id < count; id++) {
			math::matrix<float> listX0, listU;
			//Generation of Directions which can be a nested parallelism
			preLoadBalanceReachCompute(LoadBalanceDS[id],numCoreAvail);	// Step 1
			//Todo:: if newIteration is <= 1 than do not store this flow-pipe details in LoadBalanceDS
			// and accordingly reduce count or handle while computing LPsolver
		} //END of count FOR-LOOP

		parallelLoadBalance_Task(LoadBalanceDS);//So an appropriate combination of parallel with sequential GLPK object is used.
		omp_set_nested(1); //enable nested parallelism

#pragma omp parallel for // num_threads(count)
		for (unsigned int id = 0; id < count; id++) {
			S[id]->setContinuousSetptr(substitute_in_ReachAlgorithm(LoadBalanceDS[id], numCoreAvail)); // Step 3
		}

//  ********************* POST_C computation Done ********************

		bool foundUnSafe = false;
		for (unsigned int index = 0; index < count; index++) {
			if (S[index]->getContinuousSetptr()->getTotalIterations() != 0) { //computed reach_region is NOT empty
				Reachability_Region.push_back(S[index]);
			}
		}
		number_times++; //One Level or one Breadth Search over

		if (number_times > bound) {
			levelcompleted = true; //check to see how many jumps have been made(i.e., number of discrete transitions made)
		}
		if (levelcompleted || foundUnSafe) { //any true
			break; //OUT FROM WHILE LOOP 	//no need to compute rest of the locations
		}
std::vector<LoadBalanceData_PostD> loadBalPostD(count);
//  ***************** Load Balanced POST_D computation Begins::Has 3 Steps ***********************************
//Step--1 :: Pre-Load Balancing Task to populate guard, assign, etc in Data Structure
//cout<<"done 1\n";
#pragma omp parallel for // num_threads(count)
		for (int id = 0; id < count; id++) {
			if (S[id]->getContinuousSetptr()->getTotalIterations() != 0) {	//if flowpipe not Empty
				loadBalPostD[id].sfm = S[id]->getContinuousSetptr();
				int transition_size = LoadBalanceDS[id].current_location->getOut_Going_Transitions().size();
				loadBalPostD[id].trans_size = transition_size;
				//cout<<"transition_size = "<<transition_size <<"\n";
				loadBalPostD[id].guard_list.resize(transition_size);
				loadBalPostD[id].assign_list.resize(transition_size);
				loadBalPostD[id].dest_locID.resize(transition_size);
				loadBalPostD[id].trans_ID.resize(transition_size);
				loadBalPostD[id].range_list.resize(transition_size);
				loadBalPostD[id].bool_arr.resize(transition_size, 0);//initially only rows are known

				loadBalPostD[id].polys_list.resize(transition_size);
				int index=0; // for guard_list, assign_list if they are vector
				for (std::list<transition::ptr>::iterator trans =LoadBalanceDS[id].current_location->getOut_Going_Transitions().begin();
					trans!= LoadBalanceDS[id].current_location->getOut_Going_Transitions().end(); trans++) { // get each destination_location_id and push into the pwl.waiting_list
					loadBalPostD[id].guard_list[index] = (*trans)->getGaurd();
					loadBalPostD[id].assign_list[index] = (*trans)->getAssignT();
					loadBalPostD[id].dest_locID[index] = (*trans)->getDestination_Location_Id();
					index++;
				}	//each transition
			} else {//EndIf of flowpipe not Empty
				cout<<"This condition Also MET!!!\n";
			}
		}
	//	cout<<"done 2\n";
//Step--2 :: First Load Balancing computation Task
		loadBalancedPostD(loadBalPostD);	//returns array of booleans with true or FALSE based on intersected or NOT
	//	cout<<"done 3\n";
//Step--3 :: Second Task to detect sequentially
		intersectionRangeDetection(loadBalPostD);	//returns list of range pair indicating Start and End as SFM's columnIndex
	//	cout<<"done 4\n";
//Step--4 :: Third Task to perform Template_Approximation
		templateApproximation(loadBalPostD);
	//	cout<<"done 5\n";
//Step--5 :: Finally creating new symbolic states to push in the Qpw_list[1-t]
#pragma omp parallel for
		for (int id = 0; id < count; id++) {

			for (int trans=0;trans<loadBalPostD[id].trans_size;trans++){
				std::list<polytope::ptr> polys;
				polys = loadBalPostD[id].polys_list[trans];
				#pragma omp critical
				{
					iter_max += polys.size();
				}
				string locName = H.getLocation(loadBalPostD[id].dest_locID[trans])->getName();
				if ((locName.compare("BAD") == 0) || (locName.compare("GOOD") == 0) || (locName.compare("FINAL") == 0)
						|| (locName.compare("UNSAFE") == 0)) {
			/*		#pragma omp critical
					{
						iter_max = iter_max - 1;
					}*/
					continue; //do not push into the waitingList
				}
				polytope::ptr intersectedRegion;
				for (std::list<polytope::ptr>::iterator it = polys.begin(); it != polys.end(); it++) {
					intersectedRegion = (*it);
					polytope::ptr newShiftedPolytope, newPolytope; //created an object here
					newPolytope = intersectedRegion->GetPolytope_Intersection(loadBalPostD[id].guard_list[trans]); //Retuns only the intersected region as a single newpolytope. ****** with added directions

					/*newShiftedPolytope = post_assign_exact(newPolytope, loadBalPostD[id].assign_list[trans].Map,
						loadBalPostD[id].assign_list[trans].b); */

					math::matrix<double> test(
							loadBalPostD[id].assign_list[trans].Map.size1(),
							loadBalPostD[id].assign_list[trans].Map.size2());
					if (loadBalPostD[id].assign_list[trans].Map.inverse(test)) //invertible?
							{
						//std::cout << "Exact Post Assignment\n";
						newShiftedPolytope = post_assign_exact(newPolytope,
								loadBalPostD[id].assign_list[trans].Map,
								loadBalPostD[id].assign_list[trans].b);
					} else {
						//std::cout << "Approximate Post Assignment\n";
						newShiftedPolytope = post_assign_approx_deterministic( newPolytope,
								loadBalPostD[id].assign_list[trans].Map, loadBalPostD[id].assign_list[trans].b,
								reach_parameters.Directions, lp_solver_type_choosen);
					}

					initial_state::ptr newState = initial_state::ptr(new initial_state(loadBalPostD[id].dest_locID[trans], newShiftedPolytope));
					newState->setTransitionId(loadBalPostD[id].trans_ID[trans]); // keeps track of the transition_ID
					newState->setParentPtrSymbolicState(S[id]);

					Qpw_list[1 - t][id]->WaitingList_insert(newState); //False RACE CONDITION HERE
				} //end of multiple intersected region with guard
			}//end of for each guard or transitions
		}//end of for parallel

//  ************************************** POST_D computation Ends **********************************************************



		t = 1 - t; //Switching Read/Write options for Qpw_list[1-t]
		// ************************* BFS Ends *************************************
	} //end of while loop checking waiting_list != empty
	cout << "\n **********************************************************\n";
	cout << "   *** Maximum Iterations Completed = " << iter_max << "  ***\n";
	cout << "\n **********************************************************\n";

	return Reachability_Region;
}


template_polyhedra::ptr tpbfs::substitute_in_ReachAlgorithm(
		LoadBalanceData& LoadBalanceDS, int numCoreAvail) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;
	template_polyhedra::ptr reachableRegion; //template_polyhedra::ptr reachRegion;

/*	boost::timer::cpu_timer reachLoop_time;
	reachLoop_time.start();*/
	unsigned int numVectors = LoadBalanceDS.reach_param.Directions.size1(); //size2 or the dimension will be some for all sym_state
	int num_inv =
			LoadBalanceDS.current_location->getInvariant()->getColumnVector().size(); //number of Invariant's constriants
	math::matrix<double> inv_directions;
	inv_directions =
			LoadBalanceDS.current_location->getInvariant()->getCoeffMatrix();
	std::vector<double> inv_bounds(num_inv);
	inv_bounds =
			LoadBalanceDS.current_location->getInvariant()->getColumnVector();

	if (LoadBalanceDS.newIteration <= 1) {
		template_polyhedra::ptr poly_emptyp;
		return poly_emptyp;
	}
	int dimension = LoadBalanceDS.X0->getSystemDimension();
	int Min_Or_Max = 2;
	size_type row = numVectors, col = LoadBalanceDS.newIteration;
	math::matrix<double> MatrixValue(row, col);

	int cores;
	if (numVectors >= numCoreAvail)
		cores = numVectors;
	else
		cores = numCoreAvail;

	for (unsigned int eachDirection = 0; eachDirection < numVectors; eachDirection++) {
		unsigned int index_X0, index_U; //making the index suitable for parallelizing
		//unsigned int index; 	//index = eachDirection * NewTotalIteration;
		//here i have a list of result of Supp_fun_Of_UnitBall_infinity_norm
		if (eachDirection == 0) { //only starting loop begins with 0
			index_X0 = eachDirection * LoadBalanceDS.newIteration;
		} else { //
			index_X0 = eachDirection * LoadBalanceDS.newIteration + eachDirection; //only X0(list_X0) has 2 directions for first-iteration
		}
		bool U_empty =
				LoadBalanceDS.current_location->getSystem_Dynamics().U->getIsEmpty();
		if (!U_empty) {
			index_U = eachDirection * LoadBalanceDS.newIteration;
		}
		double res1, result1;
		double term1, term2, term3, term3a, term3b, res2, term3c = 0.0;
		double zIInitial = 0.0, zI = 0.0, zV = 0.0;
		double sVariable = 0.0, s1Variable; //initialize s0
		std::vector<double> rVariable(dimension), r1Variable(dimension);
		unsigned int loopIteration = 0;

		//  ************** Omega Function   ********************
		res1 = LoadBalanceDS.sf_X0[index_X0]; //X0->SF(direction)			//	0

		term3b = (double) LoadBalanceDS.sf_UnitBall[index_X0]; //  needed  0
		if (!LoadBalanceDS.current_location->getSystem_Dynamics().isEmptyC) {
			term3c = LoadBalanceDS.reach_param.time_step * LoadBalanceDS.sf_dotProduct[index_X0];
		}
		index_X0++; //	made 1
		term1 = LoadBalanceDS.sf_X0[index_X0]; //X0->SF(phi_trans_dir)		//  1
		index_X0++; //	made 2
		if (!U_empty) {
			term2 = LoadBalanceDS.reach_param.time_step * LoadBalanceDS.sf_U[index_U]; //U->SF(Btrans_dir)
			index_U++;
		} else
			term2 = 0;

		term3a = LoadBalanceDS.reach_param.result_alfa; //compute_alfa(ReachParameters.time_step,system_dynamics,Initial_X0);
		term3 = term3a * term3b;
		res2 = term1 + term2 + term3 + term3c; //term3c Added
		//zIInitial = (res1 > res2 ? res1:res2);
		if (res1 > res2)
			zIInitial = res1;
		else
			zIInitial = res2;
		//  ************** Omega Function   ********************
		MatrixValue(eachDirection, loopIteration) = zIInitial; //		index++;
		//	cout<<"zIInitial = "<< zIInitial<<std::endl;
		loopIteration++;
		for (; loopIteration < LoadBalanceDS.newIteration;) { //Now stopping condition is only "shm_NewTotalIteration"
			//  ************** W_Support Function   ********************
			double result, res_sup;
			if (!U_empty) {
				result1 = term2;
			}
			double beta = LoadBalanceDS.reach_param.result_beta;

			//double res_beta = beta * res_sup;
			double res_beta = beta * term3b;

			result = result1 + res_beta + term3c; //Added term3c
			zV = result;
			//  ************** W_Support Function   ********************
			s1Variable = sVariable + zV;
			//  ************** Omega Function   ********************
			res1 = term1; ////replace

			term3b = (double) LoadBalanceDS.sf_UnitBall[index_X0 - 1]; //Compute here	//needed 1

			if (!LoadBalanceDS.current_location->getSystem_Dynamics().isEmptyC) {
				term3c = LoadBalanceDS.reach_param.time_step * LoadBalanceDS.sf_dotProduct[index_X0 - 1];
			}
			double term3, term3a, res2;
			term1 = LoadBalanceDS.sf_X0[index_X0]; //X0->SF(phi_trans_dir)		//	2
			index_X0++; // 	made 3
			if (!U_empty) {
				term2 = LoadBalanceDS.reach_param.time_step * LoadBalanceDS.sf_U[index_U]; //U->SF(Btrans_dir)
				index_U++;
			} else {
				term2 = 0;
			}
			term3a = LoadBalanceDS.reach_param.result_alfa; //compute_alfa(ReachParameters.time_step,system_dynamics,Initial_X0);
			term3 = term3a * term3b;
			res2 = term1 + term2 + term3 + term3c;
			//zIInitial = (res1 > res2 ? res1:res2);
			if (res1 > res2)
				zI = res1;
			else
				zI = res2;
			//  ************** Omega Function   ********************
			double TempOmega;
			TempOmega = zI + s1Variable; //Y1
			//	cout<<"TempOmega = "<< TempOmega<<std::endl;
			MatrixValue(eachDirection, loopIteration) = TempOmega; //Y1
			sVariable = s1Variable; //	index++;
			loopIteration++; //for the next Omega-iteration or Time-bound
		} //end of all Iterations of each vector/direction
	} //end of for each vector/directions

	std::cout << std::fixed;
	std::cout.precision(10);
	if (LoadBalanceDS.current_location->isInvariantExists() == true) { //if invariant exist. Computing

		math::matrix<double> inv_sfm;
		inv_sfm.resize(num_inv, LoadBalanceDS.newIteration);
		for (int eachInvDirection = 0; eachInvDirection < num_inv; eachInvDirection++) {
			for (unsigned int i = 0; i < LoadBalanceDS.newIteration; i++) {
				inv_sfm(eachInvDirection, i) = inv_bounds[eachInvDirection];

			}
		}

		reachableRegion = template_polyhedra::ptr(new template_polyhedra());
		reachableRegion->setTemplateDirections(LoadBalanceDS.reach_param.Directions);
		reachableRegion->setMatrix_InvariantBound(inv_sfm);
		reachableRegion->setInvariantDirections(inv_directions);
		reachableRegion->setMatrixSupportFunction(MatrixValue);
	} else {
		reachableRegion = template_polyhedra::ptr(new template_polyhedra());
		reachableRegion->setMatrixSupportFunction(MatrixValue);
		reachableRegion->setTemplateDirections(LoadBalanceDS.reach_param.Directions);
	}

	return reachableRegion;

}


void tpbfs::parallelLoadBalance_Task(std::vector<LoadBalanceData>& LoadBalanceDS) {
//cout<<"Inside parallelLoadBalance_Task\n";
	int numCores = omp_get_num_procs(); //get the number of cores
//	cout<<"numCores = "<<numCores<<std::endl;
//	int numSymStates = LoadBalanceDS.size(); //get the number of symbolic states i.e. number of polytopes X0
//	cout<<"numSymStates = "<<numSymStates<<std::endl;
	unsigned int seq_LP;
	int size, chunk_size;
	/*if (numCores > numSymStates)
		size = numCores;
	else
		size = numSymStates;*/
	int dimension = LoadBalanceDS[0].List_dir_X0.size2();
	unsigned int countTotal_X = 0, countTotal_U = 0;
	for (int i = 0; i < LoadBalanceDS.size(); i++) { //for each symbolic-states
		countTotal_X += LoadBalanceDS[i].List_dir_X0.size1();
		//cout<<"   = "<<LoadBalanceDS[i].List_dir_X0.size1()<<std::endl;
		countTotal_U += LoadBalanceDS[i].List_dir_U.size1();
		// *********** resize all result vector  *********************
		LoadBalanceDS[i].sf_X0.resize(LoadBalanceDS[i].List_dir_X0.size1()); // resize
		LoadBalanceDS[i].sf_U.resize(LoadBalanceDS[i].List_dir_U.size1()); // resize
		LoadBalanceDS[i].sf_UnitBall.resize(LoadBalanceDS[i].List_dir_X0.size1()); // resize
		LoadBalanceDS[i].sf_dotProduct.resize(LoadBalanceDS[i].List_dir_X0.size1()); // resize
	} //getCountTotal(LoadBalanceDS, countTotal_X, countTotal_U);
// ************* Chunk_approach for polytope X ******************************
	//cout << "countTotal_X = " << countTotal_X << std::endl;
	/*if (countTotal_X <= size) {
		chunk_size = (int)countTotal_X;	//todo:: may be do it sequential
	} else
		chunk_size = compute_chunk_size(countTotal_X);//todo: find an appropriate number of chunk-size based on the size of countTotal_X*/
	if (countTotal_X < numCores)
		chunk_size = 1;
	else
		chunk_size = numCores;
	seq_LP =(unsigned int) (countTotal_X / chunk_size); //Last group/partition may not be equal if result is fraction
	//cout<<"  chunk_size = "<<chunk_size <<"  seq_LP = "<<seq_LP<<std::endl;

#pragma omp parallel for
	for (int i = 0; i < chunk_size; i++) { //number of partition for load balancing
		unsigned int j; //index for sequential inner-loop
		unsigned int lb = i * seq_LP, ub = (i * seq_LP + seq_LP);
		if (i == (chunk_size - 1))	//handling the Last group/partition
			ub = countTotal_X;
		j = lb;
		int index, oldIndex; //index of symbolic_state
		unsigned int  indexDir; //indices of symbolic_state and direction within that symbolic state
		search_SymState_dirsX0Index(j, LoadBalanceDS, index, indexDir);
		while (j < ub) {
			oldIndex = index;	//previously computed
			lp_solver lp(this->lp_solver_type_choosen); //1 for using GLPK solver
			lp.setMin_Or_Max(2); //Maximizing
			lp.setConstraints(LoadBalanceDS[index].X0->getCoeffMatrix(), LoadBalanceDS[index].X0->getColumnVector(),
					LoadBalanceDS[index].X0->getInEqualitySign());
			//for (;(j < ub) && (oldIndex == index);){
			while (j < ub){
				std::vector<double> dir(dimension);
				for (int ind = 0; ind < dimension; ind++) {
					dir[ind] = LoadBalanceDS[index].List_dir_X0(indexDir, ind);
				}
				LoadBalanceDS[index].sf_X0[indexDir] = lp.Compute_LLP(dir);	//support function of X0
			//	cout<<LoadBalanceDS[index].sf_X0[indexDir]<<"\t";
				// ******DotProduction and Support Function of UnitBall  *******
				//TODO: MAKE DECISION TO KEEP IT OUTSIDE AT ONE LOOP-VALUE LIKE BEFORE
				if (!LoadBalanceDS[index].current_location->getSystem_Dynamics().isEmptyC) {
					LoadBalanceDS[index].sf_dotProduct[indexDir] = dot_product(LoadBalanceDS[index].current_location->getSystem_Dynamics().C,dir);
				}
				LoadBalanceDS[index].sf_UnitBall[indexDir] = support_unitball_infnorm(dir);
				// ******DotProduction and Support Function of UnitBall  *******
				//cout<<j<<"\t";
				j++;//next LP problem
				search_SymState_dirsX0Index(j, LoadBalanceDS, index, indexDir);
				if (oldIndex != index){
					break;	//change of symbolic_state or polytope X0
				}
			}	//end-for	//cout<<"\n";
		}//end-while	//cout<<"\n\n";
	}//end of parallel
// ************* Chunk_approach for polytope X ******************************
	bool U_empty;
	U_empty = LoadBalanceDS[0].current_location->getSystem_Dynamics().U->getIsEmpty();//assuming all symbolic states has same setup for polytope U
	if (!U_empty) {
// ************* Chunk_approach for polytope U ******************************
	//	cout<<"polytope U is NOT empty!!!!\n";
		if (countTotal_X <= numCores)
			chunk_size = 1;
		else
			chunk_size = numCores;
		seq_LP =(unsigned int) (countTotal_U / chunk_size); //Last group/partition may not be equal if result is fraction
#pragma omp parallel for
		for (int i = 0; i < chunk_size; i++) { //number of partition for load balancing
			unsigned int j; //index for sequential inner-loop
			unsigned int lb = i * seq_LP, ub = (i * seq_LP + seq_LP);
			if (i == (chunk_size - 1))	//handling the Last group/partition
				ub = countTotal_U;
			j = lb;
			int index, oldIndex; //index of symbolic_state
			unsigned int  indexDir; //indices of symbolic_state and direction within that symbolic state
			search_SymState_dirsUIndex(j, LoadBalanceDS, index, indexDir);
			while (j < ub) {
				oldIndex = index;	//previously computed
				lp_solver lp(this->lp_solver_type_choosen); //1 for using GLPK solver
				lp.setMin_Or_Max(2); //Maximizing
				lp.setConstraints(LoadBalanceDS[index].U->getCoeffMatrix(), LoadBalanceDS[index].U->getColumnVector(),
						LoadBalanceDS[index].U->getInEqualitySign());
				while (j < ub){
					std::vector<double> dir(dimension);
					for (int ind = 0; ind < dimension; ind++) {
						dir[ind] = LoadBalanceDS[index].List_dir_U(indexDir, ind);
					}
					LoadBalanceDS[index].sf_U[indexDir] = lp.Compute_LLP(dir);	//support function of U
					//cout<<j<<"\t";
					j++;//next LP problem
					search_SymState_dirsUIndex(j, LoadBalanceDS, index, indexDir);
					if (oldIndex != index){
						break;	//change of symbolic_state or polytope X0
					}
				}	//end-for
			} //end-while
		} //end of parallel
// ************* Chunk_approach for polytope U ******************************
	} //end-if of empty check
}



void tpbfs::preLoadBalanceReachCompute(LoadBalanceData& LoadBalanceDS, int numCoresAvail){

unsigned int NewTotalIteration;

//	unsigned int NewTotalIteration = ReachParameters.Iterations;
//	cout <<"Before NewTotalIteration = " <<NewTotalIteration<<"\n";
	bool U_empty = false;

	NewTotalIteration = LoadBalanceDS.newIteration;
	if (NewTotalIteration <= 1) {
		return;
	}
	if (LoadBalanceDS.current_location->getSystem_Dynamics().U->getIsEmpty()) { //polytope U can be empty set
		U_empty = true;
	}
	int Solver = Solver_GLPK_Gurobi_GPU; //1 for CPU solver(GLPK); //2 for CPU solver(Gurobi); //3 for GPU solver(Gimplex)
//  ************* Generation of Directions Begins ***************
	unsigned int numVectors = LoadBalanceDS.reach_param.Directions.size1();
	unsigned int totalDirList1 = numVectors * (NewTotalIteration + 1); //1 extra for loop1
	math::matrix<float> List_for_X0(totalDirList1, LoadBalanceDS.reach_param.Directions.size2());
	unsigned int totalDirList2 = numVectors * NewTotalIteration; //'n' dirs for each 'n' loops
	math::matrix<float> List_for_U(totalDirList2, LoadBalanceDS.reach_param.Directions.size2());

	/*boost::timer::cpu_timer DirectionsGenerate_time;
	DirectionsGenerate_time.start();*/
	getDirectionList_X0_and_U(numCoresAvail, LoadBalanceDS.reach_param, NewTotalIteration, List_for_X0,
			List_for_U, U_empty, LoadBalanceDS.current_location->getSystem_Dynamics()); //Optimized into a single function the 2 Tasks
	/*DirectionsGenerate_time.stop();
	double wall_clock1;
	wall_clock1 = DirectionsGenerate_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
	double return_Time1 = wall_clock1 / (double) 1000;
	std::cout<< "\nDirections Generation: Boost Time taken:Wall  (in Seconds) = " << return_Time1 << std::endl;*/
//  ************* Generation of Directions Ends ***************
	/*boost::timer::cpu_timer dirCopy_time; //just to verify the Time take
	dirCopy_time.start();*/

	LoadBalanceDS.List_dir_X0 = List_for_X0;
	LoadBalanceDS.List_dir_U = List_for_U;

	/*dirCopy_time.stop();
	double wall_clock2;
	wall_clock2 = dirCopy_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
	double return_Time2 = wall_clock2 / (double) 1000;
	std::cout << "\nDirections copy: Boost Time taken:Wall  (in Seconds) = "<< return_Time2 << std::endl;*/

}



void tpbfs::loadBalancedPostD(std::vector<LoadBalanceData_PostD>& loadBalPostD){

	unsigned int Total_Length = 0;
	for (int i = 0; i < loadBalPostD.size(); i++) { //for each symbolic-states
		unsigned int sfm_len = loadBalPostD[i].sfm->getTotalIterations();	//length of each sfm
		Total_Length += sfm_len;
		loadBalPostD[i].bool_arr.resize(loadBalPostD[i].trans_size, sfm_len);	//resize(rows,col)
	}	//total_length Obtainted and bool_arrary resized
//cout<<"Total_Length = "<<Total_Length<<"\n";
#pragma omp parallel for
	for (unsigned int i = 0; i < Total_Length; i++) {
		int sfmIndex;
		unsigned int sfmColIndex;
		search_sfmIndex_colIndex(i, loadBalPostD, sfmIndex, sfmColIndex);
		polytope::ptr p;
		p = loadBalPostD[sfmIndex].sfm->getPolytope(sfmColIndex);
		std::vector<double> constraint_bound_values(loadBalPostD[sfmIndex].sfm->getInvariantDirections().size1());
		constraint_bound_values = loadBalPostD[sfmIndex].sfm->getInvariantBoundValue(sfmColIndex);
		p->setMoreConstraints(loadBalPostD[sfmIndex].sfm->getInvariantDirections(),constraint_bound_values);
		for (int trans = 0; trans < loadBalPostD[sfmIndex].trans_size; trans++) {
			loadBalPostD[sfmIndex].bool_arr(trans, sfmColIndex) = p->check_polytope_intersection(loadBalPostD[sfmIndex].guard_list[trans], lp_solver_type_choosen);//result of intersection
		}
	} //end of parallel-loop :: we have the list of intersected polys
}

/*
 * Detecting the range of sfm that are intersected with the guard
 */
void tpbfs::intersectionRangeDetection(std::vector<LoadBalanceData_PostD>& loadBalPostD) {
	int count = loadBalPostD.size();
	//cout<<"Is this same count = "<<count<<"\n";
#pragma omp parallel for
	for (int id = 0; id < count; id++) {
		for (int trans = 0; trans < loadBalPostD[id].trans_size; trans++) {
			std::list<std::pair<unsigned int, unsigned int> > intersected_range;
			bool is_intersected = false, intersection_started = false, intersection_ended = false;
			int foundIntersection = 0;
			std::pair<unsigned int, unsigned int> inte_range;
			for (unsigned int i = 0; i < loadBalPostD[id].sfm->getTotalIterations(); i++) { //sequential reading of an boolean_array that tells intersected polys
				is_intersected = loadBalPostD[id].bool_arr(trans, i);
				if (is_intersected == true) { //if intersects create a new template_polyhedra subset
					intersection_started = true;
					if (foundIntersection == 0) { //Intersection start
						foundIntersection++; //intersection_start = i;	//Intersection started at this position of 'i'
						inte_range.first = i; //	cout << "\nIntersection Found at = " << i << endl;
					}
				}
				if (intersection_started == true && is_intersected == false) {
					intersection_ended = true;
					inte_range.second = i;
					foundIntersection = 0;
				}
				if (intersection_ended == true && intersection_started == true) {
					intersected_range.push_back(inte_range);
					intersection_started = false;
					intersection_ended = false; //cout << "\nIntersection Ended at = " << i << "\n";
				}
			}
			if (intersection_started == true && intersection_ended == false) {
				inte_range.second = loadBalPostD[id].sfm->getTotalIterations() - 1;
				intersected_range.push_back(inte_range); //	cout << "\nIntersection did not End = " << i2 << "\n";
			}
			loadBalPostD[id].range_list[trans] = intersected_range;	//todo:: Verify resize problem
		}	//end of each transition or guard
	}	//end of FOR each sfm in parallel
}



void tpbfs::templateApproximation(std::vector<LoadBalanceData_PostD>& loadBalPostD){
int count = loadBalPostD.size();

#pragma omp parallel for
	for (int id=0;id<count;id++){
		for (int trans = 0; trans < loadBalPostD[id].trans_size; trans++) {
			std::list<std::pair<unsigned int, unsigned int> > range_list;
			//range_list = loadBalPostD[id].polys_list[trans];
			range_list = loadBalPostD[id].range_list[trans];
			std::list<polytope::ptr> polys;
			unsigned int poly_dir_size = loadBalPostD[id].sfm->getTotalTemplateDirections()
							+ loadBalPostD[id].sfm->getTotalInvariantDirections();
			std::vector<double> colVector(poly_dir_size);
			for (std::list<std::pair<unsigned int, unsigned int> >::iterator range_it = range_list.begin();
					range_it != range_list.end(); range_it++) {
				unsigned int start = (*range_it).first, end=(*range_it).second;
			//	cout << "first = " << start << "  second = " << end << std::endl;
				for (unsigned int eachTemplateDir=0;eachTemplateDir<loadBalPostD[id].sfm->getTotalTemplateDirections();	eachTemplateDir++){
					double Max_sf=loadBalPostD[id].sfm->getMatrixSupportFunction()(eachTemplateDir,start);
					for (int i = start+1; i <= end; i++) {
						double sf=loadBalPostD[id].sfm->getMatrixSupportFunction()(eachTemplateDir,i);
						if (sf > Max_sf)
							Max_sf = sf;
					}//end of each intersected region
					colVector[eachTemplateDir] = Max_sf;
				}//end of each template direction ALSO HAVE TO PERFORM INVARIANT DIRECTION
				unsigned int total_dir = loadBalPostD[id].sfm->getTotalTemplateDirections();
				for (unsigned int eachInvDir=0;eachInvDir<loadBalPostD[id].sfm->getTotalInvariantDirections();eachInvDir++){
					double Max_sf=loadBalPostD[id].sfm->getMatrix_InvariantBound()(eachInvDir,start);
					for (int i = start + 1; i <= end; i++) {
						double sf = loadBalPostD[id].sfm->getMatrix_InvariantBound()(eachInvDir, i);
						if (sf > Max_sf)
							Max_sf = sf;
					} //end of each intersected region
					colVector[total_dir + eachInvDir] = Max_sf;
				}
				math::matrix<double> allDirs;
				//Debug
				math::matrix<double> test;
				test =loadBalPostD[id].sfm->getTemplateDirections();
				//-----
				test.matrix_join(loadBalPostD[id].sfm->getInvariantDirections(), allDirs);
				polytope::ptr p = polytope::ptr(new polytope(allDirs,colVector,1));
				polys.push_back(p);
			}//end of multiple intersected region
			loadBalPostD[id].polys_list[trans]= polys;
		}
	}//end for each parallel

}


//Aggregrate All Flowpipe computation work into one BIG task and will run that in parallel either by multi-core CPU or GPU
//This interface is not in use now.
std::list<symbolic_states::ptr> tpbfs::computeParallelLoadBalanceReach(
		std::list<abstractCE::ptr>& ce_candidates) {

	std::list < symbolic_states::ptr > Reachability_Region; //	template_polyhedra::ptr reach_region;
	int t = 0; //0 for Read and 1 for Write

	std::vector < std::vector<pwlist::ptr> > Qpw_list(2); // QpwList[0] for read and QpwList[1] for write 	//cout << "Test 1\n";

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

	pwlist::ptr allPassedList; //so we create a permanent pwlist for storing only the passedList;
	allPassedList = pwlist::ptr(new pwlist()); //have to instantiate it
	int number_times = 0;
	bool levelcompleted = false;
	unsigned int iter_max = 1;
	int numCores = omp_get_num_procs(); //get the number of cores

	while (!isEmpty_Qpw_list(Qpw_list[t]) && (number_times <= bound)) { //	cout << "Test 5\n";
	//	cout<<"Breadth - Level === "<<number_times<<"\n";
//		boost::timer::cpu_timer t73;
//		t73.start();
		unsigned int count = getSize_Qpw_list(Qpw_list[t]); //get the size of PWList
			cout << "\nCount = " << count << "\n";
		std::vector < symbolic_states::ptr > S(count);
		vector < initial_state::ptr > list_U(count); //SubList for parallel		//	cout << "Test 6\n";
		list_U = getAllpw_list(Qpw_list, t, count, allPassedList); //All initial_state have been deleted //cout<<"Identifed pwList Done\n";
		//cout<<"list_U = "<<list_U.size();

		Qpw_list[t].resize(0);
		//cout<<"Qpw_list[t].size() = "<<Qpw_list[t].size()<<std::endl;

		Qpw_list[1 - t].resize(count); //resize to accommodate
		for (int i = 0; i < count; i++) {
			Qpw_list[1 - t][i] = pwlist::ptr(new pwlist()); //have to instantiate it
		}
//		t73.stop();
//		double clock73;
//		clock73 = t73.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
//		double return73 = clock73 / (double) 1000;
//		std::cout << "\nHow Much queue delete  Time:Wall(Seconds) = "<< return73 << std::endl;

		// ********************************* BFS Starts **********************************************************
		//	cout << "Test 7\n";
// ***************** DIRECTION and BIG_Task CREATION:: per breadth-level *****************
		std::vector < LoadBalanceData > LoadBalanceDS(count);
// ***************** DIRECTION and BIG_Task CREATION *****************

		//omp_set_nested(1);	//enables nested parallelization
#pragma omp parallel for // num_threads(count)
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
			S[id] = symbolic_states::ptr(new symbolic_states());
			S[id]->setDiscreteSet(discrete_state);
			S[id]->setParentPtrSymbolicState(U->getParentPtrSymbolicState()); //keeps track of parent pointer to symbolic_states
			S[id]->setTransitionId(U->getTransitionId()); //keeps track of originating transition_ID
			location::ptr current_location;
			current_location = H.getLocation(location_id);
			string name = current_location->getName();
			//	cout<<"Location Name = "<< name<<"\n";
			if ((name.compare("GOOD") == 0) || (name.compare("BAD") == 0)
					|| (name.compare("UNSAFE") == 0) || (name.compare("FINAL") == 0))
				continue; //do not compute the continuous reachability algorithm
			// ******************* Computing Parameters ************************ //current_location:: parameters alfa, beta and phi_trans
			//cout<<"\nBefore Compute Alfa";
			double result_alfa = compute_alfa(reach_parameter_local.time_step,
					current_location->getSystem_Dynamics(), continuous_initial_polytope, lp_solver_type_choosen); //cout<<"\nCompute Alfa Done";
			//cout<<"\nCompute alfa  = " << result_alfa;
			double result_beta = compute_beta(current_location->getSystem_Dynamics(),
					reach_parameter_local.time_step, lp_solver_type_choosen); // NO glpk object created here
			//	cout<<"\nCompute Beta  = " << result_beta;
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
			// ******************* Computing Parameters Done *******************************
// *************** POST_C computation ********** in 3 steps ***************
			bool U_empty = false;
			if (current_location->getSystem_Dynamics().U->getIsEmpty()) { //polytope U can be empty set
				U_empty = true;
			}
			LoadBalanceDS[id].X0 = continuous_initial_polytope;
			LoadBalanceDS[id].U = current_location->getSystem_Dynamics().U;
			LoadBalanceDS[id].current_location = current_location;
			LoadBalanceDS[id].symState_ID = id;
			LoadBalanceDS[id].reach_param = reach_parameter_local;
		}	//END of count FOR-LOOP  -- reach parameters computed

#pragma omp parallel for // num_threads(count)
		for (unsigned int id = 0; id < count; id++) { //separate parameters assignment with Invariant check and preLoadBalanceReachCompute task
			unsigned int NewTotalIteration;
			if (LoadBalanceDS[id].current_location->isInvariantExists()) {
				InvariantBoundaryCheck(LoadBalanceDS[id].current_location->getSystem_Dynamics(), LoadBalanceDS[id].X0,
						LoadBalanceDS[id].reach_param, LoadBalanceDS[id].current_location->getInvariant(), lp_solver_type_choosen, NewTotalIteration);
				LoadBalanceDS[id].newIteration = NewTotalIteration; //Important to take care
			}else
				LoadBalanceDS[id].newIteration = LoadBalanceDS[id].reach_param.Iterations; //Important to take care

			std::cout << "NewTotalIteration = " << NewTotalIteration << std::endl;
		}	//END of count FOR-LOOP  -- invariant boundary check done

//		double cpu_usage; //including the directionComputation + SF_computation + substitute_Algo
//		init_cpu_usage();//initializing the CPU Usage utility to start recording usages

	//	omp_set_dynamic(0);	//handles dynamic adjustment of the number of threads within a team
	//	omp_set_nested(1); //enable nested parallelism
		int numCoreAvail = 1;
		if (numCores >= count)
			numCoreAvail = numCores - count;

#pragma omp parallel for //num_threads(count)
		for (unsigned int id = 0; id < count; id++) {
			math::matrix<float> listX0, listU;
			//	cout << "Test 8\n";
			//Generation of Directions which can be a nested parallelism
			preLoadBalanceReachCompute(LoadBalanceDS[id],numCoreAvail);	// Step 1
			//	cout << "Test 9\n";
			//Todo:: if newIteration is <= 1 than do not store this flow-pipe details in LoadBalanceDS
			// and accordingly reduce count or handle while computing LPsolver
		} //END of count FOR-LOOP
		parallelLoadBalance_Task(LoadBalanceDS);
		omp_set_nested(1); //enable nested parallelism

#pragma omp parallel for // num_threads(count)
		for (unsigned int id = 0; id < count; id++) {
			S[id]->setContinuousSetptr(substitute_in_ReachAlgorithm(LoadBalanceDS[id], numCoreAvail)); // Step 3
		}
//  ********************* POST_C computation Done ********************

		bool foundUnSafe = false;
	//this being a simple task of only pushing pointer into the "Reachability_Region" so separated from inside
	//the omp parallel region to avoid having it inside the critical region
		for (unsigned int index = 0; index < count; index++) {
			if (S[index]->getContinuousSetptr()->getTotalIterations() != 0) { //computed reach_region is NOT empty
				Reachability_Region.push_back(S[index]);
			}
		}

		number_times++; //One Level or one Breadth Search over

		if (number_times > bound) {
			levelcompleted = true; //check to see how many jumps have been made(i.e., number of discrete transitions made)
		}
		if (levelcompleted || foundUnSafe) { //any true
			break; //OUT FROM WHILE LOOP 	//no need to compute rest of the locations
		}

#pragma omp parallel for // num_threads(count)
		for (int id = 0; id < count; id++) {
//  ************************************** POST_D computation Begins **********************************************************
			template_polyhedra::ptr t_poly = S[id]->getContinuousSetptr();
			if (t_poly->getTotalIterations() != 0) { //computed reach_region is empty && optimize computation
				for (std::list<transition::ptr>::iterator trans =LoadBalanceDS[id].current_location->getOut_Going_Transitions().begin();
						trans!= LoadBalanceDS[id].current_location->getOut_Going_Transitions().end(); trans++) { // get each destination_location_id and push into the pwl.waiting_list
					int transition_id = (*trans)->getTransitionId();
					location::ptr current_destination;
					Assign current_assignment;
					polytope::ptr gaurd_polytope;
					//std::list < template_polyhedra::ptr > intersected_polyhedra;
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
					std::list<polytope::ptr> polys;
					//intersected_polyhedra = t_poly->polys_intersectionParallel(gaurd_polytope, lp_solver_type_choosen); //, intersection_start_point);
					polys = t_poly->flowpipe_intersectionSequential(gaurd_polytope, lp_solver_type_choosen);
					//intersected_polyhedra =t_poly->polys_intersectionSequential(gaurd_polytope,lp_solver_type_choosen); //, intersection_start_point);
					std::cout << "Intersected = " <<polys.size() << std::endl;
					if (polys.size() > 0) { //there is intersection so new symbolic state will be inserted into the waitingList
#pragma omp critical
						{
							iter_max += polys.size();
						}
					}
					// Handle this later as In SpaceEx model we did not specified BAD or GOOD
					if ((locName.compare("BAD") == 0) || (locName.compare("GOOD") == 0) || (locName.compare("FINAL") == 0)
							|| (locName.compare("UNSAFE") == 0)) {
						continue; //do not push into the waitingList
					} //	std::cout << "Before calling getTemplate_approx\n";
					int destination_locID = (*trans)->getDestination_Location_Id();
					ds.insert_element(destination_locID);
				//	cout << "\nNumber of Intersections = "<<intersected_polyhedra.size()<<std::endl;
					for (std::list<polytope::ptr>::iterator it = polys.begin(); it != polys.end(); it++) {
						//cout << "\nNumber of Intersections #1\n";

						//intersectedRegion = (*it)->getTemplate_approx(lp_solver_type_choosen);
						intersectedRegion = (*it);
						//Returns a single over-approximated polytope from the list of intersected polytopes
						polytope::ptr newShiftedPolytope, newPolytope; //created an object here
						newPolytope = intersectedRegion->GetPolytope_Intersection(gaurd_polytope); //Retuns only the intersected region as a single newpolytope. ****** with added directions
						//std::cout << "Before calling post_assign_exact\n";
						newShiftedPolytope = post_assign_exact(newPolytope, current_assignment.Map, current_assignment.b); //initial_polytope_I = post_assign_exact(newPolytope, R, w);

						initial_state::ptr newState = initial_state::ptr(new initial_state(destination_locID, newShiftedPolytope));
						newState->setTransitionId(transition_id); // keeps track of the transition_ID
						newState->setParentPtrSymbolicState(S[id]);
//#pragma omp critical
//					{

	// Todo:: we insert all the symbolic_states into the pwlist for computing FlowPipe iff 1) and 2) holds
	// Step 1) the new current_destination.getLocID() is NOT in the Passed List and
	// Step 2) the "New Initial Polytope" or "the newShiftedPolytope" is NOT contained in the FlowPipe of LocationID of step 1
	//  as the FlowPipe in the same location will be same if the newShiftedPolytope is IN FlowPipe
	// But if Step 1 holds and Step 2 does not then it will be inserted in pwlist even if the LocationID is in Passed List
	//
						Qpw_list[1 - t][id]->WaitingList_insert(newState); //RACE CONDITION HERE
						//					}
					} //end of multiple intersected region with guard
					//cout<<"Size = "<< pwlist.getWaitingList().size()<<endl;
				} //end of multiple transaction
				  //Here we have again populated the pwlist for next-round's parallel process
			} // End-if
//  ************************************** POST_D computation Ends **********************************************************
		} //parallel for-loop
		t = 1 - t; //Switching Read/Write options for Qpw_list[1-t]
		// ************************* BFS Ends *************************************
	} //end of while loop checking waiting_list != empty
	cout << "\n **********************************************************\n";
	cout << "   *** Maximum Iterations Completed = " << iter_max << "  ***\n";
	cout << "\n **********************************************************\n";
	return Reachability_Region;
}


