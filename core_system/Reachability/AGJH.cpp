/*
 * AGJH.cpp
 *
 *  Created on: 28-Oct-2016
 *      Author: amit
 */

#include "AGJH.h"

// Holzmann algorithm adaption
std::list<symbolic_states::ptr> agjh::ParallelBFS_GH()
{
	std::list < symbolic_states::ptr > Reachability_Region; //	template_polyhedra::ptr reach_region;
	int t = 0; //0 for Read and 1 for Write
	int N = CORE; // Number of cores in the machine
	std::list<initial_state::ptr> Wlist[2][N][N];
	for(unsigned int i=0;i<2;i++){
		for(unsigned int j=0;j<N;j++){
			for(unsigned int k=0;k<N;k++){
				std::list<initial_state::ptr> l;
				Wlist[i][j][k] = l;
			}
		}
	}
	int initialState_index = 0; // first initial state
	for (std::list<initial_state::ptr>::iterator i = I.begin(); i != I.end(); i++) {

		Wlist[t][initialState_index][0].push_back(*(i));

		initialState_index++; //next initial state
	}	//Wlist[t][0][0].push_back(this->I); //old implementation


	std::list<symbolic_states::ptr> PASSED;
	unsigned int level = 0;
	bool done;
	//srand(time(NULL));
	int count=1;	//for initial symbolic state
	unsigned int previousBreadth=1;
	cout << "\nNumber of Flowpipes to be Computed (per Breadths) = " << count<< "\n";
	do {
		done = true;
#pragma omp parallel for
		for(unsigned int i=0;i<N;i++){
			for(unsigned int q=0;q<N;q++){
				for(std::list<initial_state::ptr>::iterator it = Wlist[t][i][q].begin(); it!=Wlist[t][i][q].end();it++){
					initial_state::ptr s = (*it);
					template_polyhedra::ptr R;
					R = postC(s);

					symbolic_states::ptr R1 = symbolic_states::ptr(new symbolic_states());
					R1->setContinuousSetptr(R);
					discrete_set d;
					d.insert_element(s->getLocationId());
					R1->setDiscreteSet(d);
#pragma omp critical
				{
					PASSED.push_back(R1);
				}
					//----end of postC on s

				//Currently removed the Safety Check Section from here

					std::list<initial_state::ptr> R2;
					if (level < bound){	//Check level to avoid last PostD computation
						R2 = postD(R1);
#pragma omp critical
					{
					printf("found %d new states" , R2.size());
					 count = count + R2.size();
					}
						//cout <<"postD size = " <<R2.size()<<std::endl;
						std::vector<int> ArrCores(N);	//If this is done only once outside time saved but race-condition?
						for (int id=0;id<N;id++){
							ArrCores[id] = id;	//sequential insertion into the array
						}
						int startVal=0;
						std::random_shuffle(ArrCores.begin(), ArrCores.end());
						for(std::list<initial_state::ptr>::iterator its = R2.begin();its != R2.end();its++){
							initial_state::ptr next_s = (*its);
							/* generate a random number between 1 and N: */
							int w_rand = ArrCores[startVal];	//ArrCores.pop_back();	//Removes from the back
							startVal++;	//sequentially traversing the shuffed arrary
							if (startVal==N) //traversing reached end this occurs if sym_states more than the available cores
								startVal = 0; //so start traversing again from the beginning of the shuffed arrary.
							//std::cout << "Selected Core = " << w_rand << std::endl;
							Wlist[1-t][w_rand][i].push_back(next_s);
						}
					}	//end of level check
					// call postD with res of prev step
				}
					Wlist[t][i][q].clear();
			}
		}
		// barrier synchronization
		for(unsigned int i=0;i<N;i++){
			for(unsigned int j=0;j<N;j++){
				if(!Wlist[1-t][i][j].empty()){
					done = false;
					break;
				}
			}
			if(!done)
				break;
		}
		unsigned int temp = count - previousBreadth;
		if (level < bound){
			cout << "\nNumber of Flowpipes to be Computed (per Breadths) = " << temp << "\n";
		}
		previousBreadth = count;

		level++;
		t = 1 - t;
	}while(level<=bound && !done);
cout<<"******************************************************************\n";
cout <<"Maximum number of Symbolic states Passed = " <<count<<"\n";
cout<<"******************************************************************\n";
	return PASSED;

}


template_polyhedra::ptr agjh::postC(initial_state::ptr s){
	int location_id;

	discrete_set d;
	location_id = s->getLocationId();
	d.insert_element(location_id); //creating discrete_state

	//continuous_initial_polytope = U.getInitial_ContinousSetptr();
	polytope::ptr continuous_initial_polytope;
	continuous_initial_polytope = s->getInitialSet();

	ReachabilityParameters reach_parameters; //local
	//cout<<"No error 111!!!!\n";
	reach_parameters = this->reach_parameters;
//	cout<<"No error 2222!!!!\n";
	reach_parameters.X0 = continuous_initial_polytope;
	symbolic_states::ptr S;

	location::ptr current_location;

	current_location = H.getLocation(location_id);
	string name = current_location->getName();

	double result_alfa = compute_alfa(reach_parameters.time_step,
			current_location->getSystem_Dynamics(), continuous_initial_polytope,
			lp_solver_type_choosen); //2 glpk object created here
	double result_beta = compute_beta(current_location->getSystem_Dynamics(),
			reach_parameters.time_step, lp_solver_type_choosen); // NO glpk object created here

	reach_parameters.result_alfa = result_alfa;
	reach_parameters.result_beta = result_beta;
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

	// ******************* Computing Parameters *******************************

	unsigned int NewTotalIteration = reach_parameters.Iterations;
	if (current_location->isInvariantExists()) {

		if (current_location->getSystem_Dynamics().isEmptyMatrixB==true && current_location->getSystem_Dynamics().isEmptyC==true){
			//Approach of Coarse-time-step and Fine-time-step
			jumpInvariantBoundaryCheck(current_location->getSystem_Dynamics(), continuous_initial_polytope, reach_parameters,
									current_location->getInvariant(), lp_solver_type_choosen, NewTotalIteration);
		}else{
			//Approach of Sequential invariant check will work for all case
			InvariantBoundaryCheck(current_location->getSystem_Dynamics(), continuous_initial_polytope,
					 reach_parameters, current_location->getInvariant(), lp_solver_type_choosen, NewTotalIteration);
		}

		/*jumpInvariantBoundaryCheck(current_location->getSystem_Dynamics(), continuous_initial_polytope, reach_parameters,
						current_location->getInvariant(), lp_solver_type_choosen, NewTotalIteration);*/

		/*quickInvariantBoundaryCheck(current_location->getSystem_Dynamics(),
				continuous_initial_polytope, reach_parameters,
				current_location->getInvariant(), lp_solver_type_choosen,
				NewTotalIteration);*/
	}

	template_polyhedra::ptr reach_region;
//	std::cout << "NewTotalIteration = " << NewTotalIteration << std::endl;

	//parallelReachSelection(NewTotalIteration, current_location, continuous_initial_polytope, reach_parameters, reach_region, id);
	reach_region = reachabilitySequential(NewTotalIteration,
					current_location->getSystem_Dynamics(),
					continuous_initial_polytope, reach_parameters,
					current_location->getInvariant(),
					current_location->isInvariantExists(), lp_solver_type_choosen);
	/*sequentialReachSelection(NewTotalIteration, current_location,
			continuous_initial_polytope, reach_region);*/

	return reach_region;
}
std::list<initial_state::ptr> agjh::postD(symbolic_states::ptr symb)
{
	template_polyhedra::ptr reach_region= symb->getContinuousSetptr();
	int locId = *(symb->getDiscreteSet().getDiscreteElements().begin());

	location::ptr current_location = H.getLocation(locId);
	std::list<initial_state::ptr> res;

	if (reach_region->getTotalIterations() != 0) { //computed reach_region is empty && optimize transition BreadthLevel-wise
	//	cout<<"1\n";
		for (std::list<transition::ptr>::iterator t = current_location->getOut_Going_Transitions().begin();
				t != current_location->getOut_Going_Transitions().end(); t++) { // get each destination_location_id and push into the pwl.waiting_list
			int transition_id = (*t)->getTransitionId();

			location::ptr current_destination;
			Assign current_assignment;
			polytope::ptr gaurd_polytope;
			//std::list < template_polyhedra::ptr > intersected_polyhedra;
			polytope::ptr intersectedRegion;//created two objects here
			discrete_set ds;
			current_destination = H.getLocation((*t)->getDestination_Location_Id());

			string locName = current_destination->getName();
			gaurd_polytope = (*t)->getGaurd();//	GeneratePolytopePlotter(gaurd_polytope);

		//cout<<"2\n";
			std::list<polytope::ptr> polys;
			//intersected_polyhedra = reach_region->polys_intersectionSequential(gaurd_polytope, lp_solver_type_choosen); //, intersection_start_point);
			polys = reach_region->flowpipe_intersectionSequential(gaurd_polytope, lp_solver_type_choosen);
			printf("poly size %d\n",polys.size());

		//cout<<"3\n";
			//Todo to make is even procedure with Sequential procedure.... so intersection is done first and then decide to skip this loc
			if ((locName.compare("BAD") == 0) || (locName.compare("GOOD") == 0)
					|| (locName.compare("FINAL") == 0) || (locName.compare("UNSAFE") == 0))
			continue;//do not push into the waitingList

			current_assignment = (*t)->getAssignT();
			// *** interesected_polyhedra included with invariant_directions also ******
		//	cout<<"size = "<< intersected_polyhedra.size();

			int destination_locID = (*t)->getDestination_Location_Id();
			ds.insert_element(destination_locID);

			for (std::list<polytope::ptr>::iterator i = polys.begin(); i != polys.end(); i++) {
				intersectedRegion = (*i);
				//intersectedRegion = (*i)->getTemplate_approx(lp_solver_type_choosen);
				//Returns a single over-approximated polytope from the list of intersected polytopes
				//	GeneratePolytopePlotter(intersectedRegion);
				polytope::ptr newShiftedPolytope, newPolytope;//created an object here
				newPolytope = intersectedRegion->GetPolytope_Intersection(gaurd_polytope);//Retuns the intersected region as a single newpolytope. **** with added directions
				//newShiftedPolytope = post_assign_exact(newPolytope, current_assignment.Map, current_assignment.b);//initial_polytope_I = post_assign_exact(newPolytope, R, w);

				math::matrix<double> test(current_assignment.Map.size1(),
						current_assignment.Map.size2());
				if (current_assignment.Map.inverse(test))	//invertible?
				{
					//std::cout << "Exact Post Assignment\n";
					newShiftedPolytope = post_assign_exact(newPolytope,
							current_assignment.Map, current_assignment.b);
				} else {
					//std::cout << "Approximate Post Assignment\n";
					newShiftedPolytope = post_assign_approx_deterministic(
							newPolytope, current_assignment.Map,
							current_assignment.b, reach_parameters.Directions,
							lp_solver_type_choosen);
				}

				//	newShiftedPolytope->print2file(newInitSet,0,1); //printing the New Initial Set
				initial_state::ptr newState = initial_state::ptr(new initial_state(destination_locID, newShiftedPolytope));
				newState->setTransitionId(transition_id);// keeps track of the transition_ID
				res.push_back(newState);
			} //end of multiple intersected region with guard

		} //end of multiple transaction
	} //end-if of flowpipe computation not empty
	return res;

} //end of while loop checking waiting_list != empty

