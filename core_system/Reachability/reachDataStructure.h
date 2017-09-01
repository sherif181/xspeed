/*
 * reachDataStructure.h
 *
 *  Created on: 28-Feb-2016
 *      Author: amit
 */

#ifndef REACHDATASTRUCTURE_H_
#define REACHDATASTRUCTURE_H_

#include "core_system/HybridAutomata/Location.h"

//Data structure used for POST_C computation in Load Balancing Algorithm
struct LoadBalanceData{

	unsigned int symState_ID;	// ID of the symbolic_state or flowpipe whose details are stored. #Redundant actually
	// As every symbolic_state will have different
	polytope::ptr X0;	//initial polytope for the specific symbolic_state
	polytope::ptr U;	//non-deterministic input polytope for the specific symbolic_state
	math::matrix<float> List_dir_X0;//directions for computing Support Function on X0 (or Initial) polytope
	math::matrix<float> List_dir_U;	//directions for computing Support Function on U (or non-deterministic input) polytope
	std::vector<double> sf_X0;	//Result Support Function on polytope X0 for List_dir_X0;
	std::vector<double> sf_U;	//Result Support Function on polytope U for List_dir_U;

	std::vector<double> sf_UnitBall;	//Result Support Function of UnitBall for List_dir_X0;
	std::vector<double> sf_dotProduct;	//Result dotProduct of List_dir_X0 and C (dynamic's expression + C);

	location::ptr current_location;
//location current_location;	//placed here as it is comfortable here
	unsigned int newIteration;	//placed here as it is comfortable here
	ReachabilityParameters reach_param;	//placed here as it is comfortable here
};


//Data structure used for POST_D computation in Load Balancing Algorithm
struct LoadBalanceData_PostD{
	template_polyhedra::ptr sfm;
	unsigned int trans_size;
	std::vector<polytope::ptr> guard_list;	//std::vector may be better for implementation
	std::vector<Assign> assign_list;
	std::vector<int> dest_locID;	//Destination_locationID
	std::vector<int> trans_ID;	//transition_ID

	//std::vector<std::vector<bool> > bool_arr; //Race Condition Occurred

	math::matrix<bool> bool_arr; //each row is the guard and elements are booleans indicating intersected or not

	std::vector<std::list<std::pair<unsigned int, unsigned int> > > range_list;	//each guard can generate
			//more than One intersected region. Vector of guards. with each guard generating list of intersected range

	std::vector<std::list<polytope::ptr> > polys_list;	//for each guards list of tempalate_Approximation as polytope
	//Each transition will have destination_loc, guard and assignments
};


#endif /* REACHDATASTRUCTURE_H_ */
