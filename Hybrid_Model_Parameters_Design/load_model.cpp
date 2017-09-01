/*
 * load_model.cpp
 *
 *  Created on: 09-May-2016
 *      Author: rajarshi
 */

#include "Hybrid_Model_Parameters_Design/load_model.h"

void load_model(std::list<initial_state::ptr>& init_state, hybrid_automata& ha,
		userOptions& op, ReachabilityParameters& reach_parameters,
		std::pair<int, polytope::ptr>& forbidden_set) {
	unsigned int row, col;

	reach_parameters.TimeBound = op.get_timeHorizon(); //Total Time Interval
	reach_parameters.Iterations = (unsigned int) op.get_timeHorizon()
			/ op.get_timeStep(); // number of iterations
	reach_parameters.time_step = op.get_timeStep();

//Assigning the Model of the Hybrid System
//	(1,2,3,4,5,6,7) = (BBALL, TBBALL, HELICOPTER, FIVEDIMSYS, NAVIGATION, CIRCLE, CIRCLE_FOUR_LOC)
	unsigned int HybridSystem_Model_Type = op.get_model();

//Assigning Directions
	unsigned int Directions_Type = op.get_directionTemplate(); //(1,2,>2) = (BOX, OCT, UNIFORM)

	if (HybridSystem_Model_Type == BBALL) {
		SetBouncingBall_Parameters(ha, init_state, reach_parameters);
	}
	if (HybridSystem_Model_Type == TBBALL) {

		//SetTimedBouncingBall_ParametersHystOutput(ha, init_state,	reach_parameters);
		SetTimedBouncingBall_2initSet(ha, init_state, reach_parameters);

	}
	if (HybridSystem_Model_Type == HELICOPTER) {

		SetHelicopter_Parameters3(ha, init_state, reach_parameters);
	}
	if (HybridSystem_Model_Type == FIVEDIMSYS) {

		setSysParams(ha, init_state, reach_parameters);
	}

	if (HybridSystem_Model_Type == NAVIGATION_1) {

		SetNavigationBenchMark4Var(ha, init_state, reach_parameters); //Paper's Model

	}

	if (HybridSystem_Model_Type == NAVIGATION_2) {

		SetNavigationModel2(ha, init_state, reach_parameters); //My own testing Model NAV2
	}

	if (HybridSystem_Model_Type == NAVIGATION_3) {

		SetNavigationModel4(ha, init_state, reach_parameters); //Model NAV04

	}
	if (HybridSystem_Model_Type == NAVIGATION_4) {
//		SetNavigationModel5by5(ha, init_state, reach_parameters); //My own testing Model NAV_5by5
	}

	if (HybridSystem_Model_Type == NAVIGATION_5) {
		//	SetNavigationModel9by9(ha, init_state, reach_parameters); //My own testing Model NAV_9by9
	}

	if (HybridSystem_Model_Type == CIRCLE_ONE_LOC) {
		SetRotationCircleOneLocation_Parameters(ha, init_state, reach_parameters);
	}
	if (HybridSystem_Model_Type == CIRCLE_TWO_LOC) {
		//SetRotationCircle_Parameters(ha, init_state, reach_parameters);
		SetRotationTimedCircle_Parameters(ha, init_state, reach_parameters);
	}
	if (HybridSystem_Model_Type == CIRCLE_FOUR_LOC) {
		SetRotationCircle4Location_Parameters(ha, init_state, reach_parameters);
		//SetRotation_Navtimed_Parameters(ha, init_state,reach_parameters);
	}
	if (HybridSystem_Model_Type == OSCILLATOR) {
		SetOscillatorParameters(ha, init_state, reach_parameters);
	}
//std::cout <<"HybridSystem_Model_Type = "<<HybridSystem_Model_Type <<std::endl;
	if (HybridSystem_Model_Type == 14) {
		//SetConstantMotion(ha, init_state,reach_parameters);	//Call to constant dynamic Model
		//Set_NavTimed_Parameters(ha, init_state, reach_parameters);
		//user_model(ha, init_state,reach_parameters);
		//Set_NavTimed_5by5(ha, init_state, reach_parameters);
		//Set_NavTimed_9by9(ha,init_state,reach_parameters);
		//std::cout <<"Running the Test Model "<<std::endl;
		setTTEthernetModel2(ha, init_state, reach_parameters);
		//std::cout <<"Test Model Assigned"<<std::endl;
	}
	unsigned int dims=0;
	for (std::list<initial_state::ptr>::iterator it=init_state.begin();it!=init_state.end();it++){
		dims = (*it)->getInitialSet()->getSystemDimension();
	}

//Assigning the Number of Directions and Generating the Template Directions from the above given dimension of the model
//todo:: needs to decide that is this the right place to include Invariant direction
	//and also Redundant invariant directional constraints to be removed

	math::matrix<double> Real_Directions; //List of all directions
	std::vector<std::vector<double> > newDirections;

	if (Directions_Type == BOX) {
		unsigned int dir_nums = 2 * dims; //Axis Directions
		newDirections = generate_axis_directions(dims);
		get_ublas_matrix(newDirections, Real_Directions); //it returns vector vector so need to do conversion here:: Temporary solution
		row = dir_nums;
		col = dims;
		reach_parameters.Directions.resize(row, col);
		reach_parameters.Directions = Real_Directions; //Direct Assignment
	}
	if (Directions_Type == OCT) {
		unsigned int dir_nums = 2 * dims * dims; // Octagonal directions
		newDirections = get_octagonal_directions(dims);
		get_ublas_matrix(newDirections, Real_Directions); //it returns vector vector so need to do conversion here:: Temporary solution
		row = dir_nums;
		col = dims;
		reach_parameters.Directions.resize(row, col);
		reach_parameters.Directions = Real_Directions; //Direct Assignment
	}
	if (Directions_Type > 2) {
		unsigned int dir_nums = op.get_directionTemplate(); // ASSIGN HERE Number of Vectors/Directions for UNIform spear algorithm
		newDirections = math::uni_sphere(dir_nums, dims, 100, 0.0005);

		get_ublas_matrix(newDirections, Real_Directions); //it returns vector vector so need to do conversion here:: Temporary solution
		row = dir_nums;
		col = dims;
		reach_parameters.Directions.resize(row, col);
		reach_parameters.Directions = Real_Directions; //Direct Assignment
	}
	if (!op.get_forbidden_state().empty()) {
		string_to_poly(op.get_forbidden_state(), forbidden_set);
	}
}
