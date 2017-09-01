/*
 * Rotation_Circle.cpp
 *
 *  Created on: 04-Dec-2014
 *      Author: amit
 */
//Model with Two Locations First the Upper-Half and Second the Lower-Half
#include "Hybrid_Model_Parameters_Design/Rotation_Circle.h"

void SetRotationCircle_Parameters(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;
	polytope::ptr initial_polytope_I;
	polytope::ptr invariant1, invariant2;
	polytope::ptr gaurd_polytope1, gaurd_polytope2;
	Dynamics system_dynamics;

	math::matrix<double> ConstraintsMatrixI, ConstraintsMatrixV,
			invariantConstraintsMatrix, gaurdConstraintsMatrix, Amatrix,
			Bmatrix;
	std::vector<double> boundValueI, boundValueV, invariantBoundValue,
			gaurdBoundValue;
	int boundSignI, invariantBoundSign, gaurdBoundSign, boundSignV;

	size_type row, col;

	//Polytope I Declaration in the form of Ax<=b
	//Input Polytope I as a x==5 and y==5.

	row = 4;
	col = 2;
	ConstraintsMatrixI.resize(row, col);

	/*ConstraintsMatrixI(0, 0) = 1;
	 ConstraintsMatrixI(0, 1) = 0;
	 ConstraintsMatrixI(1, 0) = -1;
	 ConstraintsMatrixI(1, 1) = 0;

	 ConstraintsMatrixI(2, 0) = 0;
	 ConstraintsMatrixI(2, 1) = 1;
	 ConstraintsMatrixI(3, 0) = 0;
	 ConstraintsMatrixI(3, 1) = -1;

	 boundValueI.resize(row); //Input Polytope I as a x==0.5 and y==0.5.
	 boundValueI[0] = 0.5;
	 boundValueI[1] = -0.5;
	 boundValueI[2] = 0.5;
	 boundValueI[3] = -0.5;
	 */
	/*
	 * Input set for Model 9
	 * //Input Polytope I: good example to show containment or inclusion of Flowpipe into another
	 **/
	ConstraintsMatrixI(0, 0) = -1;
	ConstraintsMatrixI(0, 1) = -1;
	ConstraintsMatrixI(1, 0) = 1;
	ConstraintsMatrixI(1, 1) = 1;
	ConstraintsMatrixI(2, 0) = 1;
	ConstraintsMatrixI(2, 1) = -1;
	ConstraintsMatrixI(3, 0) = -1;
	ConstraintsMatrixI(3, 1) = 1;
	boundValueI.resize(row);
	boundValueI[0] = -50; //-2;	//-5
	boundValueI[1] = 53; //5;	//13
	boundValueI[2] = 2; //5
	boundValueI[3] = 1; //3

	/*
	 * This data is used for plotting circle for Invariant Check Algorithm with Coarse and Fine time-step

	 ConstraintsMatrixI(0, 0) = 1;
	 ConstraintsMatrixI(0, 1) = 0;
	 ConstraintsMatrixI(1, 0) = -1;
	 ConstraintsMatrixI(1, 1) = 0;
	 ConstraintsMatrixI(2, 0) = 0;
	 ConstraintsMatrixI(2, 1) = 1;
	 ConstraintsMatrixI(3, 0) = 0;
	 ConstraintsMatrixI(3, 1) = -1;
	 boundValueI.resize(row);
	 boundValueI[0] = 27.5;
	 boundValueI[1] = -24.5;
	 boundValueI[2] = 27;
	 boundValueI[3] = -24;
	 */
	boundSignI = 1;

	//polytope U =0		not required
	ConstraintsMatrixV.resize(row, col);
	ConstraintsMatrixV(0, 0) = 0;
	ConstraintsMatrixV(0, 1) = -1;
	ConstraintsMatrixV(1, 0) = 0;
	ConstraintsMatrixV(1, 1) = 1;
	ConstraintsMatrixV(2, 0) = -1;
	ConstraintsMatrixV(2, 1) = 0;
	ConstraintsMatrixV(3, 0) = 1;
	ConstraintsMatrixV(3, 1) = 0;

	boundValueV.resize(row);
	boundValueV[0] = 0;
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;

	boundSignV = 1;

	//Dynamics  x' = -y and y' = x
	row = 2;
	col = 2;
	Amatrix.resize(row, col);
	Amatrix(0, 0) = 0;
	Amatrix(0, 1) = -1;
	Amatrix(1, 0) = 1;
	Amatrix(1, 1) = 0;

	Bmatrix.resize(row, col);
	Bmatrix(0, 0) = 1;
	Bmatrix(0, 1) = 0;
	Bmatrix(1, 0) = 0;
	Bmatrix(1, 1) = 1;

//Transition Dynamics  Rx + w where R is the Assignment Mapping and w is a vector
	math::matrix<double> R; //Transition Dynamics
	R.resize(row, col);
	R(0, 0) = 1;
	R(0, 1) = 0;
	R(1, 0) = 0;
	R(1, 1) = 1;

	std::vector<double> w(row);
	w[0] = 0;
	w[1] = 0;

	Assign assignment;
	assignment.Map = R;
	assignment.b = w;

//Common Parameters : initial polytope and dynamics

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = false;
	system_dynamics.MatrixB = Bmatrix;

	system_dynamics.isEmptyC = true;
	system_dynamics.U = polytope::ptr(new polytope());
	system_dynamics.U->setIsEmpty(true); //set empty

//	Dynamics Initalised ---------------------
	// Initial Polytope is initialised
	initial_polytope_I = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));
//--------------

//Location 1:: gaurd is y<=0 and No Assignment so its identity i.e., x'=x and y'=y
	row = 1;
	col = 2;
	gaurdConstraintsMatrix.resize(row, col);
	gaurdConstraintsMatrix(0, 0) = 0;
	gaurdConstraintsMatrix(0, 1) = 1;

	gaurdBoundValue.resize(row);
	gaurdBoundValue[0] = 0;

	gaurdBoundSign = 1;

	gaurd_polytope1 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t1 = transition::ptr(
			new transition(1, "T1", 1, 2, gaurd_polytope1, assignment));

//Location 1:: Invariant constraint : y >=0
	row = 1;
	col = 2;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = 0;
	invariantConstraintsMatrix(0, 1) = -1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundSign = 1;
	invariant1 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc1;
	Out_Going_Trans_fromLoc1.push_back(t1);

	location::ptr l1 = location::ptr(
			new location(1, "Loc-1", system_dynamics, invariant1, true,
					Out_Going_Trans_fromLoc1));
//	Initalised for Location 1	 ---------------------

//Location 2:: gaurd is y>=0 and No Assignment so its identity i.e., x'=x and y'=y
	row = 1;
	col = 2;
	gaurdConstraintsMatrix.resize(row, col);
	gaurdConstraintsMatrix(0, 0) = 0;
	gaurdConstraintsMatrix(0, 1) = -1;

	gaurdBoundValue.resize(row);
	gaurdBoundValue[0] = 0;
	gaurdBoundSign = 1;
	gaurd_polytope2 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t2 = transition::ptr(
			new transition(2, "T2", 2, 1, gaurd_polytope2, assignment));

//Location 2:: Invariant constraint : y <=0
	row = 1;
	col = 2;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = 0;
	invariantConstraintsMatrix(0, 1) = 1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundSign = 1;
	invariant2 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc2;
	Out_Going_Trans_fromLoc2.push_back(t2);

	location::ptr l2 = location::ptr(
			new location(2, "Loc-2", system_dynamics, invariant2, true,
					Out_Going_Trans_fromLoc2));
//Initialised Location 2	--------------------------

	int dim = initial_polytope_I->getSystemDimension();
	Hybrid_Automata.setDimension(dim);
	Hybrid_Automata.addInitial_Location(l1);
	Hybrid_Automata.addLocation(l1);
	Hybrid_Automata.addLocation(l2);

	Hybrid_Automata.insert_to_map("x", 0);
	Hybrid_Automata.insert_to_map("y", 1);

	unsigned int initial_location_id = 1; //the initial Location ID
	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));
	init_state_list.push_back(I);

}

void SetRotationTimedCircle_Parameters(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;
	polytope::ptr initial_polytope_I;
	polytope::ptr invariant1, invariant2;
	polytope::ptr gaurd_polytope1, gaurd_polytope2;
	Dynamics system_dynamics;

	math::matrix<double> ConstraintsMatrixI, ConstraintsMatrixV,
			invariantConstraintsMatrix, gaurdConstraintsMatrix, Amatrix,
			Bmatrix;
	std::vector<double> boundValueI, C0, boundValueV, invariantBoundValue,
			gaurdBoundValue;
	int boundSignI, invariantBoundSign, gaurdBoundSign, boundSignV;

	size_type row, col;

	//Polytope I Declaration in the form of Ax<=b

	row = 6;
	col = 3;
	ConstraintsMatrixI.resize(row, col);
	ConstraintsMatrixI(0, 0) = 1;
	ConstraintsMatrixI(0, 1) = 0;
	ConstraintsMatrixI(0, 2) = 0;
	ConstraintsMatrixI(1, 0) = -1;
	ConstraintsMatrixI(1, 1) = 0;
	ConstraintsMatrixI(1, 2) = 0;
	ConstraintsMatrixI(2, 0) = 0;
	ConstraintsMatrixI(2, 1) = 1;
	ConstraintsMatrixI(2, 2) = 0;
	ConstraintsMatrixI(3, 0) = 0;
	ConstraintsMatrixI(3, 1) = -1;
	ConstraintsMatrixI(3, 2) = 0;
	ConstraintsMatrixI(4, 0) = 0;
	ConstraintsMatrixI(4, 1) = 0;
	ConstraintsMatrixI(4, 2) = 1;
	ConstraintsMatrixI(5, 0) = 0;
	ConstraintsMatrixI(5, 1) = 0;
	ConstraintsMatrixI(5, 2) = -1;

	boundValueI.resize(row); //Input Polytope I as a x==0.5 and y==0.5 and t==0.
	boundValueI[0] = 0.5;
	boundValueI[1] = -0.5;
	boundValueI[2] = 0.5;
	boundValueI[3] = -0.5;
	boundValueI[4] = 0;
	boundValueI[5] = 0;

	boundSignI = 1;

	//polytope U =0		not required
	ConstraintsMatrixV.resize(row, col);
	ConstraintsMatrixV(0, 0) = 0;
	ConstraintsMatrixV(0, 1) = -1;
	ConstraintsMatrixV(1, 0) = 0;
	ConstraintsMatrixV(1, 1) = 1;
	ConstraintsMatrixV(2, 0) = -1;
	ConstraintsMatrixV(2, 1) = 0;
	ConstraintsMatrixV(3, 0) = 1;
	ConstraintsMatrixV(3, 1) = 0;

	boundValueV.resize(row);
	boundValueV[0] = 0;
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;

	boundSignV = 1;

	//Dynamics  x' = -y and y' = x and t'=1
	row = 3;
	col = 3;
	Amatrix.resize(row, col);
	Amatrix(0, 0) = 0;
	Amatrix(0, 1) = -1;
	Amatrix(0, 2) = 0;
	Amatrix(1, 0) = 1;
	Amatrix(1, 1) = 0;
	Amatrix(1, 2) = 0;
	Amatrix(2, 0) = 0;
	Amatrix(2, 1) = 0;
	Amatrix(2, 2) = 0;

	system_dynamics.isEmptyMatrixB = true;

	C0.resize(row);
	C0[0] = 0.0;
	C0[1] = 0.0;
	C0[2] = 1.0;
	system_dynamics.isEmptyC = false;
	system_dynamics.C = C0;

//Transition Dynamics  Rx + w where R is the Assignment Mapping and w is a vector
	math::matrix<double> R; //Transition Dynamics
	R.resize(row, col);
	R(0, 0) = 1;
	R(0, 1) = 0;
	R(0, 2) = 0;
	R(1, 0) = 0;
	R(1, 1) = 1;
	R(1, 2) = 0;
	R(2, 0) = 0;
	R(2, 1) = 0;
	R(2, 2) = 1;

	std::vector<double> w(row);
	w[0] = 0;
	w[1] = 0;
	w[2] = 0;

	Assign assignment;
	assignment.Map = R;
	assignment.b = w;

//Common Parameters : initial polytope and dynamics

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.U = polytope::ptr(new polytope());
	system_dynamics.U->setIsEmpty(true); //set empty

//	Dynamics Initalised ---------------------
	// Initial Polytope is initialised
	initial_polytope_I = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));
//--------------

//Location 1:: gaurd is y<=0 and No Assignment so its identity i.e., x'=x and y'=y
	row = 1;
	col = 3;
	gaurdConstraintsMatrix.resize(row, col);
	gaurdConstraintsMatrix(0, 0) = 0;
	gaurdConstraintsMatrix(0, 1) = 1;
	gaurdConstraintsMatrix(0, 2) = 0;

	gaurdBoundValue.resize(row);
	gaurdBoundValue[0] = 0;

	gaurdBoundSign = 1;

	gaurd_polytope1 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t1 = transition::ptr(
			new transition(1, "T1", 1, 2, gaurd_polytope1, assignment));

//Location 1:: Invariant constraint : y >=0
	row = 1;
	col = 3;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = 0;
	invariantConstraintsMatrix(0, 1) = -1;
	invariantConstraintsMatrix(0, 2) = 0;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundSign = 1;
	invariant1 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc1;
	Out_Going_Trans_fromLoc1.push_back(t1);

	location::ptr l1 = location::ptr(
			new location(1, "Loc-1", system_dynamics, invariant1, true,
					Out_Going_Trans_fromLoc1));
//	Initalised for Location 1	 ---------------------

//Location 2:: gaurd is y>=0 and No Assignment so its identity i.e., x'=x and y'=y
	row = 1;
	col = 3;
	gaurdConstraintsMatrix.resize(row, col);
	gaurdConstraintsMatrix(0, 0) = 0;
	gaurdConstraintsMatrix(0, 1) = -1;
	gaurdConstraintsMatrix(0, 2) = 0;

	gaurdBoundValue.resize(row);
	gaurdBoundValue[0] = 0;
	gaurdBoundSign = 1;
	gaurd_polytope2 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t2 = transition::ptr(
			new transition(2, "T2", 2, 1, gaurd_polytope2, assignment));

//Location 2:: Invariant constraint : y <=0
	row = 1;
	col = 3;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = 0;
	invariantConstraintsMatrix(0, 1) = 1;
	invariantConstraintsMatrix(0, 2) = 0;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundSign = 1;
	invariant2 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc2;
	Out_Going_Trans_fromLoc2.push_back(t2);

	location::ptr l2 = location::ptr(
			new location(2, "Loc-2", system_dynamics, invariant2, true,
					Out_Going_Trans_fromLoc2));
//Initialised Location 2	--------------------------

	int dim = initial_polytope_I->getSystemDimension();
	Hybrid_Automata.setDimension(dim);
	Hybrid_Automata.addInitial_Location(l1);
	Hybrid_Automata.addLocation(l1);
	Hybrid_Automata.addLocation(l2);
	Hybrid_Automata.insert_to_map("x", 0);
	Hybrid_Automata.insert_to_map("y", 1);
	Hybrid_Automata.insert_to_map("t", 2);

	unsigned int initial_location_id = 1; //the initial Location ID
	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));
	init_state_list.push_back(I);
}
