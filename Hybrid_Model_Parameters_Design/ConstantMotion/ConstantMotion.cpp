/*
 * ConstantMotion.cpp
 *
 *  Created on: 06-Jun-2016
 *      Author: amit
 */

#include "Hybrid_Model_Parameters_Design/ConstantMotion/ConstantMotion.h"

void SetConstantMotion(hybrid_automata& Hybrid_Automata,
		initial_state::ptr& init_state,
		ReachabilityParameters& reach_parameters) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;
	polytope::ptr initial_polytope_I, invariant, gaurd_polytope;
	Dynamics system_dynamics;

	math::matrix<double> ConstraintsMatrixI, ConstraintsMatrixV,
			invariantConstraintsMatrix, gaurdConstraintsMatrix, Amatrix,
			Bmatrix;
	std::vector<double> boundValueI, boundValueV, invariantBoundValue,
			gaurdBoundValue, C;
	int boundSignI, invariantBoundSign, gaurdBoundSign, boundSignV;

	size_type row, col;
	unsigned int initial_location_id; //the initial Location ID

	//Polytope I Declaration in the form of Ax<=b
	//Input Polytope I as a point(x1,x2,t) (0.5 <=x1<=0.5,0.2<=x2<=0.2,t==0) in the grid of cells
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

	boundValueI.resize(row);
	// ********************* start_location=1:: (0.5<=x1<=0.5, 0.2<=x2<=0.2, t==0) ************************
	initial_location_id = 1; //the initial Location ID

	boundValueI[0] = 0.5; //
	boundValueI[1] = -0.5;
	boundValueI[2] = 0.2;
	boundValueI[3] = -0.2;
	boundValueI[4] = 0;
	boundValueI[5] = 0;

	boundSignI = 1;
	initial_polytope_I = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));
	//initial_polytope_I.setPolytope(ConstraintsMatrixI, boundValueI, boundSignI);

	row = 3;
	col = 3;
	Amatrix.resize(row, col);
	Amatrix(0, 0) = 0;
	Amatrix(0, 1) = 0;
	Amatrix(0, 2) = 0;

	Amatrix(1, 0) = 0;
	Amatrix(1, 1) = 0;
	Amatrix(1, 2) = 0;

	Amatrix(2, 0) = 0;
	Amatrix(2, 1) = 0;
	Amatrix(2, 2) = 0;

	Bmatrix = Amatrix;

	math::matrix<double> R; //Transition Dynamics
	row = 3;
	col = 3;
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

// ***********************************************************

	/*	*************** Initialization of all transition *******************
	 *  List of transition are t1, t2, ... , t20 including transition towards the Locations labelled "A" and "B"
	 *  where Label "A" is the "Final location" to be reached and "B" the "Bad location" to be avoided.
	 */
	row = 6;
	col = 3;
	gaurdConstraintsMatrix.resize(row, col); //this matrix will be common for all transition except the gaurdBoundValue.
	gaurdConstraintsMatrix(0, 0) = 1;
	gaurdConstraintsMatrix(0, 1) = 0;
	gaurdConstraintsMatrix(0, 2) = 0;

	gaurdConstraintsMatrix(1, 0) = -1;
	gaurdConstraintsMatrix(1, 1) = 0;
	gaurdConstraintsMatrix(1, 2) = 0;

	gaurdConstraintsMatrix(2, 0) = 0;
	gaurdConstraintsMatrix(2, 1) = 1;
	gaurdConstraintsMatrix(2, 2) = 0;

	gaurdConstraintsMatrix(3, 0) = 0;
	gaurdConstraintsMatrix(3, 1) = -1;
	gaurdConstraintsMatrix(3, 2) = 0;

	gaurdConstraintsMatrix(4, 0) = 0;
	gaurdConstraintsMatrix(4, 1) = 0;
	gaurdConstraintsMatrix(4, 2) = 1;

	gaurdConstraintsMatrix(5, 0) = 0;
	gaurdConstraintsMatrix(5, 1) = 0;
	gaurdConstraintsMatrix(5, 2) = -1;

	gaurdBoundSign = 1;

	gaurdBoundValue.resize(row);
	gaurdBoundValue[0] = 1;	 // x1==1 and 0<=x2<=1 and  -100<=t<=100
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = 0;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;

	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));

	transition::ptr t1 = transition::ptr(
			new transition(1, "1 to 2", 1, 2, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; // x2==1 and 0<=x1<=1 and  -100<=t<=100
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;

	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t2 = transition::ptr(
			new transition(2, "1 to 6", 1, 6, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // x1==2 and 0<=x2<=1 and  -100<=t<=100
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = 0;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;

	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t3 = transition::ptr(
			new transition(3, "2 to 3", 2, 3, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // x2==1 and 1<=x1<=2 and  -100<=t<=100
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t4 = transition::ptr(
			new transition(4, "2 to 5", 2, 5, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; // x1==1 and 0<=x2<=1 and  -100<=t<=100
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = 0;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t5 = transition::ptr(
			new transition(5, "2 to 1", 2, 1, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 3; // x2==1 and 2<=x1<=3 and  -100<=t<=100
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t6 = transition::ptr(
			new transition(6, "3 to 4", 3, 4, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // x1==2 and 0<=x2<=1 and  -100<=t<=100
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = 0;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;

	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t7 = transition::ptr(
			new transition(7, "3 to 2", 3, 2, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // x1==2 and 1<=x2<=2 and  -100<=t<=100
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;
//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t8 = transition::ptr(
			new transition(8, "4 to 5", 4, 5, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 3; // x2==1 and 2<=x1<=3 and  -100<=t<=100
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t9 = transition::ptr(
			new transition(9, "4 to 3", 4, 3, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // x1==2 and 1<=x2<=2 and  -100<=t<=100
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t10 = transition::ptr(
			new transition(10, "5 to 4", 5, 4, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; // x1==1 and 1<=x2<=2 and  -100<=t<=100
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t11 = transition::ptr(
			new transition(11, "5 to 6", 5, 6, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // x2==1 and 1<=x1<=2 and  -100<=t<=100
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t12 = transition::ptr(
			new transition(12, "5 to 2", 5, 2, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; // x1==1 and 1<=x2<=2 and  -100<=t<=100
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t13 = transition::ptr(
			new transition(13, "6 to 5", 6, 5, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; // x2==1 and 0<=x1<=1 and  -100<=t<=100
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 100;
	gaurdBoundValue[5] = 100;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t14 = transition::ptr(
			new transition(14, "6 to 1", 6, 1, gaurd_polytope, assignment));
// ******************* Transition initialized **************************

	/*	*************** Initialization of all Locations *******************
	 *  List of Locations are l1, l2, ... , l6
	 */
	row = 4;
	col = 3;
	invariantConstraintsMatrix.resize(row, col); //Common for all polytope except the invariantBoundValue.
	invariantConstraintsMatrix(0, 0) = 1;
	invariantConstraintsMatrix(0, 1) = 0;
	invariantConstraintsMatrix(0, 2) = 0;

	invariantConstraintsMatrix(1, 0) = -1;
	invariantConstraintsMatrix(1, 1) = 0;
	invariantConstraintsMatrix(1, 2) = 0;

	invariantConstraintsMatrix(2, 0) = 0;
	invariantConstraintsMatrix(2, 1) = 1;
	invariantConstraintsMatrix(2, 2) = 0;

	invariantConstraintsMatrix(3, 0) = 0;
	invariantConstraintsMatrix(3, 1) = -1;
	invariantConstraintsMatrix(3, 2) = 0;

	invariantBoundSign = 1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 1; //0<=x1<=1 and 0<=x2<=1
	invariantBoundValue[1] = 0;
	invariantBoundValue[2] = 1;
	invariantBoundValue[3] = 0;

	system_dynamics.isEmptyMatrixA = true;
	system_dynamics.MatrixA = Amatrix;
	system_dynamics.isEmptyMatrixB = true;
	system_dynamics.MatrixB = Bmatrix;

	row = 3;
	C.resize(row);
	C[0] = 1;
	C[1] = 1;
	C[2] = 1;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = C;

	system_dynamics.U = polytope::ptr(new polytope(true));

	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc1;
	Out_Going_Trans_fromLoc1.push_back(t1);
	Out_Going_Trans_fromLoc1.push_back(t2);

	location::ptr l1 = location::ptr(
			new location(1, "loc1", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc1));
//  ************ Location ID=1 completed  ************
	row = 4;
	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 2; //1<=x1<=2 and 0<=x2<=1
	invariantBoundValue[1] = -1;
	invariantBoundValue[2] = 1;
	invariantBoundValue[3] = 0;

	system_dynamics.isEmptyMatrixA = true;
	system_dynamics.MatrixA = Amatrix;
	system_dynamics.isEmptyMatrixB = true;
	system_dynamics.MatrixB = Bmatrix;

	row = 3;
	C.resize(row);
	C[0] = 1;
	C[1] = 0;
	C[2] = 1;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = C;

	system_dynamics.U = polytope::ptr(new polytope(true));

	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc2;
	Out_Going_Trans_fromLoc2.push_back(t3);
	Out_Going_Trans_fromLoc2.push_back(t4);
	Out_Going_Trans_fromLoc2.push_back(t5);

	location::ptr l2 = location::ptr(
			new location(2, "loc2", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc2));
	//  ************ Location ID=2 completed  ************

	row = 4;
	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 3; //2<=x1<=3 and 0<=x2<=1
	invariantBoundValue[1] = -2;
	invariantBoundValue[2] = 1;
	invariantBoundValue[3] = 0;

	system_dynamics.isEmptyMatrixA = true;
	system_dynamics.MatrixA = Amatrix;
	system_dynamics.isEmptyMatrixB = true;
	system_dynamics.MatrixB = Bmatrix;

	row = 3;
	C.resize(row);
	C[0] = 1;
	C[1] = 1;
	C[2] = 1;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = C;

	system_dynamics.U = polytope::ptr(new polytope(true));

	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc3;
	Out_Going_Trans_fromLoc3.push_back(t6);
	Out_Going_Trans_fromLoc3.push_back(t7);

	location::ptr l3 = location::ptr(
			new location(3, "loc3", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc3));
	//  ************ Location ID=3 completed  ************

	row = 4;
	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 3; //2<=x1<=3 and 1<=x2<=2
	invariantBoundValue[1] = -2;
	invariantBoundValue[2] = 2;
	invariantBoundValue[3] = -1;

	system_dynamics.isEmptyMatrixA = true;
	system_dynamics.MatrixA = Amatrix;
	system_dynamics.isEmptyMatrixB = true;
	system_dynamics.MatrixB = Bmatrix;

	row = 3;
	C.resize(row);
	C[0] = -1;
	C[1] = 1;
	C[2] = 1;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = C;

	system_dynamics.U = polytope::ptr(new polytope(true));

	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc4;
	Out_Going_Trans_fromLoc4.push_back(t8);
	Out_Going_Trans_fromLoc4.push_back(t9);

	location::ptr l4 = location::ptr(
			new location(4, "loc4", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc4));
	//  ************ Location ID=4 completed  ************

	row = 4;
	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 2; //1<=x1<=2 and 1<=x2<=2
	invariantBoundValue[1] = -1;
	invariantBoundValue[2] = 2;
	invariantBoundValue[3] = -1;

	system_dynamics.isEmptyMatrixA = true;
	system_dynamics.MatrixA = Amatrix;
	system_dynamics.isEmptyMatrixB = true;
	system_dynamics.MatrixB = Bmatrix;

	row = 3;
	C.resize(row);
	C[0] = -1;
	C[1] = 0;
	C[2] = 1;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = C;

	system_dynamics.U = polytope::ptr(new polytope(true));

	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc5;
	Out_Going_Trans_fromLoc5.push_back(t10);
	Out_Going_Trans_fromLoc5.push_back(t11);
	Out_Going_Trans_fromLoc5.push_back(t12);

	location::ptr l5 = location::ptr(
			new location(5, "loc5", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc5));
	//  ************ Location ID=5 completed  ************

	row = 4;
	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 1; //0<=x1<=1 and 1<=x2<=2
	invariantBoundValue[1] = 0;
	invariantBoundValue[2] = 2;
	invariantBoundValue[3] = -1;

	system_dynamics.isEmptyMatrixA = true;
	system_dynamics.MatrixA = Amatrix;
	system_dynamics.isEmptyMatrixB = true;
	system_dynamics.MatrixB = Bmatrix;

	row = 3;
	C.resize(row);
	C[0] = -1;
	C[1] = 0;
	C[2] = 1;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = C;

	system_dynamics.U = polytope::ptr(new polytope(true));

	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc6;
	Out_Going_Trans_fromLoc6.push_back(t13);
	Out_Going_Trans_fromLoc6.push_back(t14);

	location::ptr l6 = location::ptr(
			new location(6, "loc6", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc6));
	//  ************ Location ID=6 completed  ************

	//	*************** Locations Initialized *******************

	//	*************** Initialization The Hybrid Automata *******************
	int dim = initial_polytope_I->getSystemDimension();

	Hybrid_Automata.setDimension(dim);
	Hybrid_Automata.addInitial_Location(l1);
	Hybrid_Automata.addLocation(l1);
	Hybrid_Automata.addLocation(l2);
	Hybrid_Automata.addLocation(l3);
	Hybrid_Automata.addLocation(l4);
	Hybrid_Automata.addLocation(l5);
	Hybrid_Automata.addLocation(l6);

	Hybrid_Automata.insert_to_map("x1", 0);
	Hybrid_Automata.insert_to_map("x2", 1);
	Hybrid_Automata.insert_to_map("t", 2);

	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));

	init_state = I;

}
