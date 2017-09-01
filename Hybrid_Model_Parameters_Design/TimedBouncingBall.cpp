/*
 * TimedBouncingBall.cpp
 *
 *  Created on: 25-Nov-2014
 *      Author: amit
 */

#include "Hybrid_Model_Parameters_Design/TimedBouncingBall.h"

void SetTimedBouncingBall_ParametersOurOutput(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;
	polytope::ptr initial_polytope_I;
	polytope::ptr invariant;
	polytope::ptr gaurd_polytope;
	Dynamics system_dynamics;

	math::matrix<double> ConstraintsMatrixI, ConstraintsMatrixV,
			invariantConstraintsMatrix, gaurdConstraintsMatrix, Amatrix,
			Bmatrix;
	std::vector<double> boundValueI, boundValueV, invariantBoundValue,
			gaurdBoundValue, C;
	int boundSignI, invariantBoundSign, gaurdBoundSign, boundSignV;

	size_type row, col;

	//Polytope I Declaration in the form of Ax<=b
	//Input Polytope I as a line(bar) 10<=x(position)<=10.2 y(velocity)== 0 and t(time)==0.

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
	boundValueI[0] = 10.2;
	boundValueI[1] = -10;
	boundValueI[2] = 0;
	boundValueI[3] = 0;
	boundValueI[4] = 0;
	boundValueI[5] = 0;

	boundSignI = 1;

	// Change in dynamics system  dx/dt=+v and dv/dt = -g
	// U polytope
	row = 6;
	col = 3;
	ConstraintsMatrixV.resize(row, col);

	ConstraintsMatrixV(0, 0) = -1;
	ConstraintsMatrixV(0, 1) = 0;
	ConstraintsMatrixV(0, 2) = 0;

	ConstraintsMatrixV(1, 0) = 1;
	ConstraintsMatrixV(1, 1) = 0;
	ConstraintsMatrixV(1, 2) = 0;

	ConstraintsMatrixV(2, 0) = 0;
	ConstraintsMatrixV(2, 1) = -1;
	ConstraintsMatrixV(2, 2) = 0;

	ConstraintsMatrixV(3, 0) = 0;
	ConstraintsMatrixV(3, 1) = 1;
	ConstraintsMatrixV(3, 2) = 0;

	ConstraintsMatrixV(4, 0) = 0;
	ConstraintsMatrixV(4, 1) = 0;
	ConstraintsMatrixV(4, 2) = 1;
	ConstraintsMatrixV(5, 0) = 0;
	ConstraintsMatrixV(5, 1) = 0;
	ConstraintsMatrixV(5, 2) = -1;

	boundValueV.resize(row);
	boundValueV[0] = 0;
	boundValueV[1] = 0;
	boundValueV[2] = 1; //10;  bound for g
	boundValueV[3] = -1; //-10; bound for g
	boundValueV[4] = 1; //1		bound for t
	boundValueV[5] = -1; //-1	bound for t

	boundSignV = 1;

	row = 3;
	col = 3;
	Amatrix.resize(row, col);
	Amatrix(0, 0) = 0;
	Amatrix(0, 1) = 1;
	Amatrix(0, 2) = 0;
	Amatrix(1, 0) = 0;
	Amatrix(1, 1) = 0;
	Amatrix(1, 2) = 0;
	Amatrix(2, 0) = 0;
	Amatrix(2, 1) = 0;
	Amatrix(2, 2) = 0;

	Bmatrix.resize(row, col);
	Bmatrix(0, 0) = 1;
	Bmatrix(0, 1) = 0;
	Bmatrix(0, 2) = 0;
	Bmatrix(1, 0) = 0;
	Bmatrix(1, 1) = 1;
	Bmatrix(1, 2) = 0;
	Bmatrix(2, 0) = 0;
	Bmatrix(2, 1) = 0;
	Bmatrix(2, 2) = 1;

	// Adding the C vector of the dynamics
	C.resize(3);
	C[0] = 0;
	C[1] = -1; // g = -1
	C[2] = 1; //  t = 1
	// TO BE CHANGED LATER  WHEN DISCRETE JUMP WILL BE CHECKED
	row = 3; //4;
	col = 3;
	gaurdConstraintsMatrix.resize(row, col);
	gaurdConstraintsMatrix(0, 0) = 1;
	gaurdConstraintsMatrix(0, 1) = 0;
	gaurdConstraintsMatrix(0, 2) = 0;

	gaurdConstraintsMatrix(1, 0) = -1;
	gaurdConstraintsMatrix(1, 1) = 0;
	gaurdConstraintsMatrix(1, 2) = 0;

	gaurdConstraintsMatrix(2, 0) = 0;
	gaurdConstraintsMatrix(2, 1) = 1;
	gaurdConstraintsMatrix(2, 2) = 0;

	/*gaurdConstraintsMatrix(3, 0) = 0;
	 gaurdConstraintsMatrix(3, 1) = 1;
	 gaurdConstraintsMatrix(3, 2) = 0;*/

	gaurdBoundValue.resize(row);
	gaurdBoundValue[0] = 0; //gaurd is:: Position==0 and velocity <=0 and time <=0        0 <=x<= 0 and y<=0 and t<=0
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 0;
	//gaurdBoundValue[3] = 0;

	gaurdBoundSign = 1;

	/*	Invariant Polytope :
	 * Invariant Constraint :: position >=0 (x >= 0) conversion into -x <= -0
	 * Invariant Directions :: L1(-1,0) as positive and L2(1,0) as its negative
	 */
	row = 1;
	col = 3;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = -1;
	invariantConstraintsMatrix(0, 1) = 0;
	invariantConstraintsMatrix(0, 2) = 0;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = -0;

	invariantBoundSign = 1;
	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

//	reach_parameters.InvariantExists = true;	//false;	//Invariant exists.
	//Invariant's Directions and  Invariant polytope Initialised above
	/*
	 ************** Transition Assignment for TimedBouncing Ball *************
	 */
	math::matrix<double> R; //Transition Dynamics
	row = 3;
	col = 3;
	R.resize(row, col);
	R(0, 0) = 1;
	R(0, 1) = 0;
	R(0, 2) = 0;
	R(1, 0) = 0;
	R(1, 1) = -0.75;
	R(1, 2) = 0;
	R(2, 0) = 0;
	R(2, 1) = 0;
	R(2, 2) = 1;
	std::vector<double> w(row);
	w[0] = 0;
	w[1] = 0;
	w[2] = 0;
// ***********************************************************

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = true;
	system_dynamics.MatrixB = Bmatrix;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = C;

	//system_dynamics.U->setPolytope(ConstraintsMatrixV, boundValueV, boundSignV);
	system_dynamics.U = polytope::ptr(
			new polytope(ConstraintsMatrixV, boundValueV, boundSignV));
//	Dynamics Initalised ---------------------

	/*
	 // ******** Testing U **********
	 std::vector<double> d(3);
	 d[0]=0;
	 d[1]=0;
	 d[2]=1;
	 glpk_lp_solver lp;
	 lp.setMin_Or_Max(2);
	 lp.setConstraints(ConstraintsMatrixV, boundValueV, boundSignV);
	 double res=system_dynamics.U.computeSupportFunction(d,lp,2);
	 cout<<"Result U = "<<res<<endl;
	 */

	// Initial Polytope is initialised
	initial_polytope_I = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));
	//--------------

	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));

	Assign assignment;
	assignment.Map = R;
	assignment.b = w;

	/*transitions t;
	 t=NULL;*/
	transition::ptr trans = transition::ptr(
			new transition(1, "hop", 1, 1, gaurd_polytope, assignment));

	location::ptr source = location::ptr(new location());
	source->setLocId(1);
	source->setName("Always");
	source->setSystem_Dynamics(system_dynamics);
	source->setInvariant(invariant);
	source->setInvariantExists(true);
	source->add_Out_Going_Transition(trans);

//transitions &trans;
//transitions trans("hop",source,destination,gaurd_polytope,assignment);

	/*
	 trans.setLabel("hop");
	 trans.setSource(source);
	 trans.setDestination(destination);
	 trans.setGaurd(gaurd_polytope);
	 trans.setAssignT(assignment);
	 source.add_Adj_Transitions(trans);

	 std::list<transitions> trans_list;
	 trans_list.push_back(trans);*/

	int dim = initial_polytope_I->getSystemDimension();

	Hybrid_Automata.addInitial_Location(source);
	Hybrid_Automata.addLocation(source);
	Hybrid_Automata.setDimension(dim);
	Hybrid_Automata.insert_to_map("x", 0);
	Hybrid_Automata.insert_to_map("v", 1);
	Hybrid_Automata.insert_to_map("t", 2);
	/*
	 Hybrid_Automata.addInitLoc(source);
	 Hybrid_Automata.addTransition(trans_list);
	 Hybrid_Automata.setDimension(dim);
	 */

	unsigned int initial_location_id = 1; //the initial Location ID
	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));

	init_state_list.push_back(I);

}

//Hyst Generated output
void SetTimedBouncingBall_ParametersHystOutput(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;

	polytope::ptr initial_polytope_I;

	polytope::ptr invariant0;

	polytope::ptr gaurd_polytope0;

	Dynamics system_dynamics0;

	math::matrix<double> ConstraintsMatrixI, ConstraintsMatrixV0,
			invariantConstraintsMatrix0, gaurdConstraintsMatrix0, A0matrix,
			Bmatrix0;

	std::vector<double> boundValueI, boundValueV0, C0, invariantBoundValue0,
			gaurdBoundValue0;

	int boundSignI, invariantBoundSign, gaurdBoundSign, boundSignV;

	size_type row, col;

// The mode name is  always_running

	row = 3;
	col = 3;
	A0matrix.resize(row, col);
	A0matrix(0, 0) = 0.0;
	A0matrix(0, 1) = 1.0;
	A0matrix(0, 2) = 0.0;
	A0matrix(1, 0) = 0.0;
	A0matrix(1, 1) = 0.0;
	A0matrix(1, 2) = 0.0;
	A0matrix(2, 0) = 0.0;
	A0matrix(2, 1) = 0.0;
	A0matrix(2, 2) = 0.0;
	system_dynamics0.isEmptyMatrixA = false;
	system_dynamics0.MatrixA = A0matrix;

	system_dynamics0.isEmptyMatrixB = true;

	C0.resize(row);
	C0[0] = 0.0;
	C0[1] = -1.0;
	C0[2] = 1.0;
	system_dynamics0.isEmptyC = false;
	system_dynamics0.C = C0;

	row = 1;
	col = 3;
	invariantConstraintsMatrix0.resize(row, col);
	invariantConstraintsMatrix0(0, 0) = -1.0;
	invariantConstraintsMatrix0(0, 1) = 0.0;
	invariantConstraintsMatrix0(0, 2) = 0.0;

	invariantBoundValue0.resize(row);
	invariantBoundValue0[0] = -0.0;
	invariantBoundSign = 1;
	invariant0 = polytope::ptr(
			new polytope(invariantConstraintsMatrix0, invariantBoundValue0,
					invariantBoundSign));

	system_dynamics0.U = polytope::ptr(new polytope(true));

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
	boundValueI[0] = 10.2;
	boundValueI[1] = -10;
	boundValueI[2] = 0;
	boundValueI[3] = 0;
	boundValueI[4] = 0;
	boundValueI[5] = 0;
	boundSignI = 1;

// The transition label ishop

// Original guard: x <= ball_eps & v < 0

	row = 2;
	col = 3;

	gaurdConstraintsMatrix0.resize(row, col);
	gaurdConstraintsMatrix0(0, 0) = 1.0;
	gaurdConstraintsMatrix0(0, 1) = 0.0;
	gaurdConstraintsMatrix0(0, 2) = 0.0;
	gaurdConstraintsMatrix0(1, 0) = 0.0;
	gaurdConstraintsMatrix0(1, 1) = 1.0;
	gaurdConstraintsMatrix0(1, 2) = 0.0;

	gaurdBoundValue0.resize(row);
	gaurdBoundValue0[0] = 0.1;
	gaurdBoundValue0[1] = 0.0;
	gaurdBoundSign = 1;
	gaurd_polytope0 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix0, gaurdBoundValue0,
					gaurdBoundSign));

// The transition label is   hop

	math::matrix<double> R0;
	row = 3;
	col = 3;
	R0.resize(row, col);
	R0(0, 0) = 1.0;
	R0(0, 1) = 0.0;
	R0(0, 2) = 0.0;

	R0(1, 0) = 0.0;
	R0(1, 1) = -0.75;
	R0(1, 2) = 0.0;

	R0(2, 0) = 0.0;
	R0(2, 1) = 0.0;
	R0(2, 2) = 1.0;
	std::vector<double> w0(row);
	w0[0] = 0.0;
	w0[1] = 0.0;
	w0[2] = 0.0;

	Assign assignment0;
	assignment0.Map = R0;
	assignment0.b = w0;

	initial_polytope_I = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));

	transition::ptr t1 = transition::ptr(
			new transition(1, "hop", 1, 1, gaurd_polytope0, assignment0));

	std::list<transition::ptr> Out_Going_Trans_fromalways_running;

	Out_Going_Trans_fromalways_running.push_back(t1);
	location::ptr l1 = location::ptr(
			new location(1, "always_running", system_dynamics0, invariant0,
					true, Out_Going_Trans_fromalways_running));

	int dim = initial_polytope_I->getSystemDimension();
	Hybrid_Automata.addInitial_Location(l1);
	Hybrid_Automata.addLocation(l1);
	Hybrid_Automata.setDimension(dim);

	Hybrid_Automata.insert_to_map("x", 0);
	Hybrid_Automata.insert_to_map("v", 1);
	Hybrid_Automata.insert_to_map("t", 2);

	unsigned int initial_location_id = 1; //the initial Location ID
	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));

	init_state_list.push_back(I);
}


/*
 * Manually Testing multiple initial set
 */

void SetTimedBouncingBall_2initSet(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;

	polytope::ptr initial_polytope_I;

	polytope::ptr invariant0;

	polytope::ptr gaurd_polytope0;

	Dynamics system_dynamics0;

	math::matrix<double> ConstraintsMatrixI, ConstraintsMatrixV0,
			invariantConstraintsMatrix0, gaurdConstraintsMatrix0, A0matrix,
			Bmatrix0;

	std::vector<double> boundValueI, boundValueV0, C0, invariantBoundValue0,
			gaurdBoundValue0;

	int boundSignI, invariantBoundSign, gaurdBoundSign, boundSignV;

	size_type row, col;

// The mode name is  always_running

	row = 3;
	col = 3;
	A0matrix.resize(row, col);
	A0matrix(0, 0) = 0.0;
	A0matrix(0, 1) = 1.0;
	A0matrix(0, 2) = 0.0;
	A0matrix(1, 0) = 0.0;
	A0matrix(1, 1) = 0.0;
	A0matrix(1, 2) = 0.0;
	A0matrix(2, 0) = 0.0;
	A0matrix(2, 1) = 0.0;
	A0matrix(2, 2) = 0.0;
	system_dynamics0.isEmptyMatrixA = false;
	system_dynamics0.MatrixA = A0matrix;

	system_dynamics0.isEmptyMatrixB = true;

	C0.resize(row);
	C0[0] = 0.0;
	C0[1] = -1.0;
	C0[2] = 1.0;
	system_dynamics0.isEmptyC = false;
	system_dynamics0.C = C0;

	row = 1;
	col = 3;
	invariantConstraintsMatrix0.resize(row, col);
	invariantConstraintsMatrix0(0, 0) = -1.0;
	invariantConstraintsMatrix0(0, 1) = 0.0;
	invariantConstraintsMatrix0(0, 2) = 0.0;

	invariantBoundValue0.resize(row);
	invariantBoundValue0[0] = -0.0;
	invariantBoundSign = 1;
	invariant0 = polytope::ptr(
			new polytope(invariantConstraintsMatrix0, invariantBoundValue0,
					invariantBoundSign));

	system_dynamics0.U = polytope::ptr(new polytope(true));

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
	boundValueI[0] = 10.2;
	boundValueI[1] = -10;
	boundValueI[2] = 0;
	boundValueI[3] = 0;
	boundValueI[4] = 0;
	boundValueI[5] = 0;
	boundSignI = 1;

// The transition label ishop

// Original guard: x <= ball_eps & v < 0

	row = 2;
	col = 3;

	gaurdConstraintsMatrix0.resize(row, col);
	gaurdConstraintsMatrix0(0, 0) = 1.0;
	gaurdConstraintsMatrix0(0, 1) = 0.0;
	gaurdConstraintsMatrix0(0, 2) = 0.0;
	gaurdConstraintsMatrix0(1, 0) = 0.0;
	gaurdConstraintsMatrix0(1, 1) = 1.0;
	gaurdConstraintsMatrix0(1, 2) = 0.0;

	gaurdBoundValue0.resize(row);
	gaurdBoundValue0[0] = 0.1;
	gaurdBoundValue0[1] = 0.0;
	gaurdBoundSign = 1;
	gaurd_polytope0 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix0, gaurdBoundValue0,
					gaurdBoundSign));

// The transition label is   hop

	math::matrix<double> R0;
	row = 3;
	col = 3;
	R0.resize(row, col);
	R0(0, 0) = 1.0;
	R0(0, 1) = 0.0;
	R0(0, 2) = 0.0;

	R0(1, 0) = 0.0;
	R0(1, 1) = -0.75;
	R0(1, 2) = 0.0;

	R0(2, 0) = 0.0;
	R0(2, 1) = 0.0;
	R0(2, 2) = 1.0;
	std::vector<double> w0(row);
	w0[0] = 0.0;
	w0[1] = 0.0;
	w0[2] = 0.0;

	Assign assignment0;
	assignment0.Map = R0;
	assignment0.b = w0;

	initial_polytope_I = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));

	transition::ptr t1 = transition::ptr(
			new transition(1, "hop", 1, 1, gaurd_polytope0, assignment0));

	std::list<transition::ptr> Out_Going_Trans_fromalways_running;

	Out_Going_Trans_fromalways_running.push_back(t1);
	location::ptr l1 = location::ptr(
			new location(1, "always_running", system_dynamics0, invariant0,
					true, Out_Going_Trans_fromalways_running));

	int dim = initial_polytope_I->getSystemDimension();
	Hybrid_Automata.addInitial_Location(l1);
	Hybrid_Automata.addLocation(l1);
	Hybrid_Automata.setDimension(dim);

	Hybrid_Automata.insert_to_map("x", 0);
	Hybrid_Automata.insert_to_map("v", 1);
	Hybrid_Automata.insert_to_map("t", 2);

	unsigned int initial_location_id = 1; //the initial Location ID
	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));

	init_state_list.push_back(I);	//pushing 1st initial set object
std::cout<<"Ist initial Set created\n";
	//2nd initial set

	symbolic_states::ptr S2; //null_pointer as there is no instantiation
	polytope::ptr initial_polytope_I2;
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
	boundValueI[0] = 15.2;
	boundValueI[1] = -15;
	boundValueI[2] = 0;
	boundValueI[3] = 0;
	boundValueI[4] = 0;
	boundValueI[5] = 0;
	boundSignI = 1;
	initial_polytope_I2 = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));

	initial_state::ptr I2 = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I2, S2,
					transition_id));
	init_state_list.push_back(I2); //pushing 2nd initial set object
	std::cout<<"2nd initial Set created\n";
}

