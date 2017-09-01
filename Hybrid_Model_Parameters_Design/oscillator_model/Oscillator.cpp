/*
 * Oscillator.cpp
 *
 *  Created on: 10-May-2016
 *      Author: amit
 *
 *
 */

#include "Hybrid_Model_Parameters_Design/oscillator_model/Oscillator.h"

/*
 * Reference for model is https://ths.rwth-aachen.de/research/projects/hypro/filtered-oscillator/
 *
 */
void SetOscillatorParameters(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;
	polytope::ptr initial_polytope_I;
	polytope::ptr invariant1, invariant2, invariant3, invariant4;
	polytope::ptr gaurd_polytope1, gaurd_polytope2, gaurd_polytope3,
			gaurd_polytope4;
	Dynamics system_dynamics;

	math::matrix<double> ConstraintsMatrixI, ConstraintsMatrixV,
			invariantConstraintsMatrix, gaurdConstraintsMatrix, Amatrix,
			Bmatrix;
	std::vector<double> boundValueI, boundValueV, invariantBoundValue,
			gaurdBoundValue;

	std::vector<double> vector_c;
	int boundSignI, invariantBoundSign, gaurdBoundSign, boundSignV;

	size_type row, col;
	// ********* constants Declaration **********
	double a1 = -2.0, a2 = -1.0, c = 0.5, x0 = 0.7, y0 = 0.7;
	// ********* constants Declaration Done **********

	unsigned int initial_location_id = 3; //the initial Location ID
// ********************* Initial Set Assignment **********************
	//Polytope I Declaration in the form of Ax<=b
	//Input Polytope I as a 0.2<=x<=0.3 and -0.1<=y<=0.1.
	row = 4;
	col = 2;
	ConstraintsMatrixI.resize(row, col);
	ConstraintsMatrixI(0, 0) = 1;
	ConstraintsMatrixI(0, 1) = 0;
	ConstraintsMatrixI(1, 0) = -1;
	ConstraintsMatrixI(1, 1) = 0;
	ConstraintsMatrixI(2, 0) = 0;
	ConstraintsMatrixI(2, 1) = 1;
	ConstraintsMatrixI(3, 0) = 0;
	ConstraintsMatrixI(3, 1) = -1;

	boundValueI.resize(row);
	boundValueI[0] = 0.3;
	boundValueI[1] = -0.2;
	boundValueI[2] = 0.1;
	boundValueI[3] = 0.1;

	boundSignI = 1;
	// Initial Polytope is initialised
	initial_polytope_I = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));
// ********************* Initial Set Assignment Done **********************
// **************************** Location ID=1 Label=np  ***************************

	//polytope U =0		not required

	//Dynamics  x' = -y and y' = x
	row = 2;
	col = 2;
	Amatrix.resize(row, col);
	Amatrix(0, 0) = -2;
	Amatrix(0, 1) = 0;
	Amatrix(1, 0) = 0;
	Amatrix(1, 1) = -1;

	vector_c.resize(row);
	vector_c[0] = 1.4;
	vector_c[1] = -0.7;

	//Common Parameters : initial polytope and dynamics
	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = true;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = vector_c;

	system_dynamics.U = polytope::ptr(new polytope(true));//system_dynamics.U->setIsEmpty(true); //set empty
//	Dynamics Initalised ---------------------
//Location 'loc1' has Transition 't1' with guard is x==0 & y + 0.714268*x >= 0 and No Assignment so its identity i.e., x'=x and y'=y
	row = 3;
	col = 2;
	gaurdConstraintsMatrix.resize(row, col);
	gaurdConstraintsMatrix(0, 0) = 1;
	gaurdConstraintsMatrix(0, 1) = 0;
	gaurdConstraintsMatrix(1, 0) = -1;
	gaurdConstraintsMatrix(1, 1) = 0;
	gaurdConstraintsMatrix(2, 0) = -0.714286;
	gaurdConstraintsMatrix(2, 1) = -1;

	gaurdBoundValue.resize(row);
	gaurdBoundValue[0] = 0;
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 0;

	gaurdBoundSign = 1;
	gaurd_polytope1 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));

//Transition Dynamics  Rx + w where R is the Assignment Mapping and w is a vector
	math::matrix<double> R;	//Transition Dynamics
	row = 2;
	col = 2;
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
//--------------
	transition::ptr t1 = transition::ptr(
			new transition(1, "t1", 1, 3, gaurd_polytope1, assignment));

//Location 1:: Invariant constraint : x<=0 &  y >= -c/x0 * x
	row = 2;
	col = 2;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = 1;
	invariantConstraintsMatrix(0, 1) = 0;
	invariantConstraintsMatrix(1, 0) = -0.714286;
	invariantConstraintsMatrix(1, 1) = -1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundValue[1] = 0;
	invariantBoundSign = 1;
	invariant1 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc1;
	Out_Going_Trans_fromLoc1.push_back(t1);

	location::ptr l1 = location::ptr(
			new location(1, "loc1", system_dynamics, invariant1, true,
					Out_Going_Trans_fromLoc1));
// ********************** Initalised for Location 1 Done **********************

// **************************** Location ID=2 Label=pp  ***************************
//Location 'loc2' has Transition 't2' with guard is x<=0 & y + 0.714268*x == 0  and No Assignment so its identity i.e., x'=x and y'=y

	//Dynamics  matrix A is common for all locations
	row = 2;
	vector_c.resize(row);
	vector_c[0] = -1.4;
	vector_c[1] = 0.7;

	//Common Parameters : initial polytope and dynamics
	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = true;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = vector_c;

	system_dynamics.U = polytope::ptr(new polytope(true));//system_dynamics.U->setIsEmpty(true); //set empty

	row = 3;
	col = 2;
	gaurdConstraintsMatrix.resize(row, col);
	gaurdConstraintsMatrix(0, 0) = 1;
	gaurdConstraintsMatrix(0, 1) = 0;
	gaurdConstraintsMatrix(1, 0) = 0.714286;
	gaurdConstraintsMatrix(1, 1) = 1;
	gaurdConstraintsMatrix(2, 0) = -0.714286;
	gaurdConstraintsMatrix(2, 1) = -1;

	gaurdBoundValue.resize(row);
	gaurdBoundValue[0] = 0;
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 0;

	gaurdBoundSign = 1;
	gaurd_polytope2 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t2 = transition::ptr(
			new transition(2, "t2", 2, 1, gaurd_polytope2, assignment));

//Location 2:: Invariant constraint : y <=0
	row = 2;
	col = 2;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = 1;
	invariantConstraintsMatrix(0, 1) = 0;
	invariantConstraintsMatrix(1, 0) = 0.714286;
	invariantConstraintsMatrix(1, 1) = 1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundValue[1] = 0;
	invariantBoundSign = 1;
	invariant2 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc2;
	Out_Going_Trans_fromLoc2.push_back(t2);

	location::ptr l2 = location::ptr(
			new location(2, "loc2", system_dynamics, invariant2, true,
					Out_Going_Trans_fromLoc2));
// ********************** Initalised for Location 2 Done **********************

// **************************** Location ID=3 Label=loc3  ***************************
//Location 'loc3' has Transition 't3' with guard is x>=0 & y + 0.714268*x == 0  and No Assignment so its identity i.e., x'=x and y'=y

	//Dynamics  x' = -y and y' = x
	row = 2;
	vector_c.resize(row);
	vector_c[0] = 1.4;
	vector_c[1] = -0.7;

	//Common Parameters : initial polytope and dynamics
	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = true;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = vector_c;

	system_dynamics.U = polytope::ptr(new polytope(true));//system_dynamics.U->setIsEmpty(true); //set empty
//	Dynamics Initalised ---------------------
			//Location 3::has transition t3::with guard is x>=0 & y<=0
	row = 3;
	col = 2;
	gaurdConstraintsMatrix.resize(row, col);
	gaurdConstraintsMatrix(0, 0) = -1;
	gaurdConstraintsMatrix(0, 1) = 0;
	gaurdConstraintsMatrix(1, 0) = 0.714286;
	gaurdConstraintsMatrix(1, 1) = 1;
	gaurdConstraintsMatrix(2, 0) = -0.714286;
	gaurdConstraintsMatrix(2, 1) = -1;

	gaurdBoundValue.resize(row);
	gaurdBoundValue[0] = 0;
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 0;
	gaurdBoundSign = 1;
	gaurd_polytope3 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t3 = transition::ptr(
			new transition(3, "t3", 3, 4, gaurd_polytope3, assignment));

	//Location 3:: Invariant constraint : x<=0 & y<=0
	row = 2;
	col = 2;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = -1;
	invariantConstraintsMatrix(0, 1) = 0;
	invariantConstraintsMatrix(1, 0) = -0.714286;
	invariantConstraintsMatrix(1, 1) = -1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundValue[1] = 0;
	invariantBoundSign = 1;
	invariant3 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc3;
	Out_Going_Trans_fromLoc3.push_back(t3);

	location::ptr l3 = location::ptr(
			new location(3, "loc3", system_dynamics, invariant3, true,
					Out_Going_Trans_fromLoc3));
// ********************** Initalised for Location 3 Done **********************

// **************************** Location ID=4 Label=loc4  ***************************
//Location 'loc4' has Transition 't4' with guard is x==0 & y + 0.714286 * x <=0 and No Assignment so its identity i.e., x'=x and y'=y
	//Dynamics
	row = 2;
	vector_c.resize(row);
	vector_c[0] = -1.4;
	vector_c[1] = 0.7;

	//Common Parameters : initial polytope and dynamics
	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = true;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = vector_c;

	system_dynamics.U = polytope::ptr(new polytope(true));//system_dynamics.U->setIsEmpty(true); //set empty
	//	Dynamics Initalised ---------------------

	row = 3;
	col = 2;
	gaurdConstraintsMatrix.resize(row, col);
	gaurdConstraintsMatrix(0, 0) = 1;
	gaurdConstraintsMatrix(0, 1) = 0;
	gaurdConstraintsMatrix(1, 0) = -1;
	gaurdConstraintsMatrix(1, 1) = 0;
	gaurdConstraintsMatrix(2, 0) = 0.714286;
	gaurdConstraintsMatrix(2, 1) = 1;

	gaurdBoundValue.resize(row);
	gaurdBoundValue[0] = 0;
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 0;
	gaurdBoundSign = 1;

	gaurd_polytope4 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t4 = transition::ptr(
			new transition(4, "t4", 4, 2, gaurd_polytope4, assignment));

	//Location 4:: Invariant constraint : x<=0 & y<=-c/x0 * x
	row = 2;
	col = 2;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = -1;
	invariantConstraintsMatrix(0, 1) = 0;
	invariantConstraintsMatrix(1, 0) = 0.714286;
	invariantConstraintsMatrix(1, 1) = 1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundValue[1] = 0;
	invariantBoundSign = 1;
	invariant4 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc4;
	Out_Going_Trans_fromLoc4.push_back(t4);

	location::ptr l4 = location::ptr(
			new location(4, "loc4", system_dynamics, invariant4, true,
					Out_Going_Trans_fromLoc4));
// ********************** Initalised for Location 3 Done **********************

	int dim = initial_polytope_I->getSystemDimension();
	Hybrid_Automata.setDimension(dim);
	Hybrid_Automata.addInitial_Location(l1);
	Hybrid_Automata.addLocation(l1);
	Hybrid_Automata.addLocation(l2);
	Hybrid_Automata.addLocation(l3);
	Hybrid_Automata.addLocation(l4);

	Hybrid_Automata.insert_to_map("x", 0);
	Hybrid_Automata.insert_to_map("y", 1);

	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));

	init_state_list.push_back(I);
}

/*This function is incomplete missing guard
 * 		The model has same dynamics matrix A for all Locations but different vector_c
 *      also with different Invariant polytope with no guard and assignment
 *      so guard polytope =
 *      and assignments  may be identity assignments
 *
 */
void SetParametersOscillator1(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;
	polytope::ptr initial_polytope_I;
	polytope::ptr invariant1, invariant2, invariant3, invariant4;
	polytope::ptr gaurd_polytope1, gaurd_polytope2, gaurd_polytope3,
			gaurd_polytope4;
	Dynamics system_dynamics;

	math::matrix<double> ConstraintsMatrixI, ConstraintsMatrixV,
			invariantConstraintsMatrix, gaurdConstraintsMatrix, Amatrix,
			Bmatrix;
	std::vector<double> boundValueI, boundValueV, invariantBoundValue,
			gaurdBoundValue;

	std::vector<double> vector_c;
	int boundSignI, invariantBoundSign, gaurdBoundSign, boundSignV;

	size_type row, col;
	// ********* constants Declaration **********
	double a1 = -2.0, a2 = -1.0, c = 0.5, x0 = 0.7, y0 = 0.7;
	// ********* constants Declaration Done **********

	unsigned int initial_location_id = 1; //the initial Location ID
// ********************* Initial Set Assignment **********************
	//Polytope I Declaration in the form of Ax<=b
	//Input Polytope I as a x==0 and -0.1<=y<=0.1.
	row = 4;
	col = 2;
	ConstraintsMatrixI.resize(row, col);
	ConstraintsMatrixI(0, 0) = 1;
	ConstraintsMatrixI(0, 1) = 0;
	ConstraintsMatrixI(1, 0) = -1;
	ConstraintsMatrixI(1, 1) = 0;
	ConstraintsMatrixI(2, 0) = 0;
	ConstraintsMatrixI(2, 1) = 1;
	ConstraintsMatrixI(3, 0) = 0;
	ConstraintsMatrixI(3, 1) = -1;

	boundValueI.resize(row);
	boundValueI[0] = 0;
	boundValueI[1] = 0;
	boundValueI[2] = 0.1;
	boundValueI[3] = 0.1;

	boundSignI = 1;
	// Initial Polytope is initialised
	initial_polytope_I = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));
// ********************* Initial Set Assignment Done **********************
// **************************** Location ID=1 Label=np  ***************************

	//polytope U =0		not required

	//Dynamics  x' = -y and y' = x
	row = 2;
	col = 2;
	Amatrix.resize(row, col);
	Amatrix(0, 0) = a1;
	Amatrix(0, 1) = 0;
	Amatrix(1, 0) = 0;
	Amatrix(1, 1) = a2;

	vector_c.resize(row);
	vector_c[0] = -1 * a1 * x0;
	vector_c[1] = a2 * y0;

	//Common Parameters : initial polytope and dynamics
	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = true;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = vector_c;

	system_dynamics.U = polytope::ptr(new polytope(true));//system_dynamics.U->setIsEmpty(true); //set empty
//	Dynamics Initalised ---------------------
//Location 'np' has Transition 'hop' with guard is ..... and No Assignment so its identity i.e., x'=x and y'=y
	row = 2;
	col = 2;
	/* Find out what is the guard for Oscillator
	 * 	gaurdConstraintsMatrix.resize(row, col);
	 gaurdConstraintsMatrix(0, 0) = 1;
	 gaurdConstraintsMatrix(0, 1) = 0;
	 gaurdConstraintsMatrix(1, 0) = 0;
	 gaurdConstraintsMatrix(1, 1) = -1;
	 gaurdBoundValue.resize(row);
	 gaurdBoundValue[0] = 0;
	 gaurdBoundValue[1] = 0;
	 gaurdBoundSign = 1;
	 gaurd_polytope1 = polytope::ptr(new polytope(gaurdConstraintsMatrix, gaurdBoundValue, gaurdBoundSign));*/

//Transition Dynamics  Rx + w where R is the Assignment Mapping and w is a vector
	math::matrix<double> R;	//Transition Dynamics
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
//--------------
	transition::ptr t1 = transition::ptr(
			new transition(1, "hop", 1, 2, gaurd_polytope1, assignment));

//Location 1:: Invariant constraint : x<=0 &  y >= -c/x0 * x
	row = 2;
	col = 2;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = 1;
	invariantConstraintsMatrix(0, 1) = 0;
	invariantConstraintsMatrix(1, 0) = -1 * c / x0;
	invariantConstraintsMatrix(1, 1) = -1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundValue[1] = 0;
	invariantBoundSign = 1;
	invariant1 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc1;
	Out_Going_Trans_fromLoc1.push_back(t1);

	location::ptr l1 = location::ptr(
			new location(1, "np", system_dynamics, invariant1, true,
					Out_Going_Trans_fromLoc1));
// ********************** Initalised for Location 1 Done **********************

// **************************** Location ID=2 Label=pp  ***************************
//Location 'pp' has Transition 'hop' with guard is ..... and No Assignment so its identity i.e., x'=x and y'=y

	//Dynamics  matrix A is common for all locations
	//vector_c is same for Locations 'np' and 'pp'

	//Common Parameters : initial polytope and dynamics
	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = true;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = vector_c;

	system_dynamics.U = polytope::ptr(new polytope(true));//system_dynamics.U->setIsEmpty(true); //set empty

	row = 2;
	col = 2;
	/*gaurdConstraintsMatrix.resize(row, col);
	 gaurdConstraintsMatrix(0, 0) = 1;
	 gaurdConstraintsMatrix(0, 1) = 0;
	 gaurdConstraintsMatrix(1, 0) = 0;
	 gaurdConstraintsMatrix(1, 1) = 1;
	 gaurdBoundValue.resize(row);
	 gaurdBoundValue[0] = 0;
	 gaurdBoundValue[1] = 0;
	 gaurdBoundSign = 1;*/
	gaurd_polytope2 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t2 = transition::ptr(
			new transition(2, "hop", 2, 3, gaurd_polytope2, assignment));

//Location 2:: Invariant constraint : y <=0
	row = 2;
	col = 2;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = 1;
	invariantConstraintsMatrix(0, 1) = 0;
	invariantConstraintsMatrix(1, 0) = -1 * c / x0;
	invariantConstraintsMatrix(1, 1) = -1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundValue[1] = 0;
	invariantBoundSign = 1;
	invariant2 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc2;
	Out_Going_Trans_fromLoc2.push_back(t2);

	location::ptr l2 = location::ptr(
			new location(2, "pp", system_dynamics, invariant2, true,
					Out_Going_Trans_fromLoc2));
// ********************** Initalised for Location 2 Done **********************

// **************************** Location ID=3 Label=pn  ***************************
//Location 'pn' has Transition 'hop' with guard is ..... and No Assignment so its identity i.e., x'=x and y'=y

	//Dynamics  x' = -y and y' = x
	row = 2;
	vector_c.resize(row);
	vector_c[0] = a1 * x0;
	vector_c[1] = -1 * a2 * y0;

	//Common Parameters : initial polytope and dynamics
	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = true;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = vector_c;

	system_dynamics.U = polytope::ptr(new polytope(true));//system_dynamics.U->setIsEmpty(true); //set empty
//	Dynamics Initalised ---------------------
			//Location 3::has transition t3::with guard is x>=0 & y<=0
	row = 2;
	col = 2;
	/*gaurdConstraintsMatrix.resize(row, col);
	 gaurdConstraintsMatrix(0, 0) = -1;
	 gaurdConstraintsMatrix(0, 1) = 0;
	 gaurdConstraintsMatrix(1, 0) = 0;
	 gaurdConstraintsMatrix(1, 1) = 1;
	 gaurdBoundValue.resize(row);
	 gaurdBoundValue[0] = 0;
	 gaurdBoundValue[1] = 0;
	 gaurdBoundSign = 1;*/
	gaurd_polytope3 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t3 = transition::ptr(
			new transition(3, "hop", 3, 4, gaurd_polytope3, assignment));

	//Location 3:: Invariant constraint : x<=0 & y<=0
	row = 2;
	col = 2;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = -1;
	invariantConstraintsMatrix(0, 1) = 0;
	invariantConstraintsMatrix(1, 0) = c / x0;
	invariantConstraintsMatrix(1, 1) = 1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundValue[1] = 0;
	invariantBoundSign = 1;
	invariant3 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc3;
	Out_Going_Trans_fromLoc3.push_back(t3);

	location::ptr l3 = location::ptr(
			new location(3, "pn", system_dynamics, invariant3, true,
					Out_Going_Trans_fromLoc3));
// ********************** Initalised for Location 3 Done **********************

// **************************** Location ID=4 Label=nn  ***************************
//Location 'nn' has Transition 'hop' with guard is ..... and No Assignment so its identity i.e., x'=x and y'=y
	//Dynamics  x' = -y and y' = x
	//vector_c is same for location 'pn' and 'nn'

	//Common Parameters : initial polytope and dynamics
	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = true;

	system_dynamics.isEmptyC = false;
	system_dynamics.C = vector_c;

	system_dynamics.U = polytope::ptr(new polytope(true));//system_dynamics.U->setIsEmpty(true); //set empty
	//	Dynamics Initalised ---------------------

	row = 2;
	col = 2;
	/*gaurdConstraintsMatrix.resize(row, col);
	 gaurdConstraintsMatrix(0, 0) = -1;
	 gaurdConstraintsMatrix(0, 1) = 0;
	 gaurdConstraintsMatrix(1, 0) = 0;
	 gaurdConstraintsMatrix(1, 1) = -1;

	 gaurdBoundValue.resize(row);
	 gaurdBoundValue[0] = 0;
	 gaurdBoundValue[1] = 0;
	 gaurdBoundSign = 1;*/
	gaurd_polytope4 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t4 = transition::ptr(
			new transition(4, "hop", 4, 1, gaurd_polytope4, assignment));

	//Location 4:: Invariant constraint : x<=0 & y<=-c/x0 * x
	row = 2;
	col = 2;
	invariantConstraintsMatrix.resize(row, col);
	invariantConstraintsMatrix(0, 0) = 1;
	invariantConstraintsMatrix(0, 1) = 0;
	invariantConstraintsMatrix(1, 0) = c / x0;
	invariantConstraintsMatrix(1, 1) = 1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 0;
	invariantBoundValue[1] = 0;
	invariantBoundSign = 1;
	invariant4 = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));

	std::list<transition::ptr> Out_Going_Trans_fromLoc4;
	Out_Going_Trans_fromLoc4.push_back(t4);

	location::ptr l4 = location::ptr(
			new location(4, "nn", system_dynamics, invariant4, true,
					Out_Going_Trans_fromLoc4));
// ********************** Initalised for Location 3 Done **********************

	int dim = initial_polytope_I->getSystemDimension();
	Hybrid_Automata.setDimension(dim);
	Hybrid_Automata.addInitial_Location(l1);
	Hybrid_Automata.addLocation(l1);
	Hybrid_Automata.addLocation(l2);
	Hybrid_Automata.addLocation(l3);
	Hybrid_Automata.addLocation(l4);

	Hybrid_Automata.insert_to_map("x", 0);
	Hybrid_Automata.insert_to_map("y", 1);

	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));

	init_state_list.push_back(I);
}

