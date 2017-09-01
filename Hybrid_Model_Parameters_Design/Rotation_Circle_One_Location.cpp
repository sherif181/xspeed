/*
 * Rotation_Circle_One_Location.cpp
 *
 *  Created on: 02-Mar-2015
 *      Author: amit
 */
/*
 * There exists only one Location with no invariant and transition.
 */

#include "Hybrid_Model_Parameters_Design/Rotation_Circle_One_Location.h"

void SetRotationCircleOneLocation_Parameters(hybrid_automata& Hybrid_Automata,
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
	/*	ConstraintsMatrixI(0, 0) = 1;
	 ConstraintsMatrixI(0, 1) = 0;
	 ConstraintsMatrixI(1, 0) = -1;
	 ConstraintsMatrixI(1, 1) = 0;
	 ConstraintsMatrixI(2, 0) = 0;
	 ConstraintsMatrixI(2, 1) = 1;
	 ConstraintsMatrixI(3, 0) = 0;
	 ConstraintsMatrixI(3, 1) = -1;*/

	ConstraintsMatrixI(0, 0) = -1;
	ConstraintsMatrixI(0, 1) = -1;
	ConstraintsMatrixI(1, 0) = 1;
	ConstraintsMatrixI(1, 1) = 1;
	ConstraintsMatrixI(2, 0) = 1;
	ConstraintsMatrixI(2, 1) = -1;
	ConstraintsMatrixI(3, 0) = -1;
	ConstraintsMatrixI(3, 1) = 1;

	boundValueI.resize(row); //Input Polytope I see my note copy
	boundValueI[0] = -50; //-2;	//-5
	boundValueI[1] = 53; //5;	//13
	boundValueI[2] = 2; //5
	boundValueI[3] = 1; //3

	/*boundValueI.resize(row);		//Input Polytope I as a x==0.5 and y==0.5.
	 boundValueI[0] = 0.5;
	 boundValueI[1] = -0.5;
	 boundValueI[2] = 0.5;
	 boundValueI[3] = -0.5;*/

	/*boundValueI.resize(row);		//Input Polytope I as a x==0.5 and y==0.5.
	 boundValueI[0] = 0.5;
	 boundValueI[1] = -0.4;
	 boundValueI[2] = 0.5;
	 boundValueI[3] = -0.4;*/

	/*boundValueI.resize(row);		//Input Polytope I as a 4<=x<=5 and 4<=y<=5
	 boundValueI[0] = 5;
	 boundValueI[1] = -4;
	 boundValueI[2] = 5;
	 boundValueI[3] = -4;*/

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

	/*
	 //Location 1:: gaurd is y<=0 and No Assignment so its identity i.e., x'=x and y'=y
	 row = 1;
	 col = 2;
	 gaurdConstraintsMatrix.resize(row, col);
	 gaurdConstraintsMatrix(0, 0) = 0;
	 gaurdConstraintsMatrix(0, 1) = 1;

	 gaurdBoundValue.resize(row);
	 gaurdBoundValue[0] = 0;

	 gaurdBoundSign = 1;

	 gaurd_polytope1 = polytope::ptr(new polytope(gaurdConstraintsMatrix, gaurdBoundValue, gaurdBoundSign));
	 transitions t1(1, "T1", 1, 2, gaurd_polytope1, assignment);
	 */

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
	//Out_Going_Trans_fromLoc1.push_back(t1);

//NO INVARIANT EXITS
	location::ptr l1 = location::ptr(
			new location(1, "Loc-1", system_dynamics, invariant1, false,
					Out_Going_Trans_fromLoc1));
//	Initalised for Location 1	 ---------------------

	int dim = initial_polytope_I->getSystemDimension();
	Hybrid_Automata.setDimension(dim);
	Hybrid_Automata.addInitial_Location(l1);
	Hybrid_Automata.addLocation(l1);
//	Hybrid_Automata.addLocation(l2);

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

