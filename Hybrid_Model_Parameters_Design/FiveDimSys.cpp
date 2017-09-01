/*
 * FiveDimSys.cpp
 *
 *  Created on: 25-Nov-2014
 *      Author: amit
 */

#include "Hybrid_Model_Parameters_Design/FiveDimSys.h"

void setSysParams(hybrid_automata& Hybrid_Automata,
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
			gaurdBoundValue;
	int boundSignI, invariantBoundSign, gaurdBoundSign, boundSignV;

	size_type row, col;

	//Polytope I Declaration in the form of Ax<=b
	//Input Polytope I as a line(bar) 10<=x(position)<=10.2 y(velocity)== 0.
	row = 10;
	col = 5;
	ConstraintsMatrixI.resize(row, col);
	ConstraintsMatrixI(0, 0) = 1;
	ConstraintsMatrixI(0, 1) = 0;
	ConstraintsMatrixI(0, 2) = 0;
	ConstraintsMatrixI(0, 3) = 0;
	ConstraintsMatrixI(0, 4) = 0;

	ConstraintsMatrixI(1, 0) = -1;
	ConstraintsMatrixI(1, 1) = 0;
	ConstraintsMatrixI(1, 2) = 0;
	ConstraintsMatrixI(1, 3) = 0;
	ConstraintsMatrixI(1, 4) = 0;

	ConstraintsMatrixI(2, 0) = 0;
	ConstraintsMatrixI(2, 1) = 1;
	ConstraintsMatrixI(2, 2) = 0;
	ConstraintsMatrixI(2, 3) = 0;
	ConstraintsMatrixI(2, 4) = 0;

	ConstraintsMatrixI(3, 0) = 0;
	ConstraintsMatrixI(3, 1) = -1;
	ConstraintsMatrixI(3, 2) = 0;
	ConstraintsMatrixI(3, 3) = 0;
	ConstraintsMatrixI(3, 4) = 0;

	ConstraintsMatrixI(4, 0) = 0;
	ConstraintsMatrixI(4, 1) = 0;
	ConstraintsMatrixI(4, 2) = 1;
	ConstraintsMatrixI(4, 3) = 0;
	ConstraintsMatrixI(4, 4) = 0;

	ConstraintsMatrixI(5, 0) = 0;
	ConstraintsMatrixI(5, 1) = 0;
	ConstraintsMatrixI(5, 2) = -1;
	ConstraintsMatrixI(5, 3) = 0;
	ConstraintsMatrixI(5, 4) = 0;

	ConstraintsMatrixI(6, 0) = 0;
	ConstraintsMatrixI(6, 1) = 0;
	ConstraintsMatrixI(6, 2) = 0;
	ConstraintsMatrixI(6, 3) = 1;
	ConstraintsMatrixI(6, 4) = 0;

	ConstraintsMatrixI(7, 0) = 0;
	ConstraintsMatrixI(7, 1) = 0;
	ConstraintsMatrixI(7, 2) = 0;
	ConstraintsMatrixI(7, 3) = -1;
	ConstraintsMatrixI(7, 4) = 0;

	ConstraintsMatrixI(8, 0) = 0;
	ConstraintsMatrixI(8, 1) = 0;
	ConstraintsMatrixI(8, 2) = 0;
	ConstraintsMatrixI(8, 3) = 0;
	ConstraintsMatrixI(8, 4) = 1;

	ConstraintsMatrixI(9, 0) = 0;
	ConstraintsMatrixI(9, 1) = 0;
	ConstraintsMatrixI(9, 2) = 0;
	ConstraintsMatrixI(9, 3) = 0;
	ConstraintsMatrixI(9, 4) = -1;

	boundValueI.resize(10);
	boundValueI[0] = 1.01;
	boundValueI[1] = -0.99;
	boundValueI[2] = 0.01;
	boundValueI[3] = 0.01;
	boundValueI[4] = 0.01;
	boundValueI[5] = 0.01;
	boundValueI[6] = 0.01;
	boundValueI[7] = 0.01;
	boundValueI[8] = 0.01;
	boundValueI[9] = 0.01;

	boundSignI = 1;

//	initial_polytope_I.setPolytope(ConstraintsMatrixI, boundValueI,boundSignI);
	initial_polytope_I = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));

	row = 5;
	col = 5;
	Amatrix.resize(row, col);
	Amatrix(0, 0) = -4.8082;
	Amatrix(0, 1) = -4.0336;
	Amatrix(0, 2) = 1.5610;
	Amatrix(0, 3) = -1.9059;
	Amatrix(0, 4) = -0.7299;

	Amatrix(1, 0) = -1.4137;
	Amatrix(1, 1) = -7.2405;
	Amatrix(1, 2) = -1.4615;
	Amatrix(1, 3) = -3.4139;
	Amatrix(1, 4) = 4.4378;

	Amatrix(2, 0) = 7.1423;
	Amatrix(2, 1) = 12.4581;
	Amatrix(2, 2) = 0.3480;
	Amatrix(2, 3) = 8.6344;
	Amatrix(2, 4) = -6.9834;

	Amatrix(3, 0) = -0.0349;
	Amatrix(3, 1) = 2.3922;
	Amatrix(3, 2) = 0.8171;
	Amatrix(3, 3) = 0.0909;
	Amatrix(3, 4) = -3.5798;

	Amatrix(4, 0) = -6.6104;
	Amatrix(4, 1) = -14.5070;
	Amatrix(4, 2) = 1.1341;
	Amatrix(4, 3) = -7.1915;
	Amatrix(4, 4) = 1.6098;

//row= ;	col=
	Bmatrix.resize(row, col);
	for (unsigned int i = 0; i < row; i++)
		for (unsigned int j = 0; j < col; j++)
			if (i == j)
				Bmatrix(i, j) = 1;
			else
				Bmatrix(i, j) = 0;

	//  * here polytope U == 0
	row = 10;
	col = 5;
	ConstraintsMatrixV.resize(row, col);
	ConstraintsMatrixV(0, 0) = 1;
	ConstraintsMatrixV(0, 1) = 0;
	ConstraintsMatrixV(0, 2) = 0;
	ConstraintsMatrixV(0, 3) = 0;
	ConstraintsMatrixV(0, 4) = 0;

	ConstraintsMatrixV(1, 0) = -1;
	ConstraintsMatrixV(1, 1) = 0;
	ConstraintsMatrixV(1, 2) = 0;
	ConstraintsMatrixV(1, 3) = 0;
	ConstraintsMatrixV(1, 4) = 0;

	ConstraintsMatrixV(2, 0) = 0;
	ConstraintsMatrixV(2, 1) = 1;
	ConstraintsMatrixV(2, 2) = 0;
	ConstraintsMatrixV(2, 3) = 0;
	ConstraintsMatrixV(2, 4) = 0;

	ConstraintsMatrixV(3, 0) = 0;
	ConstraintsMatrixV(3, 1) = -1;
	ConstraintsMatrixV(3, 2) = 0;
	ConstraintsMatrixV(3, 3) = 0;
	ConstraintsMatrixV(3, 4) = 0;

	ConstraintsMatrixV(4, 0) = 0;
	ConstraintsMatrixV(4, 1) = 0;
	ConstraintsMatrixV(4, 2) = 1;
	ConstraintsMatrixV(4, 3) = 0;
	ConstraintsMatrixV(4, 4) = 0;

	ConstraintsMatrixV(5, 0) = 0;
	ConstraintsMatrixV(5, 1) = 0;
	ConstraintsMatrixV(5, 2) = -1;
	ConstraintsMatrixV(5, 3) = 0;
	ConstraintsMatrixV(5, 4) = 0;

	ConstraintsMatrixV(6, 0) = 0;
	ConstraintsMatrixV(6, 1) = 0;
	ConstraintsMatrixV(6, 2) = 0;
	ConstraintsMatrixV(6, 3) = 1;
	ConstraintsMatrixV(6, 4) = 0;

	ConstraintsMatrixV(7, 0) = 0;
	ConstraintsMatrixV(7, 1) = 0;
	ConstraintsMatrixV(7, 2) = 0;
	ConstraintsMatrixV(7, 3) = -1;
	ConstraintsMatrixV(7, 4) = 0;

	ConstraintsMatrixV(8, 0) = 0;
	ConstraintsMatrixV(8, 1) = 0;
	ConstraintsMatrixV(8, 2) = 0;
	ConstraintsMatrixV(8, 3) = 0;
	ConstraintsMatrixV(8, 4) = 1;

	ConstraintsMatrixV(9, 0) = 0;
	ConstraintsMatrixV(9, 1) = 0;
	ConstraintsMatrixV(9, 2) = 0;
	ConstraintsMatrixV(9, 3) = 0;
	ConstraintsMatrixV(9, 4) = -1;

	boundValueV.resize(10);
	boundValueV[0] = 0.01;
	boundValueV[1] = -0.01;
	boundValueV[2] = 0.01;
	boundValueV[3] = -0.01;
	boundValueV[4] = 0.01;
	boundValueV[5] = -0.01;
	boundValueV[6] = 0.01;
	boundValueV[7] = -0.01;
	boundValueV[8] = 0.01;
	boundValueV[9] = -0.01;

	boundSignV = 1;

	invariant = polytope::ptr(new polytope()); //creating an universe polytope
	invariant->setIsEmpty(true);
	invariant->setIsUniverse(true);

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = false;
	system_dynamics.MatrixB = Bmatrix;

	system_dynamics.isEmptyC = true;

	system_dynamics.U = polytope::ptr(
			new polytope(ConstraintsMatrixV, boundValueV, boundSignV));
//	Dynamics Initalised ---------------------

	transition::ptr trans = transition::ptr(new transition()); //empty transition

	location::ptr source;
	source = location::ptr(new location());
	source->setLocId(1);
	source->setName("Round_Figure");
	source->setSystem_Dynamics(system_dynamics);
	source->setInvariant(invariant);
	source->setInvariantExists(false); //no invariant available
	source->add_Out_Going_Transition(trans);

	int dim = initial_polytope_I->getSystemDimension();

	Hybrid_Automata.addInitial_Location(source);
	Hybrid_Automata.addLocation(source);
	Hybrid_Automata.setDimension(dim);

	Hybrid_Automata.insert_to_map("x1", 0);
	Hybrid_Automata.insert_to_map("x2", 1);
	Hybrid_Automata.insert_to_map("x3", 2);
	Hybrid_Automata.insert_to_map("x4", 3);
	Hybrid_Automata.insert_to_map("x5", 4);

	unsigned int initial_location_id = 1; //the initial Location ID
	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));

	init_state_list.push_back(I);

}

