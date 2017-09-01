/*
 * NavigationModel2.cpp
 *
 *  Created on: 25-Nov-2014
 *      Author: amit
 *
 *      The Grid is labeled as below *
 *      	B 2 4
 *      	4 7 4
 *   		2 1 A
 *   	With matrix A=[-1.2 0.1; 0.1 -1.2]
 */

#include "Hybrid_Model_Parameters_Design/Navigation_Benchmark/NavigationBenchmark4Var.h"

/*
 * velocity v1 in the x-axis directions and velocity v2 in the y-coordinate directions
 * So the system has Four variables, (x,y) the positions and (v1,v2) the velocities.
 */
//converted ::  Input Polytope U into 4 Variables.
//		Balanced Flow Equations into 4x4 for A and B
//		Invariants converted to 4 Variables
//		Similarly Guard is also converted to 4 Variables
void SetNavigationModel2(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;
	polytope::ptr initial_polytope_I, invariant, gaurd_polytope;
	Dynamics system_dynamics;

	math::matrix<double> ConstraintsMatrixI, ConstraintsMatrixV,
			invariantConstraintsMatrix, gaurdConstraintsMatrix, Amatrix,
			Bmatrix;
	std::vector<double> boundValueI, boundValueV, invariantBoundValue,
			gaurdBoundValue;
	int boundSignI, invariantBoundSign, gaurdBoundSign, boundSignV;

	size_type row, col;
	unsigned int initial_location_id; //the initial Location ID

	//Polytope I Declaration in the form of Ax<=b
	//Input Polytope I as a point(x,y,x1,x2) (0.2 <=x<=0.6,0.1<=y<=0.5,x1==0,x2==0) in the grid of cells
	row = 8;
	col = 4;
	ConstraintsMatrixI.resize(row, col);
	ConstraintsMatrixI(0, 0) = 1;
	ConstraintsMatrixI(0, 1) = 0;
	ConstraintsMatrixI(0, 2) = 0;
	ConstraintsMatrixI(0, 3) = 0;

	ConstraintsMatrixI(1, 0) = -1;
	ConstraintsMatrixI(1, 1) = 0;
	ConstraintsMatrixI(1, 2) = 0;
	ConstraintsMatrixI(1, 3) = 0;

	ConstraintsMatrixI(2, 0) = 0;
	ConstraintsMatrixI(2, 1) = 1;
	ConstraintsMatrixI(2, 2) = 0;
	ConstraintsMatrixI(2, 3) = 0;

	ConstraintsMatrixI(3, 0) = 0;
	ConstraintsMatrixI(3, 1) = -1;
	ConstraintsMatrixI(3, 2) = 0;
	ConstraintsMatrixI(3, 3) = 0;

	ConstraintsMatrixI(4, 0) = 0;
	ConstraintsMatrixI(4, 1) = 0;
	ConstraintsMatrixI(4, 2) = 1;
	ConstraintsMatrixI(4, 3) = 0;

	ConstraintsMatrixI(5, 0) = 0;
	ConstraintsMatrixI(5, 1) = 0;
	ConstraintsMatrixI(5, 2) = -1;
	ConstraintsMatrixI(5, 3) = 0;

	ConstraintsMatrixI(6, 0) = 0;
	ConstraintsMatrixI(6, 1) = 0;
	ConstraintsMatrixI(6, 2) = 0;
	ConstraintsMatrixI(6, 3) = 1;

	ConstraintsMatrixI(7, 0) = 0;
	ConstraintsMatrixI(7, 1) = 0;
	ConstraintsMatrixI(7, 2) = 0;
	ConstraintsMatrixI(7, 3) = -1;

	discrete_set d_set;

	boundValueI.resize(row);
	// ********************* start_location=2:: (0.2<=x1<=0.8 & 0.2<=x2<=0.8 & 0<=v1<=0 & 0<=v2<=0) ************************
	/*	 initial_location_id = 2; //the initial Location ID = 2

	 boundValueI[0] = 0.8;
	 boundValueI[1] = -0.2;
	 boundValueI[2] = 0.8;
	 boundValueI[3] = -0.2;
	 boundValueI[4] = 0;
	 boundValueI[5] = 0;
	 boundValueI[6] = 0;
	 boundValueI[7] = 0;*/

	// ********************* start_location=3:: (1.2<=x1<=1.8 & 0.2<=x2<=0.8 & 0<=v1<=0 & 0<=v2<=0) ************************
	initial_location_id = 3; //the initial Location ID = 3

	boundValueI[0] = 1.8;
	boundValueI[1] = -1.2;
	boundValueI[2] = 0.8;
	boundValueI[3] = -0.2;
	boundValueI[4] = 0;
	boundValueI[5] = 0;
	boundValueI[6] = 0;
	boundValueI[7] = 0;

	// ********************* start_location=1:: (0.5 <=x1<=0.8, 1.5<=x2<=1.8,v1==0,v2==0) ************************
	/*
	 initial_location_id = 1; //the initial Location ID = 1

	 boundValueI[0] = 0.8; //
	 boundValueI[1] = -0.5;
	 boundValueI[2] = 1.8;
	 boundValueI[3] = -1.5;
	 boundValueI[4] = 0;
	 boundValueI[5] = 0;
	 boundValueI[6] = 0;
	 boundValueI[7] = 0;*/

	/*boundValueI[0] = 1; //
	 boundValueI[1] = 0.4;
	 boundValueI[2] = 2;
	 boundValueI[3] = -1;
	 boundValueI[4] = 0;
	 boundValueI[5] = 0;
	 boundValueI[6] = 0;
	 boundValueI[7] = 0;*/

	// ********************* start_location=5:: (1.2 <=x1<=1.4, 2.5<=x2<=2.7,v1==0,v2==0) ************************
	/*
	 initial_location_id = 5;		//the initial Location ID = 5

	 boundValueI[0] = 1.4; //
	 boundValueI[1] = -1.2;
	 boundValueI[2] = 2.7;
	 boundValueI[3] = -2.5;
	 boundValueI[4] = 0;
	 boundValueI[5] = 0;
	 boundValueI[6] = 0;
	 boundValueI[7] = 0;
	 */
	// ********************* start_location=4:: (1.5 <=x1<=1.7, 1.5<=x2<=1.7,v1==-1,v2==0) ************************
	/*initial_location_id = 4; //the initial Location ID = 4

	 boundValueI[0] = 1.7; //(1.5 <=x1<=1.7, 1.5<=x2<=1.7,v1==-1,v2==0)
	 boundValueI[1] = -1.5;
	 boundValueI[2] = 1.7;
	 boundValueI[3] = -1.5;
	 boundValueI[4] = -1;
	 boundValueI[5] = 1;
	 boundValueI[6] = 0;
	 boundValueI[7] = 0;*/

	/*	boundValueI[0] = 1.7; //(1.5 <=x1<=1.7, 1.5<=x2<=1.7,v1==-1,v2==0.5)
	 boundValueI[1] = -1.5;
	 boundValueI[2] = 1.7;
	 boundValueI[3] = -1.5;
	 boundValueI[4] =0;// -1;
	 boundValueI[5] =0;// 1;
	 boundValueI[6] = 0.5;
	 boundValueI[7] = -0.5;*/

	boundSignI = 1;
	initial_polytope_I = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));
	//initial_polytope_I.setPolytope(ConstraintsMatrixI, boundValueI, boundSignI);

	/*	*************** Common Parameter Initialization *******************
	 * Common Parameter for all Locations or transition
	 * such as Matrix A, Matrix B , Transition_dynamics such as Matrix R and vector w;
	 */
	row = 4;
	col = 4;
	Amatrix.resize(row, col);
	Amatrix(0, 0) = 0;
	Amatrix(0, 1) = 0;
	Amatrix(0, 2) = 1;
	Amatrix(0, 3) = 0;

	Amatrix(1, 0) = 0;
	Amatrix(1, 1) = 0;
	Amatrix(1, 2) = 0;
	Amatrix(1, 3) = 1;

	Amatrix(2, 0) = 0;
	Amatrix(2, 1) = 0;
	Amatrix(2, 2) = -1.2;
	Amatrix(2, 3) = 0.1;

	Amatrix(3, 0) = 0;
	Amatrix(3, 1) = 0;
	Amatrix(3, 2) = 0.1;
	Amatrix(3, 3) = -1.2;

	row = 4;
	col = 4;
	Bmatrix.resize(row, col);
	Bmatrix(0, 0) = 0;
	Bmatrix(0, 1) = 0;
	Bmatrix(0, 2) = 0;
	Bmatrix(0, 3) = 0;

	Bmatrix(1, 0) = 0;
	Bmatrix(1, 1) = 0;
	Bmatrix(1, 2) = 0;
	Bmatrix(1, 3) = 0;

	Bmatrix(2, 0) = 0;
	Bmatrix(2, 1) = 0;
	Bmatrix(2, 2) = -1.2;
	Bmatrix(2, 3) = 0.1;

	Bmatrix(3, 0) = 0;
	Bmatrix(3, 1) = 0;
	Bmatrix(3, 2) = 0.1;
	Bmatrix(3, 3) = -1.2;

	//Bmatrix = Amatrix;

	math::matrix<double> R; //Transition Dynamics
	row = 4;
	col = 4;
	R.resize(row, col);
	R(0, 0) = 1;
	R(0, 1) = 0;
	R(0, 2) = 0;
	R(0, 3) = 0;

	R(1, 0) = 0;
	R(1, 1) = 1;
	R(1, 2) = 0;
	R(1, 3) = 0;

	R(2, 0) = 0;
	R(2, 1) = 0;
	R(2, 2) = 1;
	R(2, 3) = 0;

	R(3, 0) = 0;
	R(3, 1) = 0;
	R(3, 2) = 0;
	R(3, 3) = 1;

	std::vector<double> w(row);
	w[0] = 0;
	w[1] = 0;
	w[2] = 0;
	w[3] = 0;

	Assign assignment;
	assignment.Map = R;
	assignment.b = w;

// ***********************************************************

	/*	*************** Initialization of all transition *******************
	 *  List of transition are t1, t2, ... , t20 including transition towards the Locations labelled "A" and "B"
	 *  where Label "A" is the "Final location" to be reached and "B" the "Bad location" to be avoided.
	 */
	row = 8;
	col = 4;
	gaurdConstraintsMatrix.resize(row, col); //this matrix will be common for all transition except the gaurdBoundValue.
	gaurdConstraintsMatrix(0, 0) = 1;
	gaurdConstraintsMatrix(0, 1) = 0;
	gaurdConstraintsMatrix(0, 2) = 0;
	gaurdConstraintsMatrix(0, 3) = 0;

	gaurdConstraintsMatrix(1, 0) = -1;
	gaurdConstraintsMatrix(1, 1) = 0;
	gaurdConstraintsMatrix(1, 2) = 0;
	gaurdConstraintsMatrix(1, 3) = 0;

	gaurdConstraintsMatrix(2, 0) = 0;
	gaurdConstraintsMatrix(2, 1) = 1;
	gaurdConstraintsMatrix(2, 2) = 0;
	gaurdConstraintsMatrix(2, 3) = 0;

	gaurdConstraintsMatrix(3, 0) = 0;
	gaurdConstraintsMatrix(3, 1) = -1;
	gaurdConstraintsMatrix(3, 2) = 0;
	gaurdConstraintsMatrix(3, 3) = 0;

	gaurdConstraintsMatrix(4, 0) = 0;
	gaurdConstraintsMatrix(4, 1) = 0;
	gaurdConstraintsMatrix(4, 2) = 1;
	gaurdConstraintsMatrix(4, 3) = 0;

	gaurdConstraintsMatrix(5, 0) = 0;
	gaurdConstraintsMatrix(5, 1) = 0;
	gaurdConstraintsMatrix(5, 2) = -1;
	gaurdConstraintsMatrix(5, 3) = 0;

	gaurdConstraintsMatrix(6, 0) = 0;
	gaurdConstraintsMatrix(6, 1) = 0;
	gaurdConstraintsMatrix(6, 2) = 0;
	gaurdConstraintsMatrix(6, 3) = 1;

	gaurdConstraintsMatrix(7, 0) = 0;
	gaurdConstraintsMatrix(7, 1) = 0;
	gaurdConstraintsMatrix(7, 2) = 0;
	gaurdConstraintsMatrix(7, 3) = -1;

	gaurdBoundSign = 1;

	gaurdBoundValue.resize(row); //gaurd is:: V_d[sin(loc_name * pi/4), cos(loc_name * pi/4)]
	gaurdBoundValue[0] = 1; // 0<=x<=1 and y==2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -2;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;

	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t1 = transition::ptr(
			new transition(1, "1 to Bad", 1, 9, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; // x==1 and 1<=y<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t2 = transition::ptr(
			new transition(2, "1 to 4", 1, 4, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; //  0<=x<=1 and y==1 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t3 = transition::ptr(
			new transition(3, "1 to 2", 1, 2, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; //  0<=x<=1 and y==1 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t4 = transition::ptr(
			new transition(4, "2 to 1", 2, 1, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; // x==1 and 0<=y<=1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = 0;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t5 = transition::ptr(
			new transition(5, "2 to 3", 2, 3, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // y==1 and 1<=x<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t6 = transition::ptr(
			new transition(6, "3 to 4", 3, 4, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; // x==1 and 0<=y<=1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = 0;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t7 = transition::ptr(
			new transition(7, "3 to 2", 3, 2, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // x==2 and 0<=y<=1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = 0;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t8 = transition::ptr(
			new transition(8, "3 to A", 3, 8, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; // x==1 and 1<=y<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t9 = transition::ptr(
			new transition(9, "4 to 1", 4, 1, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // y==2 and 1<=x<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -2;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t10 = transition::ptr(
			new transition(10, "4 to 5", 4, 5, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // x==2 and 1<=y<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t11 = transition::ptr(
			new transition(11, "4 to 6", 4, 6, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // y==1 and 1<=x<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t12 = transition::ptr(
			new transition(12, "4 to 3", 4, 3, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // y==2 and 1<=x<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t13 = transition::ptr(
			new transition(13, "5 to 4", 5, 4, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // x==2 and 2<=y<=3 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -2;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t14 = transition::ptr(
			new transition(14, "5 to 7", 5, 7, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 1; // x==1 and 2<=y<=3 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -2;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t15 = transition::ptr(
			new transition(15, "5 to Bad", 5, 9, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 3; // y==2 and 2<=x<=3 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -2;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t16 = transition::ptr(
			new transition(16, "6 to 7", 6, 7, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // x==2 and 1<=y<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t17 = transition::ptr(
			new transition(17, "6 to 4", 6, 4, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 3; // y==1 and 2<=x<=3 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 1;
	gaurdBoundValue[3] = -1;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t18 = transition::ptr(
			new transition(18, "6 to A", 6, 8, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 3; // y==2 and 2<=x<=3 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 2;
	gaurdBoundValue[3] = -2;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t19 = transition::ptr(
			new transition(19, "7 to 6", 7, 6, gaurd_polytope, assignment));

	gaurdBoundValue[0] = 2; // x==2 and 2<=y<=3 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -2;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t20 = transition::ptr(
			new transition(20, "7 to 5", 4, 5, gaurd_polytope, assignment));

// ******************* Transition initialized **************************

	/*	*************** Initialization of all Locations *******************
	 *  List of Locations are l1, l2, ... , l7 and Locations labelled "A" as l8("Final State") and
	 *  Locations labelled "B" as l9("Bad State")
	 *  where Label "A" is the "Final location" to be reached and "B" the "Bad location" to be avoided.
	 */// U polytope
	row = 8;
	col = 4;
	ConstraintsMatrixV.resize(row, col); //Common for all polytope u except the boundValueV.
	ConstraintsMatrixV(0, 0) = 1;
	ConstraintsMatrixV(0, 1) = 0;
	ConstraintsMatrixV(0, 2) = 0;
	ConstraintsMatrixV(0, 3) = 0;

	ConstraintsMatrixV(1, 0) = -1;
	ConstraintsMatrixV(1, 1) = 0;
	ConstraintsMatrixV(1, 2) = 0;
	ConstraintsMatrixV(1, 3) = 0;

	ConstraintsMatrixV(2, 0) = 0;
	ConstraintsMatrixV(2, 1) = 1;
	ConstraintsMatrixV(2, 2) = 0;
	ConstraintsMatrixV(2, 3) = 0;

	ConstraintsMatrixV(3, 0) = 0;
	ConstraintsMatrixV(3, 1) = -1;
	ConstraintsMatrixV(3, 2) = 0;
	ConstraintsMatrixV(3, 3) = 0;

	ConstraintsMatrixV(4, 0) = 0;
	ConstraintsMatrixV(4, 1) = 0;
	ConstraintsMatrixV(4, 2) = 1;
	ConstraintsMatrixV(4, 3) = 0;

	ConstraintsMatrixV(5, 0) = 0;
	ConstraintsMatrixV(5, 1) = 0;
	ConstraintsMatrixV(5, 2) = -1;
	ConstraintsMatrixV(5, 3) = 0;

	ConstraintsMatrixV(6, 0) = 0;
	ConstraintsMatrixV(6, 1) = 0;
	ConstraintsMatrixV(6, 2) = 0;
	ConstraintsMatrixV(6, 3) = 1;

	ConstraintsMatrixV(7, 0) = 0;
	ConstraintsMatrixV(7, 1) = 0;
	ConstraintsMatrixV(7, 2) = 0;
	ConstraintsMatrixV(7, 3) = -1;

	boundSignV = 1;
//cout<<"\nBefore Setting Invariants\n";
	row = 4;
	col = 4;
	invariantConstraintsMatrix.resize(row, col); //Common for all polytope except the invariantBoundValue.
	invariantConstraintsMatrix(0, 0) = 1;
	invariantConstraintsMatrix(0, 1) = 0;
	invariantConstraintsMatrix(0, 2) = 0;
	invariantConstraintsMatrix(0, 3) = 0;

	invariantConstraintsMatrix(1, 0) = -1;
	invariantConstraintsMatrix(1, 1) = 0;
	invariantConstraintsMatrix(1, 2) = 0;
	invariantConstraintsMatrix(1, 3) = 0;

	invariantConstraintsMatrix(2, 0) = 0;
	invariantConstraintsMatrix(2, 1) = 1;
	invariantConstraintsMatrix(2, 2) = 0;
	invariantConstraintsMatrix(2, 3) = 0;

	invariantConstraintsMatrix(3, 0) = 0;
	invariantConstraintsMatrix(3, 1) = -1;
	invariantConstraintsMatrix(3, 2) = 0;
	invariantConstraintsMatrix(3, 3) = 0;

	invariantBoundSign = 1;

	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 1; //0<=x<=1 and 1<=y<=2
	invariantBoundValue[1] = 0;
	invariantBoundValue[2] = 2;
	invariantBoundValue[3] = -1;
	row = 8;
	boundValueV.resize(row);
	boundValueV[0] = 0; //v1==0 and v2=-1 => -Vd = -(0,-1) u is (0,0,0,1)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 0;
	boundValueV[5] = 0;
	boundValueV[6] = 1;
	boundValueV[7] = -1;

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = false;
	system_dynamics.MatrixB = Bmatrix;

	system_dynamics.isEmptyC = true;

	//system_dynamics.U.setPolytope(ConstraintsMatrixV, boundValueV, boundSignV);
	system_dynamics.U = polytope::ptr(
			new polytope(ConstraintsMatrixV, boundValueV, boundSignV));

	//invariant.setPolytope(invariantConstraintsMatrix, invariantBoundValue,invariantBoundSign);
	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));
	std::list<transition::ptr> Out_Going_Trans_fromLoc1;
	Out_Going_Trans_fromLoc1.push_back(t1);
	Out_Going_Trans_fromLoc1.push_back(t2);
	Out_Going_Trans_fromLoc1.push_back(t3);

	location::ptr l1 = location::ptr(
			new location(1, "4", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc1));
//  ************ Location ID=1 completed  ************

	invariantBoundValue[0] = 1; //0<=x<=1 and 0<=y<=1
	invariantBoundValue[1] = 0;
	invariantBoundValue[2] = 1;
	invariantBoundValue[3] = 0;

	boundValueV[0] = 0; //v1==1 and v2=0 => -Vd = -(1,0) u is (0,0,-1,0)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = -1;
	boundValueV[5] = 1;
	boundValueV[6] = 0;
	boundValueV[7] = 0;

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = false;
	system_dynamics.MatrixB = Bmatrix;

	system_dynamics.isEmptyC = true;

	//system_dynamics.U.setPolytope(ConstraintsMatrixV, boundValueV, boundSignV);
	system_dynamics.U = polytope::ptr(
			new polytope(ConstraintsMatrixV, boundValueV, boundSignV));
	//invariant.setPolytope(invariantConstraintsMatrix, invariantBoundValue,invariantBoundSign);
	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));
	std::list<transition::ptr> Out_Going_Trans_fromLoc2;
	Out_Going_Trans_fromLoc2.push_back(t4);
	Out_Going_Trans_fromLoc2.push_back(t5);

	location::ptr l2 = location::ptr(
			new location(2, "2", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc2));
	//  ************ Location ID=2 completed  ************

	invariantBoundValue[0] = 2; //1<=x<=2 and 0<=y<=1
	invariantBoundValue[1] = -1;
	invariantBoundValue[2] = 1;
	invariantBoundValue[3] = 0;

	boundValueV[0] = 0; //v1==0.70711 and v2=0.70711 => -Vd = -(0.70711,0.70711) u is (0,0,-0.70711,-0.70711)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = -0.70711;
	boundValueV[5] = 0.70711;
	boundValueV[6] = -0.70711;
	boundValueV[7] = 0.70711;

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = false;
	system_dynamics.MatrixB = Bmatrix;

	system_dynamics.isEmptyC = true;

	//system_dynamics.U.setPolytope(ConstraintsMatrixV, boundValueV, boundSignV);
	system_dynamics.U = polytope::ptr(
			new polytope(ConstraintsMatrixV, boundValueV, boundSignV));
	//invariant.setPolytope(invariantConstraintsMatrix, invariantBoundValue,invariantBoundSign);
	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));
	std::list<transition::ptr> Out_Going_Trans_fromLoc3;
	Out_Going_Trans_fromLoc3.push_back(t6);
	Out_Going_Trans_fromLoc3.push_back(t7);
	Out_Going_Trans_fromLoc3.push_back(t8);

	location::ptr l3 = location::ptr(
			new location(3, "1", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc3));
	//  ************ Location ID=3 completed  ************

	invariantBoundValue[0] = 2; //1<=x<=2 and 1<=y<=2
	invariantBoundValue[1] = -1;
	invariantBoundValue[2] = 2;
	invariantBoundValue[3] = -1;

	boundValueV[0] = 0; //v1==-0.70711 and v2=0.70711 => -Vd = -(-0.70711,0.70711) u is (0,0,0.70711,-0.70711)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 0.70711;
	boundValueV[5] = -0.70711;
	boundValueV[6] = -0.70711;
	boundValueV[7] = 0.70711;

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = false;
	system_dynamics.MatrixB = Bmatrix;

	system_dynamics.isEmptyC = true;

	//system_dynamics.U.setPolytope(ConstraintsMatrixV, boundValueV, boundSignV);
	system_dynamics.U = polytope::ptr(
			new polytope(ConstraintsMatrixV, boundValueV, boundSignV));
	//invariant.setPolytope(invariantConstraintsMatrix, invariantBoundValue,invariantBoundSign);
	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));
	std::list<transition::ptr> Out_Going_Trans_fromLoc4;
	Out_Going_Trans_fromLoc4.push_back(t9);
	Out_Going_Trans_fromLoc4.push_back(t10);
	Out_Going_Trans_fromLoc4.push_back(t11);
	Out_Going_Trans_fromLoc4.push_back(t12);

	location::ptr l4 = location::ptr(
			new location(4, "7", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc4));
	//  ************ Location ID=4 completed  ************

	invariantBoundValue[0] = 2; //1<=x<=2 and 2<=y<=3
	invariantBoundValue[1] = -1;
	invariantBoundValue[2] = 3;
	invariantBoundValue[3] = -2;

	boundValueV[0] = 0; //v1==1 and v2=0 => -Vd = -(1,0) u is (0,0,-1,0)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = -1;
	boundValueV[5] = 1;
	boundValueV[6] = 0;
	boundValueV[7] = 0;

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = false;
	system_dynamics.MatrixB = Bmatrix;

	system_dynamics.isEmptyC = true;

	//system_dynamics.U.setPolytope(ConstraintsMatrixV, boundValueV, boundSignV);
	system_dynamics.U = polytope::ptr(
			new polytope(ConstraintsMatrixV, boundValueV, boundSignV));
	//invariant.setPolytope(invariantConstraintsMatrix, invariantBoundValue,invariantBoundSign);
	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));
	std::list<transition::ptr> Out_Going_Trans_fromLoc5;
	Out_Going_Trans_fromLoc5.push_back(t13);
	Out_Going_Trans_fromLoc5.push_back(t14);
	Out_Going_Trans_fromLoc5.push_back(t15);

	location::ptr l5 = location::ptr(
			new location(5, "2", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc5));
	//  ************ Location ID=5 completed  ************

	invariantBoundValue[0] = 3; //2<=x<=3 and 1<=y<=2
	invariantBoundValue[1] = -2;
	invariantBoundValue[2] = 2;
	invariantBoundValue[3] = -1;

	boundValueV[0] = 0; //v1==0 and v2=-1 => -Vd = -(0,-1) u is (0,0,0,1)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 0;
	boundValueV[5] = 0;
	boundValueV[6] = 1;
	boundValueV[7] = -1;

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = false;
	system_dynamics.MatrixB = Bmatrix;

	system_dynamics.isEmptyC = true;

	//system_dynamics.U.setPolytope(ConstraintsMatrixV, boundValueV, boundSignV);
	system_dynamics.U = polytope::ptr(
			new polytope(ConstraintsMatrixV, boundValueV, boundSignV));
	//invariant.setPolytope(invariantConstraintsMatrix, invariantBoundValue,invariantBoundSign);
	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));
	std::list<transition::ptr> Out_Going_Trans_fromLoc6;
	Out_Going_Trans_fromLoc6.push_back(t16);
	Out_Going_Trans_fromLoc6.push_back(t17);
	Out_Going_Trans_fromLoc6.push_back(t18);

	location::ptr l6 = location::ptr(
			new location(6, "4", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc6));
	//  ************ Location ID=6 completed  ************

	invariantBoundValue[0] = 3; //2<=x<=3 and 2<=y<=3
	invariantBoundValue[1] = -2;
	invariantBoundValue[2] = 3;
	invariantBoundValue[3] = -2;

	boundValueV[0] = 0; //v1==0 and v2=-1 => -Vd = -(0,-1) u is (0,0,0,1)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 0;
	boundValueV[5] = 0;
	boundValueV[6] = 1;
	boundValueV[7] = -1;

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = false;
	system_dynamics.MatrixB = Bmatrix;

	system_dynamics.isEmptyC = true;

	//system_dynamics.U.setPolytope(ConstraintsMatrixV, boundValueV, boundSignV);
	system_dynamics.U = polytope::ptr(
			new polytope(ConstraintsMatrixV, boundValueV, boundSignV));
	//invariant.setPolytope(invariantConstraintsMatrix, invariantBoundValue,invariantBoundSign);
	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));
	std::list<transition::ptr> Out_Going_Trans_fromLoc7;
	Out_Going_Trans_fromLoc7.push_back(t19);
	Out_Going_Trans_fromLoc7.push_back(t20);

	location::ptr l7 = location::ptr(
			new location(7, "4", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc7));
	//  ************ Location ID=7 completed  ************

//	************ No dynamics available for location=8/9    ************
//	either the values are zeros or NULL but they will not be processed as the Name is FINAL/GOOD/BAD/UNSAFE

	invariantBoundValue[0] = 0;
	invariantBoundValue[1] = 0;
	invariantBoundValue[2] = 0;
	invariantBoundValue[3] = 0;

	boundValueV[0] = 0; //v1==0 and v2=0 => -Vd = -(0,0) u is (0,0,0,0)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 0;
	boundValueV[5] = 0;
	boundValueV[6] = 0;
	boundValueV[7] = 0;

	system_dynamics.isEmptyMatrixA = false;
	system_dynamics.MatrixA = Amatrix;

	system_dynamics.isEmptyMatrixB = false;
	system_dynamics.MatrixB = Bmatrix;

	system_dynamics.isEmptyC = true;

	//system_dynamics.U.setPolytope(ConstraintsMatrixV, boundValueV, boundSignV);
	system_dynamics.U = polytope::ptr(
			new polytope(ConstraintsMatrixV, boundValueV, boundSignV));
	//invariant.setPolytope(invariantConstraintsMatrix, invariantBoundValue,invariantBoundSign);
	invariant = polytope::ptr(
			new polytope(invariantConstraintsMatrix, invariantBoundValue,
					invariantBoundSign));
	std::list<transition::ptr> Out_Going_Trans_fromLoc8,
			Out_Going_Trans_fromLoc9;

	location::ptr l8 = location::ptr(
			new location(8, "FINAL", system_dynamics, invariant, false,
					Out_Going_Trans_fromLoc8));
	location::ptr l9 = location::ptr(
			new location(9, "BAD", system_dynamics, invariant, false,
					Out_Going_Trans_fromLoc9));
//Location ID=8 and ID=9 completed ************

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
	Hybrid_Automata.addLocation(l7);
	Hybrid_Automata.addLocation(l8);
	Hybrid_Automata.addLocation(l9);

	Hybrid_Automata.insert_to_map("x1", 0);
	Hybrid_Automata.insert_to_map("x2", 1);
	Hybrid_Automata.insert_to_map("v1", 2);
	Hybrid_Automata.insert_to_map("v2", 3);

	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I;
	I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));

	init_state_list.push_back(I);

}
