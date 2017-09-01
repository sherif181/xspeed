/*
 * NavigationBenchmark.cpp
 *
 *  Created on: 25-Nov-2014
 *      Author: amit
 *
 *      The Grid of the Navigation model of 5 x 5 is labeled as below *
 *      	2 4 6 6 6
 *      	2 4 7 7 4
 *   		2 4 B 3 4
 *   		2 4 6 6 6
 *   		2 A 0 0 0
 */

#include "Hybrid_Model_Parameters_Design/Navigation_Benchmark/NavigationBenchmark4Var.h"

/*
 * velocity x1 in the x-axis directions and velocity x2 in the y-coordinate directions
 * So the system has Four variables, (x,y) the positions and (x1,x2) the velocities.
 */

//converted ::  Input Polytope U into 4 Variables.
//		Balanced Flow Equations into 4x4 for A and B
//		Invariants converted to 4 Variables
//		Similarly Guard is also converted to 4 Variables
void SetNavigationModel5by5OurCode(hybrid_automata& Hybrid_Automata,
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
	unsigned int initial_location_id;
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

	//discrete_set d_set;

	boundValueI.resize(row);

	// ********************* start_location=19:: (3.5 <=x1<=3.5, 3.5<=x2<=3.5, 1<=v1<=1, 1<=v2<=1) ************************
//	d_set.insert_element(19);	//the initial Location ID = 19
	initial_location_id = 19;

	//3.1<=x1<=3.4 & 3.6<=x2<=3.8 & 0.1<=v1<=0.1 & 0.1<=v2<=0.1
	boundValueI[0] = 3.4; // ************ :: (3.5<=x1<=3.5 & 3.5<=x2<=3.5 & 1<=v1<=1 & 0.1<=v2<=0.1) ************************
	boundValueI[1] = -3.1;
	boundValueI[2] = 3.8;
	boundValueI[3] = -3.6;
	boundValueI[4] = 0.1;
	boundValueI[5] = -0.1;
	boundValueI[6] = 0.1;
	boundValueI[7] = -0.1;

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
	Amatrix(2, 2) = -0.8;
	Amatrix(2, 3) = -0.2;

	Amatrix(3, 0) = 0;
	Amatrix(3, 1) = 0;
	Amatrix(3, 2) = -0.1;
	Amatrix(3, 3) = -0.8;

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
	Bmatrix(2, 2) = -0.8;
	Bmatrix(2, 3) = -0.2;

	Bmatrix(3, 0) = 0;
	Bmatrix(3, 1) = 0;
	Bmatrix(3, 2) = -0.1;
	Bmatrix(3, 3) = -0.8;

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
	 *  List of transition are t1, t2, ... , t73 including transition towards the Locations labelled "A" and "B"
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
//  ************************* Transition: t1 ***************************************************
	gaurdBoundValue.resize(row); //gaurd is:: V_d[sin(loc_name * pi/4), cos(loc_name * pi/4)]
	gaurdBoundValue[0] = 1; // 0<=x<=1 and x2==4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t1 = transition::ptr(
			new transition(1, "1 to 2", 1, 2, gaurd_polytope, assignment));
//  ************************* Transition: t1 End ************************************************
//  ************************* Transition: t2 ***************************************************
	gaurdBoundValue[0] = 1; // x1==1 and 4<=x2<=5 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 5;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t2 = transition::ptr(
			new transition(2, "1 to 10", 1, 10, gaurd_polytope, assignment));
//  ************************* Transition: t2 End ************************************************
//  ************************* Transition: t3 ***************************************************
	gaurdBoundValue[0] = 1; // 0<=x1<=1 and x2==4 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t3 = transition::ptr(
			new transition(3, "2 to 1", 2, 1, gaurd_polytope, assignment));
//  ************************* Transition: t3 End ************************************************
//  ************************* Transition: t4 ***************************************************
	gaurdBoundValue[0] = 1; // 0<=x1<=1 and x2==3 and   -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t4 = transition::ptr(
			new transition(4, "2 to 3", 2, 3, gaurd_polytope, assignment));
//  ************************* Transition: t4 End ************************************************
//  ************************* Transition: t5 ***************************************************
	gaurdBoundValue[0] = 1; // x==1 and 3<=x2<=4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1; //testing  0.95<=x<=1
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t5 = transition::ptr(
			new transition(5, "2 to 9", 2, 9, gaurd_polytope, assignment));
//  ************************* Transition: t5 End ************************************************
//  ************************* Transition: t6 ***************************************************
	gaurdBoundValue[0] = 1; // 0<=x1<=1 and x2==3 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = 0;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t6 = transition::ptr(
			new transition(6, "3 to 2", 3, 2, gaurd_polytope, assignment));
//  ************************* Transition: t6 End ************************************************
//  ************************* Transition: t7 ***************************************************
	gaurdBoundValue[0] = 1; // 0<=x1<=1 and x2==2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t7 = transition::ptr(
			new transition(7, "3 to 4", 3, 4, gaurd_polytope, assignment));
//  ************************* Transition: t7 End ************************************************
//  ************************* Transition: t8 ***************************************************
	gaurdBoundValue[0] = 1; // x1==1 and 2<=x2<=3 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t8 = transition::ptr(
			new transition(8, "3 to 8", 3, 8, gaurd_polytope, assignment));
//  ************************* Transition: t8 End ************************************************
//  ************************* Transition: t9 ***************************************************
	gaurdBoundValue[0] = 1; // 0<=x1<=1 and x2==2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t9 = transition::ptr(
			new transition(9, "4 to 3", 4, 3, gaurd_polytope, assignment));
//  ************************* Transition: t9 End ************************************************
//  ************************* Transition: t10 ***************************************************
	gaurdBoundValue[0] = 1; // 0<=x1<=1 and x2==1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t10 = transition::ptr(
			new transition(10, "4 to 5", 4, 5, gaurd_polytope, assignment));
//  ************************* Transition: t10 End ************************************************
//  ************************* Transition: t11 ***************************************************
	gaurdBoundValue[0] = 1; // x1==1 and 1<=x2<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t11 = transition::ptr(
			new transition(11, "4 to 7", 4, 7, gaurd_polytope, assignment));
//  ************************* Transition: t11 End ************************************************
//  ************************* Transition: t12 ***************************************************
	gaurdBoundValue[0] = 1; // 0<=x1<=1 and x2==1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t12 = transition::ptr(
			new transition(12, "5 to 4", 5, 4, gaurd_polytope, assignment));
//  ************************* Transition: t12 End ************************************************
//  ************************* Transition: t13 ***************************************************
	gaurdBoundValue[0] = 1; // x1==1 and 0<=x2<=1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t13 = transition::ptr(
			new transition(13, "5 to 6", 5, 6, gaurd_polytope, assignment));
//  ************************* Transition: t13 End ************************************************
//  ************************* Transition: t14 ***************************************************
	gaurdBoundValue[0] = 2; // 1<=x1<=2 and x2==2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t14 = transition::ptr(
			new transition(14, "7 to 8", 7, 8, gaurd_polytope, assignment));
//  ************************* Transition: t14 End ************************************************
//  ************************* Transition: t15 ***************************************************
	gaurdBoundValue[0] = 1; // x1==1 and 1<=x2<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t15 = transition::ptr(
			new transition(15, "7 to 4", 7, 4, gaurd_polytope, assignment));
//  ************************* Transition: t15 End ************************************************
//  ************************* Transition: t16 ***************************************************
	gaurdBoundValue[0] = 2; // 1<=x1<=2 and x2==1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t16 = transition::ptr(
			new transition(16, "7 to A", 7, 6, gaurd_polytope, assignment));
//  ************************* Transition: t16 End ************************************************
//  ************************* Transition: t17 ***************************************************
	gaurdBoundValue[0] = 2; // x1==2 and 1<=x2<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
			new transition(17, "7 to 12", 7, 12, gaurd_polytope, assignment));
//  ************************* Transition: t17 End ***********************************************
//  ************************* Transition: t18 ***************************************************
	gaurdBoundValue[0] = 2; //  1<=x1<=2 and x2==3 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t18 = transition::ptr(
			new transition(18, "8 to 9", 8, 9, gaurd_polytope, assignment));
	//  ************************* Transition: t18 End ***********************************************
	//  ************************* Transition: t19 ***************************************************
	gaurdBoundValue[0] = 1; // x1==1 and 2<=x2<=3 and -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t19 = transition::ptr(
			new transition(19, "8 to 3", 8, 3, gaurd_polytope, assignment));
//  ************************* Transition: t19 End ***********************************************
//  ************************* Transition: t20 ***************************************************
	gaurdBoundValue[0] = 2; //  1<=x1<=2 and x2==2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t20 = transition::ptr(
			new transition(20, "8 to 7", 8, 7, gaurd_polytope, assignment));
//  ************************* Transition: t20 End ***********************************************
//  ************************* Transition: t21 ***************************************************
	gaurdBoundValue[0] = 2; //  x1==2 and 2<=x2<=3 and -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t21 = transition::ptr(
			new transition(21, "8 to B", 8, 13, gaurd_polytope, assignment));
//  ************************* Transition: t21 End ***********************************************
//  ************************* Transition: t22 ***************************************************
	gaurdBoundValue[0] = 2; //  1<=x1<=2 and x2==4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t22 = transition::ptr(
			new transition(22, "9 to 10", 9, 10, gaurd_polytope, assignment));
//  ************************* Transition: t22 End ***********************************************
	//  ************************* Transition: t23 ***************************************************
	gaurdBoundValue[0] = 1; // x==1 and 3<=x2<=4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1; //testing  0.95<=x<=1
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t23 = transition::ptr(
			new transition(23, "9 to 2", 9, 2, gaurd_polytope, assignment));
	//  ************************* Transition: t23 End ************************************************
	//  ************************* Transition: t24 ***************************************************
	gaurdBoundValue[0] = 2; //  1<=x1<=2 and x2==3 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t24 = transition::ptr(
			new transition(18, "9 to 8", 9, 8, gaurd_polytope, assignment));
	//  ************************* Transition: t24 End ***********************************************
	//  ************************* Transition: t25 ***************************************************
	gaurdBoundValue[0] = 2; // x==2 and 3<=x2<=4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2; //testing  0.95<=x<=1
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t25 = transition::ptr(
			new transition(25, "9 to 14", 9, 14, gaurd_polytope, assignment));
	//  ************************* Transition: t25 End ************************************************
	//  ************************* Transition: t26 ***************************************************
	gaurdBoundValue[0] = 1; // x1==1 and 4<=x2<=5 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 5;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t26 = transition::ptr(
			new transition(26, "10 to 1", 10, 1, gaurd_polytope, assignment));
	//  ************************* Transition: t26 End ************************************************
	//  ************************* Transition: t27 ***************************************************
	gaurdBoundValue[0] = 2; //  1<=x1<=2 and x2==4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -1;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t27 = transition::ptr(
			new transition(27, "10 to 9", 10, 9, gaurd_polytope, assignment));
	//  ************************* Transition: t27 End ***********************************************
	//  ************************* Transition: t28 ***************************************************
	gaurdBoundValue[0] = 2; // x1==2 and 4<=x2<=5 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 5;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t28 = transition::ptr(
			new transition(28, "10 to 15", 10, 1, gaurd_polytope, assignment));
	//  ************************* Transition: t28 End ************************************************
	//  ************************* Transition: t29 ***************************************************
	gaurdBoundValue[0] = 3; //  2<=x1<=3 and x2==1 and -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t29 = transition::ptr(
			new transition(29, "11 to 12", 11, 12, gaurd_polytope, assignment));
	//  ************************* Transition: t29 End ***********************************************
	//  ************************* Transition: t30 ***************************************************
	gaurdBoundValue[0] = 2; // x1==2 and 0<=x2<=1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t30 = transition::ptr(
			new transition(30, "11 to A", 11, 6, gaurd_polytope, assignment));
	//  ************************* Transition: t30 End ************************************************
	//  ************************* Transition: t31 ***************************************************
	gaurdBoundValue[0] = 3; // x1==3 and 0<=x2<=1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
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
	transition::ptr t31 = transition::ptr(
			new transition(31, "11 to 16", 11, 16, gaurd_polytope, assignment));
	//  ************************* Transition: t31 End ************************************************
	//  ************************* Transition: t32 ***************************************************
	gaurdBoundValue[0] = 3; //  2<=x1<=3 and x2==2 and -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t32 = transition::ptr(
			new transition(32, "12 to B", 12, 13, gaurd_polytope, assignment));
	//  ************************* Transition: t32 End ***********************************************
	//  ************************* Transition: t33 ***************************************************
	gaurdBoundValue[0] = 2; // x1==2 and 1<=x2<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t33 = transition::ptr(
			new transition(33, "12 to 7", 12, 7, gaurd_polytope, assignment));
	//  ************************* Transition: t33 End ***********************************************
	//  ************************* Transition: t34 ***************************************************
	gaurdBoundValue[0] = 3; //  2<=x1<=3 and x2==1 and -1000<=v1<=1000 &  -1000<=v2<=1000
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
	transition::ptr t34 = transition::ptr(
			new transition(34, "12 to 11", 12, 11, gaurd_polytope, assignment));
	//  ************************* Transition: t34 End ***********************************************
	//  ************************* Transition: t35 ***************************************************
	gaurdBoundValue[0] = 3; // x1==3 and 1<=x2<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
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
	transition::ptr t35 = transition::ptr(
			new transition(35, "12 to 17", 12, 17, gaurd_polytope, assignment));
	//  ************************* Transition: t35 End ************************************************
	//  ************************* Transition: t36 ***************************************************
	gaurdBoundValue[0] = 3; //  2<=x1<=3 and x2==4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t36 = transition::ptr(
			new transition(36, "14 to 15", 14, 15, gaurd_polytope, assignment));
	//  ************************* Transition: t36 End ***********************************************
	//  ************************* Transition: t37 ***************************************************
	gaurdBoundValue[0] = 2; // x==2 and 3<=x2<=4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2; //testing  0.95<=x<=1
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t37 = transition::ptr(
			new transition(37, "14 to 9", 14, 9, gaurd_polytope, assignment));
	//  ************************* Transition: t37 End ************************************************
	//  ************************* Transition: t38 ***************************************************
	gaurdBoundValue[0] = 3; //  2=x1<=3 and x2==3 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t38 = transition::ptr(
			new transition(38, "14 to B", 14, 13, gaurd_polytope, assignment));
	//  ************************* Transition: t38 End ***********************************************
	//  ************************* Transition: t39 ***************************************************
	gaurdBoundValue[0] = 3; // x1==3 and 3<=x2<=4 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t39 = transition::ptr(
			new transition(39, "14 to 19", 14, 19, gaurd_polytope, assignment));
	//  ************************* Transition: t39 End ************************************************
	//  ************************* Transition: t40 ***************************************************
	gaurdBoundValue[0] = 2; // x1==2 and 4<=x2<=5 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 5;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t40 = transition::ptr(
			new transition(40, "15 to 10", 15, 10, gaurd_polytope, assignment));
	//  ************************* Transition: t40 End ************************************************
	//  ************************* Transition: t41 ***************************************************
	gaurdBoundValue[0] = 3; //  2<=x1<=3 and x2==4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -2;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t41 = transition::ptr(
			new transition(41, "15 to 14", 15, 14, gaurd_polytope, assignment));
	//  ************************* Transition: t41 End ***********************************************
	//  ************************* Transition: t42 ***************************************************
	gaurdBoundValue[0] = 3; // x1==3 and 4<=x2<=5 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
	gaurdBoundValue[2] = 5;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t42 = transition::ptr(
			new transition(42, "15 to 20", 15, 20, gaurd_polytope, assignment));
	//  ************************* Transition: t42 End ************************************************
	//  ************************* Transition: t43 ***************************************************
	gaurdBoundValue[0] = 4; //  3<=x1<=4 and x2==1 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
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
	transition::ptr t43 = transition::ptr(
			new transition(43, "16 to 17", 16, 17, gaurd_polytope, assignment));
	//  ************************* Transition: t43 End ***********************************************
	//  ************************* Transition: t44 ***************************************************
	gaurdBoundValue[0] = 3; // x1==3 and 0<=x2<=1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
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
	transition::ptr t44 = transition::ptr(
			new transition(44, "16 to 11", 16, 11, gaurd_polytope, assignment));
	//  ************************* Transition: t44 End ************************************************
	//  ************************* Transition: t45 ***************************************************
	gaurdBoundValue[0] = 4; // x1==4 and 0<=x2<=1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
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
	transition::ptr t45 = transition::ptr(
			new transition(45, "16 to 21", 16, 21, gaurd_polytope, assignment));
	//  ************************* Transition: t45 End ************************************************
	//  ************************* Transition: t46 ***************************************************
	gaurdBoundValue[0] = 4; //  3<=x1<=4 and x2==2 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
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
	transition::ptr t46 = transition::ptr(
			new transition(46, "17 to 18", 17, 18, gaurd_polytope, assignment));
	//  ************************* Transition: t46 End ***********************************************
	//  ************************* Transition: t47 ***************************************************
	gaurdBoundValue[0] = 3; // x1==3 and 1<=x2<=2 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
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
	transition::ptr t47 = transition::ptr(
			new transition(47, "17 to 12", 17, 12, gaurd_polytope, assignment));
	//  ************************* Transition: t47 End ************************************************
	//  ************************* Transition: t48 ***************************************************
	gaurdBoundValue[0] = 4; //  3<=x1<=4 and x2==1 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
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
	transition::ptr t48 = transition::ptr(
			new transition(48, "17 to 16", 17, 16, gaurd_polytope, assignment));
	//  ************************* Transition: t48 End ***********************************************
	//  ************************* Transition: t49 ***************************************************
	gaurdBoundValue[0] = 4; // x1==4 and 1<=x2<=2 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
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
	transition::ptr t49 = transition::ptr(
			new transition(49, "17 to 22", 17, 22, gaurd_polytope, assignment));
	//  ************************* Transition: t49 End ************************************************
	//  ************************* Transition: t50 ***************************************************
	gaurdBoundValue[0] = 4; //  3<=x1<=4 and x2==3 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t50 = transition::ptr(
			new transition(50, "18 to 19", 18, 19, gaurd_polytope, assignment));
	//  ************************* Transition: t50 End ***********************************************
	//  ************************* Transition: t51 ***************************************************
	gaurdBoundValue[0] = 3; // x1==3 and 2<=x2<=3 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
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
	transition::ptr t51 = transition::ptr(
			new transition(51, "18 to B", 18, 13, gaurd_polytope, assignment));
	//  ************************* Transition: t51 End ************************************************
	//  ************************* Transition: t52 ***************************************************
	gaurdBoundValue[0] = 4; //  3<=x1<=4 and x2==2 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
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
	transition::ptr t52 = transition::ptr(
			new transition(52, "18 to 17", 18, 17, gaurd_polytope, assignment));
	//  ************************* Transition: t52 End ***********************************************
	//  ************************* Transition: t53 ***************************************************
	gaurdBoundValue[0] = 4; // x1==4 and 2<=x2<=3 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
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
	transition::ptr t53 = transition::ptr(
			new transition(53, "18 to 23", 18, 23, gaurd_polytope, assignment));
	//  ************************* Transition: t53 End ************************************************
	//  ************************* Transition: t54 ***************************************************
	gaurdBoundValue[0] = 4; //  3<=x1<=4 and x2==4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t54 = transition::ptr(
			new transition(54, "19 to 20", 19, 20, gaurd_polytope, assignment));
	//  ************************* Transition: t54 End ***********************************************
	//  ************************* Transition: t55 ***************************************************
	gaurdBoundValue[0] = 3; // x1==3 and 3<=x2<=4 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t55 = transition::ptr(
			new transition(55, "19 to 14", 19, 14, gaurd_polytope, assignment));
	//  ************************* Transition: t55 End ************************************************
	//  ************************* Transition: t56 ***************************************************
	gaurdBoundValue[0] = 4; //  3<=x1<=4 and x2==3 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t56 = transition::ptr(
			new transition(56, "19 to 18", 19, 18, gaurd_polytope, assignment));
	//  ************************* Transition: t56 End ***********************************************
	//  ************************* Transition: t57 ***************************************************
	gaurdBoundValue[0] = 4; // x1==4 and 3<=x2<=4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t57 = transition::ptr(
			new transition(57, "19 to 24", 19, 24, gaurd_polytope, assignment));
	//  ************************* Transition: t57 End ************************************************
	//  ************************* Transition: t58 ***************************************************
	gaurdBoundValue[0] = 3; // x1==3 and 4<=x2<=5 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
	gaurdBoundValue[2] = 5;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t58 = transition::ptr(
			new transition(58, "20 to 15", 20, 15, gaurd_polytope, assignment));
	//  ************************* Transition: t58 End ************************************************
	//  ************************* Transition: t59 ***************************************************
	gaurdBoundValue[0] = 4; //  3<=x1<=4 and x2==4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -3;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t59 = transition::ptr(
			new transition(59, "20 to 19", 20, 19, gaurd_polytope, assignment));
	//  ************************* Transition: t59 End ***********************************************
	//  ************************* Transition: t60 ***************************************************
	gaurdBoundValue[0] = 4; // x1==4 and 4<=x2<=5 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
	gaurdBoundValue[2] = 5;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t60 = transition::ptr(
			new transition(60, "20 to 25", 20, 25, gaurd_polytope, assignment));
	//  ************************* Transition: t60 End ************************************************
	//  ************************* Transition: t61 ***************************************************
	gaurdBoundValue[0] = 5; //  4<=x1<=5 and x2==1 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
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
	transition::ptr t61 = transition::ptr(
			new transition(61, "21 to 22", 21, 22, gaurd_polytope, assignment));
	//  ************************* Transition: t61 End ***********************************************
	//  ************************* Transition: t62 ***************************************************
	gaurdBoundValue[0] = 4; // x1==4 and 0<=x2<=1 and  -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
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
	transition::ptr t62 = transition::ptr(
			new transition(62, "21 to 16", 21, 16, gaurd_polytope, assignment));
	//  ************************* Transition: t62 End ************************************************
	//  ************************* Transition: t63 ***************************************************
	gaurdBoundValue[0] = 5; //  4<=x1<=5 and x2==2 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
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
	transition::ptr t63 = transition::ptr(
			new transition(63, "22 to 23", 22, 23, gaurd_polytope, assignment));
	//  ************************* Transition: t63 End ***********************************************
	//  ************************* Transition: t64 ***************************************************
	gaurdBoundValue[0] = 4; // x1==4 and 1<=x2<=2 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
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
	transition::ptr t64 = transition::ptr(
			new transition(64, "22 to 17", 22, 17, gaurd_polytope, assignment));
	//  ************************* Transition: t64 End ************************************************
	//  ************************* Transition: t65 ***************************************************
	gaurdBoundValue[0] = 5; //  4<=x1<=5 and x2==1 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
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
	transition::ptr t65 = transition::ptr(
			new transition(65, "22 to 21", 22, 21, gaurd_polytope, assignment));
	//  ************************* Transition: t65 End ***********************************************
	//  ************************* Transition: t66 ***************************************************
	gaurdBoundValue[0] = 5; //  4<=x1<=5 and x2==3 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t66 = transition::ptr(
			new transition(66, "23 to 24", 23, 24, gaurd_polytope, assignment));
	//  ************************* Transition: t66 End ***********************************************
	//  ************************* Transition: t67 ***************************************************
	gaurdBoundValue[0] = 4; // x1==4 and 2<=x2<=3 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
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
	transition::ptr t67 = transition::ptr(
			new transition(67, "23 to 18", 23, 18, gaurd_polytope, assignment));
	//  ************************* Transition: t67 End ************************************************
	//  ************************* Transition: t68 ***************************************************
	gaurdBoundValue[0] = 5; //  4<=x1<=5 and x2==2 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
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
	transition::ptr t68 = transition::ptr(
			new transition(68, "23 to 23", 23, 22, gaurd_polytope, assignment));
	//  ************************* Transition: t68 End ***********************************************
	//  ************************* Transition: t69 ***************************************************
	gaurdBoundValue[0] = 5; //  4<=x1<=5 and x2==4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t69 = transition::ptr(
			new transition(69, "24 to 25", 24, 25, gaurd_polytope, assignment));
	//  ************************* Transition: t69 End ***********************************************
	//  ************************* Transition: t70 ***************************************************
	gaurdBoundValue[0] = 4; // x1==4 and 3<=x2<=4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t70 = transition::ptr(
			new transition(70, "24 to 19", 24, 19, gaurd_polytope, assignment));
	//  ************************* Transition: t70 End ************************************************
	//  ************************* Transition: t71 ***************************************************
	gaurdBoundValue[0] = 5; //  4<=x1<=5 and x2==3 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
	gaurdBoundValue[2] = 3;
	gaurdBoundValue[3] = -3;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t71 = transition::ptr(
			new transition(71, "24 to 23", 24, 23, gaurd_polytope, assignment));
	//  ************************* Transition: t71 End ***********************************************
	//  ************************* Transition: t72 ***************************************************
	gaurdBoundValue[0] = 4; // x1==4 and 4<=x2<=5 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
	gaurdBoundValue[2] = 5;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t72 = transition::ptr(
			new transition(72, "25 to 20", 25, 20, gaurd_polytope, assignment));
	//  ************************* Transition: t72 End ************************************************
	//  ************************* Transition: t73 ***************************************************
	gaurdBoundValue[0] = 5; //  4<=x1<=5 and x2==4 and -1000<=v1<=1000 &  -1000<=v2<=1000
	gaurdBoundValue[1] = -4;
	gaurdBoundValue[2] = 4;
	gaurdBoundValue[3] = -4;
	gaurdBoundValue[4] = 1000;
	gaurdBoundValue[5] = 1000;
	gaurdBoundValue[6] = 1000;
	gaurdBoundValue[7] = 1000;
	//gaurd_polytope.setPolytope(gaurdConstraintsMatrix, gaurdBoundValue,gaurdBoundSign);
	gaurd_polytope = polytope::ptr(
			new polytope(gaurdConstraintsMatrix, gaurdBoundValue,
					gaurdBoundSign));
	transition::ptr t73 = transition::ptr(
			new transition(73, "25 to 24", 25, 24, gaurd_polytope, assignment));
	//  ************************* Transition: t73 End ***********************************************
// *******************************************************************************************************
//	////// ********************* Transition Created and all Initialized **************************  //////
// *******************************************************************************************************

	/*	*************** Initialization of all Locations *******************
	 *  List of Locations are l1, l2, ... , l25 and Locations labelled "A" as l6("Final State") and
	 *  Locations labelled "B" as l13("Bad State")
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
	invariantBoundValue[0] = 1; //0<=x1<=1 and 4<=x2<=5
	invariantBoundValue[1] = 0;
	invariantBoundValue[2] = 5;
	invariantBoundValue[3] = -4;
	row = 8;
	boundValueV.resize(row);
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc1;
	Out_Going_Trans_fromLoc1.push_back(t1);
	Out_Going_Trans_fromLoc1.push_back(t2);

	location::ptr l1 = location::ptr(
			new location(1, "2", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc1));
//  ************ Location ID=1 completed  ************
	row = 4;
	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 1; //0<=x1<=1 and 3<=x2<=4
	invariantBoundValue[1] = 0;
	invariantBoundValue[2] = 4;
	invariantBoundValue[3] = -3;
	row = 8;
	boundValueV.resize(row);
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
	Out_Going_Trans_fromLoc2.push_back(t3);
	Out_Going_Trans_fromLoc2.push_back(t4);
	Out_Going_Trans_fromLoc2.push_back(t5);

	location::ptr l2 = location::ptr(
			new location(2, "2", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc2));
//  ************ Location ID=2 completed  ************
	row = 4;
	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 1; //0<=x1<=1 and 2<=x2<=3
	invariantBoundValue[1] = 0;
	invariantBoundValue[2] = 3;
	invariantBoundValue[3] = -2;
	row = 8;
	boundValueV.resize(row);
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc3;
	Out_Going_Trans_fromLoc3.push_back(t6);
	Out_Going_Trans_fromLoc3.push_back(t7);
	Out_Going_Trans_fromLoc3.push_back(t8);

	location::ptr l3 = location::ptr(
			new location(3, "2", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc3));
//  ************ Location ID=3 completed  ************
	row = 4;
	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 1; //0<=x1<=1 and 1<=x2<=2
	invariantBoundValue[1] = 0;
	invariantBoundValue[2] = 2;
	invariantBoundValue[3] = -1;
	row = 8;
	boundValueV.resize(row);
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc4;
	Out_Going_Trans_fromLoc4.push_back(t9);
	Out_Going_Trans_fromLoc4.push_back(t10);
	Out_Going_Trans_fromLoc4.push_back(t11);

	location::ptr l4 = location::ptr(
			new location(4, "2", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc4));
//  ************ Location ID=4 completed  ************
	row = 4;
	invariantBoundValue.resize(row);
	invariantBoundValue[0] = 1; //0<=x1<=1 and 0<=x2<=1
	invariantBoundValue[1] = 0;
	invariantBoundValue[2] = 1;
	invariantBoundValue[3] = 0;
	row = 8;
	boundValueV.resize(row);
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
	Out_Going_Trans_fromLoc5.push_back(t12);
	Out_Going_Trans_fromLoc5.push_back(t13);

	location::ptr l5 = location::ptr(
			new location(5, "2", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc5));
//  ************ Location ID=5 completed  ************

	//	************ No dynamics available for location=6/13    ************
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc6,
			Out_Going_Trans_fromLoc13;

	location::ptr l6 = location::ptr(
			new location(6, "FINAL", system_dynamics, invariant, false,
					Out_Going_Trans_fromLoc6));
	location::ptr l13 = location::ptr(
			new location(13, "BAD", system_dynamics, invariant, false,
					Out_Going_Trans_fromLoc13));

//  *********  Location ID=6 and ID=13 completed  *****************

	invariantBoundValue[0] = 2; //1<=x1<=2 and 1<=x2<=2
	invariantBoundValue[1] = -1;
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc7;
	Out_Going_Trans_fromLoc7.push_back(t14);
	Out_Going_Trans_fromLoc7.push_back(t15);
	Out_Going_Trans_fromLoc7.push_back(t16);
	Out_Going_Trans_fromLoc7.push_back(t17);
	location::ptr l7 = location::ptr(
			new location(7, "4", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc7));
	//  ************ Location ID=7 completed  ************

	invariantBoundValue[0] = 2; //1<=x1<=2 and 2<=x2<=3
	invariantBoundValue[1] = -1;
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc8;
	Out_Going_Trans_fromLoc8.push_back(t18);
	Out_Going_Trans_fromLoc8.push_back(t19);
	Out_Going_Trans_fromLoc8.push_back(t20);
	Out_Going_Trans_fromLoc8.push_back(t21);
	location::ptr l8 = location::ptr(
			new location(8, "4", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc8));
	//  ************ Location ID=8 completed  ************
	invariantBoundValue[0] = 2;	//1<=x1<=2 and 3<=x2<=4
	invariantBoundValue[1] = -1;
	invariantBoundValue[2] = 4;
	invariantBoundValue[3] = -3;

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
	std::list<transition::ptr> Out_Going_Trans_fromLoc9;
	Out_Going_Trans_fromLoc9.push_back(t22);
	Out_Going_Trans_fromLoc9.push_back(t23);
	Out_Going_Trans_fromLoc9.push_back(t24);
	Out_Going_Trans_fromLoc9.push_back(t25);
	location::ptr l9 = location::ptr(
			new location(9, "4", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc9));
	//  ************ Location ID=9 completed  ************

	invariantBoundValue[0] = 2; //1<=x1<=2 and 4<=x2<=5
	invariantBoundValue[1] = -1;
	invariantBoundValue[2] = 5;
	invariantBoundValue[3] = -4;

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
	std::list<transition::ptr> Out_Going_Trans_fromLoc10;
	Out_Going_Trans_fromLoc10.push_back(t26);
	Out_Going_Trans_fromLoc10.push_back(t27);
	Out_Going_Trans_fromLoc10.push_back(t28);
	location::ptr l10 = location::ptr(
			new location(10, "4", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc10));
	//  ************ Location ID=10 completed  ************

	invariantBoundValue[0] = 3; //2<=x1<=3 and 0<=x2<=1
	invariantBoundValue[1] = -2;
	invariantBoundValue[2] = 1;
	invariantBoundValue[3] = 0;

	boundValueV[0] = 0; //v1==0 and v2=1 => -Vd = -(0,1) u is (0,0,0,-1)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 0;
	boundValueV[5] = 0;
	boundValueV[6] = -1;
	boundValueV[7] = 1;

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
	std::list<transition::ptr> Out_Going_Trans_fromLoc11;
	Out_Going_Trans_fromLoc11.push_back(t29);
	Out_Going_Trans_fromLoc11.push_back(t30);
	Out_Going_Trans_fromLoc11.push_back(t31);
	location::ptr l11 = location::ptr(
			new location(11, "0", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc11));
	//  ************ Location ID=11 completed  ************

	invariantBoundValue[0] = 3; //2<=x1<=3 and 1<=x2<=2
	invariantBoundValue[1] = -2;
	invariantBoundValue[2] = 2;
	invariantBoundValue[3] = -1;

	boundValueV[0] = 0; //v1==-1 and v2=0 => -Vd = -(-1,0) u is (0,0,1,0)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 1;
	boundValueV[5] = -1;
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc12;
	Out_Going_Trans_fromLoc12.push_back(t32);
	Out_Going_Trans_fromLoc12.push_back(t33);
	Out_Going_Trans_fromLoc12.push_back(t34);
	Out_Going_Trans_fromLoc12.push_back(t35);
	location::ptr l12 = location::ptr(
			new location(12, "6", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc12));
	//  ************ Location ID=12 completed  ************
	//  ************ Location ID=13 Done above  ************

	invariantBoundValue[0] = 3; //2<=x1<=3 and 3<=x2<=4
	invariantBoundValue[1] = -2;
	invariantBoundValue[2] = 4;
	invariantBoundValue[3] = -3;

	boundValueV[0] = 0; //v1=-0.70711 and v2=0.70711 => -Vd = -(-0.70711,0.70711) u is (0,0,0.70711,-0.70711)
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc14;
	Out_Going_Trans_fromLoc14.push_back(t36);
	Out_Going_Trans_fromLoc14.push_back(t37);
	Out_Going_Trans_fromLoc14.push_back(t38);
	Out_Going_Trans_fromLoc14.push_back(t39);
	location::ptr l14 = location::ptr(
			new location(14, "7", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc14));
	//  ************ Location ID=14 completed  ************

	invariantBoundValue[0] = 3; //2<=x1<=3 and 4<=x2<=5
	invariantBoundValue[1] = -2;
	invariantBoundValue[2] = 5;
	invariantBoundValue[3] = -4;

	boundValueV[0] = 0; //v1==-1 and v2=0 => -Vd = -(-1,0) u is (0,0,1,0)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 1;
	boundValueV[5] = -1;
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc15;
	Out_Going_Trans_fromLoc15.push_back(t40);
	Out_Going_Trans_fromLoc15.push_back(t41);
	Out_Going_Trans_fromLoc15.push_back(t42);
	location::ptr l15 = location::ptr(
			new location(15, "6", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc15));
	//  ************ Location ID=15 completed  ************

	invariantBoundValue[0] = 4; //3<=x1<=4 and 0<=x2<=1
	invariantBoundValue[1] = -3;
	invariantBoundValue[2] = 1;
	invariantBoundValue[3] = 0;

	boundValueV[0] = 0; //v1==0 and v2=1 => -Vd = -(0,1) u is (0,0,0,-1)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 0;
	boundValueV[5] = 0;
	boundValueV[6] = -1;
	boundValueV[7] = 1;

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
	std::list<transition::ptr> Out_Going_Trans_fromLoc16;
	Out_Going_Trans_fromLoc16.push_back(t43);
	Out_Going_Trans_fromLoc16.push_back(t44);
	Out_Going_Trans_fromLoc16.push_back(t45);
	location::ptr l16 = location::ptr(
			new location(16, "0", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc16));
	//  ************ Location ID=16 completed  ************

	invariantBoundValue[0] = 4; //3<=x1<=4 and 1<=x2<=2
	invariantBoundValue[1] = -3;
	invariantBoundValue[2] = 2;
	invariantBoundValue[3] = -1;

	boundValueV[0] = 0; //v1==-1 and v2=0 => -Vd = -(-1,0) u is (0,0,1,0)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 1;
	boundValueV[5] = -1;
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc17;
	Out_Going_Trans_fromLoc17.push_back(t46);
	Out_Going_Trans_fromLoc17.push_back(t47);
	Out_Going_Trans_fromLoc17.push_back(t48);
	Out_Going_Trans_fromLoc17.push_back(t49);
	location::ptr l17 = location::ptr(
			new location(17, "6", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc17));
	//  ************ Location ID=17 completed  ************

	invariantBoundValue[0] = 4; //3<=x1<=4 and 2<=x2<=3
	invariantBoundValue[1] = -3;
	invariantBoundValue[2] = 3;
	invariantBoundValue[3] = -2;

	boundValueV[0] = 0; //v1=0.70711 and v2=-0.70711 => -Vd = -(0.70711, -0.70711) u is (0,0, -0.70711, 0.70711)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = -0.70711;
	boundValueV[5] = 0.70711;
	boundValueV[6] = 0.70711;
	boundValueV[7] = -0.70711;

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
	std::list<transition::ptr> Out_Going_Trans_fromLoc18;
	Out_Going_Trans_fromLoc18.push_back(t50);
	Out_Going_Trans_fromLoc18.push_back(t51);
	Out_Going_Trans_fromLoc18.push_back(t52);
	Out_Going_Trans_fromLoc18.push_back(t53);
	location::ptr l18 = location::ptr(
			new location(18, "3", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc18));
	//  ************ Location ID=18 completed  ************

	invariantBoundValue[0] = 4; //3<=x1<=4 and 3<=x2<=4
	invariantBoundValue[1] = -3;
	invariantBoundValue[2] = 4;
	invariantBoundValue[3] = -3;

	boundValueV[0] = 0; //v1=-0.70711 and v2=0.70711 => -Vd = -(-0.70711,0.70711) u is (0,0,0.70711,-0.70711)
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc19;
	Out_Going_Trans_fromLoc19.push_back(t54);
	Out_Going_Trans_fromLoc19.push_back(t55);
	Out_Going_Trans_fromLoc19.push_back(t56);
	Out_Going_Trans_fromLoc19.push_back(t57);
	location::ptr l19 = location::ptr(
			new location(19, "7", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc19));
	//  ************ Location ID=19 completed  ************

	invariantBoundValue[0] = 4; //3<=x1<=4 and 4<=x2<=5
	invariantBoundValue[1] = -3;
	invariantBoundValue[2] = 5;
	invariantBoundValue[3] = -4;

	boundValueV[0] = 0; //v1==-1 and v2=0 => -Vd = -(-1,0) u is (0,0,1,0)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 1;
	boundValueV[5] = -1;
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc20;
	Out_Going_Trans_fromLoc20.push_back(t58);
	Out_Going_Trans_fromLoc20.push_back(t59);
	Out_Going_Trans_fromLoc20.push_back(t60);
	location::ptr l20 = location::ptr(
			new location(20, "6", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc20));
	//  ************ Location ID=20 completed  ************

	invariantBoundValue[0] = 4; //3<=x1<=4 and 0<=x2<=1
	invariantBoundValue[1] = -3;
	invariantBoundValue[2] = 1;
	invariantBoundValue[3] = 0;

	boundValueV[0] = 0; //v1==0 and v2=1 => -Vd = -(0,1) u is (0,0,0,-1)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 0;
	boundValueV[5] = 0;
	boundValueV[6] = -1;
	boundValueV[7] = 1;

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
	std::list<transition::ptr> Out_Going_Trans_fromLoc21;
	Out_Going_Trans_fromLoc21.push_back(t61);
	Out_Going_Trans_fromLoc21.push_back(t62);
	location::ptr l21 = location::ptr(
			new location(21, "0", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc21));
	//  ************ Location ID=21 completed  ************

	invariantBoundValue[0] = 5; //4<=x1<=5 and 1<=x2<=2
	invariantBoundValue[1] = -4;
	invariantBoundValue[2] = 2;
	invariantBoundValue[3] = -1;

	boundValueV[0] = 0; //v1==-1 and v2=0 => -Vd = -(-1,0) u is (0,0,1,0)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 1;
	boundValueV[5] = -1;
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc22;
	Out_Going_Trans_fromLoc22.push_back(t63);
	Out_Going_Trans_fromLoc22.push_back(t64);
	Out_Going_Trans_fromLoc22.push_back(t65);
	location::ptr l22 = location::ptr(
			new location(22, "6", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc22));
	//  ************ Location ID=22 completed  ************

	invariantBoundValue[0] = 5; //4<=x1<=5 and 2<=x2<=3
	invariantBoundValue[1] = -4;
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc23;
	Out_Going_Trans_fromLoc23.push_back(t66);
	Out_Going_Trans_fromLoc23.push_back(t67);
	Out_Going_Trans_fromLoc23.push_back(t68);
	location::ptr l23 = location::ptr(
			new location(23, "4", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc23));
	//  ************ Location ID=23 completed  ************

	invariantBoundValue[0] = 5; //4<=x1<=5 and 3<=x2<=4
	invariantBoundValue[1] = -4;
	invariantBoundValue[2] = 4;
	invariantBoundValue[3] = -3;

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
	std::list<transition::ptr> Out_Going_Trans_fromLoc24;
	Out_Going_Trans_fromLoc24.push_back(t69);
	Out_Going_Trans_fromLoc24.push_back(t70);
	Out_Going_Trans_fromLoc24.push_back(t71);
	location::ptr l24 = location::ptr(
			new location(24, "4", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc24));
	//  ************ Location ID=24 completed  ************

	invariantBoundValue[0] = 5; //4<=x1<=5 and 4<=x2<=5
	invariantBoundValue[1] = -4;
	invariantBoundValue[2] = 5;
	invariantBoundValue[3] = -4;

	boundValueV[0] = 0; //v1==-1 and v2=0 => -Vd = -(-1,0) u is (0,0,1,0)
	boundValueV[1] = 0;
	boundValueV[2] = 0;
	boundValueV[3] = 0;
	boundValueV[4] = 1;
	boundValueV[5] = -1;
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
	std::list<transition::ptr> Out_Going_Trans_fromLoc25;
	Out_Going_Trans_fromLoc25.push_back(t72);
	Out_Going_Trans_fromLoc25.push_back(t73);
	location::ptr l25 = location::ptr(
			new location(25, "6", system_dynamics, invariant, true,
					Out_Going_Trans_fromLoc25));
	//  ************ Location ID=25 completed  ************
// ******************************************************************************************
//	// ************************* Locations Created and Initialized **************************
// ******************************************************************************************

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
	Hybrid_Automata.addLocation(l10);
	Hybrid_Automata.addLocation(l11);
	Hybrid_Automata.addLocation(l12);
	Hybrid_Automata.addLocation(l13);
	Hybrid_Automata.addLocation(l14);
	Hybrid_Automata.addLocation(l15);
	Hybrid_Automata.addLocation(l16);
	Hybrid_Automata.addLocation(l17);
	Hybrid_Automata.addLocation(l18);
	Hybrid_Automata.addLocation(l19);
	Hybrid_Automata.addLocation(l20);
	Hybrid_Automata.addLocation(l21);
	Hybrid_Automata.addLocation(l22);
	Hybrid_Automata.addLocation(l23);
	Hybrid_Automata.addLocation(l24);
	Hybrid_Automata.addLocation(l25);

	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));

	init_state_list.push_back(I);

}

//Hyst generated Output   HystCode
void SetNavigationModel5by5(hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state_list,
		ReachabilityParameters& reach_parameters) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;

	polytope::ptr initial_polytope_I;

	polytope::ptr invariant0, invariant1, invariant2, invariant3, invariant4,
			invariant5, invariant6, invariant7, invariant8, invariant9,
			invariant10, invariant11, invariant12, invariant13, invariant14,
			invariant15, invariant16, invariant17, invariant18, invariant19,
			invariant20, invariant21, invariant22, invariant23, invariant24;

	polytope::ptr gaurd_polytope0, gaurd_polytope1, gaurd_polytope2,
			gaurd_polytope3, gaurd_polytope4, gaurd_polytope5, gaurd_polytope6,
			gaurd_polytope7, gaurd_polytope8, gaurd_polytope9, gaurd_polytope10,
			gaurd_polytope11, gaurd_polytope12, gaurd_polytope13,
			gaurd_polytope14, gaurd_polytope15, gaurd_polytope16,
			gaurd_polytope17, gaurd_polytope18, gaurd_polytope19,
			gaurd_polytope20, gaurd_polytope21, gaurd_polytope22,
			gaurd_polytope23, gaurd_polytope24, gaurd_polytope25,
			gaurd_polytope26, gaurd_polytope27, gaurd_polytope28,
			gaurd_polytope29, gaurd_polytope30, gaurd_polytope31,
			gaurd_polytope32, gaurd_polytope33, gaurd_polytope34,
			gaurd_polytope35, gaurd_polytope36, gaurd_polytope37,
			gaurd_polytope38, gaurd_polytope39, gaurd_polytope40,
			gaurd_polytope41, gaurd_polytope42, gaurd_polytope43,
			gaurd_polytope44, gaurd_polytope45, gaurd_polytope46,
			gaurd_polytope47, gaurd_polytope48, gaurd_polytope49,
			gaurd_polytope50, gaurd_polytope51, gaurd_polytope52,
			gaurd_polytope53, gaurd_polytope54, gaurd_polytope55,
			gaurd_polytope56, gaurd_polytope57, gaurd_polytope58,
			gaurd_polytope59, gaurd_polytope60, gaurd_polytope61,
			gaurd_polytope62, gaurd_polytope63, gaurd_polytope64,
			gaurd_polytope65, gaurd_polytope66, gaurd_polytope67,
			gaurd_polytope68, gaurd_polytope69, gaurd_polytope70,
			gaurd_polytope71, gaurd_polytope72;

	Dynamics system_dynamics0, system_dynamics1, system_dynamics2,
			system_dynamics3, system_dynamics4, system_dynamics5,
			system_dynamics6, system_dynamics7, system_dynamics8,
			system_dynamics9, system_dynamics10, system_dynamics11,
			system_dynamics12, system_dynamics13, system_dynamics14,
			system_dynamics15, system_dynamics16, system_dynamics17,
			system_dynamics18, system_dynamics19, system_dynamics20,
			system_dynamics21, system_dynamics22, system_dynamics23,
			system_dynamics24;

	math::matrix<double> ConstraintsMatrixI, ConstraintsMatrixV0,
			ConstraintsMatrixV1, ConstraintsMatrixV2, ConstraintsMatrixV3,
			ConstraintsMatrixV4, ConstraintsMatrixV5, ConstraintsMatrixV6,
			ConstraintsMatrixV7, ConstraintsMatrixV8, ConstraintsMatrixV9,
			ConstraintsMatrixV10, ConstraintsMatrixV11, ConstraintsMatrixV12,
			ConstraintsMatrixV13, ConstraintsMatrixV14, ConstraintsMatrixV15,
			ConstraintsMatrixV16, ConstraintsMatrixV17, ConstraintsMatrixV18,
			ConstraintsMatrixV19, ConstraintsMatrixV20, ConstraintsMatrixV21,
			ConstraintsMatrixV22, ConstraintsMatrixV23, ConstraintsMatrixV24,
			invariantConstraintsMatrix0, invariantConstraintsMatrix1,
			invariantConstraintsMatrix2, invariantConstraintsMatrix3,
			invariantConstraintsMatrix4, invariantConstraintsMatrix5,
			invariantConstraintsMatrix6, invariantConstraintsMatrix7,
			invariantConstraintsMatrix8, invariantConstraintsMatrix9,
			invariantConstraintsMatrix10, invariantConstraintsMatrix11,
			invariantConstraintsMatrix12, invariantConstraintsMatrix13,
			invariantConstraintsMatrix14, invariantConstraintsMatrix15,
			invariantConstraintsMatrix16, invariantConstraintsMatrix17,
			invariantConstraintsMatrix18, invariantConstraintsMatrix19,
			invariantConstraintsMatrix20, invariantConstraintsMatrix21,
			invariantConstraintsMatrix22, invariantConstraintsMatrix23,
			invariantConstraintsMatrix24, gaurdConstraintsMatrix0,
			gaurdConstraintsMatrix1, gaurdConstraintsMatrix2,
			gaurdConstraintsMatrix3, gaurdConstraintsMatrix4,
			gaurdConstraintsMatrix5, gaurdConstraintsMatrix6,
			gaurdConstraintsMatrix7, gaurdConstraintsMatrix8,
			gaurdConstraintsMatrix9, gaurdConstraintsMatrix10,
			gaurdConstraintsMatrix11, gaurdConstraintsMatrix12,
			gaurdConstraintsMatrix13, gaurdConstraintsMatrix14,
			gaurdConstraintsMatrix15, gaurdConstraintsMatrix16,
			gaurdConstraintsMatrix17, gaurdConstraintsMatrix18,
			gaurdConstraintsMatrix19, gaurdConstraintsMatrix20,
			gaurdConstraintsMatrix21, gaurdConstraintsMatrix22,
			gaurdConstraintsMatrix23, gaurdConstraintsMatrix24,
			gaurdConstraintsMatrix25, gaurdConstraintsMatrix26,
			gaurdConstraintsMatrix27, gaurdConstraintsMatrix28,
			gaurdConstraintsMatrix29, gaurdConstraintsMatrix30,
			gaurdConstraintsMatrix31, gaurdConstraintsMatrix32,
			gaurdConstraintsMatrix33, gaurdConstraintsMatrix34,
			gaurdConstraintsMatrix35, gaurdConstraintsMatrix36,
			gaurdConstraintsMatrix37, gaurdConstraintsMatrix38,
			gaurdConstraintsMatrix39, gaurdConstraintsMatrix40,
			gaurdConstraintsMatrix41, gaurdConstraintsMatrix42,
			gaurdConstraintsMatrix43, gaurdConstraintsMatrix44,
			gaurdConstraintsMatrix45, gaurdConstraintsMatrix46,
			gaurdConstraintsMatrix47, gaurdConstraintsMatrix48,
			gaurdConstraintsMatrix49, gaurdConstraintsMatrix50,
			gaurdConstraintsMatrix51, gaurdConstraintsMatrix52,
			gaurdConstraintsMatrix53, gaurdConstraintsMatrix54,
			gaurdConstraintsMatrix55, gaurdConstraintsMatrix56,
			gaurdConstraintsMatrix57, gaurdConstraintsMatrix58,
			gaurdConstraintsMatrix59, gaurdConstraintsMatrix60,
			gaurdConstraintsMatrix61, gaurdConstraintsMatrix62,
			gaurdConstraintsMatrix63, gaurdConstraintsMatrix64,
			gaurdConstraintsMatrix65, gaurdConstraintsMatrix66,
			gaurdConstraintsMatrix67, gaurdConstraintsMatrix68,
			gaurdConstraintsMatrix69, gaurdConstraintsMatrix70,
			gaurdConstraintsMatrix71, gaurdConstraintsMatrix72, A0matrix,
			A1matrix, A2matrix, A3matrix, A4matrix, A5matrix, A6matrix,
			A7matrix, A8matrix, A9matrix, A10matrix, A11matrix, A12matrix,
			A13matrix, A14matrix, A15matrix, A16matrix, A17matrix, A18matrix,
			A19matrix, A20matrix, A21matrix, A22matrix, A23matrix, A24matrix,
			Bmatrix0, Bmatrix1, Bmatrix2, Bmatrix3, Bmatrix4, Bmatrix5,
			Bmatrix6, Bmatrix7, Bmatrix8, Bmatrix9, Bmatrix10, Bmatrix11,
			Bmatrix12, Bmatrix13, Bmatrix14, Bmatrix15, Bmatrix16, Bmatrix17,
			Bmatrix18, Bmatrix19, Bmatrix20, Bmatrix21, Bmatrix22, Bmatrix23,
			Bmatrix24;

	std::vector<double> boundValueI, boundValueV0, boundValueV1, boundValueV2,
			boundValueV3, boundValueV4, boundValueV5, boundValueV6,
			boundValueV7, boundValueV8, boundValueV9, boundValueV10,
			boundValueV11, boundValueV12, boundValueV13, boundValueV14,
			boundValueV15, boundValueV16, boundValueV17, boundValueV18,
			boundValueV19, boundValueV20, boundValueV21, boundValueV22,
			boundValueV23, boundValueV24, C0, C1, C2, C3, C4, C5, C6, C7, C8,
			C9, C10, C11, C12, C13, C14, C15, C16, C17, C18, C19, C20, C21, C22,
			C23, C24, invariantBoundValue0, invariantBoundValue1,
			invariantBoundValue2, invariantBoundValue3, invariantBoundValue4,
			invariantBoundValue5, invariantBoundValue6, invariantBoundValue7,
			invariantBoundValue8, invariantBoundValue9, invariantBoundValue10,
			invariantBoundValue11, invariantBoundValue12, invariantBoundValue13,
			invariantBoundValue14, invariantBoundValue15, invariantBoundValue16,
			invariantBoundValue17, invariantBoundValue18, invariantBoundValue19,
			invariantBoundValue20, invariantBoundValue21, invariantBoundValue22,
			invariantBoundValue23, invariantBoundValue24, gaurdBoundValue0,
			gaurdBoundValue1, gaurdBoundValue2, gaurdBoundValue3,
			gaurdBoundValue4, gaurdBoundValue5, gaurdBoundValue6,
			gaurdBoundValue7, gaurdBoundValue8, gaurdBoundValue9,
			gaurdBoundValue10, gaurdBoundValue11, gaurdBoundValue12,
			gaurdBoundValue13, gaurdBoundValue14, gaurdBoundValue15,
			gaurdBoundValue16, gaurdBoundValue17, gaurdBoundValue18,
			gaurdBoundValue19, gaurdBoundValue20, gaurdBoundValue21,
			gaurdBoundValue22, gaurdBoundValue23, gaurdBoundValue24,
			gaurdBoundValue25, gaurdBoundValue26, gaurdBoundValue27,
			gaurdBoundValue28, gaurdBoundValue29, gaurdBoundValue30,
			gaurdBoundValue31, gaurdBoundValue32, gaurdBoundValue33,
			gaurdBoundValue34, gaurdBoundValue35, gaurdBoundValue36,
			gaurdBoundValue37, gaurdBoundValue38, gaurdBoundValue39,
			gaurdBoundValue40, gaurdBoundValue41, gaurdBoundValue42,
			gaurdBoundValue43, gaurdBoundValue44, gaurdBoundValue45,
			gaurdBoundValue46, gaurdBoundValue47, gaurdBoundValue48,
			gaurdBoundValue49, gaurdBoundValue50, gaurdBoundValue51,
			gaurdBoundValue52, gaurdBoundValue53, gaurdBoundValue54,
			gaurdBoundValue55, gaurdBoundValue56, gaurdBoundValue57,
			gaurdBoundValue58, gaurdBoundValue59, gaurdBoundValue60,
			gaurdBoundValue61, gaurdBoundValue62, gaurdBoundValue63,
			gaurdBoundValue64, gaurdBoundValue65, gaurdBoundValue66,
			gaurdBoundValue67, gaurdBoundValue68, gaurdBoundValue69,
			gaurdBoundValue70, gaurdBoundValue71, gaurdBoundValue72;

	int boundSignI, invariantBoundSign, gaurdBoundSign, boundSignV;

	size_type row, col;

	// The mode name is  loc3

	row = 4;
	col = 4;
	A0matrix.resize(row, col);
	A0matrix(0, 0) = 0.0;
	A0matrix(0, 1) = 0.0;
	A0matrix(0, 2) = 1.0;
	A0matrix(0, 3) = 0.0;
	A0matrix(1, 0) = 0.0;
	A0matrix(1, 1) = 0.0;
	A0matrix(1, 2) = 0.0;
	A0matrix(1, 3) = 1.0;
	A0matrix(2, 0) = 0.0;
	A0matrix(2, 1) = 0.0;
	A0matrix(2, 2) = -0.8;
	A0matrix(2, 3) = -0.2;
	A0matrix(3, 0) = 0.0;
	A0matrix(3, 1) = 0.0;
	A0matrix(3, 2) = -0.1;
	A0matrix(3, 3) = -0.8;
	system_dynamics0.isEmptyMatrixA = false;
	system_dynamics0.MatrixA = A0matrix;

	system_dynamics0.isEmptyMatrixB = true;

	C0.resize(row);
	C0[0] = 0.0;
	C0[1] = 0.0;
	C0[2] = 0.8;
	C0[3] = 0.1;
	system_dynamics0.isEmptyC = false;
	system_dynamics0.C = C0;

	row = 4;
	col = 4;
	invariantConstraintsMatrix0.resize(row, col);
	invariantConstraintsMatrix0(0, 0) = -1.0;
	invariantConstraintsMatrix0(0, 1) = 0.0;
	invariantConstraintsMatrix0(0, 2) = 0.0;
	invariantConstraintsMatrix0(0, 3) = 0.0;
	invariantConstraintsMatrix0(1, 0) = 1.0;
	invariantConstraintsMatrix0(1, 1) = 0.0;
	invariantConstraintsMatrix0(1, 2) = 0.0;
	invariantConstraintsMatrix0(1, 3) = 0.0;
	invariantConstraintsMatrix0(2, 0) = 0.0;
	invariantConstraintsMatrix0(2, 1) = -1.0;
	invariantConstraintsMatrix0(2, 2) = 0.0;
	invariantConstraintsMatrix0(2, 3) = 0.0;
	invariantConstraintsMatrix0(3, 0) = 0.0;
	invariantConstraintsMatrix0(3, 1) = 1.0;
	invariantConstraintsMatrix0(3, 2) = 0.0;
	invariantConstraintsMatrix0(3, 3) = 0.0;

	invariantBoundValue0.resize(row);
	invariantBoundValue0[0] = -0.0;
	invariantBoundValue0[1] = 1.0;
	invariantBoundValue0[2] = -2.0;
	invariantBoundValue0[3] = 3.0;
	invariantBoundSign = 1;
	invariant0 = polytope::ptr(
			new polytope(invariantConstraintsMatrix0, invariantBoundValue0,
					invariantBoundSign));

	system_dynamics0.U = polytope::ptr(new polytope(true));

	// The mode name is  loc8

	row = 4;
	col = 4;
	A1matrix.resize(row, col);
	A1matrix(0, 0) = 0.0;
	A1matrix(0, 1) = 0.0;
	A1matrix(0, 2) = 1.0;
	A1matrix(0, 3) = 0.0;
	A1matrix(1, 0) = 0.0;
	A1matrix(1, 1) = 0.0;
	A1matrix(1, 2) = 0.0;
	A1matrix(1, 3) = 1.0;
	A1matrix(2, 0) = 0.0;
	A1matrix(2, 1) = 0.0;
	A1matrix(2, 2) = -0.8;
	A1matrix(2, 3) = -0.2;
	A1matrix(3, 0) = 0.0;
	A1matrix(3, 1) = 0.0;
	A1matrix(3, 2) = -0.1;
	A1matrix(3, 3) = -0.8;
	system_dynamics1.isEmptyMatrixA = false;
	system_dynamics1.MatrixA = A1matrix;

	system_dynamics1.isEmptyMatrixB = true;

	C1.resize(row);
	C1[0] = 0.0;
	C1[1] = 0.0;
	C1[2] = -0.2;
	C1[3] = -0.8;
	system_dynamics1.isEmptyC = false;
	system_dynamics1.C = C1;

	row = 4;
	col = 4;
	invariantConstraintsMatrix1.resize(row, col);
	invariantConstraintsMatrix1(0, 0) = -1.0;
	invariantConstraintsMatrix1(0, 1) = 0.0;
	invariantConstraintsMatrix1(0, 2) = 0.0;
	invariantConstraintsMatrix1(0, 3) = 0.0;
	invariantConstraintsMatrix1(1, 0) = 1.0;
	invariantConstraintsMatrix1(1, 1) = 0.0;
	invariantConstraintsMatrix1(1, 2) = 0.0;
	invariantConstraintsMatrix1(1, 3) = 0.0;
	invariantConstraintsMatrix1(2, 0) = 0.0;
	invariantConstraintsMatrix1(2, 1) = -1.0;
	invariantConstraintsMatrix1(2, 2) = 0.0;
	invariantConstraintsMatrix1(2, 3) = 0.0;
	invariantConstraintsMatrix1(3, 0) = 0.0;
	invariantConstraintsMatrix1(3, 1) = 1.0;
	invariantConstraintsMatrix1(3, 2) = 0.0;
	invariantConstraintsMatrix1(3, 3) = 0.0;

	invariantBoundValue1.resize(row);
	invariantBoundValue1[0] = -1.0;
	invariantBoundValue1[1] = 2.0;
	invariantBoundValue1[2] = -2.0;
	invariantBoundValue1[3] = 3.0;
	invariantBoundSign = 1;
	invariant1 = polytope::ptr(
			new polytope(invariantConstraintsMatrix1, invariantBoundValue1,
					invariantBoundSign));

	system_dynamics1.U = polytope::ptr(new polytope(true));

	// The mode name is  loc13

	system_dynamics2.isEmptyMatrixA = true;

	system_dynamics2.isEmptyMatrixB = true;

	system_dynamics2.isEmptyC = true;

	invariantBoundSign = 1;
	invariant2 = polytope::ptr(
			new polytope(invariantConstraintsMatrix2, invariantBoundValue2,
					invariantBoundSign));
	invariant2->setIsUniverse(true);

	row = 4;
	col = 4;
	ConstraintsMatrixV2.resize(row, col);
	ConstraintsMatrixV2(0, 0) = -1.0;
	ConstraintsMatrixV2(0, 1) = 0.0;
	ConstraintsMatrixV2(0, 2) = 0.0;
	ConstraintsMatrixV2(0, 3) = 0.0;
	ConstraintsMatrixV2(1, 0) = 1.0;
	ConstraintsMatrixV2(1, 1) = 0.0;
	ConstraintsMatrixV2(1, 2) = 0.0;
	ConstraintsMatrixV2(1, 3) = 0.0;
	ConstraintsMatrixV2(2, 0) = 0.0;
	ConstraintsMatrixV2(2, 1) = -1.0;
	ConstraintsMatrixV2(2, 2) = 0.0;
	ConstraintsMatrixV2(2, 3) = 0.0;
	ConstraintsMatrixV2(3, 0) = 0.0;
	ConstraintsMatrixV2(3, 1) = 1.0;
	ConstraintsMatrixV2(3, 2) = 0.0;
	ConstraintsMatrixV2(3, 3) = 0.0;

	boundValueV2.resize(row);
	boundValueV2[0] = -2.0;
	boundValueV2[1] = 3.0;
	boundValueV2[2] = -2.0;
	boundValueV2[3] = 3.0;
	boundSignV = 1;
	system_dynamics2.U = polytope::ptr(
			new polytope(ConstraintsMatrixV2, boundValueV2, boundSignV));

	// The mode name is  loc4

	row = 4;
	col = 4;
	A3matrix.resize(row, col);
	A3matrix(0, 0) = 0.0;
	A3matrix(0, 1) = 0.0;
	A3matrix(0, 2) = 1.0;
	A3matrix(0, 3) = 0.0;
	A3matrix(1, 0) = 0.0;
	A3matrix(1, 1) = 0.0;
	A3matrix(1, 2) = 0.0;
	A3matrix(1, 3) = 1.0;
	A3matrix(2, 0) = 0.0;
	A3matrix(2, 1) = 0.0;
	A3matrix(2, 2) = -0.8;
	A3matrix(2, 3) = -0.2;
	A3matrix(3, 0) = 0.0;
	A3matrix(3, 1) = 0.0;
	A3matrix(3, 2) = -0.1;
	A3matrix(3, 3) = -0.8;
	system_dynamics3.isEmptyMatrixA = false;
	system_dynamics3.MatrixA = A3matrix;

	system_dynamics3.isEmptyMatrixB = true;

	C3.resize(row);
	C3[0] = 0.0;
	C3[1] = 0.0;
	C3[2] = 0.8;
	C3[3] = 0.1;
	system_dynamics3.isEmptyC = false;
	system_dynamics3.C = C3;

	row = 4;
	col = 4;
	invariantConstraintsMatrix3.resize(row, col);
	invariantConstraintsMatrix3(0, 0) = -1.0;
	invariantConstraintsMatrix3(0, 1) = 0.0;
	invariantConstraintsMatrix3(0, 2) = 0.0;
	invariantConstraintsMatrix3(0, 3) = 0.0;
	invariantConstraintsMatrix3(1, 0) = 1.0;
	invariantConstraintsMatrix3(1, 1) = 0.0;
	invariantConstraintsMatrix3(1, 2) = 0.0;
	invariantConstraintsMatrix3(1, 3) = 0.0;
	invariantConstraintsMatrix3(2, 0) = 0.0;
	invariantConstraintsMatrix3(2, 1) = -1.0;
	invariantConstraintsMatrix3(2, 2) = 0.0;
	invariantConstraintsMatrix3(2, 3) = 0.0;
	invariantConstraintsMatrix3(3, 0) = 0.0;
	invariantConstraintsMatrix3(3, 1) = 1.0;
	invariantConstraintsMatrix3(3, 2) = 0.0;
	invariantConstraintsMatrix3(3, 3) = 0.0;

	invariantBoundValue3.resize(row);
	invariantBoundValue3[0] = -0.0;
	invariantBoundValue3[1] = 1.0;
	invariantBoundValue3[2] = -1.0;
	invariantBoundValue3[3] = 2.0;
	invariantBoundSign = 1;
	invariant3 = polytope::ptr(
			new polytope(invariantConstraintsMatrix3, invariantBoundValue3,
					invariantBoundSign));

	system_dynamics3.U = polytope::ptr(new polytope(true));

	// The mode name is  loc7

	row = 4;
	col = 4;
	A4matrix.resize(row, col);
	A4matrix(0, 0) = 0.0;
	A4matrix(0, 1) = 0.0;
	A4matrix(0, 2) = 1.0;
	A4matrix(0, 3) = 0.0;
	A4matrix(1, 0) = 0.0;
	A4matrix(1, 1) = 0.0;
	A4matrix(1, 2) = 0.0;
	A4matrix(1, 3) = 1.0;
	A4matrix(2, 0) = 0.0;
	A4matrix(2, 1) = 0.0;
	A4matrix(2, 2) = -0.8;
	A4matrix(2, 3) = -0.2;
	A4matrix(3, 0) = 0.0;
	A4matrix(3, 1) = 0.0;
	A4matrix(3, 2) = -0.1;
	A4matrix(3, 3) = -0.8;
	system_dynamics4.isEmptyMatrixA = false;
	system_dynamics4.MatrixA = A4matrix;

	system_dynamics4.isEmptyMatrixB = true;

	C4.resize(row);
	C4[0] = 0.0;
	C4[1] = 0.0;
	C4[2] = -0.2;
	C4[3] = -0.8;
	system_dynamics4.isEmptyC = false;
	system_dynamics4.C = C4;

	row = 4;
	col = 4;
	invariantConstraintsMatrix4.resize(row, col);
	invariantConstraintsMatrix4(0, 0) = -1.0;
	invariantConstraintsMatrix4(0, 1) = 0.0;
	invariantConstraintsMatrix4(0, 2) = 0.0;
	invariantConstraintsMatrix4(0, 3) = 0.0;
	invariantConstraintsMatrix4(1, 0) = 1.0;
	invariantConstraintsMatrix4(1, 1) = 0.0;
	invariantConstraintsMatrix4(1, 2) = 0.0;
	invariantConstraintsMatrix4(1, 3) = 0.0;
	invariantConstraintsMatrix4(2, 0) = 0.0;
	invariantConstraintsMatrix4(2, 1) = -1.0;
	invariantConstraintsMatrix4(2, 2) = 0.0;
	invariantConstraintsMatrix4(2, 3) = 0.0;
	invariantConstraintsMatrix4(3, 0) = 0.0;
	invariantConstraintsMatrix4(3, 1) = 1.0;
	invariantConstraintsMatrix4(3, 2) = 0.0;
	invariantConstraintsMatrix4(3, 3) = 0.0;

	invariantBoundValue4.resize(row);
	invariantBoundValue4[0] = -1.0;
	invariantBoundValue4[1] = 2.0;
	invariantBoundValue4[2] = -1.0;
	invariantBoundValue4[3] = 2.0;
	invariantBoundSign = 1;
	invariant4 = polytope::ptr(
			new polytope(invariantConstraintsMatrix4, invariantBoundValue4,
					invariantBoundSign));

	system_dynamics4.U = polytope::ptr(new polytope(true));

	// The mode name is  loc12

	row = 4;
	col = 4;
	A5matrix.resize(row, col);
	A5matrix(0, 0) = 0.0;
	A5matrix(0, 1) = 0.0;
	A5matrix(0, 2) = 1.0;
	A5matrix(0, 3) = 0.0;
	A5matrix(1, 0) = 0.0;
	A5matrix(1, 1) = 0.0;
	A5matrix(1, 2) = 0.0;
	A5matrix(1, 3) = 1.0;
	A5matrix(2, 0) = 0.0;
	A5matrix(2, 1) = 0.0;
	A5matrix(2, 2) = -0.8;
	A5matrix(2, 3) = -0.2;
	A5matrix(3, 0) = 0.0;
	A5matrix(3, 1) = 0.0;
	A5matrix(3, 2) = -0.1;
	A5matrix(3, 3) = -0.8;
	system_dynamics5.isEmptyMatrixA = false;
	system_dynamics5.MatrixA = A5matrix;

	system_dynamics5.isEmptyMatrixB = true;

	C5.resize(row);
	C5[0] = 0.0;
	C5[1] = 0.0;
	C5[2] = -0.8;
	C5[3] = -0.1;
	system_dynamics5.isEmptyC = false;
	system_dynamics5.C = C5;

	row = 4;
	col = 4;
	invariantConstraintsMatrix5.resize(row, col);
	invariantConstraintsMatrix5(0, 0) = -1.0;
	invariantConstraintsMatrix5(0, 1) = 0.0;
	invariantConstraintsMatrix5(0, 2) = 0.0;
	invariantConstraintsMatrix5(0, 3) = 0.0;
	invariantConstraintsMatrix5(1, 0) = 1.0;
	invariantConstraintsMatrix5(1, 1) = 0.0;
	invariantConstraintsMatrix5(1, 2) = 0.0;
	invariantConstraintsMatrix5(1, 3) = 0.0;
	invariantConstraintsMatrix5(2, 0) = 0.0;
	invariantConstraintsMatrix5(2, 1) = -1.0;
	invariantConstraintsMatrix5(2, 2) = 0.0;
	invariantConstraintsMatrix5(2, 3) = 0.0;
	invariantConstraintsMatrix5(3, 0) = 0.0;
	invariantConstraintsMatrix5(3, 1) = 1.0;
	invariantConstraintsMatrix5(3, 2) = 0.0;
	invariantConstraintsMatrix5(3, 3) = 0.0;

	invariantBoundValue5.resize(row);
	invariantBoundValue5[0] = -2.0;
	invariantBoundValue5[1] = 3.0;
	invariantBoundValue5[2] = -1.0;
	invariantBoundValue5[3] = 2.0;
	invariantBoundSign = 1;
	invariant5 = polytope::ptr(
			new polytope(invariantConstraintsMatrix5, invariantBoundValue5,
					invariantBoundSign));

	system_dynamics5.U = polytope::ptr(new polytope(true));

	// The mode name is  loc6

	system_dynamics6.isEmptyMatrixA = true;

	system_dynamics6.isEmptyMatrixB = true;

	system_dynamics6.isEmptyC = true;

	invariantBoundSign = 1;
	invariant6 = polytope::ptr(
			new polytope(invariantConstraintsMatrix6, invariantBoundValue6,
					invariantBoundSign));
	invariant6->setIsUniverse(true);

	row = 4;
	col = 4;
	ConstraintsMatrixV6.resize(row, col);
	ConstraintsMatrixV6(0, 0) = -1.0;
	ConstraintsMatrixV6(0, 1) = 0.0;
	ConstraintsMatrixV6(0, 2) = 0.0;
	ConstraintsMatrixV6(0, 3) = 0.0;
	ConstraintsMatrixV6(1, 0) = 1.0;
	ConstraintsMatrixV6(1, 1) = 0.0;
	ConstraintsMatrixV6(1, 2) = 0.0;
	ConstraintsMatrixV6(1, 3) = 0.0;
	ConstraintsMatrixV6(2, 0) = 0.0;
	ConstraintsMatrixV6(2, 1) = -1.0;
	ConstraintsMatrixV6(2, 2) = 0.0;
	ConstraintsMatrixV6(2, 3) = 0.0;
	ConstraintsMatrixV6(3, 0) = 0.0;
	ConstraintsMatrixV6(3, 1) = 1.0;
	ConstraintsMatrixV6(3, 2) = 0.0;
	ConstraintsMatrixV6(3, 3) = 0.0;

	boundValueV6.resize(row);
	boundValueV6[0] = -1.0;
	boundValueV6[1] = 2.0;
	boundValueV6[2] = -0.0;
	boundValueV6[3] = 1.0;
	boundSignV = 1;
	system_dynamics6.U = polytope::ptr(
			new polytope(ConstraintsMatrixV6, boundValueV6, boundSignV));

	// The mode name is  loc11

	row = 4;
	col = 4;
	A7matrix.resize(row, col);
	A7matrix(0, 0) = 0.0;
	A7matrix(0, 1) = 0.0;
	A7matrix(0, 2) = 1.0;
	A7matrix(0, 3) = 0.0;
	A7matrix(1, 0) = 0.0;
	A7matrix(1, 1) = 0.0;
	A7matrix(1, 2) = 0.0;
	A7matrix(1, 3) = 1.0;
	A7matrix(2, 0) = 0.0;
	A7matrix(2, 1) = 0.0;
	A7matrix(2, 2) = -0.8;
	A7matrix(2, 3) = -0.2;
	A7matrix(3, 0) = 0.0;
	A7matrix(3, 1) = 0.0;
	A7matrix(3, 2) = -0.1;
	A7matrix(3, 3) = -0.8;
	system_dynamics7.isEmptyMatrixA = false;
	system_dynamics7.MatrixA = A7matrix;

	system_dynamics7.isEmptyMatrixB = true;

	C7.resize(row);
	C7[0] = 0.0;
	C7[1] = 0.0;
	C7[2] = 0.2;
	C7[3] = 0.8;
	system_dynamics7.isEmptyC = false;
	system_dynamics7.C = C7;

	row = 4;
	col = 4;
	invariantConstraintsMatrix7.resize(row, col);
	invariantConstraintsMatrix7(0, 0) = -1.0;
	invariantConstraintsMatrix7(0, 1) = 0.0;
	invariantConstraintsMatrix7(0, 2) = 0.0;
	invariantConstraintsMatrix7(0, 3) = 0.0;
	invariantConstraintsMatrix7(1, 0) = 1.0;
	invariantConstraintsMatrix7(1, 1) = 0.0;
	invariantConstraintsMatrix7(1, 2) = 0.0;
	invariantConstraintsMatrix7(1, 3) = 0.0;
	invariantConstraintsMatrix7(2, 0) = 0.0;
	invariantConstraintsMatrix7(2, 1) = -1.0;
	invariantConstraintsMatrix7(2, 2) = 0.0;
	invariantConstraintsMatrix7(2, 3) = 0.0;
	invariantConstraintsMatrix7(3, 0) = 0.0;
	invariantConstraintsMatrix7(3, 1) = 1.0;
	invariantConstraintsMatrix7(3, 2) = 0.0;
	invariantConstraintsMatrix7(3, 3) = 0.0;

	invariantBoundValue7.resize(row);
	invariantBoundValue7[0] = -2.0;
	invariantBoundValue7[1] = 3.0;
	invariantBoundValue7[2] = -0.0;
	invariantBoundValue7[3] = 1.0;
	invariantBoundSign = 1;
	invariant7 = polytope::ptr(
			new polytope(invariantConstraintsMatrix7, invariantBoundValue7,
					invariantBoundSign));

	system_dynamics7.U = polytope::ptr(new polytope(true));

	// The mode name is  loc5

	row = 4;
	col = 4;
	A8matrix.resize(row, col);
	A8matrix(0, 0) = 0.0;
	A8matrix(0, 1) = 0.0;
	A8matrix(0, 2) = 1.0;
	A8matrix(0, 3) = 0.0;
	A8matrix(1, 0) = 0.0;
	A8matrix(1, 1) = 0.0;
	A8matrix(1, 2) = 0.0;
	A8matrix(1, 3) = 1.0;
	A8matrix(2, 0) = 0.0;
	A8matrix(2, 1) = 0.0;
	A8matrix(2, 2) = -0.8;
	A8matrix(2, 3) = -0.2;
	A8matrix(3, 0) = 0.0;
	A8matrix(3, 1) = 0.0;
	A8matrix(3, 2) = -0.1;
	A8matrix(3, 3) = -0.8;
	system_dynamics8.isEmptyMatrixA = false;
	system_dynamics8.MatrixA = A8matrix;

	system_dynamics8.isEmptyMatrixB = true;

	C8.resize(row);
	C8[0] = 0.0;
	C8[1] = 0.0;
	C8[2] = 0.8;
	C8[3] = 0.1;
	system_dynamics8.isEmptyC = false;
	system_dynamics8.C = C8;

	row = 4;
	col = 4;
	invariantConstraintsMatrix8.resize(row, col);
	invariantConstraintsMatrix8(0, 0) = -1.0;
	invariantConstraintsMatrix8(0, 1) = 0.0;
	invariantConstraintsMatrix8(0, 2) = 0.0;
	invariantConstraintsMatrix8(0, 3) = 0.0;
	invariantConstraintsMatrix8(1, 0) = 1.0;
	invariantConstraintsMatrix8(1, 1) = 0.0;
	invariantConstraintsMatrix8(1, 2) = 0.0;
	invariantConstraintsMatrix8(1, 3) = 0.0;
	invariantConstraintsMatrix8(2, 0) = 0.0;
	invariantConstraintsMatrix8(2, 1) = -1.0;
	invariantConstraintsMatrix8(2, 2) = 0.0;
	invariantConstraintsMatrix8(2, 3) = 0.0;
	invariantConstraintsMatrix8(3, 0) = 0.0;
	invariantConstraintsMatrix8(3, 1) = 1.0;
	invariantConstraintsMatrix8(3, 2) = 0.0;
	invariantConstraintsMatrix8(3, 3) = 0.0;

	invariantBoundValue8.resize(row);
	invariantBoundValue8[0] = -0.0;
	invariantBoundValue8[1] = 1.0;
	invariantBoundValue8[2] = -0.0;
	invariantBoundValue8[3] = 1.0;
	invariantBoundSign = 1;
	invariant8 = polytope::ptr(
			new polytope(invariantConstraintsMatrix8, invariantBoundValue8,
					invariantBoundSign));

	system_dynamics8.U = polytope::ptr(new polytope(true));

	// The mode name is  loc2

	row = 4;
	col = 4;
	A9matrix.resize(row, col);
	A9matrix(0, 0) = 0.0;
	A9matrix(0, 1) = 0.0;
	A9matrix(0, 2) = 1.0;
	A9matrix(0, 3) = 0.0;
	A9matrix(1, 0) = 0.0;
	A9matrix(1, 1) = 0.0;
	A9matrix(1, 2) = 0.0;
	A9matrix(1, 3) = 1.0;
	A9matrix(2, 0) = 0.0;
	A9matrix(2, 1) = 0.0;
	A9matrix(2, 2) = -0.8;
	A9matrix(2, 3) = -0.2;
	A9matrix(3, 0) = 0.0;
	A9matrix(3, 1) = 0.0;
	A9matrix(3, 2) = -0.1;
	A9matrix(3, 3) = -0.8;
	system_dynamics9.isEmptyMatrixA = false;
	system_dynamics9.MatrixA = A9matrix;

	system_dynamics9.isEmptyMatrixB = true;

	C9.resize(row);
	C9[0] = 0.0;
	C9[1] = 0.0;
	C9[2] = 0.8;
	C9[3] = 0.1;
	system_dynamics9.isEmptyC = false;
	system_dynamics9.C = C9;

	row = 4;
	col = 4;
	invariantConstraintsMatrix9.resize(row, col);
	invariantConstraintsMatrix9(0, 0) = -1.0;
	invariantConstraintsMatrix9(0, 1) = 0.0;
	invariantConstraintsMatrix9(0, 2) = 0.0;
	invariantConstraintsMatrix9(0, 3) = 0.0;
	invariantConstraintsMatrix9(1, 0) = 1.0;
	invariantConstraintsMatrix9(1, 1) = 0.0;
	invariantConstraintsMatrix9(1, 2) = 0.0;
	invariantConstraintsMatrix9(1, 3) = 0.0;
	invariantConstraintsMatrix9(2, 0) = 0.0;
	invariantConstraintsMatrix9(2, 1) = -1.0;
	invariantConstraintsMatrix9(2, 2) = 0.0;
	invariantConstraintsMatrix9(2, 3) = 0.0;
	invariantConstraintsMatrix9(3, 0) = 0.0;
	invariantConstraintsMatrix9(3, 1) = 1.0;
	invariantConstraintsMatrix9(3, 2) = 0.0;
	invariantConstraintsMatrix9(3, 3) = 0.0;

	invariantBoundValue9.resize(row);
	invariantBoundValue9[0] = -0.0;
	invariantBoundValue9[1] = 1.0;
	invariantBoundValue9[2] = -3.0;
	invariantBoundValue9[3] = 4.0;
	invariantBoundSign = 1;
	invariant9 = polytope::ptr(
			new polytope(invariantConstraintsMatrix9, invariantBoundValue9,
					invariantBoundSign));

	system_dynamics9.U = polytope::ptr(new polytope(true));

	// The mode name is  loc1

	row = 4;
	col = 4;
	A10matrix.resize(row, col);
	A10matrix(0, 0) = 0.0;
	A10matrix(0, 1) = 0.0;
	A10matrix(0, 2) = 1.0;
	A10matrix(0, 3) = 0.0;
	A10matrix(1, 0) = 0.0;
	A10matrix(1, 1) = 0.0;
	A10matrix(1, 2) = 0.0;
	A10matrix(1, 3) = 1.0;
	A10matrix(2, 0) = 0.0;
	A10matrix(2, 1) = 0.0;
	A10matrix(2, 2) = -0.8;
	A10matrix(2, 3) = -0.2;
	A10matrix(3, 0) = 0.0;
	A10matrix(3, 1) = 0.0;
	A10matrix(3, 2) = -0.1;
	A10matrix(3, 3) = -0.8;
	system_dynamics10.isEmptyMatrixA = false;
	system_dynamics10.MatrixA = A10matrix;

	system_dynamics10.isEmptyMatrixB = true;

	C10.resize(row);
	C10[0] = 0.0;
	C10[1] = 0.0;
	C10[2] = 0.8;
	C10[3] = 0.1;
	system_dynamics10.isEmptyC = false;
	system_dynamics10.C = C10;

	row = 4;
	col = 4;
	invariantConstraintsMatrix10.resize(row, col);
	invariantConstraintsMatrix10(0, 0) = -1.0;
	invariantConstraintsMatrix10(0, 1) = 0.0;
	invariantConstraintsMatrix10(0, 2) = 0.0;
	invariantConstraintsMatrix10(0, 3) = 0.0;
	invariantConstraintsMatrix10(1, 0) = 1.0;
	invariantConstraintsMatrix10(1, 1) = 0.0;
	invariantConstraintsMatrix10(1, 2) = 0.0;
	invariantConstraintsMatrix10(1, 3) = 0.0;
	invariantConstraintsMatrix10(2, 0) = 0.0;
	invariantConstraintsMatrix10(2, 1) = -1.0;
	invariantConstraintsMatrix10(2, 2) = 0.0;
	invariantConstraintsMatrix10(2, 3) = 0.0;
	invariantConstraintsMatrix10(3, 0) = 0.0;
	invariantConstraintsMatrix10(3, 1) = 1.0;
	invariantConstraintsMatrix10(3, 2) = 0.0;
	invariantConstraintsMatrix10(3, 3) = 0.0;

	invariantBoundValue10.resize(row);
	invariantBoundValue10[0] = -0.0;
	invariantBoundValue10[1] = 1.0;
	invariantBoundValue10[2] = -4.0;
	invariantBoundValue10[3] = 5.0;
	invariantBoundSign = 1;
	invariant10 = polytope::ptr(
			new polytope(invariantConstraintsMatrix10, invariantBoundValue10,
					invariantBoundSign));

	system_dynamics10.U = polytope::ptr(new polytope(true));

	// The mode name is  loc10

	row = 4;
	col = 4;
	A11matrix.resize(row, col);
	A11matrix(0, 0) = 0.0;
	A11matrix(0, 1) = 0.0;
	A11matrix(0, 2) = 1.0;
	A11matrix(0, 3) = 0.0;
	A11matrix(1, 0) = 0.0;
	A11matrix(1, 1) = 0.0;
	A11matrix(1, 2) = 0.0;
	A11matrix(1, 3) = 1.0;
	A11matrix(2, 0) = 0.0;
	A11matrix(2, 1) = 0.0;
	A11matrix(2, 2) = -0.8;
	A11matrix(2, 3) = -0.2;
	A11matrix(3, 0) = 0.0;
	A11matrix(3, 1) = 0.0;
	A11matrix(3, 2) = -0.1;
	A11matrix(3, 3) = -0.8;
	system_dynamics11.isEmptyMatrixA = false;
	system_dynamics11.MatrixA = A11matrix;

	system_dynamics11.isEmptyMatrixB = true;

	C11.resize(row);
	C11[0] = 0.0;
	C11[1] = 0.0;
	C11[2] = -0.2;
	C11[3] = -0.8;
	system_dynamics11.isEmptyC = false;
	system_dynamics11.C = C11;

	row = 4;
	col = 4;
	invariantConstraintsMatrix11.resize(row, col);
	invariantConstraintsMatrix11(0, 0) = -1.0;
	invariantConstraintsMatrix11(0, 1) = 0.0;
	invariantConstraintsMatrix11(0, 2) = 0.0;
	invariantConstraintsMatrix11(0, 3) = 0.0;
	invariantConstraintsMatrix11(1, 0) = 1.0;
	invariantConstraintsMatrix11(1, 1) = 0.0;
	invariantConstraintsMatrix11(1, 2) = 0.0;
	invariantConstraintsMatrix11(1, 3) = 0.0;
	invariantConstraintsMatrix11(2, 0) = 0.0;
	invariantConstraintsMatrix11(2, 1) = -1.0;
	invariantConstraintsMatrix11(2, 2) = 0.0;
	invariantConstraintsMatrix11(2, 3) = 0.0;
	invariantConstraintsMatrix11(3, 0) = 0.0;
	invariantConstraintsMatrix11(3, 1) = 1.0;
	invariantConstraintsMatrix11(3, 2) = 0.0;
	invariantConstraintsMatrix11(3, 3) = 0.0;

	invariantBoundValue11.resize(row);
	invariantBoundValue11[0] = -1.0;
	invariantBoundValue11[1] = 2.0;
	invariantBoundValue11[2] = -4.0;
	invariantBoundValue11[3] = 5.0;
	invariantBoundSign = 1;
	invariant11 = polytope::ptr(
			new polytope(invariantConstraintsMatrix11, invariantBoundValue11,
					invariantBoundSign));

	system_dynamics11.U = polytope::ptr(new polytope(true));

	// The mode name is  loc9

	row = 4;
	col = 4;
	A12matrix.resize(row, col);
	A12matrix(0, 0) = 0.0;
	A12matrix(0, 1) = 0.0;
	A12matrix(0, 2) = 1.0;
	A12matrix(0, 3) = 0.0;
	A12matrix(1, 0) = 0.0;
	A12matrix(1, 1) = 0.0;
	A12matrix(1, 2) = 0.0;
	A12matrix(1, 3) = 1.0;
	A12matrix(2, 0) = 0.0;
	A12matrix(2, 1) = 0.0;
	A12matrix(2, 2) = -0.8;
	A12matrix(2, 3) = -0.2;
	A12matrix(3, 0) = 0.0;
	A12matrix(3, 1) = 0.0;
	A12matrix(3, 2) = -0.1;
	A12matrix(3, 3) = -0.8;
	system_dynamics12.isEmptyMatrixA = false;
	system_dynamics12.MatrixA = A12matrix;

	system_dynamics12.isEmptyMatrixB = true;

	C12.resize(row);
	C12[0] = 0.0;
	C12[1] = 0.0;
	C12[2] = -0.2;
	C12[3] = -0.8;
	system_dynamics12.isEmptyC = false;
	system_dynamics12.C = C12;

	row = 4;
	col = 4;
	invariantConstraintsMatrix12.resize(row, col);
	invariantConstraintsMatrix12(0, 0) = -1.0;
	invariantConstraintsMatrix12(0, 1) = 0.0;
	invariantConstraintsMatrix12(0, 2) = 0.0;
	invariantConstraintsMatrix12(0, 3) = 0.0;
	invariantConstraintsMatrix12(1, 0) = 1.0;
	invariantConstraintsMatrix12(1, 1) = 0.0;
	invariantConstraintsMatrix12(1, 2) = 0.0;
	invariantConstraintsMatrix12(1, 3) = 0.0;
	invariantConstraintsMatrix12(2, 0) = 0.0;
	invariantConstraintsMatrix12(2, 1) = -1.0;
	invariantConstraintsMatrix12(2, 2) = 0.0;
	invariantConstraintsMatrix12(2, 3) = 0.0;
	invariantConstraintsMatrix12(3, 0) = 0.0;
	invariantConstraintsMatrix12(3, 1) = 1.0;
	invariantConstraintsMatrix12(3, 2) = 0.0;
	invariantConstraintsMatrix12(3, 3) = 0.0;

	invariantBoundValue12.resize(row);
	invariantBoundValue12[0] = -1.0;
	invariantBoundValue12[1] = 2.0;
	invariantBoundValue12[2] = -3.0;
	invariantBoundValue12[3] = 4.0;
	invariantBoundSign = 1;
	invariant12 = polytope::ptr(
			new polytope(invariantConstraintsMatrix12, invariantBoundValue12,
					invariantBoundSign));

	system_dynamics12.U = polytope::ptr(new polytope(true));

	// The mode name is  loc15

	row = 4;
	col = 4;
	A13matrix.resize(row, col);
	A13matrix(0, 0) = 0.0;
	A13matrix(0, 1) = 0.0;
	A13matrix(0, 2) = 1.0;
	A13matrix(0, 3) = 0.0;
	A13matrix(1, 0) = 0.0;
	A13matrix(1, 1) = 0.0;
	A13matrix(1, 2) = 0.0;
	A13matrix(1, 3) = 1.0;
	A13matrix(2, 0) = 0.0;
	A13matrix(2, 1) = 0.0;
	A13matrix(2, 2) = -0.8;
	A13matrix(2, 3) = -0.2;
	A13matrix(3, 0) = 0.0;
	A13matrix(3, 1) = 0.0;
	A13matrix(3, 2) = -0.1;
	A13matrix(3, 3) = -0.8;
	system_dynamics13.isEmptyMatrixA = false;
	system_dynamics13.MatrixA = A13matrix;

	system_dynamics13.isEmptyMatrixB = true;

	C13.resize(row);
	C13[0] = 0.0;
	C13[1] = 0.0;
	C13[2] = -0.8;
	C13[3] = -0.1;
	system_dynamics13.isEmptyC = false;
	system_dynamics13.C = C13;

	row = 4;
	col = 4;
	invariantConstraintsMatrix13.resize(row, col);
	invariantConstraintsMatrix13(0, 0) = -1.0;
	invariantConstraintsMatrix13(0, 1) = 0.0;
	invariantConstraintsMatrix13(0, 2) = 0.0;
	invariantConstraintsMatrix13(0, 3) = 0.0;
	invariantConstraintsMatrix13(1, 0) = 1.0;
	invariantConstraintsMatrix13(1, 1) = 0.0;
	invariantConstraintsMatrix13(1, 2) = 0.0;
	invariantConstraintsMatrix13(1, 3) = 0.0;
	invariantConstraintsMatrix13(2, 0) = 0.0;
	invariantConstraintsMatrix13(2, 1) = -1.0;
	invariantConstraintsMatrix13(2, 2) = 0.0;
	invariantConstraintsMatrix13(2, 3) = 0.0;
	invariantConstraintsMatrix13(3, 0) = 0.0;
	invariantConstraintsMatrix13(3, 1) = 1.0;
	invariantConstraintsMatrix13(3, 2) = 0.0;
	invariantConstraintsMatrix13(3, 3) = 0.0;

	invariantBoundValue13.resize(row);
	invariantBoundValue13[0] = -2.0;
	invariantBoundValue13[1] = 3.0;
	invariantBoundValue13[2] = -4.0;
	invariantBoundValue13[3] = 5.0;
	invariantBoundSign = 1;
	invariant13 = polytope::ptr(
			new polytope(invariantConstraintsMatrix13, invariantBoundValue13,
					invariantBoundSign));

	system_dynamics13.U = polytope::ptr(new polytope(true));

	// The mode name is  loc14

	row = 4;
	col = 4;
	A14matrix.resize(row, col);
	A14matrix(0, 0) = 0.0;
	A14matrix(0, 1) = 0.0;
	A14matrix(0, 2) = 1.0;
	A14matrix(0, 3) = 0.0;
	A14matrix(1, 0) = 0.0;
	A14matrix(1, 1) = 0.0;
	A14matrix(1, 2) = 0.0;
	A14matrix(1, 3) = 1.0;
	A14matrix(2, 0) = 0.0;
	A14matrix(2, 1) = 0.0;
	A14matrix(2, 2) = -0.8;
	A14matrix(2, 3) = -0.2;
	A14matrix(3, 0) = 0.0;
	A14matrix(3, 1) = 0.0;
	A14matrix(3, 2) = -0.1;
	A14matrix(3, 3) = -0.8;
	system_dynamics14.isEmptyMatrixA = false;
	system_dynamics14.MatrixA = A14matrix;

	system_dynamics14.isEmptyMatrixB = true;

	C14.resize(row);
	C14[0] = 0.0;
	C14[1] = 0.0;
	C14[2] = -0.424266;
	C14[3] = 0.494977;
	system_dynamics14.isEmptyC = false;
	system_dynamics14.C = C14;

	row = 4;
	col = 4;
	invariantConstraintsMatrix14.resize(row, col);
	invariantConstraintsMatrix14(0, 0) = -1.0;
	invariantConstraintsMatrix14(0, 1) = 0.0;
	invariantConstraintsMatrix14(0, 2) = 0.0;
	invariantConstraintsMatrix14(0, 3) = 0.0;
	invariantConstraintsMatrix14(1, 0) = 1.0;
	invariantConstraintsMatrix14(1, 1) = 0.0;
	invariantConstraintsMatrix14(1, 2) = 0.0;
	invariantConstraintsMatrix14(1, 3) = 0.0;
	invariantConstraintsMatrix14(2, 0) = 0.0;
	invariantConstraintsMatrix14(2, 1) = -1.0;
	invariantConstraintsMatrix14(2, 2) = 0.0;
	invariantConstraintsMatrix14(2, 3) = 0.0;
	invariantConstraintsMatrix14(3, 0) = 0.0;
	invariantConstraintsMatrix14(3, 1) = 1.0;
	invariantConstraintsMatrix14(3, 2) = 0.0;
	invariantConstraintsMatrix14(3, 3) = 0.0;

	invariantBoundValue14.resize(row);
	invariantBoundValue14[0] = -2.0;
	invariantBoundValue14[1] = 3.0;
	invariantBoundValue14[2] = -3.0;
	invariantBoundValue14[3] = 4.0;
	invariantBoundSign = 1;
	invariant14 = polytope::ptr(
			new polytope(invariantConstraintsMatrix14, invariantBoundValue14,
					invariantBoundSign));

	system_dynamics14.U = polytope::ptr(new polytope(true));

	// The mode name is  loc20

	row = 4;
	col = 4;
	A15matrix.resize(row, col);
	A15matrix(0, 0) = 0.0;
	A15matrix(0, 1) = 0.0;
	A15matrix(0, 2) = 1.0;
	A15matrix(0, 3) = 0.0;
	A15matrix(1, 0) = 0.0;
	A15matrix(1, 1) = 0.0;
	A15matrix(1, 2) = 0.0;
	A15matrix(1, 3) = 1.0;
	A15matrix(2, 0) = 0.0;
	A15matrix(2, 1) = 0.0;
	A15matrix(2, 2) = -0.8;
	A15matrix(2, 3) = -0.2;
	A15matrix(3, 0) = 0.0;
	A15matrix(3, 1) = 0.0;
	A15matrix(3, 2) = -0.1;
	A15matrix(3, 3) = -0.8;
	system_dynamics15.isEmptyMatrixA = false;
	system_dynamics15.MatrixA = A15matrix;

	system_dynamics15.isEmptyMatrixB = true;

	C15.resize(row);
	C15[0] = 0.0;
	C15[1] = 0.0;
	C15[2] = -0.8;
	C15[3] = -0.1;
	system_dynamics15.isEmptyC = false;
	system_dynamics15.C = C15;

	row = 4;
	col = 4;
	invariantConstraintsMatrix15.resize(row, col);
	invariantConstraintsMatrix15(0, 0) = -1.0;
	invariantConstraintsMatrix15(0, 1) = 0.0;
	invariantConstraintsMatrix15(0, 2) = 0.0;
	invariantConstraintsMatrix15(0, 3) = 0.0;
	invariantConstraintsMatrix15(1, 0) = 1.0;
	invariantConstraintsMatrix15(1, 1) = 0.0;
	invariantConstraintsMatrix15(1, 2) = 0.0;
	invariantConstraintsMatrix15(1, 3) = 0.0;
	invariantConstraintsMatrix15(2, 0) = 0.0;
	invariantConstraintsMatrix15(2, 1) = -1.0;
	invariantConstraintsMatrix15(2, 2) = 0.0;
	invariantConstraintsMatrix15(2, 3) = 0.0;
	invariantConstraintsMatrix15(3, 0) = 0.0;
	invariantConstraintsMatrix15(3, 1) = 1.0;
	invariantConstraintsMatrix15(3, 2) = 0.0;
	invariantConstraintsMatrix15(3, 3) = 0.0;

	invariantBoundValue15.resize(row);
	invariantBoundValue15[0] = -3.0;
	invariantBoundValue15[1] = 4.0;
	invariantBoundValue15[2] = -4.0;
	invariantBoundValue15[3] = 5.0;
	invariantBoundSign = 1;
	invariant15 = polytope::ptr(
			new polytope(invariantConstraintsMatrix15, invariantBoundValue15,
					invariantBoundSign));

	system_dynamics15.U = polytope::ptr(new polytope(true));

	// The mode name is  loc19

	row = 4;
	col = 4;
	A16matrix.resize(row, col);
	A16matrix(0, 0) = 0.0;
	A16matrix(0, 1) = 0.0;
	A16matrix(0, 2) = 1.0;
	A16matrix(0, 3) = 0.0;
	A16matrix(1, 0) = 0.0;
	A16matrix(1, 1) = 0.0;
	A16matrix(1, 2) = 0.0;
	A16matrix(1, 3) = 1.0;
	A16matrix(2, 0) = 0.0;
	A16matrix(2, 1) = 0.0;
	A16matrix(2, 2) = -0.8;
	A16matrix(2, 3) = -0.2;
	A16matrix(3, 0) = 0.0;
	A16matrix(3, 1) = 0.0;
	A16matrix(3, 2) = -0.1;
	A16matrix(3, 3) = -0.8;
	system_dynamics16.isEmptyMatrixA = false;
	system_dynamics16.MatrixA = A16matrix;

	system_dynamics16.isEmptyMatrixB = true;

	C16.resize(row);
	C16[0] = 0.0;
	C16[1] = 0.0;
	C16[2] = -0.424266;
	C16[3] = 0.494977;
	system_dynamics16.isEmptyC = false;
	system_dynamics16.C = C16;

	row = 4;
	col = 4;
	invariantConstraintsMatrix16.resize(row, col);
	invariantConstraintsMatrix16(0, 0) = -1.0;
	invariantConstraintsMatrix16(0, 1) = 0.0;
	invariantConstraintsMatrix16(0, 2) = 0.0;
	invariantConstraintsMatrix16(0, 3) = 0.0;
	invariantConstraintsMatrix16(1, 0) = 1.0;
	invariantConstraintsMatrix16(1, 1) = 0.0;
	invariantConstraintsMatrix16(1, 2) = 0.0;
	invariantConstraintsMatrix16(1, 3) = 0.0;
	invariantConstraintsMatrix16(2, 0) = 0.0;
	invariantConstraintsMatrix16(2, 1) = -1.0;
	invariantConstraintsMatrix16(2, 2) = 0.0;
	invariantConstraintsMatrix16(2, 3) = 0.0;
	invariantConstraintsMatrix16(3, 0) = 0.0;
	invariantConstraintsMatrix16(3, 1) = 1.0;
	invariantConstraintsMatrix16(3, 2) = 0.0;
	invariantConstraintsMatrix16(3, 3) = 0.0;

	invariantBoundValue16.resize(row);
	invariantBoundValue16[0] = -3.0;
	invariantBoundValue16[1] = 4.0;
	invariantBoundValue16[2] = -3.0;
	invariantBoundValue16[3] = 4.0;
	invariantBoundSign = 1;
	invariant16 = polytope::ptr(
			new polytope(invariantConstraintsMatrix16, invariantBoundValue16,
					invariantBoundSign));

	system_dynamics16.U = polytope::ptr(new polytope(true));

	// The mode name is  loc17

	row = 4;
	col = 4;
	A17matrix.resize(row, col);
	A17matrix(0, 0) = 0.0;
	A17matrix(0, 1) = 0.0;
	A17matrix(0, 2) = 1.0;
	A17matrix(0, 3) = 0.0;
	A17matrix(1, 0) = 0.0;
	A17matrix(1, 1) = 0.0;
	A17matrix(1, 2) = 0.0;
	A17matrix(1, 3) = 1.0;
	A17matrix(2, 0) = 0.0;
	A17matrix(2, 1) = 0.0;
	A17matrix(2, 2) = -0.8;
	A17matrix(2, 3) = -0.2;
	A17matrix(3, 0) = 0.0;
	A17matrix(3, 1) = 0.0;
	A17matrix(3, 2) = -0.1;
	A17matrix(3, 3) = -0.8;
	system_dynamics17.isEmptyMatrixA = false;
	system_dynamics17.MatrixA = A17matrix;

	system_dynamics17.isEmptyMatrixB = true;

	C17.resize(row);
	C17[0] = 0.0;
	C17[1] = 0.0;
	C17[2] = -0.8;
	C17[3] = -0.1;
	system_dynamics17.isEmptyC = false;
	system_dynamics17.C = C17;

	row = 4;
	col = 4;
	invariantConstraintsMatrix17.resize(row, col);
	invariantConstraintsMatrix17(0, 0) = -1.0;
	invariantConstraintsMatrix17(0, 1) = 0.0;
	invariantConstraintsMatrix17(0, 2) = 0.0;
	invariantConstraintsMatrix17(0, 3) = 0.0;
	invariantConstraintsMatrix17(1, 0) = 1.0;
	invariantConstraintsMatrix17(1, 1) = 0.0;
	invariantConstraintsMatrix17(1, 2) = 0.0;
	invariantConstraintsMatrix17(1, 3) = 0.0;
	invariantConstraintsMatrix17(2, 0) = 0.0;
	invariantConstraintsMatrix17(2, 1) = -1.0;
	invariantConstraintsMatrix17(2, 2) = 0.0;
	invariantConstraintsMatrix17(2, 3) = 0.0;
	invariantConstraintsMatrix17(3, 0) = 0.0;
	invariantConstraintsMatrix17(3, 1) = 1.0;
	invariantConstraintsMatrix17(3, 2) = 0.0;
	invariantConstraintsMatrix17(3, 3) = 0.0;

	invariantBoundValue17.resize(row);
	invariantBoundValue17[0] = -3.0;
	invariantBoundValue17[1] = 4.0;
	invariantBoundValue17[2] = -1.0;
	invariantBoundValue17[3] = 2.0;
	invariantBoundSign = 1;
	invariant17 = polytope::ptr(
			new polytope(invariantConstraintsMatrix17, invariantBoundValue17,
					invariantBoundSign));

	system_dynamics17.U = polytope::ptr(new polytope(true));

	// The mode name is  loc18

	row = 4;
	col = 4;
	A18matrix.resize(row, col);
	A18matrix(0, 0) = 0.0;
	A18matrix(0, 1) = 0.0;
	A18matrix(0, 2) = 1.0;
	A18matrix(0, 3) = 0.0;
	A18matrix(1, 0) = 0.0;
	A18matrix(1, 1) = 0.0;
	A18matrix(1, 2) = 0.0;
	A18matrix(1, 3) = 1.0;
	A18matrix(2, 0) = 0.0;
	A18matrix(2, 1) = 0.0;
	A18matrix(2, 2) = -0.8;
	A18matrix(2, 3) = -0.2;
	A18matrix(3, 0) = 0.0;
	A18matrix(3, 1) = 0.0;
	A18matrix(3, 2) = -0.1;
	A18matrix(3, 3) = -0.8;
	system_dynamics18.isEmptyMatrixA = false;
	system_dynamics18.MatrixA = A18matrix;

	system_dynamics18.isEmptyMatrixB = true;

	C18.resize(row);
	C18[0] = 0.0;
	C18[1] = 0.0;
	C18[2] = 0.424266;
	C18[3] = -0.494977;
	system_dynamics18.isEmptyC = false;
	system_dynamics18.C = C18;

	row = 4;
	col = 4;
	invariantConstraintsMatrix18.resize(row, col);
	invariantConstraintsMatrix18(0, 0) = -1.0;
	invariantConstraintsMatrix18(0, 1) = 0.0;
	invariantConstraintsMatrix18(0, 2) = 0.0;
	invariantConstraintsMatrix18(0, 3) = 0.0;
	invariantConstraintsMatrix18(1, 0) = 1.0;
	invariantConstraintsMatrix18(1, 1) = 0.0;
	invariantConstraintsMatrix18(1, 2) = 0.0;
	invariantConstraintsMatrix18(1, 3) = 0.0;
	invariantConstraintsMatrix18(2, 0) = 0.0;
	invariantConstraintsMatrix18(2, 1) = -1.0;
	invariantConstraintsMatrix18(2, 2) = 0.0;
	invariantConstraintsMatrix18(2, 3) = 0.0;
	invariantConstraintsMatrix18(3, 0) = 0.0;
	invariantConstraintsMatrix18(3, 1) = 1.0;
	invariantConstraintsMatrix18(3, 2) = 0.0;
	invariantConstraintsMatrix18(3, 3) = 0.0;

	invariantBoundValue18.resize(row);
	invariantBoundValue18[0] = -3.0;
	invariantBoundValue18[1] = 4.0;
	invariantBoundValue18[2] = -2.0;
	invariantBoundValue18[3] = 3.0;
	invariantBoundSign = 1;
	invariant18 = polytope::ptr(
			new polytope(invariantConstraintsMatrix18, invariantBoundValue18,
					invariantBoundSign));

	system_dynamics18.U = polytope::ptr(new polytope(true));

	// The mode name is  loc16

	row = 4;
	col = 4;
	A19matrix.resize(row, col);
	A19matrix(0, 0) = 0.0;
	A19matrix(0, 1) = 0.0;
	A19matrix(0, 2) = 1.0;
	A19matrix(0, 3) = 0.0;
	A19matrix(1, 0) = 0.0;
	A19matrix(1, 1) = 0.0;
	A19matrix(1, 2) = 0.0;
	A19matrix(1, 3) = 1.0;
	A19matrix(2, 0) = 0.0;
	A19matrix(2, 1) = 0.0;
	A19matrix(2, 2) = -0.8;
	A19matrix(2, 3) = -0.2;
	A19matrix(3, 0) = 0.0;
	A19matrix(3, 1) = 0.0;
	A19matrix(3, 2) = -0.1;
	A19matrix(3, 3) = -0.8;
	system_dynamics19.isEmptyMatrixA = false;
	system_dynamics19.MatrixA = A19matrix;

	system_dynamics19.isEmptyMatrixB = true;

	C19.resize(row);
	C19[0] = 0.0;
	C19[1] = 0.0;
	C19[2] = 0.2;
	C19[3] = 0.8;
	system_dynamics19.isEmptyC = false;
	system_dynamics19.C = C19;

	row = 4;
	col = 4;
	invariantConstraintsMatrix19.resize(row, col);
	invariantConstraintsMatrix19(0, 0) = -1.0;
	invariantConstraintsMatrix19(0, 1) = 0.0;
	invariantConstraintsMatrix19(0, 2) = 0.0;
	invariantConstraintsMatrix19(0, 3) = 0.0;
	invariantConstraintsMatrix19(1, 0) = 1.0;
	invariantConstraintsMatrix19(1, 1) = 0.0;
	invariantConstraintsMatrix19(1, 2) = 0.0;
	invariantConstraintsMatrix19(1, 3) = 0.0;
	invariantConstraintsMatrix19(2, 0) = 0.0;
	invariantConstraintsMatrix19(2, 1) = -1.0;
	invariantConstraintsMatrix19(2, 2) = 0.0;
	invariantConstraintsMatrix19(2, 3) = 0.0;
	invariantConstraintsMatrix19(3, 0) = 0.0;
	invariantConstraintsMatrix19(3, 1) = 1.0;
	invariantConstraintsMatrix19(3, 2) = 0.0;
	invariantConstraintsMatrix19(3, 3) = 0.0;

	invariantBoundValue19.resize(row);
	invariantBoundValue19[0] = -3.0;
	invariantBoundValue19[1] = 4.0;
	invariantBoundValue19[2] = -0.0;
	invariantBoundValue19[3] = 1.0;
	invariantBoundSign = 1;
	invariant19 = polytope::ptr(
			new polytope(invariantConstraintsMatrix19, invariantBoundValue19,
					invariantBoundSign));

	system_dynamics19.U = polytope::ptr(new polytope(true));

	// The mode name is  loc25

	row = 4;
	col = 4;
	A20matrix.resize(row, col);
	A20matrix(0, 0) = 0.0;
	A20matrix(0, 1) = 0.0;
	A20matrix(0, 2) = 1.0;
	A20matrix(0, 3) = 0.0;
	A20matrix(1, 0) = 0.0;
	A20matrix(1, 1) = 0.0;
	A20matrix(1, 2) = 0.0;
	A20matrix(1, 3) = 1.0;
	A20matrix(2, 0) = 0.0;
	A20matrix(2, 1) = 0.0;
	A20matrix(2, 2) = -0.8;
	A20matrix(2, 3) = -0.2;
	A20matrix(3, 0) = 0.0;
	A20matrix(3, 1) = 0.0;
	A20matrix(3, 2) = -0.1;
	A20matrix(3, 3) = -0.8;
	system_dynamics20.isEmptyMatrixA = false;
	system_dynamics20.MatrixA = A20matrix;

	system_dynamics20.isEmptyMatrixB = true;

	C20.resize(row);
	C20[0] = 0.0;
	C20[1] = 0.0;
	C20[2] = -0.8;
	C20[3] = -0.1;
	system_dynamics20.isEmptyC = false;
	system_dynamics20.C = C20;

	row = 4;
	col = 4;
	invariantConstraintsMatrix20.resize(row, col);
	invariantConstraintsMatrix20(0, 0) = -1.0;
	invariantConstraintsMatrix20(0, 1) = 0.0;
	invariantConstraintsMatrix20(0, 2) = 0.0;
	invariantConstraintsMatrix20(0, 3) = 0.0;
	invariantConstraintsMatrix20(1, 0) = 1.0;
	invariantConstraintsMatrix20(1, 1) = 0.0;
	invariantConstraintsMatrix20(1, 2) = 0.0;
	invariantConstraintsMatrix20(1, 3) = 0.0;
	invariantConstraintsMatrix20(2, 0) = 0.0;
	invariantConstraintsMatrix20(2, 1) = -1.0;
	invariantConstraintsMatrix20(2, 2) = 0.0;
	invariantConstraintsMatrix20(2, 3) = 0.0;
	invariantConstraintsMatrix20(3, 0) = 0.0;
	invariantConstraintsMatrix20(3, 1) = 1.0;
	invariantConstraintsMatrix20(3, 2) = 0.0;
	invariantConstraintsMatrix20(3, 3) = 0.0;

	invariantBoundValue20.resize(row);
	invariantBoundValue20[0] = -4.0;
	invariantBoundValue20[1] = 5.0;
	invariantBoundValue20[2] = -4.0;
	invariantBoundValue20[3] = 5.0;
	invariantBoundSign = 1;
	invariant20 = polytope::ptr(
			new polytope(invariantConstraintsMatrix20, invariantBoundValue20,
					invariantBoundSign));

	system_dynamics20.U = polytope::ptr(new polytope(true));

	// The mode name is  loc24

	row = 4;
	col = 4;
	A21matrix.resize(row, col);
	A21matrix(0, 0) = 0.0;
	A21matrix(0, 1) = 0.0;
	A21matrix(0, 2) = 1.0;
	A21matrix(0, 3) = 0.0;
	A21matrix(1, 0) = 0.0;
	A21matrix(1, 1) = 0.0;
	A21matrix(1, 2) = 0.0;
	A21matrix(1, 3) = 1.0;
	A21matrix(2, 0) = 0.0;
	A21matrix(2, 1) = 0.0;
	A21matrix(2, 2) = -0.8;
	A21matrix(2, 3) = -0.2;
	A21matrix(3, 0) = 0.0;
	A21matrix(3, 1) = 0.0;
	A21matrix(3, 2) = -0.1;
	A21matrix(3, 3) = -0.8;
	system_dynamics21.isEmptyMatrixA = false;
	system_dynamics21.MatrixA = A21matrix;

	system_dynamics21.isEmptyMatrixB = true;

	C21.resize(row);
	C21[0] = 0.0;
	C21[1] = 0.0;
	C21[2] = -0.2;
	C21[3] = -0.8;
	system_dynamics21.isEmptyC = false;
	system_dynamics21.C = C21;

	row = 4;
	col = 4;
	invariantConstraintsMatrix21.resize(row, col);
	invariantConstraintsMatrix21(0, 0) = -1.0;
	invariantConstraintsMatrix21(0, 1) = 0.0;
	invariantConstraintsMatrix21(0, 2) = 0.0;
	invariantConstraintsMatrix21(0, 3) = 0.0;
	invariantConstraintsMatrix21(1, 0) = 1.0;
	invariantConstraintsMatrix21(1, 1) = 0.0;
	invariantConstraintsMatrix21(1, 2) = 0.0;
	invariantConstraintsMatrix21(1, 3) = 0.0;
	invariantConstraintsMatrix21(2, 0) = 0.0;
	invariantConstraintsMatrix21(2, 1) = -1.0;
	invariantConstraintsMatrix21(2, 2) = 0.0;
	invariantConstraintsMatrix21(2, 3) = 0.0;
	invariantConstraintsMatrix21(3, 0) = 0.0;
	invariantConstraintsMatrix21(3, 1) = 1.0;
	invariantConstraintsMatrix21(3, 2) = 0.0;
	invariantConstraintsMatrix21(3, 3) = 0.0;

	invariantBoundValue21.resize(row);
	invariantBoundValue21[0] = -4.0;
	invariantBoundValue21[1] = 5.0;
	invariantBoundValue21[2] = -3.0;
	invariantBoundValue21[3] = 4.0;
	invariantBoundSign = 1;
	invariant21 = polytope::ptr(
			new polytope(invariantConstraintsMatrix21, invariantBoundValue21,
					invariantBoundSign));

	system_dynamics21.U = polytope::ptr(new polytope(true));

	// The mode name is  loc23

	row = 4;
	col = 4;
	A22matrix.resize(row, col);
	A22matrix(0, 0) = 0.0;
	A22matrix(0, 1) = 0.0;
	A22matrix(0, 2) = 1.0;
	A22matrix(0, 3) = 0.0;
	A22matrix(1, 0) = 0.0;
	A22matrix(1, 1) = 0.0;
	A22matrix(1, 2) = 0.0;
	A22matrix(1, 3) = 1.0;
	A22matrix(2, 0) = 0.0;
	A22matrix(2, 1) = 0.0;
	A22matrix(2, 2) = -0.8;
	A22matrix(2, 3) = -0.2;
	A22matrix(3, 0) = 0.0;
	A22matrix(3, 1) = 0.0;
	A22matrix(3, 2) = -0.1;
	A22matrix(3, 3) = -0.8;
	system_dynamics22.isEmptyMatrixA = false;
	system_dynamics22.MatrixA = A22matrix;

	system_dynamics22.isEmptyMatrixB = true;

	C22.resize(row);
	C22[0] = 0.0;
	C22[1] = 0.0;
	C22[2] = -0.2;
	C22[3] = -0.8;
	system_dynamics22.isEmptyC = false;
	system_dynamics22.C = C22;

	row = 4;
	col = 4;
	invariantConstraintsMatrix22.resize(row, col);
	invariantConstraintsMatrix22(0, 0) = -1.0;
	invariantConstraintsMatrix22(0, 1) = 0.0;
	invariantConstraintsMatrix22(0, 2) = 0.0;
	invariantConstraintsMatrix22(0, 3) = 0.0;
	invariantConstraintsMatrix22(1, 0) = 1.0;
	invariantConstraintsMatrix22(1, 1) = 0.0;
	invariantConstraintsMatrix22(1, 2) = 0.0;
	invariantConstraintsMatrix22(1, 3) = 0.0;
	invariantConstraintsMatrix22(2, 0) = 0.0;
	invariantConstraintsMatrix22(2, 1) = -1.0;
	invariantConstraintsMatrix22(2, 2) = 0.0;
	invariantConstraintsMatrix22(2, 3) = 0.0;
	invariantConstraintsMatrix22(3, 0) = 0.0;
	invariantConstraintsMatrix22(3, 1) = 1.0;
	invariantConstraintsMatrix22(3, 2) = 0.0;
	invariantConstraintsMatrix22(3, 3) = 0.0;

	invariantBoundValue22.resize(row);
	invariantBoundValue22[0] = -4.0;
	invariantBoundValue22[1] = 5.0;
	invariantBoundValue22[2] = -2.0;
	invariantBoundValue22[3] = 3.0;
	invariantBoundSign = 1;
	invariant22 = polytope::ptr(
			new polytope(invariantConstraintsMatrix22, invariantBoundValue22,
					invariantBoundSign));

	system_dynamics22.U = polytope::ptr(new polytope(true));

	// The mode name is  loc22

	row = 4;
	col = 4;
	A23matrix.resize(row, col);
	A23matrix(0, 0) = 0.0;
	A23matrix(0, 1) = 0.0;
	A23matrix(0, 2) = 1.0;
	A23matrix(0, 3) = 0.0;
	A23matrix(1, 0) = 0.0;
	A23matrix(1, 1) = 0.0;
	A23matrix(1, 2) = 0.0;
	A23matrix(1, 3) = 1.0;
	A23matrix(2, 0) = 0.0;
	A23matrix(2, 1) = 0.0;
	A23matrix(2, 2) = -0.8;
	A23matrix(2, 3) = -0.2;
	A23matrix(3, 0) = 0.0;
	A23matrix(3, 1) = 0.0;
	A23matrix(3, 2) = -0.1;
	A23matrix(3, 3) = -0.8;
	system_dynamics23.isEmptyMatrixA = false;
	system_dynamics23.MatrixA = A23matrix;

	system_dynamics23.isEmptyMatrixB = true;

	C23.resize(row);
	C23[0] = 0.0;
	C23[1] = 0.0;
	C23[2] = -0.8;
	C23[3] = -0.1;
	system_dynamics23.isEmptyC = false;
	system_dynamics23.C = C23;

	row = 4;
	col = 4;
	invariantConstraintsMatrix23.resize(row, col);
	invariantConstraintsMatrix23(0, 0) = -1.0;
	invariantConstraintsMatrix23(0, 1) = 0.0;
	invariantConstraintsMatrix23(0, 2) = 0.0;
	invariantConstraintsMatrix23(0, 3) = 0.0;
	invariantConstraintsMatrix23(1, 0) = 1.0;
	invariantConstraintsMatrix23(1, 1) = 0.0;
	invariantConstraintsMatrix23(1, 2) = 0.0;
	invariantConstraintsMatrix23(1, 3) = 0.0;
	invariantConstraintsMatrix23(2, 0) = 0.0;
	invariantConstraintsMatrix23(2, 1) = -1.0;
	invariantConstraintsMatrix23(2, 2) = 0.0;
	invariantConstraintsMatrix23(2, 3) = 0.0;
	invariantConstraintsMatrix23(3, 0) = 0.0;
	invariantConstraintsMatrix23(3, 1) = 1.0;
	invariantConstraintsMatrix23(3, 2) = 0.0;
	invariantConstraintsMatrix23(3, 3) = 0.0;

	invariantBoundValue23.resize(row);
	invariantBoundValue23[0] = -4.0;
	invariantBoundValue23[1] = 5.0;
	invariantBoundValue23[2] = -1.0;
	invariantBoundValue23[3] = 2.0;
	invariantBoundSign = 1;
	invariant23 = polytope::ptr(
			new polytope(invariantConstraintsMatrix23, invariantBoundValue23,
					invariantBoundSign));

	system_dynamics23.U = polytope::ptr(new polytope(true));

	// The mode name is  loc21

	row = 4;
	col = 4;
	A24matrix.resize(row, col);
	A24matrix(0, 0) = 0.0;
	A24matrix(0, 1) = 0.0;
	A24matrix(0, 2) = 1.0;
	A24matrix(0, 3) = 0.0;
	A24matrix(1, 0) = 0.0;
	A24matrix(1, 1) = 0.0;
	A24matrix(1, 2) = 0.0;
	A24matrix(1, 3) = 1.0;
	A24matrix(2, 0) = 0.0;
	A24matrix(2, 1) = 0.0;
	A24matrix(2, 2) = -0.8;
	A24matrix(2, 3) = -0.2;
	A24matrix(3, 0) = 0.0;
	A24matrix(3, 1) = 0.0;
	A24matrix(3, 2) = -0.1;
	A24matrix(3, 3) = -0.8;
	system_dynamics24.isEmptyMatrixA = false;
	system_dynamics24.MatrixA = A24matrix;

	system_dynamics24.isEmptyMatrixB = true;

	C24.resize(row);
	C24[0] = 0.0;
	C24[1] = 0.0;
	C24[2] = 0.2;
	C24[3] = 0.8;
	system_dynamics24.isEmptyC = false;
	system_dynamics24.C = C24;

	row = 4;
	col = 4;
	invariantConstraintsMatrix24.resize(row, col);
	invariantConstraintsMatrix24(0, 0) = -1.0;
	invariantConstraintsMatrix24(0, 1) = 0.0;
	invariantConstraintsMatrix24(0, 2) = 0.0;
	invariantConstraintsMatrix24(0, 3) = 0.0;
	invariantConstraintsMatrix24(1, 0) = 1.0;
	invariantConstraintsMatrix24(1, 1) = 0.0;
	invariantConstraintsMatrix24(1, 2) = 0.0;
	invariantConstraintsMatrix24(1, 3) = 0.0;
	invariantConstraintsMatrix24(2, 0) = 0.0;
	invariantConstraintsMatrix24(2, 1) = -1.0;
	invariantConstraintsMatrix24(2, 2) = 0.0;
	invariantConstraintsMatrix24(2, 3) = 0.0;
	invariantConstraintsMatrix24(3, 0) = 0.0;
	invariantConstraintsMatrix24(3, 1) = 1.0;
	invariantConstraintsMatrix24(3, 2) = 0.0;
	invariantConstraintsMatrix24(3, 3) = 0.0;

	invariantBoundValue24.resize(row);
	invariantBoundValue24[0] = -4.0;
	invariantBoundValue24[1] = 5.0;
	invariantBoundValue24[2] = -0.0;
	invariantBoundValue24[3] = 1.0;
	invariantBoundSign = 1;
	invariant24 = polytope::ptr(
			new polytope(invariantConstraintsMatrix24, invariantBoundValue24,
					invariantBoundSign));

	system_dynamics24.U = polytope::ptr(new polytope(true));

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

	boundValueI.resize(row);
	boundValueI[0] = 3.4;
	boundValueI[1] = -3.1;
	boundValueI[2] = 3.8;
	boundValueI[3] = -3.6;
	boundValueI[4] = 0.1;
	boundValueI[5] = -0.1;
	boundValueI[6] = 0.1;
	boundValueI[7] = -0.1;

	/*	boundValueI[0] = 3.5;
	 boundValueI[1] = -3.5;
	 boundValueI[2] = 3.5;
	 boundValueI[3] = -3.5;
	 boundValueI[4] = 0.1;
	 boundValueI[5] = -0.1;
	 boundValueI[6] = 0.1;
	 boundValueI[7] = -0.1;*/

	boundSignI = 1;

	// The transition label ist6

	// Original guard: 0 <= x1 & x1 <= 1 & x2 = 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix0.resize(row, col);
	gaurdConstraintsMatrix0(0, 0) = -1.0;
	gaurdConstraintsMatrix0(0, 1) = 0.0;
	gaurdConstraintsMatrix0(0, 2) = 0.0;
	gaurdConstraintsMatrix0(0, 3) = 0.0;
	gaurdConstraintsMatrix0(1, 0) = 1.0;
	gaurdConstraintsMatrix0(1, 1) = 0.0;
	gaurdConstraintsMatrix0(1, 2) = 0.0;
	gaurdConstraintsMatrix0(1, 3) = 0.0;
	gaurdConstraintsMatrix0(2, 0) = 0.0;
	gaurdConstraintsMatrix0(2, 1) = -1.0;
	gaurdConstraintsMatrix0(2, 2) = 0.0;
	gaurdConstraintsMatrix0(2, 3) = 0.0;
	gaurdConstraintsMatrix0(3, 0) = 0.0;
	gaurdConstraintsMatrix0(3, 1) = 1.0;
	gaurdConstraintsMatrix0(3, 2) = 0.0;
	gaurdConstraintsMatrix0(3, 3) = 0.0;
	gaurdConstraintsMatrix0(4, 0) = 0.0;
	gaurdConstraintsMatrix0(4, 1) = 0.0;
	gaurdConstraintsMatrix0(4, 2) = -1.0;
	gaurdConstraintsMatrix0(4, 3) = 0.0;
	gaurdConstraintsMatrix0(5, 0) = 0.0;
	gaurdConstraintsMatrix0(5, 1) = 0.0;
	gaurdConstraintsMatrix0(5, 2) = 1.0;
	gaurdConstraintsMatrix0(5, 3) = 0.0;
	gaurdConstraintsMatrix0(6, 0) = 0.0;
	gaurdConstraintsMatrix0(6, 1) = 0.0;
	gaurdConstraintsMatrix0(6, 2) = 0.0;
	gaurdConstraintsMatrix0(6, 3) = -1.0;
	gaurdConstraintsMatrix0(7, 0) = 0.0;
	gaurdConstraintsMatrix0(7, 1) = 0.0;
	gaurdConstraintsMatrix0(7, 2) = 0.0;
	gaurdConstraintsMatrix0(7, 3) = 1.0;

	gaurdBoundValue0.resize(row);
	gaurdBoundValue0[0] = -0.0;
	gaurdBoundValue0[1] = 1.0;
	gaurdBoundValue0[2] = -3.0;
	gaurdBoundValue0[3] = 3.0;
	gaurdBoundValue0[4] = 1000.0;
	gaurdBoundValue0[5] = 1000.0;
	gaurdBoundValue0[6] = 1000.0;
	gaurdBoundValue0[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope0 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix0, gaurdBoundValue0,
					gaurdBoundSign));

	// The transition label ist8

	// Original guard: x1 = 1 & 2 <= x2 & x2 <= 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix1.resize(row, col);
	gaurdConstraintsMatrix1(0, 0) = -1.0;
	gaurdConstraintsMatrix1(0, 1) = 0.0;
	gaurdConstraintsMatrix1(0, 2) = 0.0;
	gaurdConstraintsMatrix1(0, 3) = 0.0;
	gaurdConstraintsMatrix1(1, 0) = 1.0;
	gaurdConstraintsMatrix1(1, 1) = 0.0;
	gaurdConstraintsMatrix1(1, 2) = 0.0;
	gaurdConstraintsMatrix1(1, 3) = 0.0;
	gaurdConstraintsMatrix1(2, 0) = 0.0;
	gaurdConstraintsMatrix1(2, 1) = -1.0;
	gaurdConstraintsMatrix1(2, 2) = 0.0;
	gaurdConstraintsMatrix1(2, 3) = 0.0;
	gaurdConstraintsMatrix1(3, 0) = 0.0;
	gaurdConstraintsMatrix1(3, 1) = 1.0;
	gaurdConstraintsMatrix1(3, 2) = 0.0;
	gaurdConstraintsMatrix1(3, 3) = 0.0;
	gaurdConstraintsMatrix1(4, 0) = 0.0;
	gaurdConstraintsMatrix1(4, 1) = 0.0;
	gaurdConstraintsMatrix1(4, 2) = -1.0;
	gaurdConstraintsMatrix1(4, 3) = 0.0;
	gaurdConstraintsMatrix1(5, 0) = 0.0;
	gaurdConstraintsMatrix1(5, 1) = 0.0;
	gaurdConstraintsMatrix1(5, 2) = 1.0;
	gaurdConstraintsMatrix1(5, 3) = 0.0;
	gaurdConstraintsMatrix1(6, 0) = 0.0;
	gaurdConstraintsMatrix1(6, 1) = 0.0;
	gaurdConstraintsMatrix1(6, 2) = 0.0;
	gaurdConstraintsMatrix1(6, 3) = -1.0;
	gaurdConstraintsMatrix1(7, 0) = 0.0;
	gaurdConstraintsMatrix1(7, 1) = 0.0;
	gaurdConstraintsMatrix1(7, 2) = 0.0;
	gaurdConstraintsMatrix1(7, 3) = 1.0;

	gaurdBoundValue1.resize(row);
	gaurdBoundValue1[0] = -1.0;
	gaurdBoundValue1[1] = 1.0;
	gaurdBoundValue1[2] = -2.0;
	gaurdBoundValue1[3] = 3.0;
	gaurdBoundValue1[4] = 1000.0;
	gaurdBoundValue1[5] = 1000.0;
	gaurdBoundValue1[6] = 1000.0;
	gaurdBoundValue1[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope1 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix1, gaurdBoundValue1,
					gaurdBoundSign));

	// The transition label ist7

	// Original guard: 0 <= x1 & x1 <= 1 & x2 = 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix2.resize(row, col);
	gaurdConstraintsMatrix2(0, 0) = -1.0;
	gaurdConstraintsMatrix2(0, 1) = 0.0;
	gaurdConstraintsMatrix2(0, 2) = 0.0;
	gaurdConstraintsMatrix2(0, 3) = 0.0;
	gaurdConstraintsMatrix2(1, 0) = 1.0;
	gaurdConstraintsMatrix2(1, 1) = 0.0;
	gaurdConstraintsMatrix2(1, 2) = 0.0;
	gaurdConstraintsMatrix2(1, 3) = 0.0;
	gaurdConstraintsMatrix2(2, 0) = 0.0;
	gaurdConstraintsMatrix2(2, 1) = -1.0;
	gaurdConstraintsMatrix2(2, 2) = 0.0;
	gaurdConstraintsMatrix2(2, 3) = 0.0;
	gaurdConstraintsMatrix2(3, 0) = 0.0;
	gaurdConstraintsMatrix2(3, 1) = 1.0;
	gaurdConstraintsMatrix2(3, 2) = 0.0;
	gaurdConstraintsMatrix2(3, 3) = 0.0;
	gaurdConstraintsMatrix2(4, 0) = 0.0;
	gaurdConstraintsMatrix2(4, 1) = 0.0;
	gaurdConstraintsMatrix2(4, 2) = -1.0;
	gaurdConstraintsMatrix2(4, 3) = 0.0;
	gaurdConstraintsMatrix2(5, 0) = 0.0;
	gaurdConstraintsMatrix2(5, 1) = 0.0;
	gaurdConstraintsMatrix2(5, 2) = 1.0;
	gaurdConstraintsMatrix2(5, 3) = 0.0;
	gaurdConstraintsMatrix2(6, 0) = 0.0;
	gaurdConstraintsMatrix2(6, 1) = 0.0;
	gaurdConstraintsMatrix2(6, 2) = 0.0;
	gaurdConstraintsMatrix2(6, 3) = -1.0;
	gaurdConstraintsMatrix2(7, 0) = 0.0;
	gaurdConstraintsMatrix2(7, 1) = 0.0;
	gaurdConstraintsMatrix2(7, 2) = 0.0;
	gaurdConstraintsMatrix2(7, 3) = 1.0;

	gaurdBoundValue2.resize(row);
	gaurdBoundValue2[0] = -0.0;
	gaurdBoundValue2[1] = 1.0;
	gaurdBoundValue2[2] = -2.0;
	gaurdBoundValue2[3] = 2.0;
	gaurdBoundValue2[4] = 1000.0;
	gaurdBoundValue2[5] = 1000.0;
	gaurdBoundValue2[6] = 1000.0;
	gaurdBoundValue2[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope2 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix2, gaurdBoundValue2,
					gaurdBoundSign));

	// The transition label ist20

	// Original guard: 1 <= x1 & x1 <= 2 & x2 = 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix3.resize(row, col);
	gaurdConstraintsMatrix3(0, 0) = -1.0;
	gaurdConstraintsMatrix3(0, 1) = 0.0;
	gaurdConstraintsMatrix3(0, 2) = 0.0;
	gaurdConstraintsMatrix3(0, 3) = 0.0;
	gaurdConstraintsMatrix3(1, 0) = 1.0;
	gaurdConstraintsMatrix3(1, 1) = 0.0;
	gaurdConstraintsMatrix3(1, 2) = 0.0;
	gaurdConstraintsMatrix3(1, 3) = 0.0;
	gaurdConstraintsMatrix3(2, 0) = 0.0;
	gaurdConstraintsMatrix3(2, 1) = -1.0;
	gaurdConstraintsMatrix3(2, 2) = 0.0;
	gaurdConstraintsMatrix3(2, 3) = 0.0;
	gaurdConstraintsMatrix3(3, 0) = 0.0;
	gaurdConstraintsMatrix3(3, 1) = 1.0;
	gaurdConstraintsMatrix3(3, 2) = 0.0;
	gaurdConstraintsMatrix3(3, 3) = 0.0;
	gaurdConstraintsMatrix3(4, 0) = 0.0;
	gaurdConstraintsMatrix3(4, 1) = 0.0;
	gaurdConstraintsMatrix3(4, 2) = -1.0;
	gaurdConstraintsMatrix3(4, 3) = 0.0;
	gaurdConstraintsMatrix3(5, 0) = 0.0;
	gaurdConstraintsMatrix3(5, 1) = 0.0;
	gaurdConstraintsMatrix3(5, 2) = 1.0;
	gaurdConstraintsMatrix3(5, 3) = 0.0;
	gaurdConstraintsMatrix3(6, 0) = 0.0;
	gaurdConstraintsMatrix3(6, 1) = 0.0;
	gaurdConstraintsMatrix3(6, 2) = 0.0;
	gaurdConstraintsMatrix3(6, 3) = -1.0;
	gaurdConstraintsMatrix3(7, 0) = 0.0;
	gaurdConstraintsMatrix3(7, 1) = 0.0;
	gaurdConstraintsMatrix3(7, 2) = 0.0;
	gaurdConstraintsMatrix3(7, 3) = 1.0;

	gaurdBoundValue3.resize(row);
	gaurdBoundValue3[0] = -1.0;
	gaurdBoundValue3[1] = 2.0;
	gaurdBoundValue3[2] = -2.0;
	gaurdBoundValue3[3] = 2.0;
	gaurdBoundValue3[4] = 1000.0;
	gaurdBoundValue3[5] = 1000.0;
	gaurdBoundValue3[6] = 1000.0;
	gaurdBoundValue3[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope3 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix3, gaurdBoundValue3,
					gaurdBoundSign));

	// The transition label ist21

	// Original guard: x1 = 2 & 2 <= x2 & x2 <= 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix4.resize(row, col);
	gaurdConstraintsMatrix4(0, 0) = -1.0;
	gaurdConstraintsMatrix4(0, 1) = 0.0;
	gaurdConstraintsMatrix4(0, 2) = 0.0;
	gaurdConstraintsMatrix4(0, 3) = 0.0;
	gaurdConstraintsMatrix4(1, 0) = 1.0;
	gaurdConstraintsMatrix4(1, 1) = 0.0;
	gaurdConstraintsMatrix4(1, 2) = 0.0;
	gaurdConstraintsMatrix4(1, 3) = 0.0;
	gaurdConstraintsMatrix4(2, 0) = 0.0;
	gaurdConstraintsMatrix4(2, 1) = -1.0;
	gaurdConstraintsMatrix4(2, 2) = 0.0;
	gaurdConstraintsMatrix4(2, 3) = 0.0;
	gaurdConstraintsMatrix4(3, 0) = 0.0;
	gaurdConstraintsMatrix4(3, 1) = 1.0;
	gaurdConstraintsMatrix4(3, 2) = 0.0;
	gaurdConstraintsMatrix4(3, 3) = 0.0;
	gaurdConstraintsMatrix4(4, 0) = 0.0;
	gaurdConstraintsMatrix4(4, 1) = 0.0;
	gaurdConstraintsMatrix4(4, 2) = -1.0;
	gaurdConstraintsMatrix4(4, 3) = 0.0;
	gaurdConstraintsMatrix4(5, 0) = 0.0;
	gaurdConstraintsMatrix4(5, 1) = 0.0;
	gaurdConstraintsMatrix4(5, 2) = 1.0;
	gaurdConstraintsMatrix4(5, 3) = 0.0;
	gaurdConstraintsMatrix4(6, 0) = 0.0;
	gaurdConstraintsMatrix4(6, 1) = 0.0;
	gaurdConstraintsMatrix4(6, 2) = 0.0;
	gaurdConstraintsMatrix4(6, 3) = -1.0;
	gaurdConstraintsMatrix4(7, 0) = 0.0;
	gaurdConstraintsMatrix4(7, 1) = 0.0;
	gaurdConstraintsMatrix4(7, 2) = 0.0;
	gaurdConstraintsMatrix4(7, 3) = 1.0;

	gaurdBoundValue4.resize(row);
	gaurdBoundValue4[0] = -2.0;
	gaurdBoundValue4[1] = 2.0;
	gaurdBoundValue4[2] = -2.0;
	gaurdBoundValue4[3] = 3.0;
	gaurdBoundValue4[4] = 1000.0;
	gaurdBoundValue4[5] = 1000.0;
	gaurdBoundValue4[6] = 1000.0;
	gaurdBoundValue4[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope4 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix4, gaurdBoundValue4,
					gaurdBoundSign));

	// The transition label ist19

	// Original guard: x1 = 1 & 2 <= x2 & x2 <= 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix5.resize(row, col);
	gaurdConstraintsMatrix5(0, 0) = -1.0;
	gaurdConstraintsMatrix5(0, 1) = 0.0;
	gaurdConstraintsMatrix5(0, 2) = 0.0;
	gaurdConstraintsMatrix5(0, 3) = 0.0;
	gaurdConstraintsMatrix5(1, 0) = 1.0;
	gaurdConstraintsMatrix5(1, 1) = 0.0;
	gaurdConstraintsMatrix5(1, 2) = 0.0;
	gaurdConstraintsMatrix5(1, 3) = 0.0;
	gaurdConstraintsMatrix5(2, 0) = 0.0;
	gaurdConstraintsMatrix5(2, 1) = -1.0;
	gaurdConstraintsMatrix5(2, 2) = 0.0;
	gaurdConstraintsMatrix5(2, 3) = 0.0;
	gaurdConstraintsMatrix5(3, 0) = 0.0;
	gaurdConstraintsMatrix5(3, 1) = 1.0;
	gaurdConstraintsMatrix5(3, 2) = 0.0;
	gaurdConstraintsMatrix5(3, 3) = 0.0;
	gaurdConstraintsMatrix5(4, 0) = 0.0;
	gaurdConstraintsMatrix5(4, 1) = 0.0;
	gaurdConstraintsMatrix5(4, 2) = -1.0;
	gaurdConstraintsMatrix5(4, 3) = 0.0;
	gaurdConstraintsMatrix5(5, 0) = 0.0;
	gaurdConstraintsMatrix5(5, 1) = 0.0;
	gaurdConstraintsMatrix5(5, 2) = 1.0;
	gaurdConstraintsMatrix5(5, 3) = 0.0;
	gaurdConstraintsMatrix5(6, 0) = 0.0;
	gaurdConstraintsMatrix5(6, 1) = 0.0;
	gaurdConstraintsMatrix5(6, 2) = 0.0;
	gaurdConstraintsMatrix5(6, 3) = -1.0;
	gaurdConstraintsMatrix5(7, 0) = 0.0;
	gaurdConstraintsMatrix5(7, 1) = 0.0;
	gaurdConstraintsMatrix5(7, 2) = 0.0;
	gaurdConstraintsMatrix5(7, 3) = 1.0;

	gaurdBoundValue5.resize(row);
	gaurdBoundValue5[0] = -1.0;
	gaurdBoundValue5[1] = 1.0;
	gaurdBoundValue5[2] = -2.0;
	gaurdBoundValue5[3] = 3.0;
	gaurdBoundValue5[4] = 1000.0;
	gaurdBoundValue5[5] = 1000.0;
	gaurdBoundValue5[6] = 1000.0;
	gaurdBoundValue5[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope5 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix5, gaurdBoundValue5,
					gaurdBoundSign));

	// The transition label ist18

	// Original guard: 1 <= x1 & x1 <= 2 & x2 = 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix6.resize(row, col);
	gaurdConstraintsMatrix6(0, 0) = -1.0;
	gaurdConstraintsMatrix6(0, 1) = 0.0;
	gaurdConstraintsMatrix6(0, 2) = 0.0;
	gaurdConstraintsMatrix6(0, 3) = 0.0;
	gaurdConstraintsMatrix6(1, 0) = 1.0;
	gaurdConstraintsMatrix6(1, 1) = 0.0;
	gaurdConstraintsMatrix6(1, 2) = 0.0;
	gaurdConstraintsMatrix6(1, 3) = 0.0;
	gaurdConstraintsMatrix6(2, 0) = 0.0;
	gaurdConstraintsMatrix6(2, 1) = -1.0;
	gaurdConstraintsMatrix6(2, 2) = 0.0;
	gaurdConstraintsMatrix6(2, 3) = 0.0;
	gaurdConstraintsMatrix6(3, 0) = 0.0;
	gaurdConstraintsMatrix6(3, 1) = 1.0;
	gaurdConstraintsMatrix6(3, 2) = 0.0;
	gaurdConstraintsMatrix6(3, 3) = 0.0;
	gaurdConstraintsMatrix6(4, 0) = 0.0;
	gaurdConstraintsMatrix6(4, 1) = 0.0;
	gaurdConstraintsMatrix6(4, 2) = -1.0;
	gaurdConstraintsMatrix6(4, 3) = 0.0;
	gaurdConstraintsMatrix6(5, 0) = 0.0;
	gaurdConstraintsMatrix6(5, 1) = 0.0;
	gaurdConstraintsMatrix6(5, 2) = 1.0;
	gaurdConstraintsMatrix6(5, 3) = 0.0;
	gaurdConstraintsMatrix6(6, 0) = 0.0;
	gaurdConstraintsMatrix6(6, 1) = 0.0;
	gaurdConstraintsMatrix6(6, 2) = 0.0;
	gaurdConstraintsMatrix6(6, 3) = -1.0;
	gaurdConstraintsMatrix6(7, 0) = 0.0;
	gaurdConstraintsMatrix6(7, 1) = 0.0;
	gaurdConstraintsMatrix6(7, 2) = 0.0;
	gaurdConstraintsMatrix6(7, 3) = 1.0;

	gaurdBoundValue6.resize(row);
	gaurdBoundValue6[0] = -1.0;
	gaurdBoundValue6[1] = 2.0;
	gaurdBoundValue6[2] = -3.0;
	gaurdBoundValue6[3] = 3.0;
	gaurdBoundValue6[4] = 1000.0;
	gaurdBoundValue6[5] = 1000.0;
	gaurdBoundValue6[6] = 1000.0;
	gaurdBoundValue6[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope6 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix6, gaurdBoundValue6,
					gaurdBoundSign));

	// The transition label ist9

	// Original guard: 0 <= x1 & x1 <= 1 & x2 = 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix7.resize(row, col);
	gaurdConstraintsMatrix7(0, 0) = -1.0;
	gaurdConstraintsMatrix7(0, 1) = 0.0;
	gaurdConstraintsMatrix7(0, 2) = 0.0;
	gaurdConstraintsMatrix7(0, 3) = 0.0;
	gaurdConstraintsMatrix7(1, 0) = 1.0;
	gaurdConstraintsMatrix7(1, 1) = 0.0;
	gaurdConstraintsMatrix7(1, 2) = 0.0;
	gaurdConstraintsMatrix7(1, 3) = 0.0;
	gaurdConstraintsMatrix7(2, 0) = 0.0;
	gaurdConstraintsMatrix7(2, 1) = -1.0;
	gaurdConstraintsMatrix7(2, 2) = 0.0;
	gaurdConstraintsMatrix7(2, 3) = 0.0;
	gaurdConstraintsMatrix7(3, 0) = 0.0;
	gaurdConstraintsMatrix7(3, 1) = 1.0;
	gaurdConstraintsMatrix7(3, 2) = 0.0;
	gaurdConstraintsMatrix7(3, 3) = 0.0;
	gaurdConstraintsMatrix7(4, 0) = 0.0;
	gaurdConstraintsMatrix7(4, 1) = 0.0;
	gaurdConstraintsMatrix7(4, 2) = -1.0;
	gaurdConstraintsMatrix7(4, 3) = 0.0;
	gaurdConstraintsMatrix7(5, 0) = 0.0;
	gaurdConstraintsMatrix7(5, 1) = 0.0;
	gaurdConstraintsMatrix7(5, 2) = 1.0;
	gaurdConstraintsMatrix7(5, 3) = 0.0;
	gaurdConstraintsMatrix7(6, 0) = 0.0;
	gaurdConstraintsMatrix7(6, 1) = 0.0;
	gaurdConstraintsMatrix7(6, 2) = 0.0;
	gaurdConstraintsMatrix7(6, 3) = -1.0;
	gaurdConstraintsMatrix7(7, 0) = 0.0;
	gaurdConstraintsMatrix7(7, 1) = 0.0;
	gaurdConstraintsMatrix7(7, 2) = 0.0;
	gaurdConstraintsMatrix7(7, 3) = 1.0;

	gaurdBoundValue7.resize(row);
	gaurdBoundValue7[0] = -0.0;
	gaurdBoundValue7[1] = 1.0;
	gaurdBoundValue7[2] = -2.0;
	gaurdBoundValue7[3] = 2.0;
	gaurdBoundValue7[4] = 1000.0;
	gaurdBoundValue7[5] = 1000.0;
	gaurdBoundValue7[6] = 1000.0;
	gaurdBoundValue7[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope7 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix7, gaurdBoundValue7,
					gaurdBoundSign));

	// The transition label ist11

	// Original guard: x1 = 1 & 1 <= x2 & x2 <= 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix8.resize(row, col);
	gaurdConstraintsMatrix8(0, 0) = -1.0;
	gaurdConstraintsMatrix8(0, 1) = 0.0;
	gaurdConstraintsMatrix8(0, 2) = 0.0;
	gaurdConstraintsMatrix8(0, 3) = 0.0;
	gaurdConstraintsMatrix8(1, 0) = 1.0;
	gaurdConstraintsMatrix8(1, 1) = 0.0;
	gaurdConstraintsMatrix8(1, 2) = 0.0;
	gaurdConstraintsMatrix8(1, 3) = 0.0;
	gaurdConstraintsMatrix8(2, 0) = 0.0;
	gaurdConstraintsMatrix8(2, 1) = -1.0;
	gaurdConstraintsMatrix8(2, 2) = 0.0;
	gaurdConstraintsMatrix8(2, 3) = 0.0;
	gaurdConstraintsMatrix8(3, 0) = 0.0;
	gaurdConstraintsMatrix8(3, 1) = 1.0;
	gaurdConstraintsMatrix8(3, 2) = 0.0;
	gaurdConstraintsMatrix8(3, 3) = 0.0;
	gaurdConstraintsMatrix8(4, 0) = 0.0;
	gaurdConstraintsMatrix8(4, 1) = 0.0;
	gaurdConstraintsMatrix8(4, 2) = -1.0;
	gaurdConstraintsMatrix8(4, 3) = 0.0;
	gaurdConstraintsMatrix8(5, 0) = 0.0;
	gaurdConstraintsMatrix8(5, 1) = 0.0;
	gaurdConstraintsMatrix8(5, 2) = 1.0;
	gaurdConstraintsMatrix8(5, 3) = 0.0;
	gaurdConstraintsMatrix8(6, 0) = 0.0;
	gaurdConstraintsMatrix8(6, 1) = 0.0;
	gaurdConstraintsMatrix8(6, 2) = 0.0;
	gaurdConstraintsMatrix8(6, 3) = -1.0;
	gaurdConstraintsMatrix8(7, 0) = 0.0;
	gaurdConstraintsMatrix8(7, 1) = 0.0;
	gaurdConstraintsMatrix8(7, 2) = 0.0;
	gaurdConstraintsMatrix8(7, 3) = 1.0;

	gaurdBoundValue8.resize(row);
	gaurdBoundValue8[0] = -1.0;
	gaurdBoundValue8[1] = 1.0;
	gaurdBoundValue8[2] = -1.0;
	gaurdBoundValue8[3] = 2.0;
	gaurdBoundValue8[4] = 1000.0;
	gaurdBoundValue8[5] = 1000.0;
	gaurdBoundValue8[6] = 1000.0;
	gaurdBoundValue8[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope8 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix8, gaurdBoundValue8,
					gaurdBoundSign));

	// The transition label ist10

	// Original guard: 0 <= x1 & x1 <= 1 & x2 = 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix9.resize(row, col);
	gaurdConstraintsMatrix9(0, 0) = -1.0;
	gaurdConstraintsMatrix9(0, 1) = 0.0;
	gaurdConstraintsMatrix9(0, 2) = 0.0;
	gaurdConstraintsMatrix9(0, 3) = 0.0;
	gaurdConstraintsMatrix9(1, 0) = 1.0;
	gaurdConstraintsMatrix9(1, 1) = 0.0;
	gaurdConstraintsMatrix9(1, 2) = 0.0;
	gaurdConstraintsMatrix9(1, 3) = 0.0;
	gaurdConstraintsMatrix9(2, 0) = 0.0;
	gaurdConstraintsMatrix9(2, 1) = -1.0;
	gaurdConstraintsMatrix9(2, 2) = 0.0;
	gaurdConstraintsMatrix9(2, 3) = 0.0;
	gaurdConstraintsMatrix9(3, 0) = 0.0;
	gaurdConstraintsMatrix9(3, 1) = 1.0;
	gaurdConstraintsMatrix9(3, 2) = 0.0;
	gaurdConstraintsMatrix9(3, 3) = 0.0;
	gaurdConstraintsMatrix9(4, 0) = 0.0;
	gaurdConstraintsMatrix9(4, 1) = 0.0;
	gaurdConstraintsMatrix9(4, 2) = -1.0;
	gaurdConstraintsMatrix9(4, 3) = 0.0;
	gaurdConstraintsMatrix9(5, 0) = 0.0;
	gaurdConstraintsMatrix9(5, 1) = 0.0;
	gaurdConstraintsMatrix9(5, 2) = 1.0;
	gaurdConstraintsMatrix9(5, 3) = 0.0;
	gaurdConstraintsMatrix9(6, 0) = 0.0;
	gaurdConstraintsMatrix9(6, 1) = 0.0;
	gaurdConstraintsMatrix9(6, 2) = 0.0;
	gaurdConstraintsMatrix9(6, 3) = -1.0;
	gaurdConstraintsMatrix9(7, 0) = 0.0;
	gaurdConstraintsMatrix9(7, 1) = 0.0;
	gaurdConstraintsMatrix9(7, 2) = 0.0;
	gaurdConstraintsMatrix9(7, 3) = 1.0;

	gaurdBoundValue9.resize(row);
	gaurdBoundValue9[0] = -0.0;
	gaurdBoundValue9[1] = 1.0;
	gaurdBoundValue9[2] = -1.0;
	gaurdBoundValue9[3] = 1.0;
	gaurdBoundValue9[4] = 1000.0;
	gaurdBoundValue9[5] = 1000.0;
	gaurdBoundValue9[6] = 1000.0;
	gaurdBoundValue9[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope9 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix9, gaurdBoundValue9,
					gaurdBoundSign));

	// The transition label ist15

	// Original guard: x1 = 1 & 1 <= x2 & x2 <= 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix10.resize(row, col);
	gaurdConstraintsMatrix10(0, 0) = -1.0;
	gaurdConstraintsMatrix10(0, 1) = 0.0;
	gaurdConstraintsMatrix10(0, 2) = 0.0;
	gaurdConstraintsMatrix10(0, 3) = 0.0;
	gaurdConstraintsMatrix10(1, 0) = 1.0;
	gaurdConstraintsMatrix10(1, 1) = 0.0;
	gaurdConstraintsMatrix10(1, 2) = 0.0;
	gaurdConstraintsMatrix10(1, 3) = 0.0;
	gaurdConstraintsMatrix10(2, 0) = 0.0;
	gaurdConstraintsMatrix10(2, 1) = -1.0;
	gaurdConstraintsMatrix10(2, 2) = 0.0;
	gaurdConstraintsMatrix10(2, 3) = 0.0;
	gaurdConstraintsMatrix10(3, 0) = 0.0;
	gaurdConstraintsMatrix10(3, 1) = 1.0;
	gaurdConstraintsMatrix10(3, 2) = 0.0;
	gaurdConstraintsMatrix10(3, 3) = 0.0;
	gaurdConstraintsMatrix10(4, 0) = 0.0;
	gaurdConstraintsMatrix10(4, 1) = 0.0;
	gaurdConstraintsMatrix10(4, 2) = -1.0;
	gaurdConstraintsMatrix10(4, 3) = 0.0;
	gaurdConstraintsMatrix10(5, 0) = 0.0;
	gaurdConstraintsMatrix10(5, 1) = 0.0;
	gaurdConstraintsMatrix10(5, 2) = 1.0;
	gaurdConstraintsMatrix10(5, 3) = 0.0;
	gaurdConstraintsMatrix10(6, 0) = 0.0;
	gaurdConstraintsMatrix10(6, 1) = 0.0;
	gaurdConstraintsMatrix10(6, 2) = 0.0;
	gaurdConstraintsMatrix10(6, 3) = -1.0;
	gaurdConstraintsMatrix10(7, 0) = 0.0;
	gaurdConstraintsMatrix10(7, 1) = 0.0;
	gaurdConstraintsMatrix10(7, 2) = 0.0;
	gaurdConstraintsMatrix10(7, 3) = 1.0;

	gaurdBoundValue10.resize(row);
	gaurdBoundValue10[0] = -1.0;
	gaurdBoundValue10[1] = 1.0;
	gaurdBoundValue10[2] = -1.0;
	gaurdBoundValue10[3] = 2.0;
	gaurdBoundValue10[4] = 1000.0;
	gaurdBoundValue10[5] = 1000.0;
	gaurdBoundValue10[6] = 1000.0;
	gaurdBoundValue10[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope10 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix10, gaurdBoundValue10,
					gaurdBoundSign));

	// The transition label ist14

	// Original guard: 1 <= x1 & x1 <= 2 & x2 = 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix11.resize(row, col);
	gaurdConstraintsMatrix11(0, 0) = -1.0;
	gaurdConstraintsMatrix11(0, 1) = 0.0;
	gaurdConstraintsMatrix11(0, 2) = 0.0;
	gaurdConstraintsMatrix11(0, 3) = 0.0;
	gaurdConstraintsMatrix11(1, 0) = 1.0;
	gaurdConstraintsMatrix11(1, 1) = 0.0;
	gaurdConstraintsMatrix11(1, 2) = 0.0;
	gaurdConstraintsMatrix11(1, 3) = 0.0;
	gaurdConstraintsMatrix11(2, 0) = 0.0;
	gaurdConstraintsMatrix11(2, 1) = -1.0;
	gaurdConstraintsMatrix11(2, 2) = 0.0;
	gaurdConstraintsMatrix11(2, 3) = 0.0;
	gaurdConstraintsMatrix11(3, 0) = 0.0;
	gaurdConstraintsMatrix11(3, 1) = 1.0;
	gaurdConstraintsMatrix11(3, 2) = 0.0;
	gaurdConstraintsMatrix11(3, 3) = 0.0;
	gaurdConstraintsMatrix11(4, 0) = 0.0;
	gaurdConstraintsMatrix11(4, 1) = 0.0;
	gaurdConstraintsMatrix11(4, 2) = -1.0;
	gaurdConstraintsMatrix11(4, 3) = 0.0;
	gaurdConstraintsMatrix11(5, 0) = 0.0;
	gaurdConstraintsMatrix11(5, 1) = 0.0;
	gaurdConstraintsMatrix11(5, 2) = 1.0;
	gaurdConstraintsMatrix11(5, 3) = 0.0;
	gaurdConstraintsMatrix11(6, 0) = 0.0;
	gaurdConstraintsMatrix11(6, 1) = 0.0;
	gaurdConstraintsMatrix11(6, 2) = 0.0;
	gaurdConstraintsMatrix11(6, 3) = -1.0;
	gaurdConstraintsMatrix11(7, 0) = 0.0;
	gaurdConstraintsMatrix11(7, 1) = 0.0;
	gaurdConstraintsMatrix11(7, 2) = 0.0;
	gaurdConstraintsMatrix11(7, 3) = 1.0;

	gaurdBoundValue11.resize(row);
	gaurdBoundValue11[0] = -1.0;
	gaurdBoundValue11[1] = 2.0;
	gaurdBoundValue11[2] = -2.0;
	gaurdBoundValue11[3] = 2.0;
	gaurdBoundValue11[4] = 1000.0;
	gaurdBoundValue11[5] = 1000.0;
	gaurdBoundValue11[6] = 1000.0;
	gaurdBoundValue11[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope11 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix11, gaurdBoundValue11,
					gaurdBoundSign));

	// The transition label ist17

	// Original guard: x1 = 2 & 1 <= x2 & x2 <= 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix12.resize(row, col);
	gaurdConstraintsMatrix12(0, 0) = -1.0;
	gaurdConstraintsMatrix12(0, 1) = 0.0;
	gaurdConstraintsMatrix12(0, 2) = 0.0;
	gaurdConstraintsMatrix12(0, 3) = 0.0;
	gaurdConstraintsMatrix12(1, 0) = 1.0;
	gaurdConstraintsMatrix12(1, 1) = 0.0;
	gaurdConstraintsMatrix12(1, 2) = 0.0;
	gaurdConstraintsMatrix12(1, 3) = 0.0;
	gaurdConstraintsMatrix12(2, 0) = 0.0;
	gaurdConstraintsMatrix12(2, 1) = -1.0;
	gaurdConstraintsMatrix12(2, 2) = 0.0;
	gaurdConstraintsMatrix12(2, 3) = 0.0;
	gaurdConstraintsMatrix12(3, 0) = 0.0;
	gaurdConstraintsMatrix12(3, 1) = 1.0;
	gaurdConstraintsMatrix12(3, 2) = 0.0;
	gaurdConstraintsMatrix12(3, 3) = 0.0;
	gaurdConstraintsMatrix12(4, 0) = 0.0;
	gaurdConstraintsMatrix12(4, 1) = 0.0;
	gaurdConstraintsMatrix12(4, 2) = -1.0;
	gaurdConstraintsMatrix12(4, 3) = 0.0;
	gaurdConstraintsMatrix12(5, 0) = 0.0;
	gaurdConstraintsMatrix12(5, 1) = 0.0;
	gaurdConstraintsMatrix12(5, 2) = 1.0;
	gaurdConstraintsMatrix12(5, 3) = 0.0;
	gaurdConstraintsMatrix12(6, 0) = 0.0;
	gaurdConstraintsMatrix12(6, 1) = 0.0;
	gaurdConstraintsMatrix12(6, 2) = 0.0;
	gaurdConstraintsMatrix12(6, 3) = -1.0;
	gaurdConstraintsMatrix12(7, 0) = 0.0;
	gaurdConstraintsMatrix12(7, 1) = 0.0;
	gaurdConstraintsMatrix12(7, 2) = 0.0;
	gaurdConstraintsMatrix12(7, 3) = 1.0;

	gaurdBoundValue12.resize(row);
	gaurdBoundValue12[0] = -2.0;
	gaurdBoundValue12[1] = 2.0;
	gaurdBoundValue12[2] = -1.0;
	gaurdBoundValue12[3] = 2.0;
	gaurdBoundValue12[4] = 1000.0;
	gaurdBoundValue12[5] = 1000.0;
	gaurdBoundValue12[6] = 1000.0;
	gaurdBoundValue12[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope12 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix12, gaurdBoundValue12,
					gaurdBoundSign));

	// The transition label ist16

	// Original guard: 1 <= x1 & x1 <= 2 & x2 = 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix13.resize(row, col);
	gaurdConstraintsMatrix13(0, 0) = -1.0;
	gaurdConstraintsMatrix13(0, 1) = 0.0;
	gaurdConstraintsMatrix13(0, 2) = 0.0;
	gaurdConstraintsMatrix13(0, 3) = 0.0;
	gaurdConstraintsMatrix13(1, 0) = 1.0;
	gaurdConstraintsMatrix13(1, 1) = 0.0;
	gaurdConstraintsMatrix13(1, 2) = 0.0;
	gaurdConstraintsMatrix13(1, 3) = 0.0;
	gaurdConstraintsMatrix13(2, 0) = 0.0;
	gaurdConstraintsMatrix13(2, 1) = -1.0;
	gaurdConstraintsMatrix13(2, 2) = 0.0;
	gaurdConstraintsMatrix13(2, 3) = 0.0;
	gaurdConstraintsMatrix13(3, 0) = 0.0;
	gaurdConstraintsMatrix13(3, 1) = 1.0;
	gaurdConstraintsMatrix13(3, 2) = 0.0;
	gaurdConstraintsMatrix13(3, 3) = 0.0;
	gaurdConstraintsMatrix13(4, 0) = 0.0;
	gaurdConstraintsMatrix13(4, 1) = 0.0;
	gaurdConstraintsMatrix13(4, 2) = -1.0;
	gaurdConstraintsMatrix13(4, 3) = 0.0;
	gaurdConstraintsMatrix13(5, 0) = 0.0;
	gaurdConstraintsMatrix13(5, 1) = 0.0;
	gaurdConstraintsMatrix13(5, 2) = 1.0;
	gaurdConstraintsMatrix13(5, 3) = 0.0;
	gaurdConstraintsMatrix13(6, 0) = 0.0;
	gaurdConstraintsMatrix13(6, 1) = 0.0;
	gaurdConstraintsMatrix13(6, 2) = 0.0;
	gaurdConstraintsMatrix13(6, 3) = -1.0;
	gaurdConstraintsMatrix13(7, 0) = 0.0;
	gaurdConstraintsMatrix13(7, 1) = 0.0;
	gaurdConstraintsMatrix13(7, 2) = 0.0;
	gaurdConstraintsMatrix13(7, 3) = 1.0;

	gaurdBoundValue13.resize(row);
	gaurdBoundValue13[0] = -1.0;
	gaurdBoundValue13[1] = 2.0;
	gaurdBoundValue13[2] = -1.0;
	gaurdBoundValue13[3] = 1.0;
	gaurdBoundValue13[4] = 1000.0;
	gaurdBoundValue13[5] = 1000.0;
	gaurdBoundValue13[6] = 1000.0;
	gaurdBoundValue13[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope13 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix13, gaurdBoundValue13,
					gaurdBoundSign));

	// The transition label ist32

	// Original guard: 2 <= x1 & x1 <= 3 & x2 = 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix14.resize(row, col);
	gaurdConstraintsMatrix14(0, 0) = -1.0;
	gaurdConstraintsMatrix14(0, 1) = 0.0;
	gaurdConstraintsMatrix14(0, 2) = 0.0;
	gaurdConstraintsMatrix14(0, 3) = 0.0;
	gaurdConstraintsMatrix14(1, 0) = 1.0;
	gaurdConstraintsMatrix14(1, 1) = 0.0;
	gaurdConstraintsMatrix14(1, 2) = 0.0;
	gaurdConstraintsMatrix14(1, 3) = 0.0;
	gaurdConstraintsMatrix14(2, 0) = 0.0;
	gaurdConstraintsMatrix14(2, 1) = -1.0;
	gaurdConstraintsMatrix14(2, 2) = 0.0;
	gaurdConstraintsMatrix14(2, 3) = 0.0;
	gaurdConstraintsMatrix14(3, 0) = 0.0;
	gaurdConstraintsMatrix14(3, 1) = 1.0;
	gaurdConstraintsMatrix14(3, 2) = 0.0;
	gaurdConstraintsMatrix14(3, 3) = 0.0;
	gaurdConstraintsMatrix14(4, 0) = 0.0;
	gaurdConstraintsMatrix14(4, 1) = 0.0;
	gaurdConstraintsMatrix14(4, 2) = -1.0;
	gaurdConstraintsMatrix14(4, 3) = 0.0;
	gaurdConstraintsMatrix14(5, 0) = 0.0;
	gaurdConstraintsMatrix14(5, 1) = 0.0;
	gaurdConstraintsMatrix14(5, 2) = 1.0;
	gaurdConstraintsMatrix14(5, 3) = 0.0;
	gaurdConstraintsMatrix14(6, 0) = 0.0;
	gaurdConstraintsMatrix14(6, 1) = 0.0;
	gaurdConstraintsMatrix14(6, 2) = 0.0;
	gaurdConstraintsMatrix14(6, 3) = -1.0;
	gaurdConstraintsMatrix14(7, 0) = 0.0;
	gaurdConstraintsMatrix14(7, 1) = 0.0;
	gaurdConstraintsMatrix14(7, 2) = 0.0;
	gaurdConstraintsMatrix14(7, 3) = 1.0;

	gaurdBoundValue14.resize(row);
	gaurdBoundValue14[0] = -2.0;
	gaurdBoundValue14[1] = 3.0;
	gaurdBoundValue14[2] = -2.0;
	gaurdBoundValue14[3] = 2.0;
	gaurdBoundValue14[4] = 1000.0;
	gaurdBoundValue14[5] = 1000.0;
	gaurdBoundValue14[6] = 1000.0;
	gaurdBoundValue14[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope14 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix14, gaurdBoundValue14,
					gaurdBoundSign));

	// The transition label ist33

	// Original guard: x1 = 2 & 1 <= x2 & x2 <= 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix15.resize(row, col);
	gaurdConstraintsMatrix15(0, 0) = -1.0;
	gaurdConstraintsMatrix15(0, 1) = 0.0;
	gaurdConstraintsMatrix15(0, 2) = 0.0;
	gaurdConstraintsMatrix15(0, 3) = 0.0;
	gaurdConstraintsMatrix15(1, 0) = 1.0;
	gaurdConstraintsMatrix15(1, 1) = 0.0;
	gaurdConstraintsMatrix15(1, 2) = 0.0;
	gaurdConstraintsMatrix15(1, 3) = 0.0;
	gaurdConstraintsMatrix15(2, 0) = 0.0;
	gaurdConstraintsMatrix15(2, 1) = -1.0;
	gaurdConstraintsMatrix15(2, 2) = 0.0;
	gaurdConstraintsMatrix15(2, 3) = 0.0;
	gaurdConstraintsMatrix15(3, 0) = 0.0;
	gaurdConstraintsMatrix15(3, 1) = 1.0;
	gaurdConstraintsMatrix15(3, 2) = 0.0;
	gaurdConstraintsMatrix15(3, 3) = 0.0;
	gaurdConstraintsMatrix15(4, 0) = 0.0;
	gaurdConstraintsMatrix15(4, 1) = 0.0;
	gaurdConstraintsMatrix15(4, 2) = -1.0;
	gaurdConstraintsMatrix15(4, 3) = 0.0;
	gaurdConstraintsMatrix15(5, 0) = 0.0;
	gaurdConstraintsMatrix15(5, 1) = 0.0;
	gaurdConstraintsMatrix15(5, 2) = 1.0;
	gaurdConstraintsMatrix15(5, 3) = 0.0;
	gaurdConstraintsMatrix15(6, 0) = 0.0;
	gaurdConstraintsMatrix15(6, 1) = 0.0;
	gaurdConstraintsMatrix15(6, 2) = 0.0;
	gaurdConstraintsMatrix15(6, 3) = -1.0;
	gaurdConstraintsMatrix15(7, 0) = 0.0;
	gaurdConstraintsMatrix15(7, 1) = 0.0;
	gaurdConstraintsMatrix15(7, 2) = 0.0;
	gaurdConstraintsMatrix15(7, 3) = 1.0;

	gaurdBoundValue15.resize(row);
	gaurdBoundValue15[0] = -2.0;
	gaurdBoundValue15[1] = 2.0;
	gaurdBoundValue15[2] = -1.0;
	gaurdBoundValue15[3] = 2.0;
	gaurdBoundValue15[4] = 1000.0;
	gaurdBoundValue15[5] = 1000.0;
	gaurdBoundValue15[6] = 1000.0;
	gaurdBoundValue15[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope15 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix15, gaurdBoundValue15,
					gaurdBoundSign));

	// The transition label ist34

	// Original guard: 2 <= x1 & x1 <= 3 & x2 = 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix16.resize(row, col);
	gaurdConstraintsMatrix16(0, 0) = -1.0;
	gaurdConstraintsMatrix16(0, 1) = 0.0;
	gaurdConstraintsMatrix16(0, 2) = 0.0;
	gaurdConstraintsMatrix16(0, 3) = 0.0;
	gaurdConstraintsMatrix16(1, 0) = 1.0;
	gaurdConstraintsMatrix16(1, 1) = 0.0;
	gaurdConstraintsMatrix16(1, 2) = 0.0;
	gaurdConstraintsMatrix16(1, 3) = 0.0;
	gaurdConstraintsMatrix16(2, 0) = 0.0;
	gaurdConstraintsMatrix16(2, 1) = -1.0;
	gaurdConstraintsMatrix16(2, 2) = 0.0;
	gaurdConstraintsMatrix16(2, 3) = 0.0;
	gaurdConstraintsMatrix16(3, 0) = 0.0;
	gaurdConstraintsMatrix16(3, 1) = 1.0;
	gaurdConstraintsMatrix16(3, 2) = 0.0;
	gaurdConstraintsMatrix16(3, 3) = 0.0;
	gaurdConstraintsMatrix16(4, 0) = 0.0;
	gaurdConstraintsMatrix16(4, 1) = 0.0;
	gaurdConstraintsMatrix16(4, 2) = -1.0;
	gaurdConstraintsMatrix16(4, 3) = 0.0;
	gaurdConstraintsMatrix16(5, 0) = 0.0;
	gaurdConstraintsMatrix16(5, 1) = 0.0;
	gaurdConstraintsMatrix16(5, 2) = 1.0;
	gaurdConstraintsMatrix16(5, 3) = 0.0;
	gaurdConstraintsMatrix16(6, 0) = 0.0;
	gaurdConstraintsMatrix16(6, 1) = 0.0;
	gaurdConstraintsMatrix16(6, 2) = 0.0;
	gaurdConstraintsMatrix16(6, 3) = -1.0;
	gaurdConstraintsMatrix16(7, 0) = 0.0;
	gaurdConstraintsMatrix16(7, 1) = 0.0;
	gaurdConstraintsMatrix16(7, 2) = 0.0;
	gaurdConstraintsMatrix16(7, 3) = 1.0;

	gaurdBoundValue16.resize(row);
	gaurdBoundValue16[0] = -2.0;
	gaurdBoundValue16[1] = 3.0;
	gaurdBoundValue16[2] = -1.0;
	gaurdBoundValue16[3] = 1.0;
	gaurdBoundValue16[4] = 1000.0;
	gaurdBoundValue16[5] = 1000.0;
	gaurdBoundValue16[6] = 1000.0;
	gaurdBoundValue16[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope16 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix16, gaurdBoundValue16,
					gaurdBoundSign));

	// The transition label ist35

	// Original guard: x1 = 3 & 1 <= x2 & x2 <= 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix17.resize(row, col);
	gaurdConstraintsMatrix17(0, 0) = -1.0;
	gaurdConstraintsMatrix17(0, 1) = 0.0;
	gaurdConstraintsMatrix17(0, 2) = 0.0;
	gaurdConstraintsMatrix17(0, 3) = 0.0;
	gaurdConstraintsMatrix17(1, 0) = 1.0;
	gaurdConstraintsMatrix17(1, 1) = 0.0;
	gaurdConstraintsMatrix17(1, 2) = 0.0;
	gaurdConstraintsMatrix17(1, 3) = 0.0;
	gaurdConstraintsMatrix17(2, 0) = 0.0;
	gaurdConstraintsMatrix17(2, 1) = -1.0;
	gaurdConstraintsMatrix17(2, 2) = 0.0;
	gaurdConstraintsMatrix17(2, 3) = 0.0;
	gaurdConstraintsMatrix17(3, 0) = 0.0;
	gaurdConstraintsMatrix17(3, 1) = 1.0;
	gaurdConstraintsMatrix17(3, 2) = 0.0;
	gaurdConstraintsMatrix17(3, 3) = 0.0;
	gaurdConstraintsMatrix17(4, 0) = 0.0;
	gaurdConstraintsMatrix17(4, 1) = 0.0;
	gaurdConstraintsMatrix17(4, 2) = -1.0;
	gaurdConstraintsMatrix17(4, 3) = 0.0;
	gaurdConstraintsMatrix17(5, 0) = 0.0;
	gaurdConstraintsMatrix17(5, 1) = 0.0;
	gaurdConstraintsMatrix17(5, 2) = 1.0;
	gaurdConstraintsMatrix17(5, 3) = 0.0;
	gaurdConstraintsMatrix17(6, 0) = 0.0;
	gaurdConstraintsMatrix17(6, 1) = 0.0;
	gaurdConstraintsMatrix17(6, 2) = 0.0;
	gaurdConstraintsMatrix17(6, 3) = -1.0;
	gaurdConstraintsMatrix17(7, 0) = 0.0;
	gaurdConstraintsMatrix17(7, 1) = 0.0;
	gaurdConstraintsMatrix17(7, 2) = 0.0;
	gaurdConstraintsMatrix17(7, 3) = 1.0;

	gaurdBoundValue17.resize(row);
	gaurdBoundValue17[0] = -3.0;
	gaurdBoundValue17[1] = 3.0;
	gaurdBoundValue17[2] = -1.0;
	gaurdBoundValue17[3] = 2.0;
	gaurdBoundValue17[4] = 1000.0;
	gaurdBoundValue17[5] = 1000.0;
	gaurdBoundValue17[6] = 1000.0;
	gaurdBoundValue17[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope17 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix17, gaurdBoundValue17,
					gaurdBoundSign));

	// The transition label ist31

	// Original guard: x1 = 3 & 0 <= x2 & x2 <= 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix18.resize(row, col);
	gaurdConstraintsMatrix18(0, 0) = -1.0;
	gaurdConstraintsMatrix18(0, 1) = 0.0;
	gaurdConstraintsMatrix18(0, 2) = 0.0;
	gaurdConstraintsMatrix18(0, 3) = 0.0;
	gaurdConstraintsMatrix18(1, 0) = 1.0;
	gaurdConstraintsMatrix18(1, 1) = 0.0;
	gaurdConstraintsMatrix18(1, 2) = 0.0;
	gaurdConstraintsMatrix18(1, 3) = 0.0;
	gaurdConstraintsMatrix18(2, 0) = 0.0;
	gaurdConstraintsMatrix18(2, 1) = -1.0;
	gaurdConstraintsMatrix18(2, 2) = 0.0;
	gaurdConstraintsMatrix18(2, 3) = 0.0;
	gaurdConstraintsMatrix18(3, 0) = 0.0;
	gaurdConstraintsMatrix18(3, 1) = 1.0;
	gaurdConstraintsMatrix18(3, 2) = 0.0;
	gaurdConstraintsMatrix18(3, 3) = 0.0;
	gaurdConstraintsMatrix18(4, 0) = 0.0;
	gaurdConstraintsMatrix18(4, 1) = 0.0;
	gaurdConstraintsMatrix18(4, 2) = -1.0;
	gaurdConstraintsMatrix18(4, 3) = 0.0;
	gaurdConstraintsMatrix18(5, 0) = 0.0;
	gaurdConstraintsMatrix18(5, 1) = 0.0;
	gaurdConstraintsMatrix18(5, 2) = 1.0;
	gaurdConstraintsMatrix18(5, 3) = 0.0;
	gaurdConstraintsMatrix18(6, 0) = 0.0;
	gaurdConstraintsMatrix18(6, 1) = 0.0;
	gaurdConstraintsMatrix18(6, 2) = 0.0;
	gaurdConstraintsMatrix18(6, 3) = -1.0;
	gaurdConstraintsMatrix18(7, 0) = 0.0;
	gaurdConstraintsMatrix18(7, 1) = 0.0;
	gaurdConstraintsMatrix18(7, 2) = 0.0;
	gaurdConstraintsMatrix18(7, 3) = 1.0;

	gaurdBoundValue18.resize(row);
	gaurdBoundValue18[0] = -3.0;
	gaurdBoundValue18[1] = 3.0;
	gaurdBoundValue18[2] = -0.0;
	gaurdBoundValue18[3] = 1.0;
	gaurdBoundValue18[4] = 1000.0;
	gaurdBoundValue18[5] = 1000.0;
	gaurdBoundValue18[6] = 1000.0;
	gaurdBoundValue18[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope18 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix18, gaurdBoundValue18,
					gaurdBoundSign));

	// The transition label ist30

	// Original guard: x1 = 2 & 0 <= x2 & x2 <= 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix19.resize(row, col);
	gaurdConstraintsMatrix19(0, 0) = -1.0;
	gaurdConstraintsMatrix19(0, 1) = 0.0;
	gaurdConstraintsMatrix19(0, 2) = 0.0;
	gaurdConstraintsMatrix19(0, 3) = 0.0;
	gaurdConstraintsMatrix19(1, 0) = 1.0;
	gaurdConstraintsMatrix19(1, 1) = 0.0;
	gaurdConstraintsMatrix19(1, 2) = 0.0;
	gaurdConstraintsMatrix19(1, 3) = 0.0;
	gaurdConstraintsMatrix19(2, 0) = 0.0;
	gaurdConstraintsMatrix19(2, 1) = -1.0;
	gaurdConstraintsMatrix19(2, 2) = 0.0;
	gaurdConstraintsMatrix19(2, 3) = 0.0;
	gaurdConstraintsMatrix19(3, 0) = 0.0;
	gaurdConstraintsMatrix19(3, 1) = 1.0;
	gaurdConstraintsMatrix19(3, 2) = 0.0;
	gaurdConstraintsMatrix19(3, 3) = 0.0;
	gaurdConstraintsMatrix19(4, 0) = 0.0;
	gaurdConstraintsMatrix19(4, 1) = 0.0;
	gaurdConstraintsMatrix19(4, 2) = -1.0;
	gaurdConstraintsMatrix19(4, 3) = 0.0;
	gaurdConstraintsMatrix19(5, 0) = 0.0;
	gaurdConstraintsMatrix19(5, 1) = 0.0;
	gaurdConstraintsMatrix19(5, 2) = 1.0;
	gaurdConstraintsMatrix19(5, 3) = 0.0;
	gaurdConstraintsMatrix19(6, 0) = 0.0;
	gaurdConstraintsMatrix19(6, 1) = 0.0;
	gaurdConstraintsMatrix19(6, 2) = 0.0;
	gaurdConstraintsMatrix19(6, 3) = -1.0;
	gaurdConstraintsMatrix19(7, 0) = 0.0;
	gaurdConstraintsMatrix19(7, 1) = 0.0;
	gaurdConstraintsMatrix19(7, 2) = 0.0;
	gaurdConstraintsMatrix19(7, 3) = 1.0;

	gaurdBoundValue19.resize(row);
	gaurdBoundValue19[0] = -2.0;
	gaurdBoundValue19[1] = 2.0;
	gaurdBoundValue19[2] = -0.0;
	gaurdBoundValue19[3] = 1.0;
	gaurdBoundValue19[4] = 1000.0;
	gaurdBoundValue19[5] = 1000.0;
	gaurdBoundValue19[6] = 1000.0;
	gaurdBoundValue19[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope19 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix19, gaurdBoundValue19,
					gaurdBoundSign));

	// The transition label ist29

	// Original guard: 2 <= x1 & x1 <= 3 & x2 = 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix20.resize(row, col);
	gaurdConstraintsMatrix20(0, 0) = -1.0;
	gaurdConstraintsMatrix20(0, 1) = 0.0;
	gaurdConstraintsMatrix20(0, 2) = 0.0;
	gaurdConstraintsMatrix20(0, 3) = 0.0;
	gaurdConstraintsMatrix20(1, 0) = 1.0;
	gaurdConstraintsMatrix20(1, 1) = 0.0;
	gaurdConstraintsMatrix20(1, 2) = 0.0;
	gaurdConstraintsMatrix20(1, 3) = 0.0;
	gaurdConstraintsMatrix20(2, 0) = 0.0;
	gaurdConstraintsMatrix20(2, 1) = -1.0;
	gaurdConstraintsMatrix20(2, 2) = 0.0;
	gaurdConstraintsMatrix20(2, 3) = 0.0;
	gaurdConstraintsMatrix20(3, 0) = 0.0;
	gaurdConstraintsMatrix20(3, 1) = 1.0;
	gaurdConstraintsMatrix20(3, 2) = 0.0;
	gaurdConstraintsMatrix20(3, 3) = 0.0;
	gaurdConstraintsMatrix20(4, 0) = 0.0;
	gaurdConstraintsMatrix20(4, 1) = 0.0;
	gaurdConstraintsMatrix20(4, 2) = -1.0;
	gaurdConstraintsMatrix20(4, 3) = 0.0;
	gaurdConstraintsMatrix20(5, 0) = 0.0;
	gaurdConstraintsMatrix20(5, 1) = 0.0;
	gaurdConstraintsMatrix20(5, 2) = 1.0;
	gaurdConstraintsMatrix20(5, 3) = 0.0;
	gaurdConstraintsMatrix20(6, 0) = 0.0;
	gaurdConstraintsMatrix20(6, 1) = 0.0;
	gaurdConstraintsMatrix20(6, 2) = 0.0;
	gaurdConstraintsMatrix20(6, 3) = -1.0;
	gaurdConstraintsMatrix20(7, 0) = 0.0;
	gaurdConstraintsMatrix20(7, 1) = 0.0;
	gaurdConstraintsMatrix20(7, 2) = 0.0;
	gaurdConstraintsMatrix20(7, 3) = 1.0;

	gaurdBoundValue20.resize(row);
	gaurdBoundValue20[0] = -2.0;
	gaurdBoundValue20[1] = 3.0;
	gaurdBoundValue20[2] = -1.0;
	gaurdBoundValue20[3] = 1.0;
	gaurdBoundValue20[4] = 1000.0;
	gaurdBoundValue20[5] = 1000.0;
	gaurdBoundValue20[6] = 1000.0;
	gaurdBoundValue20[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope20 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix20, gaurdBoundValue20,
					gaurdBoundSign));

	// The transition label ist13

	// Original guard: x1 = 1 & 0 <= x2 & x2 <= 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix21.resize(row, col);
	gaurdConstraintsMatrix21(0, 0) = -1.0;
	gaurdConstraintsMatrix21(0, 1) = 0.0;
	gaurdConstraintsMatrix21(0, 2) = 0.0;
	gaurdConstraintsMatrix21(0, 3) = 0.0;
	gaurdConstraintsMatrix21(1, 0) = 1.0;
	gaurdConstraintsMatrix21(1, 1) = 0.0;
	gaurdConstraintsMatrix21(1, 2) = 0.0;
	gaurdConstraintsMatrix21(1, 3) = 0.0;
	gaurdConstraintsMatrix21(2, 0) = 0.0;
	gaurdConstraintsMatrix21(2, 1) = -1.0;
	gaurdConstraintsMatrix21(2, 2) = 0.0;
	gaurdConstraintsMatrix21(2, 3) = 0.0;
	gaurdConstraintsMatrix21(3, 0) = 0.0;
	gaurdConstraintsMatrix21(3, 1) = 1.0;
	gaurdConstraintsMatrix21(3, 2) = 0.0;
	gaurdConstraintsMatrix21(3, 3) = 0.0;
	gaurdConstraintsMatrix21(4, 0) = 0.0;
	gaurdConstraintsMatrix21(4, 1) = 0.0;
	gaurdConstraintsMatrix21(4, 2) = -1.0;
	gaurdConstraintsMatrix21(4, 3) = 0.0;
	gaurdConstraintsMatrix21(5, 0) = 0.0;
	gaurdConstraintsMatrix21(5, 1) = 0.0;
	gaurdConstraintsMatrix21(5, 2) = 1.0;
	gaurdConstraintsMatrix21(5, 3) = 0.0;
	gaurdConstraintsMatrix21(6, 0) = 0.0;
	gaurdConstraintsMatrix21(6, 1) = 0.0;
	gaurdConstraintsMatrix21(6, 2) = 0.0;
	gaurdConstraintsMatrix21(6, 3) = -1.0;
	gaurdConstraintsMatrix21(7, 0) = 0.0;
	gaurdConstraintsMatrix21(7, 1) = 0.0;
	gaurdConstraintsMatrix21(7, 2) = 0.0;
	gaurdConstraintsMatrix21(7, 3) = 1.0;

	gaurdBoundValue21.resize(row);
	gaurdBoundValue21[0] = -1.0;
	gaurdBoundValue21[1] = 1.0;
	gaurdBoundValue21[2] = -0.0;
	gaurdBoundValue21[3] = 1.0;
	gaurdBoundValue21[4] = 1000.0;
	gaurdBoundValue21[5] = 1000.0;
	gaurdBoundValue21[6] = 1000.0;
	gaurdBoundValue21[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope21 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix21, gaurdBoundValue21,
					gaurdBoundSign));

	// The transition label ist12

	// Original guard: 0 <= x1 & x1 <= 1 & x2 = 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix22.resize(row, col);
	gaurdConstraintsMatrix22(0, 0) = -1.0;
	gaurdConstraintsMatrix22(0, 1) = 0.0;
	gaurdConstraintsMatrix22(0, 2) = 0.0;
	gaurdConstraintsMatrix22(0, 3) = 0.0;
	gaurdConstraintsMatrix22(1, 0) = 1.0;
	gaurdConstraintsMatrix22(1, 1) = 0.0;
	gaurdConstraintsMatrix22(1, 2) = 0.0;
	gaurdConstraintsMatrix22(1, 3) = 0.0;
	gaurdConstraintsMatrix22(2, 0) = 0.0;
	gaurdConstraintsMatrix22(2, 1) = -1.0;
	gaurdConstraintsMatrix22(2, 2) = 0.0;
	gaurdConstraintsMatrix22(2, 3) = 0.0;
	gaurdConstraintsMatrix22(3, 0) = 0.0;
	gaurdConstraintsMatrix22(3, 1) = 1.0;
	gaurdConstraintsMatrix22(3, 2) = 0.0;
	gaurdConstraintsMatrix22(3, 3) = 0.0;
	gaurdConstraintsMatrix22(4, 0) = 0.0;
	gaurdConstraintsMatrix22(4, 1) = 0.0;
	gaurdConstraintsMatrix22(4, 2) = -1.0;
	gaurdConstraintsMatrix22(4, 3) = 0.0;
	gaurdConstraintsMatrix22(5, 0) = 0.0;
	gaurdConstraintsMatrix22(5, 1) = 0.0;
	gaurdConstraintsMatrix22(5, 2) = 1.0;
	gaurdConstraintsMatrix22(5, 3) = 0.0;
	gaurdConstraintsMatrix22(6, 0) = 0.0;
	gaurdConstraintsMatrix22(6, 1) = 0.0;
	gaurdConstraintsMatrix22(6, 2) = 0.0;
	gaurdConstraintsMatrix22(6, 3) = -1.0;
	gaurdConstraintsMatrix22(7, 0) = 0.0;
	gaurdConstraintsMatrix22(7, 1) = 0.0;
	gaurdConstraintsMatrix22(7, 2) = 0.0;
	gaurdConstraintsMatrix22(7, 3) = 1.0;

	gaurdBoundValue22.resize(row);
	gaurdBoundValue22[0] = -0.0;
	gaurdBoundValue22[1] = 1.0;
	gaurdBoundValue22[2] = -1.0;
	gaurdBoundValue22[3] = 1.0;
	gaurdBoundValue22[4] = 1000.0;
	gaurdBoundValue22[5] = 1000.0;
	gaurdBoundValue22[6] = 1000.0;
	gaurdBoundValue22[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope22 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix22, gaurdBoundValue22,
					gaurdBoundSign));

	// The transition label ist4

	// Original guard: 0 <= x1 & x1 <= 1 & x2 = 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix23.resize(row, col);
	gaurdConstraintsMatrix23(0, 0) = -1.0;
	gaurdConstraintsMatrix23(0, 1) = 0.0;
	gaurdConstraintsMatrix23(0, 2) = 0.0;
	gaurdConstraintsMatrix23(0, 3) = 0.0;
	gaurdConstraintsMatrix23(1, 0) = 1.0;
	gaurdConstraintsMatrix23(1, 1) = 0.0;
	gaurdConstraintsMatrix23(1, 2) = 0.0;
	gaurdConstraintsMatrix23(1, 3) = 0.0;
	gaurdConstraintsMatrix23(2, 0) = 0.0;
	gaurdConstraintsMatrix23(2, 1) = -1.0;
	gaurdConstraintsMatrix23(2, 2) = 0.0;
	gaurdConstraintsMatrix23(2, 3) = 0.0;
	gaurdConstraintsMatrix23(3, 0) = 0.0;
	gaurdConstraintsMatrix23(3, 1) = 1.0;
	gaurdConstraintsMatrix23(3, 2) = 0.0;
	gaurdConstraintsMatrix23(3, 3) = 0.0;
	gaurdConstraintsMatrix23(4, 0) = 0.0;
	gaurdConstraintsMatrix23(4, 1) = 0.0;
	gaurdConstraintsMatrix23(4, 2) = -1.0;
	gaurdConstraintsMatrix23(4, 3) = 0.0;
	gaurdConstraintsMatrix23(5, 0) = 0.0;
	gaurdConstraintsMatrix23(5, 1) = 0.0;
	gaurdConstraintsMatrix23(5, 2) = 1.0;
	gaurdConstraintsMatrix23(5, 3) = 0.0;
	gaurdConstraintsMatrix23(6, 0) = 0.0;
	gaurdConstraintsMatrix23(6, 1) = 0.0;
	gaurdConstraintsMatrix23(6, 2) = 0.0;
	gaurdConstraintsMatrix23(6, 3) = -1.0;
	gaurdConstraintsMatrix23(7, 0) = 0.0;
	gaurdConstraintsMatrix23(7, 1) = 0.0;
	gaurdConstraintsMatrix23(7, 2) = 0.0;
	gaurdConstraintsMatrix23(7, 3) = 1.0;

	gaurdBoundValue23.resize(row);
	gaurdBoundValue23[0] = -0.0;
	gaurdBoundValue23[1] = 1.0;
	gaurdBoundValue23[2] = -3.0;
	gaurdBoundValue23[3] = 3.0;
	gaurdBoundValue23[4] = 1000.0;
	gaurdBoundValue23[5] = 1000.0;
	gaurdBoundValue23[6] = 1000.0;
	gaurdBoundValue23[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope23 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix23, gaurdBoundValue23,
					gaurdBoundSign));

	// The transition label ist3

	// Original guard: 0 <= x1 & x1 <= 1 & x2 = 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix24.resize(row, col);
	gaurdConstraintsMatrix24(0, 0) = -1.0;
	gaurdConstraintsMatrix24(0, 1) = 0.0;
	gaurdConstraintsMatrix24(0, 2) = 0.0;
	gaurdConstraintsMatrix24(0, 3) = 0.0;
	gaurdConstraintsMatrix24(1, 0) = 1.0;
	gaurdConstraintsMatrix24(1, 1) = 0.0;
	gaurdConstraintsMatrix24(1, 2) = 0.0;
	gaurdConstraintsMatrix24(1, 3) = 0.0;
	gaurdConstraintsMatrix24(2, 0) = 0.0;
	gaurdConstraintsMatrix24(2, 1) = -1.0;
	gaurdConstraintsMatrix24(2, 2) = 0.0;
	gaurdConstraintsMatrix24(2, 3) = 0.0;
	gaurdConstraintsMatrix24(3, 0) = 0.0;
	gaurdConstraintsMatrix24(3, 1) = 1.0;
	gaurdConstraintsMatrix24(3, 2) = 0.0;
	gaurdConstraintsMatrix24(3, 3) = 0.0;
	gaurdConstraintsMatrix24(4, 0) = 0.0;
	gaurdConstraintsMatrix24(4, 1) = 0.0;
	gaurdConstraintsMatrix24(4, 2) = -1.0;
	gaurdConstraintsMatrix24(4, 3) = 0.0;
	gaurdConstraintsMatrix24(5, 0) = 0.0;
	gaurdConstraintsMatrix24(5, 1) = 0.0;
	gaurdConstraintsMatrix24(5, 2) = 1.0;
	gaurdConstraintsMatrix24(5, 3) = 0.0;
	gaurdConstraintsMatrix24(6, 0) = 0.0;
	gaurdConstraintsMatrix24(6, 1) = 0.0;
	gaurdConstraintsMatrix24(6, 2) = 0.0;
	gaurdConstraintsMatrix24(6, 3) = -1.0;
	gaurdConstraintsMatrix24(7, 0) = 0.0;
	gaurdConstraintsMatrix24(7, 1) = 0.0;
	gaurdConstraintsMatrix24(7, 2) = 0.0;
	gaurdConstraintsMatrix24(7, 3) = 1.0;

	gaurdBoundValue24.resize(row);
	gaurdBoundValue24[0] = -0.0;
	gaurdBoundValue24[1] = 1.0;
	gaurdBoundValue24[2] = -4.0;
	gaurdBoundValue24[3] = 4.0;
	gaurdBoundValue24[4] = 1000.0;
	gaurdBoundValue24[5] = 1000.0;
	gaurdBoundValue24[6] = 1000.0;
	gaurdBoundValue24[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope24 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix24, gaurdBoundValue24,
					gaurdBoundSign));

	// The transition label ist5

	// Original guard: x1 = 1 & 3 <= x2 & x2 <= 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix25.resize(row, col);
	gaurdConstraintsMatrix25(0, 0) = -1.0;
	gaurdConstraintsMatrix25(0, 1) = 0.0;
	gaurdConstraintsMatrix25(0, 2) = 0.0;
	gaurdConstraintsMatrix25(0, 3) = 0.0;
	gaurdConstraintsMatrix25(1, 0) = 1.0;
	gaurdConstraintsMatrix25(1, 1) = 0.0;
	gaurdConstraintsMatrix25(1, 2) = 0.0;
	gaurdConstraintsMatrix25(1, 3) = 0.0;
	gaurdConstraintsMatrix25(2, 0) = 0.0;
	gaurdConstraintsMatrix25(2, 1) = -1.0;
	gaurdConstraintsMatrix25(2, 2) = 0.0;
	gaurdConstraintsMatrix25(2, 3) = 0.0;
	gaurdConstraintsMatrix25(3, 0) = 0.0;
	gaurdConstraintsMatrix25(3, 1) = 1.0;
	gaurdConstraintsMatrix25(3, 2) = 0.0;
	gaurdConstraintsMatrix25(3, 3) = 0.0;
	gaurdConstraintsMatrix25(4, 0) = 0.0;
	gaurdConstraintsMatrix25(4, 1) = 0.0;
	gaurdConstraintsMatrix25(4, 2) = -1.0;
	gaurdConstraintsMatrix25(4, 3) = 0.0;
	gaurdConstraintsMatrix25(5, 0) = 0.0;
	gaurdConstraintsMatrix25(5, 1) = 0.0;
	gaurdConstraintsMatrix25(5, 2) = 1.0;
	gaurdConstraintsMatrix25(5, 3) = 0.0;
	gaurdConstraintsMatrix25(6, 0) = 0.0;
	gaurdConstraintsMatrix25(6, 1) = 0.0;
	gaurdConstraintsMatrix25(6, 2) = 0.0;
	gaurdConstraintsMatrix25(6, 3) = -1.0;
	gaurdConstraintsMatrix25(7, 0) = 0.0;
	gaurdConstraintsMatrix25(7, 1) = 0.0;
	gaurdConstraintsMatrix25(7, 2) = 0.0;
	gaurdConstraintsMatrix25(7, 3) = 1.0;

	gaurdBoundValue25.resize(row);
	gaurdBoundValue25[0] = -1.0;
	gaurdBoundValue25[1] = 1.0;
	gaurdBoundValue25[2] = -3.0;
	gaurdBoundValue25[3] = 4.0;
	gaurdBoundValue25[4] = 1000.0;
	gaurdBoundValue25[5] = 1000.0;
	gaurdBoundValue25[6] = 1000.0;
	gaurdBoundValue25[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope25 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix25, gaurdBoundValue25,
					gaurdBoundSign));

	// The transition label ist1

	// Original guard: 0 <= x1 & x1 <= 1 & x2 = 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix26.resize(row, col);
	gaurdConstraintsMatrix26(0, 0) = -1.0;
	gaurdConstraintsMatrix26(0, 1) = 0.0;
	gaurdConstraintsMatrix26(0, 2) = 0.0;
	gaurdConstraintsMatrix26(0, 3) = 0.0;
	gaurdConstraintsMatrix26(1, 0) = 1.0;
	gaurdConstraintsMatrix26(1, 1) = 0.0;
	gaurdConstraintsMatrix26(1, 2) = 0.0;
	gaurdConstraintsMatrix26(1, 3) = 0.0;
	gaurdConstraintsMatrix26(2, 0) = 0.0;
	gaurdConstraintsMatrix26(2, 1) = -1.0;
	gaurdConstraintsMatrix26(2, 2) = 0.0;
	gaurdConstraintsMatrix26(2, 3) = 0.0;
	gaurdConstraintsMatrix26(3, 0) = 0.0;
	gaurdConstraintsMatrix26(3, 1) = 1.0;
	gaurdConstraintsMatrix26(3, 2) = 0.0;
	gaurdConstraintsMatrix26(3, 3) = 0.0;
	gaurdConstraintsMatrix26(4, 0) = 0.0;
	gaurdConstraintsMatrix26(4, 1) = 0.0;
	gaurdConstraintsMatrix26(4, 2) = -1.0;
	gaurdConstraintsMatrix26(4, 3) = 0.0;
	gaurdConstraintsMatrix26(5, 0) = 0.0;
	gaurdConstraintsMatrix26(5, 1) = 0.0;
	gaurdConstraintsMatrix26(5, 2) = 1.0;
	gaurdConstraintsMatrix26(5, 3) = 0.0;
	gaurdConstraintsMatrix26(6, 0) = 0.0;
	gaurdConstraintsMatrix26(6, 1) = 0.0;
	gaurdConstraintsMatrix26(6, 2) = 0.0;
	gaurdConstraintsMatrix26(6, 3) = -1.0;
	gaurdConstraintsMatrix26(7, 0) = 0.0;
	gaurdConstraintsMatrix26(7, 1) = 0.0;
	gaurdConstraintsMatrix26(7, 2) = 0.0;
	gaurdConstraintsMatrix26(7, 3) = 1.0;

	gaurdBoundValue26.resize(row);
	gaurdBoundValue26[0] = -0.0;
	gaurdBoundValue26[1] = 1.0;
	gaurdBoundValue26[2] = -4.0;
	gaurdBoundValue26[3] = 4.0;
	gaurdBoundValue26[4] = 1000.0;
	gaurdBoundValue26[5] = 1000.0;
	gaurdBoundValue26[6] = 1000.0;
	gaurdBoundValue26[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope26 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix26, gaurdBoundValue26,
					gaurdBoundSign));

	// The transition label ist2

	// Original guard: x1 = 1 & 4 <= x2 & x2 <= 5 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix27.resize(row, col);
	gaurdConstraintsMatrix27(0, 0) = -1.0;
	gaurdConstraintsMatrix27(0, 1) = 0.0;
	gaurdConstraintsMatrix27(0, 2) = 0.0;
	gaurdConstraintsMatrix27(0, 3) = 0.0;
	gaurdConstraintsMatrix27(1, 0) = 1.0;
	gaurdConstraintsMatrix27(1, 1) = 0.0;
	gaurdConstraintsMatrix27(1, 2) = 0.0;
	gaurdConstraintsMatrix27(1, 3) = 0.0;
	gaurdConstraintsMatrix27(2, 0) = 0.0;
	gaurdConstraintsMatrix27(2, 1) = -1.0;
	gaurdConstraintsMatrix27(2, 2) = 0.0;
	gaurdConstraintsMatrix27(2, 3) = 0.0;
	gaurdConstraintsMatrix27(3, 0) = 0.0;
	gaurdConstraintsMatrix27(3, 1) = 1.0;
	gaurdConstraintsMatrix27(3, 2) = 0.0;
	gaurdConstraintsMatrix27(3, 3) = 0.0;
	gaurdConstraintsMatrix27(4, 0) = 0.0;
	gaurdConstraintsMatrix27(4, 1) = 0.0;
	gaurdConstraintsMatrix27(4, 2) = -1.0;
	gaurdConstraintsMatrix27(4, 3) = 0.0;
	gaurdConstraintsMatrix27(5, 0) = 0.0;
	gaurdConstraintsMatrix27(5, 1) = 0.0;
	gaurdConstraintsMatrix27(5, 2) = 1.0;
	gaurdConstraintsMatrix27(5, 3) = 0.0;
	gaurdConstraintsMatrix27(6, 0) = 0.0;
	gaurdConstraintsMatrix27(6, 1) = 0.0;
	gaurdConstraintsMatrix27(6, 2) = 0.0;
	gaurdConstraintsMatrix27(6, 3) = -1.0;
	gaurdConstraintsMatrix27(7, 0) = 0.0;
	gaurdConstraintsMatrix27(7, 1) = 0.0;
	gaurdConstraintsMatrix27(7, 2) = 0.0;
	gaurdConstraintsMatrix27(7, 3) = 1.0;

	gaurdBoundValue27.resize(row);
	gaurdBoundValue27[0] = -1.0;
	gaurdBoundValue27[1] = 1.0;
	gaurdBoundValue27[2] = -4.0;
	gaurdBoundValue27[3] = 5.0;
	gaurdBoundValue27[4] = 1000.0;
	gaurdBoundValue27[5] = 1000.0;
	gaurdBoundValue27[6] = 1000.0;
	gaurdBoundValue27[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope27 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix27, gaurdBoundValue27,
					gaurdBoundSign));

	// The transition label ist26

	// Original guard: x1 = 1 & 4 <= x2 & x2 <= 5 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix28.resize(row, col);
	gaurdConstraintsMatrix28(0, 0) = -1.0;
	gaurdConstraintsMatrix28(0, 1) = 0.0;
	gaurdConstraintsMatrix28(0, 2) = 0.0;
	gaurdConstraintsMatrix28(0, 3) = 0.0;
	gaurdConstraintsMatrix28(1, 0) = 1.0;
	gaurdConstraintsMatrix28(1, 1) = 0.0;
	gaurdConstraintsMatrix28(1, 2) = 0.0;
	gaurdConstraintsMatrix28(1, 3) = 0.0;
	gaurdConstraintsMatrix28(2, 0) = 0.0;
	gaurdConstraintsMatrix28(2, 1) = -1.0;
	gaurdConstraintsMatrix28(2, 2) = 0.0;
	gaurdConstraintsMatrix28(2, 3) = 0.0;
	gaurdConstraintsMatrix28(3, 0) = 0.0;
	gaurdConstraintsMatrix28(3, 1) = 1.0;
	gaurdConstraintsMatrix28(3, 2) = 0.0;
	gaurdConstraintsMatrix28(3, 3) = 0.0;
	gaurdConstraintsMatrix28(4, 0) = 0.0;
	gaurdConstraintsMatrix28(4, 1) = 0.0;
	gaurdConstraintsMatrix28(4, 2) = -1.0;
	gaurdConstraintsMatrix28(4, 3) = 0.0;
	gaurdConstraintsMatrix28(5, 0) = 0.0;
	gaurdConstraintsMatrix28(5, 1) = 0.0;
	gaurdConstraintsMatrix28(5, 2) = 1.0;
	gaurdConstraintsMatrix28(5, 3) = 0.0;
	gaurdConstraintsMatrix28(6, 0) = 0.0;
	gaurdConstraintsMatrix28(6, 1) = 0.0;
	gaurdConstraintsMatrix28(6, 2) = 0.0;
	gaurdConstraintsMatrix28(6, 3) = -1.0;
	gaurdConstraintsMatrix28(7, 0) = 0.0;
	gaurdConstraintsMatrix28(7, 1) = 0.0;
	gaurdConstraintsMatrix28(7, 2) = 0.0;
	gaurdConstraintsMatrix28(7, 3) = 1.0;

	gaurdBoundValue28.resize(row);
	gaurdBoundValue28[0] = -1.0;
	gaurdBoundValue28[1] = 1.0;
	gaurdBoundValue28[2] = -4.0;
	gaurdBoundValue28[3] = 5.0;
	gaurdBoundValue28[4] = 1000.0;
	gaurdBoundValue28[5] = 1000.0;
	gaurdBoundValue28[6] = 1000.0;
	gaurdBoundValue28[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope28 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix28, gaurdBoundValue28,
					gaurdBoundSign));

	// The transition label ist28

	// Original guard: x1 = 2 & 4 <= x2 & x2 <= 5 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix29.resize(row, col);
	gaurdConstraintsMatrix29(0, 0) = -1.0;
	gaurdConstraintsMatrix29(0, 1) = 0.0;
	gaurdConstraintsMatrix29(0, 2) = 0.0;
	gaurdConstraintsMatrix29(0, 3) = 0.0;
	gaurdConstraintsMatrix29(1, 0) = 1.0;
	gaurdConstraintsMatrix29(1, 1) = 0.0;
	gaurdConstraintsMatrix29(1, 2) = 0.0;
	gaurdConstraintsMatrix29(1, 3) = 0.0;
	gaurdConstraintsMatrix29(2, 0) = 0.0;
	gaurdConstraintsMatrix29(2, 1) = -1.0;
	gaurdConstraintsMatrix29(2, 2) = 0.0;
	gaurdConstraintsMatrix29(2, 3) = 0.0;
	gaurdConstraintsMatrix29(3, 0) = 0.0;
	gaurdConstraintsMatrix29(3, 1) = 1.0;
	gaurdConstraintsMatrix29(3, 2) = 0.0;
	gaurdConstraintsMatrix29(3, 3) = 0.0;
	gaurdConstraintsMatrix29(4, 0) = 0.0;
	gaurdConstraintsMatrix29(4, 1) = 0.0;
	gaurdConstraintsMatrix29(4, 2) = -1.0;
	gaurdConstraintsMatrix29(4, 3) = 0.0;
	gaurdConstraintsMatrix29(5, 0) = 0.0;
	gaurdConstraintsMatrix29(5, 1) = 0.0;
	gaurdConstraintsMatrix29(5, 2) = 1.0;
	gaurdConstraintsMatrix29(5, 3) = 0.0;
	gaurdConstraintsMatrix29(6, 0) = 0.0;
	gaurdConstraintsMatrix29(6, 1) = 0.0;
	gaurdConstraintsMatrix29(6, 2) = 0.0;
	gaurdConstraintsMatrix29(6, 3) = -1.0;
	gaurdConstraintsMatrix29(7, 0) = 0.0;
	gaurdConstraintsMatrix29(7, 1) = 0.0;
	gaurdConstraintsMatrix29(7, 2) = 0.0;
	gaurdConstraintsMatrix29(7, 3) = 1.0;

	gaurdBoundValue29.resize(row);
	gaurdBoundValue29[0] = -2.0;
	gaurdBoundValue29[1] = 2.0;
	gaurdBoundValue29[2] = -4.0;
	gaurdBoundValue29[3] = 5.0;
	gaurdBoundValue29[4] = 1000.0;
	gaurdBoundValue29[5] = 1000.0;
	gaurdBoundValue29[6] = 1000.0;
	gaurdBoundValue29[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope29 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix29, gaurdBoundValue29,
					gaurdBoundSign));

	// The transition label ist27

	// Original guard: 1 <= x1 & x1 <= 2 & x2 = 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix30.resize(row, col);
	gaurdConstraintsMatrix30(0, 0) = -1.0;
	gaurdConstraintsMatrix30(0, 1) = 0.0;
	gaurdConstraintsMatrix30(0, 2) = 0.0;
	gaurdConstraintsMatrix30(0, 3) = 0.0;
	gaurdConstraintsMatrix30(1, 0) = 1.0;
	gaurdConstraintsMatrix30(1, 1) = 0.0;
	gaurdConstraintsMatrix30(1, 2) = 0.0;
	gaurdConstraintsMatrix30(1, 3) = 0.0;
	gaurdConstraintsMatrix30(2, 0) = 0.0;
	gaurdConstraintsMatrix30(2, 1) = -1.0;
	gaurdConstraintsMatrix30(2, 2) = 0.0;
	gaurdConstraintsMatrix30(2, 3) = 0.0;
	gaurdConstraintsMatrix30(3, 0) = 0.0;
	gaurdConstraintsMatrix30(3, 1) = 1.0;
	gaurdConstraintsMatrix30(3, 2) = 0.0;
	gaurdConstraintsMatrix30(3, 3) = 0.0;
	gaurdConstraintsMatrix30(4, 0) = 0.0;
	gaurdConstraintsMatrix30(4, 1) = 0.0;
	gaurdConstraintsMatrix30(4, 2) = -1.0;
	gaurdConstraintsMatrix30(4, 3) = 0.0;
	gaurdConstraintsMatrix30(5, 0) = 0.0;
	gaurdConstraintsMatrix30(5, 1) = 0.0;
	gaurdConstraintsMatrix30(5, 2) = 1.0;
	gaurdConstraintsMatrix30(5, 3) = 0.0;
	gaurdConstraintsMatrix30(6, 0) = 0.0;
	gaurdConstraintsMatrix30(6, 1) = 0.0;
	gaurdConstraintsMatrix30(6, 2) = 0.0;
	gaurdConstraintsMatrix30(6, 3) = -1.0;
	gaurdConstraintsMatrix30(7, 0) = 0.0;
	gaurdConstraintsMatrix30(7, 1) = 0.0;
	gaurdConstraintsMatrix30(7, 2) = 0.0;
	gaurdConstraintsMatrix30(7, 3) = 1.0;

	gaurdBoundValue30.resize(row);
	gaurdBoundValue30[0] = -1.0;
	gaurdBoundValue30[1] = 2.0;
	gaurdBoundValue30[2] = -4.0;
	gaurdBoundValue30[3] = 4.0;
	gaurdBoundValue30[4] = 1000.0;
	gaurdBoundValue30[5] = 1000.0;
	gaurdBoundValue30[6] = 1000.0;
	gaurdBoundValue30[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope30 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix30, gaurdBoundValue30,
					gaurdBoundSign));

	// The transition label ist23

	// Original guard: x1 = 1 & 3 <= x2 & x2 <= 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix31.resize(row, col);
	gaurdConstraintsMatrix31(0, 0) = -1.0;
	gaurdConstraintsMatrix31(0, 1) = 0.0;
	gaurdConstraintsMatrix31(0, 2) = 0.0;
	gaurdConstraintsMatrix31(0, 3) = 0.0;
	gaurdConstraintsMatrix31(1, 0) = 1.0;
	gaurdConstraintsMatrix31(1, 1) = 0.0;
	gaurdConstraintsMatrix31(1, 2) = 0.0;
	gaurdConstraintsMatrix31(1, 3) = 0.0;
	gaurdConstraintsMatrix31(2, 0) = 0.0;
	gaurdConstraintsMatrix31(2, 1) = -1.0;
	gaurdConstraintsMatrix31(2, 2) = 0.0;
	gaurdConstraintsMatrix31(2, 3) = 0.0;
	gaurdConstraintsMatrix31(3, 0) = 0.0;
	gaurdConstraintsMatrix31(3, 1) = 1.0;
	gaurdConstraintsMatrix31(3, 2) = 0.0;
	gaurdConstraintsMatrix31(3, 3) = 0.0;
	gaurdConstraintsMatrix31(4, 0) = 0.0;
	gaurdConstraintsMatrix31(4, 1) = 0.0;
	gaurdConstraintsMatrix31(4, 2) = -1.0;
	gaurdConstraintsMatrix31(4, 3) = 0.0;
	gaurdConstraintsMatrix31(5, 0) = 0.0;
	gaurdConstraintsMatrix31(5, 1) = 0.0;
	gaurdConstraintsMatrix31(5, 2) = 1.0;
	gaurdConstraintsMatrix31(5, 3) = 0.0;
	gaurdConstraintsMatrix31(6, 0) = 0.0;
	gaurdConstraintsMatrix31(6, 1) = 0.0;
	gaurdConstraintsMatrix31(6, 2) = 0.0;
	gaurdConstraintsMatrix31(6, 3) = -1.0;
	gaurdConstraintsMatrix31(7, 0) = 0.0;
	gaurdConstraintsMatrix31(7, 1) = 0.0;
	gaurdConstraintsMatrix31(7, 2) = 0.0;
	gaurdConstraintsMatrix31(7, 3) = 1.0;

	gaurdBoundValue31.resize(row);
	gaurdBoundValue31[0] = -1.0;
	gaurdBoundValue31[1] = 1.0;
	gaurdBoundValue31[2] = -3.0;
	gaurdBoundValue31[3] = 4.0;
	gaurdBoundValue31[4] = 1000.0;
	gaurdBoundValue31[5] = 1000.0;
	gaurdBoundValue31[6] = 1000.0;
	gaurdBoundValue31[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope31 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix31, gaurdBoundValue31,
					gaurdBoundSign));

	// The transition label ist22

	// Original guard: 1 <= x1 & x1 <= 2 & x2 = 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix32.resize(row, col);
	gaurdConstraintsMatrix32(0, 0) = -1.0;
	gaurdConstraintsMatrix32(0, 1) = 0.0;
	gaurdConstraintsMatrix32(0, 2) = 0.0;
	gaurdConstraintsMatrix32(0, 3) = 0.0;
	gaurdConstraintsMatrix32(1, 0) = 1.0;
	gaurdConstraintsMatrix32(1, 1) = 0.0;
	gaurdConstraintsMatrix32(1, 2) = 0.0;
	gaurdConstraintsMatrix32(1, 3) = 0.0;
	gaurdConstraintsMatrix32(2, 0) = 0.0;
	gaurdConstraintsMatrix32(2, 1) = -1.0;
	gaurdConstraintsMatrix32(2, 2) = 0.0;
	gaurdConstraintsMatrix32(2, 3) = 0.0;
	gaurdConstraintsMatrix32(3, 0) = 0.0;
	gaurdConstraintsMatrix32(3, 1) = 1.0;
	gaurdConstraintsMatrix32(3, 2) = 0.0;
	gaurdConstraintsMatrix32(3, 3) = 0.0;
	gaurdConstraintsMatrix32(4, 0) = 0.0;
	gaurdConstraintsMatrix32(4, 1) = 0.0;
	gaurdConstraintsMatrix32(4, 2) = -1.0;
	gaurdConstraintsMatrix32(4, 3) = 0.0;
	gaurdConstraintsMatrix32(5, 0) = 0.0;
	gaurdConstraintsMatrix32(5, 1) = 0.0;
	gaurdConstraintsMatrix32(5, 2) = 1.0;
	gaurdConstraintsMatrix32(5, 3) = 0.0;
	gaurdConstraintsMatrix32(6, 0) = 0.0;
	gaurdConstraintsMatrix32(6, 1) = 0.0;
	gaurdConstraintsMatrix32(6, 2) = 0.0;
	gaurdConstraintsMatrix32(6, 3) = -1.0;
	gaurdConstraintsMatrix32(7, 0) = 0.0;
	gaurdConstraintsMatrix32(7, 1) = 0.0;
	gaurdConstraintsMatrix32(7, 2) = 0.0;
	gaurdConstraintsMatrix32(7, 3) = 1.0;

	gaurdBoundValue32.resize(row);
	gaurdBoundValue32[0] = -1.0;
	gaurdBoundValue32[1] = 2.0;
	gaurdBoundValue32[2] = -4.0;
	gaurdBoundValue32[3] = 4.0;
	gaurdBoundValue32[4] = 1000.0;
	gaurdBoundValue32[5] = 1000.0;
	gaurdBoundValue32[6] = 1000.0;
	gaurdBoundValue32[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope32 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix32, gaurdBoundValue32,
					gaurdBoundSign));

	// The transition label ist25

	// Original guard: x1 = 2 & 3 <= x2 & x2 <= 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix33.resize(row, col);
	gaurdConstraintsMatrix33(0, 0) = -1.0;
	gaurdConstraintsMatrix33(0, 1) = 0.0;
	gaurdConstraintsMatrix33(0, 2) = 0.0;
	gaurdConstraintsMatrix33(0, 3) = 0.0;
	gaurdConstraintsMatrix33(1, 0) = 1.0;
	gaurdConstraintsMatrix33(1, 1) = 0.0;
	gaurdConstraintsMatrix33(1, 2) = 0.0;
	gaurdConstraintsMatrix33(1, 3) = 0.0;
	gaurdConstraintsMatrix33(2, 0) = 0.0;
	gaurdConstraintsMatrix33(2, 1) = -1.0;
	gaurdConstraintsMatrix33(2, 2) = 0.0;
	gaurdConstraintsMatrix33(2, 3) = 0.0;
	gaurdConstraintsMatrix33(3, 0) = 0.0;
	gaurdConstraintsMatrix33(3, 1) = 1.0;
	gaurdConstraintsMatrix33(3, 2) = 0.0;
	gaurdConstraintsMatrix33(3, 3) = 0.0;
	gaurdConstraintsMatrix33(4, 0) = 0.0;
	gaurdConstraintsMatrix33(4, 1) = 0.0;
	gaurdConstraintsMatrix33(4, 2) = -1.0;
	gaurdConstraintsMatrix33(4, 3) = 0.0;
	gaurdConstraintsMatrix33(5, 0) = 0.0;
	gaurdConstraintsMatrix33(5, 1) = 0.0;
	gaurdConstraintsMatrix33(5, 2) = 1.0;
	gaurdConstraintsMatrix33(5, 3) = 0.0;
	gaurdConstraintsMatrix33(6, 0) = 0.0;
	gaurdConstraintsMatrix33(6, 1) = 0.0;
	gaurdConstraintsMatrix33(6, 2) = 0.0;
	gaurdConstraintsMatrix33(6, 3) = -1.0;
	gaurdConstraintsMatrix33(7, 0) = 0.0;
	gaurdConstraintsMatrix33(7, 1) = 0.0;
	gaurdConstraintsMatrix33(7, 2) = 0.0;
	gaurdConstraintsMatrix33(7, 3) = 1.0;

	gaurdBoundValue33.resize(row);
	gaurdBoundValue33[0] = -2.0;
	gaurdBoundValue33[1] = 2.0;
	gaurdBoundValue33[2] = -3.0;
	gaurdBoundValue33[3] = 4.0;
	gaurdBoundValue33[4] = 1000.0;
	gaurdBoundValue33[5] = 1000.0;
	gaurdBoundValue33[6] = 1000.0;
	gaurdBoundValue33[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope33 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix33, gaurdBoundValue33,
					gaurdBoundSign));

	// The transition label ist24

	// Original guard: 1 <= x1 & x1 <= 2 & x2 = 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix34.resize(row, col);
	gaurdConstraintsMatrix34(0, 0) = -1.0;
	gaurdConstraintsMatrix34(0, 1) = 0.0;
	gaurdConstraintsMatrix34(0, 2) = 0.0;
	gaurdConstraintsMatrix34(0, 3) = 0.0;
	gaurdConstraintsMatrix34(1, 0) = 1.0;
	gaurdConstraintsMatrix34(1, 1) = 0.0;
	gaurdConstraintsMatrix34(1, 2) = 0.0;
	gaurdConstraintsMatrix34(1, 3) = 0.0;
	gaurdConstraintsMatrix34(2, 0) = 0.0;
	gaurdConstraintsMatrix34(2, 1) = -1.0;
	gaurdConstraintsMatrix34(2, 2) = 0.0;
	gaurdConstraintsMatrix34(2, 3) = 0.0;
	gaurdConstraintsMatrix34(3, 0) = 0.0;
	gaurdConstraintsMatrix34(3, 1) = 1.0;
	gaurdConstraintsMatrix34(3, 2) = 0.0;
	gaurdConstraintsMatrix34(3, 3) = 0.0;
	gaurdConstraintsMatrix34(4, 0) = 0.0;
	gaurdConstraintsMatrix34(4, 1) = 0.0;
	gaurdConstraintsMatrix34(4, 2) = -1.0;
	gaurdConstraintsMatrix34(4, 3) = 0.0;
	gaurdConstraintsMatrix34(5, 0) = 0.0;
	gaurdConstraintsMatrix34(5, 1) = 0.0;
	gaurdConstraintsMatrix34(5, 2) = 1.0;
	gaurdConstraintsMatrix34(5, 3) = 0.0;
	gaurdConstraintsMatrix34(6, 0) = 0.0;
	gaurdConstraintsMatrix34(6, 1) = 0.0;
	gaurdConstraintsMatrix34(6, 2) = 0.0;
	gaurdConstraintsMatrix34(6, 3) = -1.0;
	gaurdConstraintsMatrix34(7, 0) = 0.0;
	gaurdConstraintsMatrix34(7, 1) = 0.0;
	gaurdConstraintsMatrix34(7, 2) = 0.0;
	gaurdConstraintsMatrix34(7, 3) = 1.0;

	gaurdBoundValue34.resize(row);
	gaurdBoundValue34[0] = -1.0;
	gaurdBoundValue34[1] = 2.0;
	gaurdBoundValue34[2] = -3.0;
	gaurdBoundValue34[3] = 3.0;
	gaurdBoundValue34[4] = 1000.0;
	gaurdBoundValue34[5] = 1000.0;
	gaurdBoundValue34[6] = 1000.0;
	gaurdBoundValue34[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope34 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix34, gaurdBoundValue34,
					gaurdBoundSign));

	// The transition label ist40

	// Original guard: x1 = 2 & 4 <= x2 & x2 <= 5 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix35.resize(row, col);
	gaurdConstraintsMatrix35(0, 0) = -1.0;
	gaurdConstraintsMatrix35(0, 1) = 0.0;
	gaurdConstraintsMatrix35(0, 2) = 0.0;
	gaurdConstraintsMatrix35(0, 3) = 0.0;
	gaurdConstraintsMatrix35(1, 0) = 1.0;
	gaurdConstraintsMatrix35(1, 1) = 0.0;
	gaurdConstraintsMatrix35(1, 2) = 0.0;
	gaurdConstraintsMatrix35(1, 3) = 0.0;
	gaurdConstraintsMatrix35(2, 0) = 0.0;
	gaurdConstraintsMatrix35(2, 1) = -1.0;
	gaurdConstraintsMatrix35(2, 2) = 0.0;
	gaurdConstraintsMatrix35(2, 3) = 0.0;
	gaurdConstraintsMatrix35(3, 0) = 0.0;
	gaurdConstraintsMatrix35(3, 1) = 1.0;
	gaurdConstraintsMatrix35(3, 2) = 0.0;
	gaurdConstraintsMatrix35(3, 3) = 0.0;
	gaurdConstraintsMatrix35(4, 0) = 0.0;
	gaurdConstraintsMatrix35(4, 1) = 0.0;
	gaurdConstraintsMatrix35(4, 2) = -1.0;
	gaurdConstraintsMatrix35(4, 3) = 0.0;
	gaurdConstraintsMatrix35(5, 0) = 0.0;
	gaurdConstraintsMatrix35(5, 1) = 0.0;
	gaurdConstraintsMatrix35(5, 2) = 1.0;
	gaurdConstraintsMatrix35(5, 3) = 0.0;
	gaurdConstraintsMatrix35(6, 0) = 0.0;
	gaurdConstraintsMatrix35(6, 1) = 0.0;
	gaurdConstraintsMatrix35(6, 2) = 0.0;
	gaurdConstraintsMatrix35(6, 3) = -1.0;
	gaurdConstraintsMatrix35(7, 0) = 0.0;
	gaurdConstraintsMatrix35(7, 1) = 0.0;
	gaurdConstraintsMatrix35(7, 2) = 0.0;
	gaurdConstraintsMatrix35(7, 3) = 1.0;

	gaurdBoundValue35.resize(row);
	gaurdBoundValue35[0] = -2.0;
	gaurdBoundValue35[1] = 2.0;
	gaurdBoundValue35[2] = -4.0;
	gaurdBoundValue35[3] = 5.0;
	gaurdBoundValue35[4] = 1000.0;
	gaurdBoundValue35[5] = 1000.0;
	gaurdBoundValue35[6] = 1000.0;
	gaurdBoundValue35[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope35 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix35, gaurdBoundValue35,
					gaurdBoundSign));

	// The transition label ist42

	// Original guard: x1 = 3 & 4 <= x2 & x2 <= 5 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix36.resize(row, col);
	gaurdConstraintsMatrix36(0, 0) = -1.0;
	gaurdConstraintsMatrix36(0, 1) = 0.0;
	gaurdConstraintsMatrix36(0, 2) = 0.0;
	gaurdConstraintsMatrix36(0, 3) = 0.0;
	gaurdConstraintsMatrix36(1, 0) = 1.0;
	gaurdConstraintsMatrix36(1, 1) = 0.0;
	gaurdConstraintsMatrix36(1, 2) = 0.0;
	gaurdConstraintsMatrix36(1, 3) = 0.0;
	gaurdConstraintsMatrix36(2, 0) = 0.0;
	gaurdConstraintsMatrix36(2, 1) = -1.0;
	gaurdConstraintsMatrix36(2, 2) = 0.0;
	gaurdConstraintsMatrix36(2, 3) = 0.0;
	gaurdConstraintsMatrix36(3, 0) = 0.0;
	gaurdConstraintsMatrix36(3, 1) = 1.0;
	gaurdConstraintsMatrix36(3, 2) = 0.0;
	gaurdConstraintsMatrix36(3, 3) = 0.0;
	gaurdConstraintsMatrix36(4, 0) = 0.0;
	gaurdConstraintsMatrix36(4, 1) = 0.0;
	gaurdConstraintsMatrix36(4, 2) = -1.0;
	gaurdConstraintsMatrix36(4, 3) = 0.0;
	gaurdConstraintsMatrix36(5, 0) = 0.0;
	gaurdConstraintsMatrix36(5, 1) = 0.0;
	gaurdConstraintsMatrix36(5, 2) = 1.0;
	gaurdConstraintsMatrix36(5, 3) = 0.0;
	gaurdConstraintsMatrix36(6, 0) = 0.0;
	gaurdConstraintsMatrix36(6, 1) = 0.0;
	gaurdConstraintsMatrix36(6, 2) = 0.0;
	gaurdConstraintsMatrix36(6, 3) = -1.0;
	gaurdConstraintsMatrix36(7, 0) = 0.0;
	gaurdConstraintsMatrix36(7, 1) = 0.0;
	gaurdConstraintsMatrix36(7, 2) = 0.0;
	gaurdConstraintsMatrix36(7, 3) = 1.0;

	gaurdBoundValue36.resize(row);
	gaurdBoundValue36[0] = -3.0;
	gaurdBoundValue36[1] = 3.0;
	gaurdBoundValue36[2] = -4.0;
	gaurdBoundValue36[3] = 5.0;
	gaurdBoundValue36[4] = 1000.0;
	gaurdBoundValue36[5] = 1000.0;
	gaurdBoundValue36[6] = 1000.0;
	gaurdBoundValue36[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope36 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix36, gaurdBoundValue36,
					gaurdBoundSign));

	// The transition label ist41

	// Original guard: 2 <= x1 & x1 <= 3 & x2 = 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix37.resize(row, col);
	gaurdConstraintsMatrix37(0, 0) = -1.0;
	gaurdConstraintsMatrix37(0, 1) = 0.0;
	gaurdConstraintsMatrix37(0, 2) = 0.0;
	gaurdConstraintsMatrix37(0, 3) = 0.0;
	gaurdConstraintsMatrix37(1, 0) = 1.0;
	gaurdConstraintsMatrix37(1, 1) = 0.0;
	gaurdConstraintsMatrix37(1, 2) = 0.0;
	gaurdConstraintsMatrix37(1, 3) = 0.0;
	gaurdConstraintsMatrix37(2, 0) = 0.0;
	gaurdConstraintsMatrix37(2, 1) = -1.0;
	gaurdConstraintsMatrix37(2, 2) = 0.0;
	gaurdConstraintsMatrix37(2, 3) = 0.0;
	gaurdConstraintsMatrix37(3, 0) = 0.0;
	gaurdConstraintsMatrix37(3, 1) = 1.0;
	gaurdConstraintsMatrix37(3, 2) = 0.0;
	gaurdConstraintsMatrix37(3, 3) = 0.0;
	gaurdConstraintsMatrix37(4, 0) = 0.0;
	gaurdConstraintsMatrix37(4, 1) = 0.0;
	gaurdConstraintsMatrix37(4, 2) = -1.0;
	gaurdConstraintsMatrix37(4, 3) = 0.0;
	gaurdConstraintsMatrix37(5, 0) = 0.0;
	gaurdConstraintsMatrix37(5, 1) = 0.0;
	gaurdConstraintsMatrix37(5, 2) = 1.0;
	gaurdConstraintsMatrix37(5, 3) = 0.0;
	gaurdConstraintsMatrix37(6, 0) = 0.0;
	gaurdConstraintsMatrix37(6, 1) = 0.0;
	gaurdConstraintsMatrix37(6, 2) = 0.0;
	gaurdConstraintsMatrix37(6, 3) = -1.0;
	gaurdConstraintsMatrix37(7, 0) = 0.0;
	gaurdConstraintsMatrix37(7, 1) = 0.0;
	gaurdConstraintsMatrix37(7, 2) = 0.0;
	gaurdConstraintsMatrix37(7, 3) = 1.0;

	gaurdBoundValue37.resize(row);
	gaurdBoundValue37[0] = -2.0;
	gaurdBoundValue37[1] = 3.0;
	gaurdBoundValue37[2] = -4.0;
	gaurdBoundValue37[3] = 4.0;
	gaurdBoundValue37[4] = 1000.0;
	gaurdBoundValue37[5] = 1000.0;
	gaurdBoundValue37[6] = 1000.0;
	gaurdBoundValue37[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope37 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix37, gaurdBoundValue37,
					gaurdBoundSign));

	// The transition label ist36

	// Original guard: 2 <= x1 & x1 <= 3 & x2 = 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix38.resize(row, col);
	gaurdConstraintsMatrix38(0, 0) = -1.0;
	gaurdConstraintsMatrix38(0, 1) = 0.0;
	gaurdConstraintsMatrix38(0, 2) = 0.0;
	gaurdConstraintsMatrix38(0, 3) = 0.0;
	gaurdConstraintsMatrix38(1, 0) = 1.0;
	gaurdConstraintsMatrix38(1, 1) = 0.0;
	gaurdConstraintsMatrix38(1, 2) = 0.0;
	gaurdConstraintsMatrix38(1, 3) = 0.0;
	gaurdConstraintsMatrix38(2, 0) = 0.0;
	gaurdConstraintsMatrix38(2, 1) = -1.0;
	gaurdConstraintsMatrix38(2, 2) = 0.0;
	gaurdConstraintsMatrix38(2, 3) = 0.0;
	gaurdConstraintsMatrix38(3, 0) = 0.0;
	gaurdConstraintsMatrix38(3, 1) = 1.0;
	gaurdConstraintsMatrix38(3, 2) = 0.0;
	gaurdConstraintsMatrix38(3, 3) = 0.0;
	gaurdConstraintsMatrix38(4, 0) = 0.0;
	gaurdConstraintsMatrix38(4, 1) = 0.0;
	gaurdConstraintsMatrix38(4, 2) = -1.0;
	gaurdConstraintsMatrix38(4, 3) = 0.0;
	gaurdConstraintsMatrix38(5, 0) = 0.0;
	gaurdConstraintsMatrix38(5, 1) = 0.0;
	gaurdConstraintsMatrix38(5, 2) = 1.0;
	gaurdConstraintsMatrix38(5, 3) = 0.0;
	gaurdConstraintsMatrix38(6, 0) = 0.0;
	gaurdConstraintsMatrix38(6, 1) = 0.0;
	gaurdConstraintsMatrix38(6, 2) = 0.0;
	gaurdConstraintsMatrix38(6, 3) = -1.0;
	gaurdConstraintsMatrix38(7, 0) = 0.0;
	gaurdConstraintsMatrix38(7, 1) = 0.0;
	gaurdConstraintsMatrix38(7, 2) = 0.0;
	gaurdConstraintsMatrix38(7, 3) = 1.0;

	gaurdBoundValue38.resize(row);
	gaurdBoundValue38[0] = -2.0;
	gaurdBoundValue38[1] = 3.0;
	gaurdBoundValue38[2] = -4.0;
	gaurdBoundValue38[3] = 4.0;
	gaurdBoundValue38[4] = 1000.0;
	gaurdBoundValue38[5] = 1000.0;
	gaurdBoundValue38[6] = 1000.0;
	gaurdBoundValue38[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope38 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix38, gaurdBoundValue38,
					gaurdBoundSign));

	// The transition label ist37

	// Original guard: x1 = 2 & 3 <= x2 & x2 <= 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix39.resize(row, col);
	gaurdConstraintsMatrix39(0, 0) = -1.0;
	gaurdConstraintsMatrix39(0, 1) = 0.0;
	gaurdConstraintsMatrix39(0, 2) = 0.0;
	gaurdConstraintsMatrix39(0, 3) = 0.0;
	gaurdConstraintsMatrix39(1, 0) = 1.0;
	gaurdConstraintsMatrix39(1, 1) = 0.0;
	gaurdConstraintsMatrix39(1, 2) = 0.0;
	gaurdConstraintsMatrix39(1, 3) = 0.0;
	gaurdConstraintsMatrix39(2, 0) = 0.0;
	gaurdConstraintsMatrix39(2, 1) = -1.0;
	gaurdConstraintsMatrix39(2, 2) = 0.0;
	gaurdConstraintsMatrix39(2, 3) = 0.0;
	gaurdConstraintsMatrix39(3, 0) = 0.0;
	gaurdConstraintsMatrix39(3, 1) = 1.0;
	gaurdConstraintsMatrix39(3, 2) = 0.0;
	gaurdConstraintsMatrix39(3, 3) = 0.0;
	gaurdConstraintsMatrix39(4, 0) = 0.0;
	gaurdConstraintsMatrix39(4, 1) = 0.0;
	gaurdConstraintsMatrix39(4, 2) = -1.0;
	gaurdConstraintsMatrix39(4, 3) = 0.0;
	gaurdConstraintsMatrix39(5, 0) = 0.0;
	gaurdConstraintsMatrix39(5, 1) = 0.0;
	gaurdConstraintsMatrix39(5, 2) = 1.0;
	gaurdConstraintsMatrix39(5, 3) = 0.0;
	gaurdConstraintsMatrix39(6, 0) = 0.0;
	gaurdConstraintsMatrix39(6, 1) = 0.0;
	gaurdConstraintsMatrix39(6, 2) = 0.0;
	gaurdConstraintsMatrix39(6, 3) = -1.0;
	gaurdConstraintsMatrix39(7, 0) = 0.0;
	gaurdConstraintsMatrix39(7, 1) = 0.0;
	gaurdConstraintsMatrix39(7, 2) = 0.0;
	gaurdConstraintsMatrix39(7, 3) = 1.0;

	gaurdBoundValue39.resize(row);
	gaurdBoundValue39[0] = -2.0;
	gaurdBoundValue39[1] = 2.0;
	gaurdBoundValue39[2] = -3.0;
	gaurdBoundValue39[3] = 4.0;
	gaurdBoundValue39[4] = 1000.0;
	gaurdBoundValue39[5] = 1000.0;
	gaurdBoundValue39[6] = 1000.0;
	gaurdBoundValue39[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope39 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix39, gaurdBoundValue39,
					gaurdBoundSign));

	// The transition label ist39

	// Original guard: x1 = 3 & 3 <= x2 & x2 <= 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix40.resize(row, col);
	gaurdConstraintsMatrix40(0, 0) = -1.0;
	gaurdConstraintsMatrix40(0, 1) = 0.0;
	gaurdConstraintsMatrix40(0, 2) = 0.0;
	gaurdConstraintsMatrix40(0, 3) = 0.0;
	gaurdConstraintsMatrix40(1, 0) = 1.0;
	gaurdConstraintsMatrix40(1, 1) = 0.0;
	gaurdConstraintsMatrix40(1, 2) = 0.0;
	gaurdConstraintsMatrix40(1, 3) = 0.0;
	gaurdConstraintsMatrix40(2, 0) = 0.0;
	gaurdConstraintsMatrix40(2, 1) = -1.0;
	gaurdConstraintsMatrix40(2, 2) = 0.0;
	gaurdConstraintsMatrix40(2, 3) = 0.0;
	gaurdConstraintsMatrix40(3, 0) = 0.0;
	gaurdConstraintsMatrix40(3, 1) = 1.0;
	gaurdConstraintsMatrix40(3, 2) = 0.0;
	gaurdConstraintsMatrix40(3, 3) = 0.0;
	gaurdConstraintsMatrix40(4, 0) = 0.0;
	gaurdConstraintsMatrix40(4, 1) = 0.0;
	gaurdConstraintsMatrix40(4, 2) = -1.0;
	gaurdConstraintsMatrix40(4, 3) = 0.0;
	gaurdConstraintsMatrix40(5, 0) = 0.0;
	gaurdConstraintsMatrix40(5, 1) = 0.0;
	gaurdConstraintsMatrix40(5, 2) = 1.0;
	gaurdConstraintsMatrix40(5, 3) = 0.0;
	gaurdConstraintsMatrix40(6, 0) = 0.0;
	gaurdConstraintsMatrix40(6, 1) = 0.0;
	gaurdConstraintsMatrix40(6, 2) = 0.0;
	gaurdConstraintsMatrix40(6, 3) = -1.0;
	gaurdConstraintsMatrix40(7, 0) = 0.0;
	gaurdConstraintsMatrix40(7, 1) = 0.0;
	gaurdConstraintsMatrix40(7, 2) = 0.0;
	gaurdConstraintsMatrix40(7, 3) = 1.0;

	gaurdBoundValue40.resize(row);
	gaurdBoundValue40[0] = -3.0;
	gaurdBoundValue40[1] = 3.0;
	gaurdBoundValue40[2] = -3.0;
	gaurdBoundValue40[3] = 4.0;
	gaurdBoundValue40[4] = 1000.0;
	gaurdBoundValue40[5] = 1000.0;
	gaurdBoundValue40[6] = 1000.0;
	gaurdBoundValue40[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope40 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix40, gaurdBoundValue40,
					gaurdBoundSign));

	// The transition label ist38

	// Original guard: 2 <= x1 & x1 <= 3 & x2 = 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix41.resize(row, col);
	gaurdConstraintsMatrix41(0, 0) = -1.0;
	gaurdConstraintsMatrix41(0, 1) = 0.0;
	gaurdConstraintsMatrix41(0, 2) = 0.0;
	gaurdConstraintsMatrix41(0, 3) = 0.0;
	gaurdConstraintsMatrix41(1, 0) = 1.0;
	gaurdConstraintsMatrix41(1, 1) = 0.0;
	gaurdConstraintsMatrix41(1, 2) = 0.0;
	gaurdConstraintsMatrix41(1, 3) = 0.0;
	gaurdConstraintsMatrix41(2, 0) = 0.0;
	gaurdConstraintsMatrix41(2, 1) = -1.0;
	gaurdConstraintsMatrix41(2, 2) = 0.0;
	gaurdConstraintsMatrix41(2, 3) = 0.0;
	gaurdConstraintsMatrix41(3, 0) = 0.0;
	gaurdConstraintsMatrix41(3, 1) = 1.0;
	gaurdConstraintsMatrix41(3, 2) = 0.0;
	gaurdConstraintsMatrix41(3, 3) = 0.0;
	gaurdConstraintsMatrix41(4, 0) = 0.0;
	gaurdConstraintsMatrix41(4, 1) = 0.0;
	gaurdConstraintsMatrix41(4, 2) = -1.0;
	gaurdConstraintsMatrix41(4, 3) = 0.0;
	gaurdConstraintsMatrix41(5, 0) = 0.0;
	gaurdConstraintsMatrix41(5, 1) = 0.0;
	gaurdConstraintsMatrix41(5, 2) = 1.0;
	gaurdConstraintsMatrix41(5, 3) = 0.0;
	gaurdConstraintsMatrix41(6, 0) = 0.0;
	gaurdConstraintsMatrix41(6, 1) = 0.0;
	gaurdConstraintsMatrix41(6, 2) = 0.0;
	gaurdConstraintsMatrix41(6, 3) = -1.0;
	gaurdConstraintsMatrix41(7, 0) = 0.0;
	gaurdConstraintsMatrix41(7, 1) = 0.0;
	gaurdConstraintsMatrix41(7, 2) = 0.0;
	gaurdConstraintsMatrix41(7, 3) = 1.0;

	gaurdBoundValue41.resize(row);
	gaurdBoundValue41[0] = -2.0;
	gaurdBoundValue41[1] = 3.0;
	gaurdBoundValue41[2] = -3.0;
	gaurdBoundValue41[3] = 3.0;
	gaurdBoundValue41[4] = 1000.0;
	gaurdBoundValue41[5] = 1000.0;
	gaurdBoundValue41[6] = 1000.0;
	gaurdBoundValue41[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope41 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix41, gaurdBoundValue41,
					gaurdBoundSign));

	// The transition label ist58

	// Original guard: x1 = 3 & 4 <= x2 & x2 <= 5 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix42.resize(row, col);
	gaurdConstraintsMatrix42(0, 0) = -1.0;
	gaurdConstraintsMatrix42(0, 1) = 0.0;
	gaurdConstraintsMatrix42(0, 2) = 0.0;
	gaurdConstraintsMatrix42(0, 3) = 0.0;
	gaurdConstraintsMatrix42(1, 0) = 1.0;
	gaurdConstraintsMatrix42(1, 1) = 0.0;
	gaurdConstraintsMatrix42(1, 2) = 0.0;
	gaurdConstraintsMatrix42(1, 3) = 0.0;
	gaurdConstraintsMatrix42(2, 0) = 0.0;
	gaurdConstraintsMatrix42(2, 1) = -1.0;
	gaurdConstraintsMatrix42(2, 2) = 0.0;
	gaurdConstraintsMatrix42(2, 3) = 0.0;
	gaurdConstraintsMatrix42(3, 0) = 0.0;
	gaurdConstraintsMatrix42(3, 1) = 1.0;
	gaurdConstraintsMatrix42(3, 2) = 0.0;
	gaurdConstraintsMatrix42(3, 3) = 0.0;
	gaurdConstraintsMatrix42(4, 0) = 0.0;
	gaurdConstraintsMatrix42(4, 1) = 0.0;
	gaurdConstraintsMatrix42(4, 2) = -1.0;
	gaurdConstraintsMatrix42(4, 3) = 0.0;
	gaurdConstraintsMatrix42(5, 0) = 0.0;
	gaurdConstraintsMatrix42(5, 1) = 0.0;
	gaurdConstraintsMatrix42(5, 2) = 1.0;
	gaurdConstraintsMatrix42(5, 3) = 0.0;
	gaurdConstraintsMatrix42(6, 0) = 0.0;
	gaurdConstraintsMatrix42(6, 1) = 0.0;
	gaurdConstraintsMatrix42(6, 2) = 0.0;
	gaurdConstraintsMatrix42(6, 3) = -1.0;
	gaurdConstraintsMatrix42(7, 0) = 0.0;
	gaurdConstraintsMatrix42(7, 1) = 0.0;
	gaurdConstraintsMatrix42(7, 2) = 0.0;
	gaurdConstraintsMatrix42(7, 3) = 1.0;

	gaurdBoundValue42.resize(row);
	gaurdBoundValue42[0] = -3.0;
	gaurdBoundValue42[1] = 3.0;
	gaurdBoundValue42[2] = -4.0;
	gaurdBoundValue42[3] = 5.0;
	gaurdBoundValue42[4] = 1000.0;
	gaurdBoundValue42[5] = 1000.0;
	gaurdBoundValue42[6] = 1000.0;
	gaurdBoundValue42[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope42 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix42, gaurdBoundValue42,
					gaurdBoundSign));

	// The transition label ist60

	// Original guard: x1 = 4 & 4 <= x2 & x2 <= 5 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix43.resize(row, col);
	gaurdConstraintsMatrix43(0, 0) = -1.0;
	gaurdConstraintsMatrix43(0, 1) = 0.0;
	gaurdConstraintsMatrix43(0, 2) = 0.0;
	gaurdConstraintsMatrix43(0, 3) = 0.0;
	gaurdConstraintsMatrix43(1, 0) = 1.0;
	gaurdConstraintsMatrix43(1, 1) = 0.0;
	gaurdConstraintsMatrix43(1, 2) = 0.0;
	gaurdConstraintsMatrix43(1, 3) = 0.0;
	gaurdConstraintsMatrix43(2, 0) = 0.0;
	gaurdConstraintsMatrix43(2, 1) = -1.0;
	gaurdConstraintsMatrix43(2, 2) = 0.0;
	gaurdConstraintsMatrix43(2, 3) = 0.0;
	gaurdConstraintsMatrix43(3, 0) = 0.0;
	gaurdConstraintsMatrix43(3, 1) = 1.0;
	gaurdConstraintsMatrix43(3, 2) = 0.0;
	gaurdConstraintsMatrix43(3, 3) = 0.0;
	gaurdConstraintsMatrix43(4, 0) = 0.0;
	gaurdConstraintsMatrix43(4, 1) = 0.0;
	gaurdConstraintsMatrix43(4, 2) = -1.0;
	gaurdConstraintsMatrix43(4, 3) = 0.0;
	gaurdConstraintsMatrix43(5, 0) = 0.0;
	gaurdConstraintsMatrix43(5, 1) = 0.0;
	gaurdConstraintsMatrix43(5, 2) = 1.0;
	gaurdConstraintsMatrix43(5, 3) = 0.0;
	gaurdConstraintsMatrix43(6, 0) = 0.0;
	gaurdConstraintsMatrix43(6, 1) = 0.0;
	gaurdConstraintsMatrix43(6, 2) = 0.0;
	gaurdConstraintsMatrix43(6, 3) = -1.0;
	gaurdConstraintsMatrix43(7, 0) = 0.0;
	gaurdConstraintsMatrix43(7, 1) = 0.0;
	gaurdConstraintsMatrix43(7, 2) = 0.0;
	gaurdConstraintsMatrix43(7, 3) = 1.0;

	gaurdBoundValue43.resize(row);
	gaurdBoundValue43[0] = -4.0;
	gaurdBoundValue43[1] = 4.0;
	gaurdBoundValue43[2] = -4.0;
	gaurdBoundValue43[3] = 5.0;
	gaurdBoundValue43[4] = 1000.0;
	gaurdBoundValue43[5] = 1000.0;
	gaurdBoundValue43[6] = 1000.0;
	gaurdBoundValue43[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope43 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix43, gaurdBoundValue43,
					gaurdBoundSign));

	// The transition label ist59

	// Original guard: 3 <= x1 & x1 <= 4 & x2 = 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix44.resize(row, col);
	gaurdConstraintsMatrix44(0, 0) = -1.0;
	gaurdConstraintsMatrix44(0, 1) = 0.0;
	gaurdConstraintsMatrix44(0, 2) = 0.0;
	gaurdConstraintsMatrix44(0, 3) = 0.0;
	gaurdConstraintsMatrix44(1, 0) = 1.0;
	gaurdConstraintsMatrix44(1, 1) = 0.0;
	gaurdConstraintsMatrix44(1, 2) = 0.0;
	gaurdConstraintsMatrix44(1, 3) = 0.0;
	gaurdConstraintsMatrix44(2, 0) = 0.0;
	gaurdConstraintsMatrix44(2, 1) = -1.0;
	gaurdConstraintsMatrix44(2, 2) = 0.0;
	gaurdConstraintsMatrix44(2, 3) = 0.0;
	gaurdConstraintsMatrix44(3, 0) = 0.0;
	gaurdConstraintsMatrix44(3, 1) = 1.0;
	gaurdConstraintsMatrix44(3, 2) = 0.0;
	gaurdConstraintsMatrix44(3, 3) = 0.0;
	gaurdConstraintsMatrix44(4, 0) = 0.0;
	gaurdConstraintsMatrix44(4, 1) = 0.0;
	gaurdConstraintsMatrix44(4, 2) = -1.0;
	gaurdConstraintsMatrix44(4, 3) = 0.0;
	gaurdConstraintsMatrix44(5, 0) = 0.0;
	gaurdConstraintsMatrix44(5, 1) = 0.0;
	gaurdConstraintsMatrix44(5, 2) = 1.0;
	gaurdConstraintsMatrix44(5, 3) = 0.0;
	gaurdConstraintsMatrix44(6, 0) = 0.0;
	gaurdConstraintsMatrix44(6, 1) = 0.0;
	gaurdConstraintsMatrix44(6, 2) = 0.0;
	gaurdConstraintsMatrix44(6, 3) = -1.0;
	gaurdConstraintsMatrix44(7, 0) = 0.0;
	gaurdConstraintsMatrix44(7, 1) = 0.0;
	gaurdConstraintsMatrix44(7, 2) = 0.0;
	gaurdConstraintsMatrix44(7, 3) = 1.0;

	gaurdBoundValue44.resize(row);
	gaurdBoundValue44[0] = -3.0;
	gaurdBoundValue44[1] = 4.0;
	gaurdBoundValue44[2] = -4.0;
	gaurdBoundValue44[3] = 4.0;
	gaurdBoundValue44[4] = 1000.0;
	gaurdBoundValue44[5] = 1000.0;
	gaurdBoundValue44[6] = 1000.0;
	gaurdBoundValue44[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope44 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix44, gaurdBoundValue44,
					gaurdBoundSign));

	// The transition label ist55

	// Original guard: x1 = 3 & 3 <= x2 & x2 <= 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix45.resize(row, col);
	gaurdConstraintsMatrix45(0, 0) = -1.0;
	gaurdConstraintsMatrix45(0, 1) = 0.0;
	gaurdConstraintsMatrix45(0, 2) = 0.0;
	gaurdConstraintsMatrix45(0, 3) = 0.0;
	gaurdConstraintsMatrix45(1, 0) = 1.0;
	gaurdConstraintsMatrix45(1, 1) = 0.0;
	gaurdConstraintsMatrix45(1, 2) = 0.0;
	gaurdConstraintsMatrix45(1, 3) = 0.0;
	gaurdConstraintsMatrix45(2, 0) = 0.0;
	gaurdConstraintsMatrix45(2, 1) = -1.0;
	gaurdConstraintsMatrix45(2, 2) = 0.0;
	gaurdConstraintsMatrix45(2, 3) = 0.0;
	gaurdConstraintsMatrix45(3, 0) = 0.0;
	gaurdConstraintsMatrix45(3, 1) = 1.0;
	gaurdConstraintsMatrix45(3, 2) = 0.0;
	gaurdConstraintsMatrix45(3, 3) = 0.0;
	gaurdConstraintsMatrix45(4, 0) = 0.0;
	gaurdConstraintsMatrix45(4, 1) = 0.0;
	gaurdConstraintsMatrix45(4, 2) = -1.0;
	gaurdConstraintsMatrix45(4, 3) = 0.0;
	gaurdConstraintsMatrix45(5, 0) = 0.0;
	gaurdConstraintsMatrix45(5, 1) = 0.0;
	gaurdConstraintsMatrix45(5, 2) = 1.0;
	gaurdConstraintsMatrix45(5, 3) = 0.0;
	gaurdConstraintsMatrix45(6, 0) = 0.0;
	gaurdConstraintsMatrix45(6, 1) = 0.0;
	gaurdConstraintsMatrix45(6, 2) = 0.0;
	gaurdConstraintsMatrix45(6, 3) = -1.0;
	gaurdConstraintsMatrix45(7, 0) = 0.0;
	gaurdConstraintsMatrix45(7, 1) = 0.0;
	gaurdConstraintsMatrix45(7, 2) = 0.0;
	gaurdConstraintsMatrix45(7, 3) = 1.0;

	gaurdBoundValue45.resize(row);
	gaurdBoundValue45[0] = -3.0;
	gaurdBoundValue45[1] = 3.0;
	gaurdBoundValue45[2] = -3.0;
	gaurdBoundValue45[3] = 4.0;
	gaurdBoundValue45[4] = 1000.0;
	gaurdBoundValue45[5] = 1000.0;
	gaurdBoundValue45[6] = 1000.0;
	gaurdBoundValue45[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope45 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix45, gaurdBoundValue45,
					gaurdBoundSign));

	// The transition label ist57

	// Original guard: x1 = 4 & 3 <= x2 & x2 <= 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix46.resize(row, col);
	gaurdConstraintsMatrix46(0, 0) = -1.0;
	gaurdConstraintsMatrix46(0, 1) = 0.0;
	gaurdConstraintsMatrix46(0, 2) = 0.0;
	gaurdConstraintsMatrix46(0, 3) = 0.0;
	gaurdConstraintsMatrix46(1, 0) = 1.0;
	gaurdConstraintsMatrix46(1, 1) = 0.0;
	gaurdConstraintsMatrix46(1, 2) = 0.0;
	gaurdConstraintsMatrix46(1, 3) = 0.0;
	gaurdConstraintsMatrix46(2, 0) = 0.0;
	gaurdConstraintsMatrix46(2, 1) = -1.0;
	gaurdConstraintsMatrix46(2, 2) = 0.0;
	gaurdConstraintsMatrix46(2, 3) = 0.0;
	gaurdConstraintsMatrix46(3, 0) = 0.0;
	gaurdConstraintsMatrix46(3, 1) = 1.0;
	gaurdConstraintsMatrix46(3, 2) = 0.0;
	gaurdConstraintsMatrix46(3, 3) = 0.0;
	gaurdConstraintsMatrix46(4, 0) = 0.0;
	gaurdConstraintsMatrix46(4, 1) = 0.0;
	gaurdConstraintsMatrix46(4, 2) = -1.0;
	gaurdConstraintsMatrix46(4, 3) = 0.0;
	gaurdConstraintsMatrix46(5, 0) = 0.0;
	gaurdConstraintsMatrix46(5, 1) = 0.0;
	gaurdConstraintsMatrix46(5, 2) = 1.0;
	gaurdConstraintsMatrix46(5, 3) = 0.0;
	gaurdConstraintsMatrix46(6, 0) = 0.0;
	gaurdConstraintsMatrix46(6, 1) = 0.0;
	gaurdConstraintsMatrix46(6, 2) = 0.0;
	gaurdConstraintsMatrix46(6, 3) = -1.0;
	gaurdConstraintsMatrix46(7, 0) = 0.0;
	gaurdConstraintsMatrix46(7, 1) = 0.0;
	gaurdConstraintsMatrix46(7, 2) = 0.0;
	gaurdConstraintsMatrix46(7, 3) = 1.0;

	gaurdBoundValue46.resize(row);
	gaurdBoundValue46[0] = -4.0;
	gaurdBoundValue46[1] = 4.0;
	gaurdBoundValue46[2] = -3.0;
	gaurdBoundValue46[3] = 4.0;
	gaurdBoundValue46[4] = 1000.0;
	gaurdBoundValue46[5] = 1000.0;
	gaurdBoundValue46[6] = 1000.0;
	gaurdBoundValue46[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope46 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix46, gaurdBoundValue46,
					gaurdBoundSign));

	// The transition label ist54

	// Original guard: 3 <= x1 & x1 <= 4 & x2 = 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix47.resize(row, col);
	gaurdConstraintsMatrix47(0, 0) = -1.0;
	gaurdConstraintsMatrix47(0, 1) = 0.0;
	gaurdConstraintsMatrix47(0, 2) = 0.0;
	gaurdConstraintsMatrix47(0, 3) = 0.0;
	gaurdConstraintsMatrix47(1, 0) = 1.0;
	gaurdConstraintsMatrix47(1, 1) = 0.0;
	gaurdConstraintsMatrix47(1, 2) = 0.0;
	gaurdConstraintsMatrix47(1, 3) = 0.0;
	gaurdConstraintsMatrix47(2, 0) = 0.0;
	gaurdConstraintsMatrix47(2, 1) = -1.0;
	gaurdConstraintsMatrix47(2, 2) = 0.0;
	gaurdConstraintsMatrix47(2, 3) = 0.0;
	gaurdConstraintsMatrix47(3, 0) = 0.0;
	gaurdConstraintsMatrix47(3, 1) = 1.0;
	gaurdConstraintsMatrix47(3, 2) = 0.0;
	gaurdConstraintsMatrix47(3, 3) = 0.0;
	gaurdConstraintsMatrix47(4, 0) = 0.0;
	gaurdConstraintsMatrix47(4, 1) = 0.0;
	gaurdConstraintsMatrix47(4, 2) = -1.0;
	gaurdConstraintsMatrix47(4, 3) = 0.0;
	gaurdConstraintsMatrix47(5, 0) = 0.0;
	gaurdConstraintsMatrix47(5, 1) = 0.0;
	gaurdConstraintsMatrix47(5, 2) = 1.0;
	gaurdConstraintsMatrix47(5, 3) = 0.0;
	gaurdConstraintsMatrix47(6, 0) = 0.0;
	gaurdConstraintsMatrix47(6, 1) = 0.0;
	gaurdConstraintsMatrix47(6, 2) = 0.0;
	gaurdConstraintsMatrix47(6, 3) = -1.0;
	gaurdConstraintsMatrix47(7, 0) = 0.0;
	gaurdConstraintsMatrix47(7, 1) = 0.0;
	gaurdConstraintsMatrix47(7, 2) = 0.0;
	gaurdConstraintsMatrix47(7, 3) = 1.0;

	gaurdBoundValue47.resize(row);
	gaurdBoundValue47[0] = -3.0;
	gaurdBoundValue47[1] = 4.0;
	gaurdBoundValue47[2] = -4.0;
	gaurdBoundValue47[3] = 4.0;
	gaurdBoundValue47[4] = 1000.0;
	gaurdBoundValue47[5] = 1000.0;
	gaurdBoundValue47[6] = 1000.0;
	gaurdBoundValue47[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope47 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix47, gaurdBoundValue47,
					gaurdBoundSign));

	// The transition label ist56

	// Original guard: 3 <= x1 & x1 <= 4 & x2 = 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix48.resize(row, col);
	gaurdConstraintsMatrix48(0, 0) = -1.0;
	gaurdConstraintsMatrix48(0, 1) = 0.0;
	gaurdConstraintsMatrix48(0, 2) = 0.0;
	gaurdConstraintsMatrix48(0, 3) = 0.0;
	gaurdConstraintsMatrix48(1, 0) = 1.0;
	gaurdConstraintsMatrix48(1, 1) = 0.0;
	gaurdConstraintsMatrix48(1, 2) = 0.0;
	gaurdConstraintsMatrix48(1, 3) = 0.0;
	gaurdConstraintsMatrix48(2, 0) = 0.0;
	gaurdConstraintsMatrix48(2, 1) = -1.0;
	gaurdConstraintsMatrix48(2, 2) = 0.0;
	gaurdConstraintsMatrix48(2, 3) = 0.0;
	gaurdConstraintsMatrix48(3, 0) = 0.0;
	gaurdConstraintsMatrix48(3, 1) = 1.0;
	gaurdConstraintsMatrix48(3, 2) = 0.0;
	gaurdConstraintsMatrix48(3, 3) = 0.0;
	gaurdConstraintsMatrix48(4, 0) = 0.0;
	gaurdConstraintsMatrix48(4, 1) = 0.0;
	gaurdConstraintsMatrix48(4, 2) = -1.0;
	gaurdConstraintsMatrix48(4, 3) = 0.0;
	gaurdConstraintsMatrix48(5, 0) = 0.0;
	gaurdConstraintsMatrix48(5, 1) = 0.0;
	gaurdConstraintsMatrix48(5, 2) = 1.0;
	gaurdConstraintsMatrix48(5, 3) = 0.0;
	gaurdConstraintsMatrix48(6, 0) = 0.0;
	gaurdConstraintsMatrix48(6, 1) = 0.0;
	gaurdConstraintsMatrix48(6, 2) = 0.0;
	gaurdConstraintsMatrix48(6, 3) = -1.0;
	gaurdConstraintsMatrix48(7, 0) = 0.0;
	gaurdConstraintsMatrix48(7, 1) = 0.0;
	gaurdConstraintsMatrix48(7, 2) = 0.0;
	gaurdConstraintsMatrix48(7, 3) = 1.0;

	gaurdBoundValue48.resize(row);
	gaurdBoundValue48[0] = -3.0;
	gaurdBoundValue48[1] = 4.0;
	gaurdBoundValue48[2] = -3.0;
	gaurdBoundValue48[3] = 3.0;
	gaurdBoundValue48[4] = 1000.0;
	gaurdBoundValue48[5] = 1000.0;
	gaurdBoundValue48[6] = 1000.0;
	gaurdBoundValue48[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope48 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix48, gaurdBoundValue48,
					gaurdBoundSign));

	// The transition label ist46

	// Original guard: 3 <= x1 & x1 <= 4 & x2 = 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix49.resize(row, col);
	gaurdConstraintsMatrix49(0, 0) = -1.0;
	gaurdConstraintsMatrix49(0, 1) = 0.0;
	gaurdConstraintsMatrix49(0, 2) = 0.0;
	gaurdConstraintsMatrix49(0, 3) = 0.0;
	gaurdConstraintsMatrix49(1, 0) = 1.0;
	gaurdConstraintsMatrix49(1, 1) = 0.0;
	gaurdConstraintsMatrix49(1, 2) = 0.0;
	gaurdConstraintsMatrix49(1, 3) = 0.0;
	gaurdConstraintsMatrix49(2, 0) = 0.0;
	gaurdConstraintsMatrix49(2, 1) = -1.0;
	gaurdConstraintsMatrix49(2, 2) = 0.0;
	gaurdConstraintsMatrix49(2, 3) = 0.0;
	gaurdConstraintsMatrix49(3, 0) = 0.0;
	gaurdConstraintsMatrix49(3, 1) = 1.0;
	gaurdConstraintsMatrix49(3, 2) = 0.0;
	gaurdConstraintsMatrix49(3, 3) = 0.0;
	gaurdConstraintsMatrix49(4, 0) = 0.0;
	gaurdConstraintsMatrix49(4, 1) = 0.0;
	gaurdConstraintsMatrix49(4, 2) = -1.0;
	gaurdConstraintsMatrix49(4, 3) = 0.0;
	gaurdConstraintsMatrix49(5, 0) = 0.0;
	gaurdConstraintsMatrix49(5, 1) = 0.0;
	gaurdConstraintsMatrix49(5, 2) = 1.0;
	gaurdConstraintsMatrix49(5, 3) = 0.0;
	gaurdConstraintsMatrix49(6, 0) = 0.0;
	gaurdConstraintsMatrix49(6, 1) = 0.0;
	gaurdConstraintsMatrix49(6, 2) = 0.0;
	gaurdConstraintsMatrix49(6, 3) = -1.0;
	gaurdConstraintsMatrix49(7, 0) = 0.0;
	gaurdConstraintsMatrix49(7, 1) = 0.0;
	gaurdConstraintsMatrix49(7, 2) = 0.0;
	gaurdConstraintsMatrix49(7, 3) = 1.0;

	gaurdBoundValue49.resize(row);
	gaurdBoundValue49[0] = -3.0;
	gaurdBoundValue49[1] = 4.0;
	gaurdBoundValue49[2] = -2.0;
	gaurdBoundValue49[3] = 2.0;
	gaurdBoundValue49[4] = 1000.0;
	gaurdBoundValue49[5] = 1000.0;
	gaurdBoundValue49[6] = 1000.0;
	gaurdBoundValue49[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope49 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix49, gaurdBoundValue49,
					gaurdBoundSign));

	// The transition label ist49

	// Original guard: x1 = 4 & 1 <= x2 & x2 <= 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix50.resize(row, col);
	gaurdConstraintsMatrix50(0, 0) = -1.0;
	gaurdConstraintsMatrix50(0, 1) = 0.0;
	gaurdConstraintsMatrix50(0, 2) = 0.0;
	gaurdConstraintsMatrix50(0, 3) = 0.0;
	gaurdConstraintsMatrix50(1, 0) = 1.0;
	gaurdConstraintsMatrix50(1, 1) = 0.0;
	gaurdConstraintsMatrix50(1, 2) = 0.0;
	gaurdConstraintsMatrix50(1, 3) = 0.0;
	gaurdConstraintsMatrix50(2, 0) = 0.0;
	gaurdConstraintsMatrix50(2, 1) = -1.0;
	gaurdConstraintsMatrix50(2, 2) = 0.0;
	gaurdConstraintsMatrix50(2, 3) = 0.0;
	gaurdConstraintsMatrix50(3, 0) = 0.0;
	gaurdConstraintsMatrix50(3, 1) = 1.0;
	gaurdConstraintsMatrix50(3, 2) = 0.0;
	gaurdConstraintsMatrix50(3, 3) = 0.0;
	gaurdConstraintsMatrix50(4, 0) = 0.0;
	gaurdConstraintsMatrix50(4, 1) = 0.0;
	gaurdConstraintsMatrix50(4, 2) = -1.0;
	gaurdConstraintsMatrix50(4, 3) = 0.0;
	gaurdConstraintsMatrix50(5, 0) = 0.0;
	gaurdConstraintsMatrix50(5, 1) = 0.0;
	gaurdConstraintsMatrix50(5, 2) = 1.0;
	gaurdConstraintsMatrix50(5, 3) = 0.0;
	gaurdConstraintsMatrix50(6, 0) = 0.0;
	gaurdConstraintsMatrix50(6, 1) = 0.0;
	gaurdConstraintsMatrix50(6, 2) = 0.0;
	gaurdConstraintsMatrix50(6, 3) = -1.0;
	gaurdConstraintsMatrix50(7, 0) = 0.0;
	gaurdConstraintsMatrix50(7, 1) = 0.0;
	gaurdConstraintsMatrix50(7, 2) = 0.0;
	gaurdConstraintsMatrix50(7, 3) = 1.0;

	gaurdBoundValue50.resize(row);
	gaurdBoundValue50[0] = -4.0;
	gaurdBoundValue50[1] = 4.0;
	gaurdBoundValue50[2] = -1.0;
	gaurdBoundValue50[3] = 2.0;
	gaurdBoundValue50[4] = 1000.0;
	gaurdBoundValue50[5] = 1000.0;
	gaurdBoundValue50[6] = 1000.0;
	gaurdBoundValue50[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope50 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix50, gaurdBoundValue50,
					gaurdBoundSign));

	// The transition label ist48

	// Original guard: 3 <= x1 & x1 <= 4 & x2 = 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix51.resize(row, col);
	gaurdConstraintsMatrix51(0, 0) = -1.0;
	gaurdConstraintsMatrix51(0, 1) = 0.0;
	gaurdConstraintsMatrix51(0, 2) = 0.0;
	gaurdConstraintsMatrix51(0, 3) = 0.0;
	gaurdConstraintsMatrix51(1, 0) = 1.0;
	gaurdConstraintsMatrix51(1, 1) = 0.0;
	gaurdConstraintsMatrix51(1, 2) = 0.0;
	gaurdConstraintsMatrix51(1, 3) = 0.0;
	gaurdConstraintsMatrix51(2, 0) = 0.0;
	gaurdConstraintsMatrix51(2, 1) = -1.0;
	gaurdConstraintsMatrix51(2, 2) = 0.0;
	gaurdConstraintsMatrix51(2, 3) = 0.0;
	gaurdConstraintsMatrix51(3, 0) = 0.0;
	gaurdConstraintsMatrix51(3, 1) = 1.0;
	gaurdConstraintsMatrix51(3, 2) = 0.0;
	gaurdConstraintsMatrix51(3, 3) = 0.0;
	gaurdConstraintsMatrix51(4, 0) = 0.0;
	gaurdConstraintsMatrix51(4, 1) = 0.0;
	gaurdConstraintsMatrix51(4, 2) = -1.0;
	gaurdConstraintsMatrix51(4, 3) = 0.0;
	gaurdConstraintsMatrix51(5, 0) = 0.0;
	gaurdConstraintsMatrix51(5, 1) = 0.0;
	gaurdConstraintsMatrix51(5, 2) = 1.0;
	gaurdConstraintsMatrix51(5, 3) = 0.0;
	gaurdConstraintsMatrix51(6, 0) = 0.0;
	gaurdConstraintsMatrix51(6, 1) = 0.0;
	gaurdConstraintsMatrix51(6, 2) = 0.0;
	gaurdConstraintsMatrix51(6, 3) = -1.0;
	gaurdConstraintsMatrix51(7, 0) = 0.0;
	gaurdConstraintsMatrix51(7, 1) = 0.0;
	gaurdConstraintsMatrix51(7, 2) = 0.0;
	gaurdConstraintsMatrix51(7, 3) = 1.0;

	gaurdBoundValue51.resize(row);
	gaurdBoundValue51[0] = -3.0;
	gaurdBoundValue51[1] = 4.0;
	gaurdBoundValue51[2] = -1.0;
	gaurdBoundValue51[3] = 1.0;
	gaurdBoundValue51[4] = 1000.0;
	gaurdBoundValue51[5] = 1000.0;
	gaurdBoundValue51[6] = 1000.0;
	gaurdBoundValue51[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope51 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix51, gaurdBoundValue51,
					gaurdBoundSign));

	// The transition label ist47

	// Original guard: x1 = 3 & 1 <= x2 & x2 <= 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix52.resize(row, col);
	gaurdConstraintsMatrix52(0, 0) = -1.0;
	gaurdConstraintsMatrix52(0, 1) = 0.0;
	gaurdConstraintsMatrix52(0, 2) = 0.0;
	gaurdConstraintsMatrix52(0, 3) = 0.0;
	gaurdConstraintsMatrix52(1, 0) = 1.0;
	gaurdConstraintsMatrix52(1, 1) = 0.0;
	gaurdConstraintsMatrix52(1, 2) = 0.0;
	gaurdConstraintsMatrix52(1, 3) = 0.0;
	gaurdConstraintsMatrix52(2, 0) = 0.0;
	gaurdConstraintsMatrix52(2, 1) = -1.0;
	gaurdConstraintsMatrix52(2, 2) = 0.0;
	gaurdConstraintsMatrix52(2, 3) = 0.0;
	gaurdConstraintsMatrix52(3, 0) = 0.0;
	gaurdConstraintsMatrix52(3, 1) = 1.0;
	gaurdConstraintsMatrix52(3, 2) = 0.0;
	gaurdConstraintsMatrix52(3, 3) = 0.0;
	gaurdConstraintsMatrix52(4, 0) = 0.0;
	gaurdConstraintsMatrix52(4, 1) = 0.0;
	gaurdConstraintsMatrix52(4, 2) = -1.0;
	gaurdConstraintsMatrix52(4, 3) = 0.0;
	gaurdConstraintsMatrix52(5, 0) = 0.0;
	gaurdConstraintsMatrix52(5, 1) = 0.0;
	gaurdConstraintsMatrix52(5, 2) = 1.0;
	gaurdConstraintsMatrix52(5, 3) = 0.0;
	gaurdConstraintsMatrix52(6, 0) = 0.0;
	gaurdConstraintsMatrix52(6, 1) = 0.0;
	gaurdConstraintsMatrix52(6, 2) = 0.0;
	gaurdConstraintsMatrix52(6, 3) = -1.0;
	gaurdConstraintsMatrix52(7, 0) = 0.0;
	gaurdConstraintsMatrix52(7, 1) = 0.0;
	gaurdConstraintsMatrix52(7, 2) = 0.0;
	gaurdConstraintsMatrix52(7, 3) = 1.0;

	gaurdBoundValue52.resize(row);
	gaurdBoundValue52[0] = -3.0;
	gaurdBoundValue52[1] = 3.0;
	gaurdBoundValue52[2] = -1.0;
	gaurdBoundValue52[3] = 2.0;
	gaurdBoundValue52[4] = 1000.0;
	gaurdBoundValue52[5] = 1000.0;
	gaurdBoundValue52[6] = 1000.0;
	gaurdBoundValue52[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope52 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix52, gaurdBoundValue52,
					gaurdBoundSign));

	// The transition label ist50

	// Original guard: 3 <= x1 & x1 <= 4 & x2 = 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix53.resize(row, col);
	gaurdConstraintsMatrix53(0, 0) = -1.0;
	gaurdConstraintsMatrix53(0, 1) = 0.0;
	gaurdConstraintsMatrix53(0, 2) = 0.0;
	gaurdConstraintsMatrix53(0, 3) = 0.0;
	gaurdConstraintsMatrix53(1, 0) = 1.0;
	gaurdConstraintsMatrix53(1, 1) = 0.0;
	gaurdConstraintsMatrix53(1, 2) = 0.0;
	gaurdConstraintsMatrix53(1, 3) = 0.0;
	gaurdConstraintsMatrix53(2, 0) = 0.0;
	gaurdConstraintsMatrix53(2, 1) = -1.0;
	gaurdConstraintsMatrix53(2, 2) = 0.0;
	gaurdConstraintsMatrix53(2, 3) = 0.0;
	gaurdConstraintsMatrix53(3, 0) = 0.0;
	gaurdConstraintsMatrix53(3, 1) = 1.0;
	gaurdConstraintsMatrix53(3, 2) = 0.0;
	gaurdConstraintsMatrix53(3, 3) = 0.0;
	gaurdConstraintsMatrix53(4, 0) = 0.0;
	gaurdConstraintsMatrix53(4, 1) = 0.0;
	gaurdConstraintsMatrix53(4, 2) = -1.0;
	gaurdConstraintsMatrix53(4, 3) = 0.0;
	gaurdConstraintsMatrix53(5, 0) = 0.0;
	gaurdConstraintsMatrix53(5, 1) = 0.0;
	gaurdConstraintsMatrix53(5, 2) = 1.0;
	gaurdConstraintsMatrix53(5, 3) = 0.0;
	gaurdConstraintsMatrix53(6, 0) = 0.0;
	gaurdConstraintsMatrix53(6, 1) = 0.0;
	gaurdConstraintsMatrix53(6, 2) = 0.0;
	gaurdConstraintsMatrix53(6, 3) = -1.0;
	gaurdConstraintsMatrix53(7, 0) = 0.0;
	gaurdConstraintsMatrix53(7, 1) = 0.0;
	gaurdConstraintsMatrix53(7, 2) = 0.0;
	gaurdConstraintsMatrix53(7, 3) = 1.0;

	gaurdBoundValue53.resize(row);
	gaurdBoundValue53[0] = -3.0;
	gaurdBoundValue53[1] = 4.0;
	gaurdBoundValue53[2] = -3.0;
	gaurdBoundValue53[3] = 3.0;
	gaurdBoundValue53[4] = 1000.0;
	gaurdBoundValue53[5] = 1000.0;
	gaurdBoundValue53[6] = 1000.0;
	gaurdBoundValue53[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope53 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix53, gaurdBoundValue53,
					gaurdBoundSign));

	// The transition label ist51

	// Original guard: x1 = 3 & 2 <= x2 & x2 <= 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix54.resize(row, col);
	gaurdConstraintsMatrix54(0, 0) = -1.0;
	gaurdConstraintsMatrix54(0, 1) = 0.0;
	gaurdConstraintsMatrix54(0, 2) = 0.0;
	gaurdConstraintsMatrix54(0, 3) = 0.0;
	gaurdConstraintsMatrix54(1, 0) = 1.0;
	gaurdConstraintsMatrix54(1, 1) = 0.0;
	gaurdConstraintsMatrix54(1, 2) = 0.0;
	gaurdConstraintsMatrix54(1, 3) = 0.0;
	gaurdConstraintsMatrix54(2, 0) = 0.0;
	gaurdConstraintsMatrix54(2, 1) = -1.0;
	gaurdConstraintsMatrix54(2, 2) = 0.0;
	gaurdConstraintsMatrix54(2, 3) = 0.0;
	gaurdConstraintsMatrix54(3, 0) = 0.0;
	gaurdConstraintsMatrix54(3, 1) = 1.0;
	gaurdConstraintsMatrix54(3, 2) = 0.0;
	gaurdConstraintsMatrix54(3, 3) = 0.0;
	gaurdConstraintsMatrix54(4, 0) = 0.0;
	gaurdConstraintsMatrix54(4, 1) = 0.0;
	gaurdConstraintsMatrix54(4, 2) = -1.0;
	gaurdConstraintsMatrix54(4, 3) = 0.0;
	gaurdConstraintsMatrix54(5, 0) = 0.0;
	gaurdConstraintsMatrix54(5, 1) = 0.0;
	gaurdConstraintsMatrix54(5, 2) = 1.0;
	gaurdConstraintsMatrix54(5, 3) = 0.0;
	gaurdConstraintsMatrix54(6, 0) = 0.0;
	gaurdConstraintsMatrix54(6, 1) = 0.0;
	gaurdConstraintsMatrix54(6, 2) = 0.0;
	gaurdConstraintsMatrix54(6, 3) = -1.0;
	gaurdConstraintsMatrix54(7, 0) = 0.0;
	gaurdConstraintsMatrix54(7, 1) = 0.0;
	gaurdConstraintsMatrix54(7, 2) = 0.0;
	gaurdConstraintsMatrix54(7, 3) = 1.0;

	gaurdBoundValue54.resize(row);
	gaurdBoundValue54[0] = -3.0;
	gaurdBoundValue54[1] = 3.0;
	gaurdBoundValue54[2] = -2.0;
	gaurdBoundValue54[3] = 3.0;
	gaurdBoundValue54[4] = 1000.0;
	gaurdBoundValue54[5] = 1000.0;
	gaurdBoundValue54[6] = 1000.0;
	gaurdBoundValue54[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope54 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix54, gaurdBoundValue54,
					gaurdBoundSign));

	// The transition label ist53

	// Original guard: x1 = 4 & 2 <= x2 & x2 <= 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix55.resize(row, col);
	gaurdConstraintsMatrix55(0, 0) = -1.0;
	gaurdConstraintsMatrix55(0, 1) = 0.0;
	gaurdConstraintsMatrix55(0, 2) = 0.0;
	gaurdConstraintsMatrix55(0, 3) = 0.0;
	gaurdConstraintsMatrix55(1, 0) = 1.0;
	gaurdConstraintsMatrix55(1, 1) = 0.0;
	gaurdConstraintsMatrix55(1, 2) = 0.0;
	gaurdConstraintsMatrix55(1, 3) = 0.0;
	gaurdConstraintsMatrix55(2, 0) = 0.0;
	gaurdConstraintsMatrix55(2, 1) = -1.0;
	gaurdConstraintsMatrix55(2, 2) = 0.0;
	gaurdConstraintsMatrix55(2, 3) = 0.0;
	gaurdConstraintsMatrix55(3, 0) = 0.0;
	gaurdConstraintsMatrix55(3, 1) = 1.0;
	gaurdConstraintsMatrix55(3, 2) = 0.0;
	gaurdConstraintsMatrix55(3, 3) = 0.0;
	gaurdConstraintsMatrix55(4, 0) = 0.0;
	gaurdConstraintsMatrix55(4, 1) = 0.0;
	gaurdConstraintsMatrix55(4, 2) = -1.0;
	gaurdConstraintsMatrix55(4, 3) = 0.0;
	gaurdConstraintsMatrix55(5, 0) = 0.0;
	gaurdConstraintsMatrix55(5, 1) = 0.0;
	gaurdConstraintsMatrix55(5, 2) = 1.0;
	gaurdConstraintsMatrix55(5, 3) = 0.0;
	gaurdConstraintsMatrix55(6, 0) = 0.0;
	gaurdConstraintsMatrix55(6, 1) = 0.0;
	gaurdConstraintsMatrix55(6, 2) = 0.0;
	gaurdConstraintsMatrix55(6, 3) = -1.0;
	gaurdConstraintsMatrix55(7, 0) = 0.0;
	gaurdConstraintsMatrix55(7, 1) = 0.0;
	gaurdConstraintsMatrix55(7, 2) = 0.0;
	gaurdConstraintsMatrix55(7, 3) = 1.0;

	gaurdBoundValue55.resize(row);
	gaurdBoundValue55[0] = -4.0;
	gaurdBoundValue55[1] = 4.0;
	gaurdBoundValue55[2] = -2.0;
	gaurdBoundValue55[3] = 3.0;
	gaurdBoundValue55[4] = 1000.0;
	gaurdBoundValue55[5] = 1000.0;
	gaurdBoundValue55[6] = 1000.0;
	gaurdBoundValue55[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope55 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix55, gaurdBoundValue55,
					gaurdBoundSign));

	// The transition label ist52

	// Original guard: 3 <= x1 & x1 <= 4 & x2 = 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix56.resize(row, col);
	gaurdConstraintsMatrix56(0, 0) = -1.0;
	gaurdConstraintsMatrix56(0, 1) = 0.0;
	gaurdConstraintsMatrix56(0, 2) = 0.0;
	gaurdConstraintsMatrix56(0, 3) = 0.0;
	gaurdConstraintsMatrix56(1, 0) = 1.0;
	gaurdConstraintsMatrix56(1, 1) = 0.0;
	gaurdConstraintsMatrix56(1, 2) = 0.0;
	gaurdConstraintsMatrix56(1, 3) = 0.0;
	gaurdConstraintsMatrix56(2, 0) = 0.0;
	gaurdConstraintsMatrix56(2, 1) = -1.0;
	gaurdConstraintsMatrix56(2, 2) = 0.0;
	gaurdConstraintsMatrix56(2, 3) = 0.0;
	gaurdConstraintsMatrix56(3, 0) = 0.0;
	gaurdConstraintsMatrix56(3, 1) = 1.0;
	gaurdConstraintsMatrix56(3, 2) = 0.0;
	gaurdConstraintsMatrix56(3, 3) = 0.0;
	gaurdConstraintsMatrix56(4, 0) = 0.0;
	gaurdConstraintsMatrix56(4, 1) = 0.0;
	gaurdConstraintsMatrix56(4, 2) = -1.0;
	gaurdConstraintsMatrix56(4, 3) = 0.0;
	gaurdConstraintsMatrix56(5, 0) = 0.0;
	gaurdConstraintsMatrix56(5, 1) = 0.0;
	gaurdConstraintsMatrix56(5, 2) = 1.0;
	gaurdConstraintsMatrix56(5, 3) = 0.0;
	gaurdConstraintsMatrix56(6, 0) = 0.0;
	gaurdConstraintsMatrix56(6, 1) = 0.0;
	gaurdConstraintsMatrix56(6, 2) = 0.0;
	gaurdConstraintsMatrix56(6, 3) = -1.0;
	gaurdConstraintsMatrix56(7, 0) = 0.0;
	gaurdConstraintsMatrix56(7, 1) = 0.0;
	gaurdConstraintsMatrix56(7, 2) = 0.0;
	gaurdConstraintsMatrix56(7, 3) = 1.0;

	gaurdBoundValue56.resize(row);
	gaurdBoundValue56[0] = -3.0;
	gaurdBoundValue56[1] = 4.0;
	gaurdBoundValue56[2] = -2.0;
	gaurdBoundValue56[3] = 2.0;
	gaurdBoundValue56[4] = 1000.0;
	gaurdBoundValue56[5] = 1000.0;
	gaurdBoundValue56[6] = 1000.0;
	gaurdBoundValue56[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope56 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix56, gaurdBoundValue56,
					gaurdBoundSign));

	// The transition label ist43

	// Original guard: 3 <= x1 & x1 <= 4 & x2 = 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix57.resize(row, col);
	gaurdConstraintsMatrix57(0, 0) = -1.0;
	gaurdConstraintsMatrix57(0, 1) = 0.0;
	gaurdConstraintsMatrix57(0, 2) = 0.0;
	gaurdConstraintsMatrix57(0, 3) = 0.0;
	gaurdConstraintsMatrix57(1, 0) = 1.0;
	gaurdConstraintsMatrix57(1, 1) = 0.0;
	gaurdConstraintsMatrix57(1, 2) = 0.0;
	gaurdConstraintsMatrix57(1, 3) = 0.0;
	gaurdConstraintsMatrix57(2, 0) = 0.0;
	gaurdConstraintsMatrix57(2, 1) = -1.0;
	gaurdConstraintsMatrix57(2, 2) = 0.0;
	gaurdConstraintsMatrix57(2, 3) = 0.0;
	gaurdConstraintsMatrix57(3, 0) = 0.0;
	gaurdConstraintsMatrix57(3, 1) = 1.0;
	gaurdConstraintsMatrix57(3, 2) = 0.0;
	gaurdConstraintsMatrix57(3, 3) = 0.0;
	gaurdConstraintsMatrix57(4, 0) = 0.0;
	gaurdConstraintsMatrix57(4, 1) = 0.0;
	gaurdConstraintsMatrix57(4, 2) = -1.0;
	gaurdConstraintsMatrix57(4, 3) = 0.0;
	gaurdConstraintsMatrix57(5, 0) = 0.0;
	gaurdConstraintsMatrix57(5, 1) = 0.0;
	gaurdConstraintsMatrix57(5, 2) = 1.0;
	gaurdConstraintsMatrix57(5, 3) = 0.0;
	gaurdConstraintsMatrix57(6, 0) = 0.0;
	gaurdConstraintsMatrix57(6, 1) = 0.0;
	gaurdConstraintsMatrix57(6, 2) = 0.0;
	gaurdConstraintsMatrix57(6, 3) = -1.0;
	gaurdConstraintsMatrix57(7, 0) = 0.0;
	gaurdConstraintsMatrix57(7, 1) = 0.0;
	gaurdConstraintsMatrix57(7, 2) = 0.0;
	gaurdConstraintsMatrix57(7, 3) = 1.0;

	gaurdBoundValue57.resize(row);
	gaurdBoundValue57[0] = -3.0;
	gaurdBoundValue57[1] = 4.0;
	gaurdBoundValue57[2] = -1.0;
	gaurdBoundValue57[3] = 1.0;
	gaurdBoundValue57[4] = 1000.0;
	gaurdBoundValue57[5] = 1000.0;
	gaurdBoundValue57[6] = 1000.0;
	gaurdBoundValue57[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope57 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix57, gaurdBoundValue57,
					gaurdBoundSign));

	// The transition label ist45

	// Original guard: x1 = 4 & 0 <= x2 & x2 <= 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix58.resize(row, col);
	gaurdConstraintsMatrix58(0, 0) = -1.0;
	gaurdConstraintsMatrix58(0, 1) = 0.0;
	gaurdConstraintsMatrix58(0, 2) = 0.0;
	gaurdConstraintsMatrix58(0, 3) = 0.0;
	gaurdConstraintsMatrix58(1, 0) = 1.0;
	gaurdConstraintsMatrix58(1, 1) = 0.0;
	gaurdConstraintsMatrix58(1, 2) = 0.0;
	gaurdConstraintsMatrix58(1, 3) = 0.0;
	gaurdConstraintsMatrix58(2, 0) = 0.0;
	gaurdConstraintsMatrix58(2, 1) = -1.0;
	gaurdConstraintsMatrix58(2, 2) = 0.0;
	gaurdConstraintsMatrix58(2, 3) = 0.0;
	gaurdConstraintsMatrix58(3, 0) = 0.0;
	gaurdConstraintsMatrix58(3, 1) = 1.0;
	gaurdConstraintsMatrix58(3, 2) = 0.0;
	gaurdConstraintsMatrix58(3, 3) = 0.0;
	gaurdConstraintsMatrix58(4, 0) = 0.0;
	gaurdConstraintsMatrix58(4, 1) = 0.0;
	gaurdConstraintsMatrix58(4, 2) = -1.0;
	gaurdConstraintsMatrix58(4, 3) = 0.0;
	gaurdConstraintsMatrix58(5, 0) = 0.0;
	gaurdConstraintsMatrix58(5, 1) = 0.0;
	gaurdConstraintsMatrix58(5, 2) = 1.0;
	gaurdConstraintsMatrix58(5, 3) = 0.0;
	gaurdConstraintsMatrix58(6, 0) = 0.0;
	gaurdConstraintsMatrix58(6, 1) = 0.0;
	gaurdConstraintsMatrix58(6, 2) = 0.0;
	gaurdConstraintsMatrix58(6, 3) = -1.0;
	gaurdConstraintsMatrix58(7, 0) = 0.0;
	gaurdConstraintsMatrix58(7, 1) = 0.0;
	gaurdConstraintsMatrix58(7, 2) = 0.0;
	gaurdConstraintsMatrix58(7, 3) = 1.0;

	gaurdBoundValue58.resize(row);
	gaurdBoundValue58[0] = -4.0;
	gaurdBoundValue58[1] = 4.0;
	gaurdBoundValue58[2] = -0.0;
	gaurdBoundValue58[3] = 1.0;
	gaurdBoundValue58[4] = 1000.0;
	gaurdBoundValue58[5] = 1000.0;
	gaurdBoundValue58[6] = 1000.0;
	gaurdBoundValue58[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope58 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix58, gaurdBoundValue58,
					gaurdBoundSign));

	// The transition label ist44

	// Original guard: x1 = 3 & 0 <= x2 & x2 <= 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix59.resize(row, col);
	gaurdConstraintsMatrix59(0, 0) = -1.0;
	gaurdConstraintsMatrix59(0, 1) = 0.0;
	gaurdConstraintsMatrix59(0, 2) = 0.0;
	gaurdConstraintsMatrix59(0, 3) = 0.0;
	gaurdConstraintsMatrix59(1, 0) = 1.0;
	gaurdConstraintsMatrix59(1, 1) = 0.0;
	gaurdConstraintsMatrix59(1, 2) = 0.0;
	gaurdConstraintsMatrix59(1, 3) = 0.0;
	gaurdConstraintsMatrix59(2, 0) = 0.0;
	gaurdConstraintsMatrix59(2, 1) = -1.0;
	gaurdConstraintsMatrix59(2, 2) = 0.0;
	gaurdConstraintsMatrix59(2, 3) = 0.0;
	gaurdConstraintsMatrix59(3, 0) = 0.0;
	gaurdConstraintsMatrix59(3, 1) = 1.0;
	gaurdConstraintsMatrix59(3, 2) = 0.0;
	gaurdConstraintsMatrix59(3, 3) = 0.0;
	gaurdConstraintsMatrix59(4, 0) = 0.0;
	gaurdConstraintsMatrix59(4, 1) = 0.0;
	gaurdConstraintsMatrix59(4, 2) = -1.0;
	gaurdConstraintsMatrix59(4, 3) = 0.0;
	gaurdConstraintsMatrix59(5, 0) = 0.0;
	gaurdConstraintsMatrix59(5, 1) = 0.0;
	gaurdConstraintsMatrix59(5, 2) = 1.0;
	gaurdConstraintsMatrix59(5, 3) = 0.0;
	gaurdConstraintsMatrix59(6, 0) = 0.0;
	gaurdConstraintsMatrix59(6, 1) = 0.0;
	gaurdConstraintsMatrix59(6, 2) = 0.0;
	gaurdConstraintsMatrix59(6, 3) = -1.0;
	gaurdConstraintsMatrix59(7, 0) = 0.0;
	gaurdConstraintsMatrix59(7, 1) = 0.0;
	gaurdConstraintsMatrix59(7, 2) = 0.0;
	gaurdConstraintsMatrix59(7, 3) = 1.0;

	gaurdBoundValue59.resize(row);
	gaurdBoundValue59[0] = -3.0;
	gaurdBoundValue59[1] = 3.0;
	gaurdBoundValue59[2] = -0.0;
	gaurdBoundValue59[3] = 1.0;
	gaurdBoundValue59[4] = 1000.0;
	gaurdBoundValue59[5] = 1000.0;
	gaurdBoundValue59[6] = 1000.0;
	gaurdBoundValue59[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope59 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix59, gaurdBoundValue59,
					gaurdBoundSign));

	// The transition label ist72

	// Original guard: x1 = 4 & 4 <= x2 & x2 <= 5 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix60.resize(row, col);
	gaurdConstraintsMatrix60(0, 0) = -1.0;
	gaurdConstraintsMatrix60(0, 1) = 0.0;
	gaurdConstraintsMatrix60(0, 2) = 0.0;
	gaurdConstraintsMatrix60(0, 3) = 0.0;
	gaurdConstraintsMatrix60(1, 0) = 1.0;
	gaurdConstraintsMatrix60(1, 1) = 0.0;
	gaurdConstraintsMatrix60(1, 2) = 0.0;
	gaurdConstraintsMatrix60(1, 3) = 0.0;
	gaurdConstraintsMatrix60(2, 0) = 0.0;
	gaurdConstraintsMatrix60(2, 1) = -1.0;
	gaurdConstraintsMatrix60(2, 2) = 0.0;
	gaurdConstraintsMatrix60(2, 3) = 0.0;
	gaurdConstraintsMatrix60(3, 0) = 0.0;
	gaurdConstraintsMatrix60(3, 1) = 1.0;
	gaurdConstraintsMatrix60(3, 2) = 0.0;
	gaurdConstraintsMatrix60(3, 3) = 0.0;
	gaurdConstraintsMatrix60(4, 0) = 0.0;
	gaurdConstraintsMatrix60(4, 1) = 0.0;
	gaurdConstraintsMatrix60(4, 2) = -1.0;
	gaurdConstraintsMatrix60(4, 3) = 0.0;
	gaurdConstraintsMatrix60(5, 0) = 0.0;
	gaurdConstraintsMatrix60(5, 1) = 0.0;
	gaurdConstraintsMatrix60(5, 2) = 1.0;
	gaurdConstraintsMatrix60(5, 3) = 0.0;
	gaurdConstraintsMatrix60(6, 0) = 0.0;
	gaurdConstraintsMatrix60(6, 1) = 0.0;
	gaurdConstraintsMatrix60(6, 2) = 0.0;
	gaurdConstraintsMatrix60(6, 3) = -1.0;
	gaurdConstraintsMatrix60(7, 0) = 0.0;
	gaurdConstraintsMatrix60(7, 1) = 0.0;
	gaurdConstraintsMatrix60(7, 2) = 0.0;
	gaurdConstraintsMatrix60(7, 3) = 1.0;

	gaurdBoundValue60.resize(row);
	gaurdBoundValue60[0] = -4.0;
	gaurdBoundValue60[1] = 4.0;
	gaurdBoundValue60[2] = -4.0;
	gaurdBoundValue60[3] = 5.0;
	gaurdBoundValue60[4] = 1000.0;
	gaurdBoundValue60[5] = 1000.0;
	gaurdBoundValue60[6] = 1000.0;
	gaurdBoundValue60[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope60 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix60, gaurdBoundValue60,
					gaurdBoundSign));

	// The transition label ist73

	// Original guard: 4 <= x1 & x1 <= 5 & x2 = 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix61.resize(row, col);
	gaurdConstraintsMatrix61(0, 0) = -1.0;
	gaurdConstraintsMatrix61(0, 1) = 0.0;
	gaurdConstraintsMatrix61(0, 2) = 0.0;
	gaurdConstraintsMatrix61(0, 3) = 0.0;
	gaurdConstraintsMatrix61(1, 0) = 1.0;
	gaurdConstraintsMatrix61(1, 1) = 0.0;
	gaurdConstraintsMatrix61(1, 2) = 0.0;
	gaurdConstraintsMatrix61(1, 3) = 0.0;
	gaurdConstraintsMatrix61(2, 0) = 0.0;
	gaurdConstraintsMatrix61(2, 1) = -1.0;
	gaurdConstraintsMatrix61(2, 2) = 0.0;
	gaurdConstraintsMatrix61(2, 3) = 0.0;
	gaurdConstraintsMatrix61(3, 0) = 0.0;
	gaurdConstraintsMatrix61(3, 1) = 1.0;
	gaurdConstraintsMatrix61(3, 2) = 0.0;
	gaurdConstraintsMatrix61(3, 3) = 0.0;
	gaurdConstraintsMatrix61(4, 0) = 0.0;
	gaurdConstraintsMatrix61(4, 1) = 0.0;
	gaurdConstraintsMatrix61(4, 2) = -1.0;
	gaurdConstraintsMatrix61(4, 3) = 0.0;
	gaurdConstraintsMatrix61(5, 0) = 0.0;
	gaurdConstraintsMatrix61(5, 1) = 0.0;
	gaurdConstraintsMatrix61(5, 2) = 1.0;
	gaurdConstraintsMatrix61(5, 3) = 0.0;
	gaurdConstraintsMatrix61(6, 0) = 0.0;
	gaurdConstraintsMatrix61(6, 1) = 0.0;
	gaurdConstraintsMatrix61(6, 2) = 0.0;
	gaurdConstraintsMatrix61(6, 3) = -1.0;
	gaurdConstraintsMatrix61(7, 0) = 0.0;
	gaurdConstraintsMatrix61(7, 1) = 0.0;
	gaurdConstraintsMatrix61(7, 2) = 0.0;
	gaurdConstraintsMatrix61(7, 3) = 1.0;

	gaurdBoundValue61.resize(row);
	gaurdBoundValue61[0] = -4.0;
	gaurdBoundValue61[1] = 5.0;
	gaurdBoundValue61[2] = -4.0;
	gaurdBoundValue61[3] = 4.0;
	gaurdBoundValue61[4] = 1000.0;
	gaurdBoundValue61[5] = 1000.0;
	gaurdBoundValue61[6] = 1000.0;
	gaurdBoundValue61[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope61 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix61, gaurdBoundValue61,
					gaurdBoundSign));

	// The transition label ist70

	// Original guard: x1 = 4 & 3 <= x2 & x2 <= 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix62.resize(row, col);
	gaurdConstraintsMatrix62(0, 0) = -1.0;
	gaurdConstraintsMatrix62(0, 1) = 0.0;
	gaurdConstraintsMatrix62(0, 2) = 0.0;
	gaurdConstraintsMatrix62(0, 3) = 0.0;
	gaurdConstraintsMatrix62(1, 0) = 1.0;
	gaurdConstraintsMatrix62(1, 1) = 0.0;
	gaurdConstraintsMatrix62(1, 2) = 0.0;
	gaurdConstraintsMatrix62(1, 3) = 0.0;
	gaurdConstraintsMatrix62(2, 0) = 0.0;
	gaurdConstraintsMatrix62(2, 1) = -1.0;
	gaurdConstraintsMatrix62(2, 2) = 0.0;
	gaurdConstraintsMatrix62(2, 3) = 0.0;
	gaurdConstraintsMatrix62(3, 0) = 0.0;
	gaurdConstraintsMatrix62(3, 1) = 1.0;
	gaurdConstraintsMatrix62(3, 2) = 0.0;
	gaurdConstraintsMatrix62(3, 3) = 0.0;
	gaurdConstraintsMatrix62(4, 0) = 0.0;
	gaurdConstraintsMatrix62(4, 1) = 0.0;
	gaurdConstraintsMatrix62(4, 2) = -1.0;
	gaurdConstraintsMatrix62(4, 3) = 0.0;
	gaurdConstraintsMatrix62(5, 0) = 0.0;
	gaurdConstraintsMatrix62(5, 1) = 0.0;
	gaurdConstraintsMatrix62(5, 2) = 1.0;
	gaurdConstraintsMatrix62(5, 3) = 0.0;
	gaurdConstraintsMatrix62(6, 0) = 0.0;
	gaurdConstraintsMatrix62(6, 1) = 0.0;
	gaurdConstraintsMatrix62(6, 2) = 0.0;
	gaurdConstraintsMatrix62(6, 3) = -1.0;
	gaurdConstraintsMatrix62(7, 0) = 0.0;
	gaurdConstraintsMatrix62(7, 1) = 0.0;
	gaurdConstraintsMatrix62(7, 2) = 0.0;
	gaurdConstraintsMatrix62(7, 3) = 1.0;

	gaurdBoundValue62.resize(row);
	gaurdBoundValue62[0] = -4.0;
	gaurdBoundValue62[1] = 4.0;
	gaurdBoundValue62[2] = -3.0;
	gaurdBoundValue62[3] = 4.0;
	gaurdBoundValue62[4] = 1000.0;
	gaurdBoundValue62[5] = 1000.0;
	gaurdBoundValue62[6] = 1000.0;
	gaurdBoundValue62[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope62 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix62, gaurdBoundValue62,
					gaurdBoundSign));

	// The transition label ist69

	// Original guard: 4 <= x1 & x1 <= 5 & x2 = 4 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix63.resize(row, col);
	gaurdConstraintsMatrix63(0, 0) = -1.0;
	gaurdConstraintsMatrix63(0, 1) = 0.0;
	gaurdConstraintsMatrix63(0, 2) = 0.0;
	gaurdConstraintsMatrix63(0, 3) = 0.0;
	gaurdConstraintsMatrix63(1, 0) = 1.0;
	gaurdConstraintsMatrix63(1, 1) = 0.0;
	gaurdConstraintsMatrix63(1, 2) = 0.0;
	gaurdConstraintsMatrix63(1, 3) = 0.0;
	gaurdConstraintsMatrix63(2, 0) = 0.0;
	gaurdConstraintsMatrix63(2, 1) = -1.0;
	gaurdConstraintsMatrix63(2, 2) = 0.0;
	gaurdConstraintsMatrix63(2, 3) = 0.0;
	gaurdConstraintsMatrix63(3, 0) = 0.0;
	gaurdConstraintsMatrix63(3, 1) = 1.0;
	gaurdConstraintsMatrix63(3, 2) = 0.0;
	gaurdConstraintsMatrix63(3, 3) = 0.0;
	gaurdConstraintsMatrix63(4, 0) = 0.0;
	gaurdConstraintsMatrix63(4, 1) = 0.0;
	gaurdConstraintsMatrix63(4, 2) = -1.0;
	gaurdConstraintsMatrix63(4, 3) = 0.0;
	gaurdConstraintsMatrix63(5, 0) = 0.0;
	gaurdConstraintsMatrix63(5, 1) = 0.0;
	gaurdConstraintsMatrix63(5, 2) = 1.0;
	gaurdConstraintsMatrix63(5, 3) = 0.0;
	gaurdConstraintsMatrix63(6, 0) = 0.0;
	gaurdConstraintsMatrix63(6, 1) = 0.0;
	gaurdConstraintsMatrix63(6, 2) = 0.0;
	gaurdConstraintsMatrix63(6, 3) = -1.0;
	gaurdConstraintsMatrix63(7, 0) = 0.0;
	gaurdConstraintsMatrix63(7, 1) = 0.0;
	gaurdConstraintsMatrix63(7, 2) = 0.0;
	gaurdConstraintsMatrix63(7, 3) = 1.0;

	gaurdBoundValue63.resize(row);
	gaurdBoundValue63[0] = -4.0;
	gaurdBoundValue63[1] = 5.0;
	gaurdBoundValue63[2] = -4.0;
	gaurdBoundValue63[3] = 4.0;
	gaurdBoundValue63[4] = 1000.0;
	gaurdBoundValue63[5] = 1000.0;
	gaurdBoundValue63[6] = 1000.0;
	gaurdBoundValue63[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope63 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix63, gaurdBoundValue63,
					gaurdBoundSign));

	// The transition label ist71

	// Original guard: 4 <= x1 & x1 <= 5 & x2 = 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix64.resize(row, col);
	gaurdConstraintsMatrix64(0, 0) = -1.0;
	gaurdConstraintsMatrix64(0, 1) = 0.0;
	gaurdConstraintsMatrix64(0, 2) = 0.0;
	gaurdConstraintsMatrix64(0, 3) = 0.0;
	gaurdConstraintsMatrix64(1, 0) = 1.0;
	gaurdConstraintsMatrix64(1, 1) = 0.0;
	gaurdConstraintsMatrix64(1, 2) = 0.0;
	gaurdConstraintsMatrix64(1, 3) = 0.0;
	gaurdConstraintsMatrix64(2, 0) = 0.0;
	gaurdConstraintsMatrix64(2, 1) = -1.0;
	gaurdConstraintsMatrix64(2, 2) = 0.0;
	gaurdConstraintsMatrix64(2, 3) = 0.0;
	gaurdConstraintsMatrix64(3, 0) = 0.0;
	gaurdConstraintsMatrix64(3, 1) = 1.0;
	gaurdConstraintsMatrix64(3, 2) = 0.0;
	gaurdConstraintsMatrix64(3, 3) = 0.0;
	gaurdConstraintsMatrix64(4, 0) = 0.0;
	gaurdConstraintsMatrix64(4, 1) = 0.0;
	gaurdConstraintsMatrix64(4, 2) = -1.0;
	gaurdConstraintsMatrix64(4, 3) = 0.0;
	gaurdConstraintsMatrix64(5, 0) = 0.0;
	gaurdConstraintsMatrix64(5, 1) = 0.0;
	gaurdConstraintsMatrix64(5, 2) = 1.0;
	gaurdConstraintsMatrix64(5, 3) = 0.0;
	gaurdConstraintsMatrix64(6, 0) = 0.0;
	gaurdConstraintsMatrix64(6, 1) = 0.0;
	gaurdConstraintsMatrix64(6, 2) = 0.0;
	gaurdConstraintsMatrix64(6, 3) = -1.0;
	gaurdConstraintsMatrix64(7, 0) = 0.0;
	gaurdConstraintsMatrix64(7, 1) = 0.0;
	gaurdConstraintsMatrix64(7, 2) = 0.0;
	gaurdConstraintsMatrix64(7, 3) = 1.0;

	gaurdBoundValue64.resize(row);
	gaurdBoundValue64[0] = -4.0;
	gaurdBoundValue64[1] = 5.0;
	gaurdBoundValue64[2] = -3.0;
	gaurdBoundValue64[3] = 3.0;
	gaurdBoundValue64[4] = 1000.0;
	gaurdBoundValue64[5] = 1000.0;
	gaurdBoundValue64[6] = 1000.0;
	gaurdBoundValue64[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope64 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix64, gaurdBoundValue64,
					gaurdBoundSign));

	// The transition label ist66

	// Original guard: 4 <= x1 & x1 <= 5 & x2 = 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix65.resize(row, col);
	gaurdConstraintsMatrix65(0, 0) = -1.0;
	gaurdConstraintsMatrix65(0, 1) = 0.0;
	gaurdConstraintsMatrix65(0, 2) = 0.0;
	gaurdConstraintsMatrix65(0, 3) = 0.0;
	gaurdConstraintsMatrix65(1, 0) = 1.0;
	gaurdConstraintsMatrix65(1, 1) = 0.0;
	gaurdConstraintsMatrix65(1, 2) = 0.0;
	gaurdConstraintsMatrix65(1, 3) = 0.0;
	gaurdConstraintsMatrix65(2, 0) = 0.0;
	gaurdConstraintsMatrix65(2, 1) = -1.0;
	gaurdConstraintsMatrix65(2, 2) = 0.0;
	gaurdConstraintsMatrix65(2, 3) = 0.0;
	gaurdConstraintsMatrix65(3, 0) = 0.0;
	gaurdConstraintsMatrix65(3, 1) = 1.0;
	gaurdConstraintsMatrix65(3, 2) = 0.0;
	gaurdConstraintsMatrix65(3, 3) = 0.0;
	gaurdConstraintsMatrix65(4, 0) = 0.0;
	gaurdConstraintsMatrix65(4, 1) = 0.0;
	gaurdConstraintsMatrix65(4, 2) = -1.0;
	gaurdConstraintsMatrix65(4, 3) = 0.0;
	gaurdConstraintsMatrix65(5, 0) = 0.0;
	gaurdConstraintsMatrix65(5, 1) = 0.0;
	gaurdConstraintsMatrix65(5, 2) = 1.0;
	gaurdConstraintsMatrix65(5, 3) = 0.0;
	gaurdConstraintsMatrix65(6, 0) = 0.0;
	gaurdConstraintsMatrix65(6, 1) = 0.0;
	gaurdConstraintsMatrix65(6, 2) = 0.0;
	gaurdConstraintsMatrix65(6, 3) = -1.0;
	gaurdConstraintsMatrix65(7, 0) = 0.0;
	gaurdConstraintsMatrix65(7, 1) = 0.0;
	gaurdConstraintsMatrix65(7, 2) = 0.0;
	gaurdConstraintsMatrix65(7, 3) = 1.0;

	gaurdBoundValue65.resize(row);
	gaurdBoundValue65[0] = -4.0;
	gaurdBoundValue65[1] = 5.0;
	gaurdBoundValue65[2] = -3.0;
	gaurdBoundValue65[3] = 3.0;
	gaurdBoundValue65[4] = 1000.0;
	gaurdBoundValue65[5] = 1000.0;
	gaurdBoundValue65[6] = 1000.0;
	gaurdBoundValue65[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope65 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix65, gaurdBoundValue65,
					gaurdBoundSign));

	// The transition label ist67

	// Original guard: x1 = 4 & 2 <= x2 & x2 <= 3 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix66.resize(row, col);
	gaurdConstraintsMatrix66(0, 0) = -1.0;
	gaurdConstraintsMatrix66(0, 1) = 0.0;
	gaurdConstraintsMatrix66(0, 2) = 0.0;
	gaurdConstraintsMatrix66(0, 3) = 0.0;
	gaurdConstraintsMatrix66(1, 0) = 1.0;
	gaurdConstraintsMatrix66(1, 1) = 0.0;
	gaurdConstraintsMatrix66(1, 2) = 0.0;
	gaurdConstraintsMatrix66(1, 3) = 0.0;
	gaurdConstraintsMatrix66(2, 0) = 0.0;
	gaurdConstraintsMatrix66(2, 1) = -1.0;
	gaurdConstraintsMatrix66(2, 2) = 0.0;
	gaurdConstraintsMatrix66(2, 3) = 0.0;
	gaurdConstraintsMatrix66(3, 0) = 0.0;
	gaurdConstraintsMatrix66(3, 1) = 1.0;
	gaurdConstraintsMatrix66(3, 2) = 0.0;
	gaurdConstraintsMatrix66(3, 3) = 0.0;
	gaurdConstraintsMatrix66(4, 0) = 0.0;
	gaurdConstraintsMatrix66(4, 1) = 0.0;
	gaurdConstraintsMatrix66(4, 2) = -1.0;
	gaurdConstraintsMatrix66(4, 3) = 0.0;
	gaurdConstraintsMatrix66(5, 0) = 0.0;
	gaurdConstraintsMatrix66(5, 1) = 0.0;
	gaurdConstraintsMatrix66(5, 2) = 1.0;
	gaurdConstraintsMatrix66(5, 3) = 0.0;
	gaurdConstraintsMatrix66(6, 0) = 0.0;
	gaurdConstraintsMatrix66(6, 1) = 0.0;
	gaurdConstraintsMatrix66(6, 2) = 0.0;
	gaurdConstraintsMatrix66(6, 3) = -1.0;
	gaurdConstraintsMatrix66(7, 0) = 0.0;
	gaurdConstraintsMatrix66(7, 1) = 0.0;
	gaurdConstraintsMatrix66(7, 2) = 0.0;
	gaurdConstraintsMatrix66(7, 3) = 1.0;

	gaurdBoundValue66.resize(row);
	gaurdBoundValue66[0] = -4.0;
	gaurdBoundValue66[1] = 4.0;
	gaurdBoundValue66[2] = -2.0;
	gaurdBoundValue66[3] = 3.0;
	gaurdBoundValue66[4] = 1000.0;
	gaurdBoundValue66[5] = 1000.0;
	gaurdBoundValue66[6] = 1000.0;
	gaurdBoundValue66[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope66 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix66, gaurdBoundValue66,
					gaurdBoundSign));

	// The transition label ist68

	// Original guard: 4 <= x1 & x1 <= 5 & x2 = 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix67.resize(row, col);
	gaurdConstraintsMatrix67(0, 0) = -1.0;
	gaurdConstraintsMatrix67(0, 1) = 0.0;
	gaurdConstraintsMatrix67(0, 2) = 0.0;
	gaurdConstraintsMatrix67(0, 3) = 0.0;
	gaurdConstraintsMatrix67(1, 0) = 1.0;
	gaurdConstraintsMatrix67(1, 1) = 0.0;
	gaurdConstraintsMatrix67(1, 2) = 0.0;
	gaurdConstraintsMatrix67(1, 3) = 0.0;
	gaurdConstraintsMatrix67(2, 0) = 0.0;
	gaurdConstraintsMatrix67(2, 1) = -1.0;
	gaurdConstraintsMatrix67(2, 2) = 0.0;
	gaurdConstraintsMatrix67(2, 3) = 0.0;
	gaurdConstraintsMatrix67(3, 0) = 0.0;
	gaurdConstraintsMatrix67(3, 1) = 1.0;
	gaurdConstraintsMatrix67(3, 2) = 0.0;
	gaurdConstraintsMatrix67(3, 3) = 0.0;
	gaurdConstraintsMatrix67(4, 0) = 0.0;
	gaurdConstraintsMatrix67(4, 1) = 0.0;
	gaurdConstraintsMatrix67(4, 2) = -1.0;
	gaurdConstraintsMatrix67(4, 3) = 0.0;
	gaurdConstraintsMatrix67(5, 0) = 0.0;
	gaurdConstraintsMatrix67(5, 1) = 0.0;
	gaurdConstraintsMatrix67(5, 2) = 1.0;
	gaurdConstraintsMatrix67(5, 3) = 0.0;
	gaurdConstraintsMatrix67(6, 0) = 0.0;
	gaurdConstraintsMatrix67(6, 1) = 0.0;
	gaurdConstraintsMatrix67(6, 2) = 0.0;
	gaurdConstraintsMatrix67(6, 3) = -1.0;
	gaurdConstraintsMatrix67(7, 0) = 0.0;
	gaurdConstraintsMatrix67(7, 1) = 0.0;
	gaurdConstraintsMatrix67(7, 2) = 0.0;
	gaurdConstraintsMatrix67(7, 3) = 1.0;

	gaurdBoundValue67.resize(row);
	gaurdBoundValue67[0] = -4.0;
	gaurdBoundValue67[1] = 5.0;
	gaurdBoundValue67[2] = -2.0;
	gaurdBoundValue67[3] = 2.0;
	gaurdBoundValue67[4] = 1000.0;
	gaurdBoundValue67[5] = 1000.0;
	gaurdBoundValue67[6] = 1000.0;
	gaurdBoundValue67[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope67 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix67, gaurdBoundValue67,
					gaurdBoundSign));

	// The transition label ist63

	// Original guard: 4 <= x1 & x1 <= 5 & x2 = 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix68.resize(row, col);
	gaurdConstraintsMatrix68(0, 0) = -1.0;
	gaurdConstraintsMatrix68(0, 1) = 0.0;
	gaurdConstraintsMatrix68(0, 2) = 0.0;
	gaurdConstraintsMatrix68(0, 3) = 0.0;
	gaurdConstraintsMatrix68(1, 0) = 1.0;
	gaurdConstraintsMatrix68(1, 1) = 0.0;
	gaurdConstraintsMatrix68(1, 2) = 0.0;
	gaurdConstraintsMatrix68(1, 3) = 0.0;
	gaurdConstraintsMatrix68(2, 0) = 0.0;
	gaurdConstraintsMatrix68(2, 1) = -1.0;
	gaurdConstraintsMatrix68(2, 2) = 0.0;
	gaurdConstraintsMatrix68(2, 3) = 0.0;
	gaurdConstraintsMatrix68(3, 0) = 0.0;
	gaurdConstraintsMatrix68(3, 1) = 1.0;
	gaurdConstraintsMatrix68(3, 2) = 0.0;
	gaurdConstraintsMatrix68(3, 3) = 0.0;
	gaurdConstraintsMatrix68(4, 0) = 0.0;
	gaurdConstraintsMatrix68(4, 1) = 0.0;
	gaurdConstraintsMatrix68(4, 2) = -1.0;
	gaurdConstraintsMatrix68(4, 3) = 0.0;
	gaurdConstraintsMatrix68(5, 0) = 0.0;
	gaurdConstraintsMatrix68(5, 1) = 0.0;
	gaurdConstraintsMatrix68(5, 2) = 1.0;
	gaurdConstraintsMatrix68(5, 3) = 0.0;
	gaurdConstraintsMatrix68(6, 0) = 0.0;
	gaurdConstraintsMatrix68(6, 1) = 0.0;
	gaurdConstraintsMatrix68(6, 2) = 0.0;
	gaurdConstraintsMatrix68(6, 3) = -1.0;
	gaurdConstraintsMatrix68(7, 0) = 0.0;
	gaurdConstraintsMatrix68(7, 1) = 0.0;
	gaurdConstraintsMatrix68(7, 2) = 0.0;
	gaurdConstraintsMatrix68(7, 3) = 1.0;

	gaurdBoundValue68.resize(row);
	gaurdBoundValue68[0] = -4.0;
	gaurdBoundValue68[1] = 5.0;
	gaurdBoundValue68[2] = -2.0;
	gaurdBoundValue68[3] = 2.0;
	gaurdBoundValue68[4] = 1000.0;
	gaurdBoundValue68[5] = 1000.0;
	gaurdBoundValue68[6] = 1000.0;
	gaurdBoundValue68[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope68 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix68, gaurdBoundValue68,
					gaurdBoundSign));

	// The transition label ist65

	// Original guard: 4 <= x1 & x1 <= 5 & x2 = 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix69.resize(row, col);
	gaurdConstraintsMatrix69(0, 0) = -1.0;
	gaurdConstraintsMatrix69(0, 1) = 0.0;
	gaurdConstraintsMatrix69(0, 2) = 0.0;
	gaurdConstraintsMatrix69(0, 3) = 0.0;
	gaurdConstraintsMatrix69(1, 0) = 1.0;
	gaurdConstraintsMatrix69(1, 1) = 0.0;
	gaurdConstraintsMatrix69(1, 2) = 0.0;
	gaurdConstraintsMatrix69(1, 3) = 0.0;
	gaurdConstraintsMatrix69(2, 0) = 0.0;
	gaurdConstraintsMatrix69(2, 1) = -1.0;
	gaurdConstraintsMatrix69(2, 2) = 0.0;
	gaurdConstraintsMatrix69(2, 3) = 0.0;
	gaurdConstraintsMatrix69(3, 0) = 0.0;
	gaurdConstraintsMatrix69(3, 1) = 1.0;
	gaurdConstraintsMatrix69(3, 2) = 0.0;
	gaurdConstraintsMatrix69(3, 3) = 0.0;
	gaurdConstraintsMatrix69(4, 0) = 0.0;
	gaurdConstraintsMatrix69(4, 1) = 0.0;
	gaurdConstraintsMatrix69(4, 2) = -1.0;
	gaurdConstraintsMatrix69(4, 3) = 0.0;
	gaurdConstraintsMatrix69(5, 0) = 0.0;
	gaurdConstraintsMatrix69(5, 1) = 0.0;
	gaurdConstraintsMatrix69(5, 2) = 1.0;
	gaurdConstraintsMatrix69(5, 3) = 0.0;
	gaurdConstraintsMatrix69(6, 0) = 0.0;
	gaurdConstraintsMatrix69(6, 1) = 0.0;
	gaurdConstraintsMatrix69(6, 2) = 0.0;
	gaurdConstraintsMatrix69(6, 3) = -1.0;
	gaurdConstraintsMatrix69(7, 0) = 0.0;
	gaurdConstraintsMatrix69(7, 1) = 0.0;
	gaurdConstraintsMatrix69(7, 2) = 0.0;
	gaurdConstraintsMatrix69(7, 3) = 1.0;

	gaurdBoundValue69.resize(row);
	gaurdBoundValue69[0] = -4.0;
	gaurdBoundValue69[1] = 5.0;
	gaurdBoundValue69[2] = -1.0;
	gaurdBoundValue69[3] = 1.0;
	gaurdBoundValue69[4] = 1000.0;
	gaurdBoundValue69[5] = 1000.0;
	gaurdBoundValue69[6] = 1000.0;
	gaurdBoundValue69[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope69 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix69, gaurdBoundValue69,
					gaurdBoundSign));

	// The transition label ist64

	// Original guard: x1 = 4 & 1 <= x2 & x2 <= 2 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix70.resize(row, col);
	gaurdConstraintsMatrix70(0, 0) = -1.0;
	gaurdConstraintsMatrix70(0, 1) = 0.0;
	gaurdConstraintsMatrix70(0, 2) = 0.0;
	gaurdConstraintsMatrix70(0, 3) = 0.0;
	gaurdConstraintsMatrix70(1, 0) = 1.0;
	gaurdConstraintsMatrix70(1, 1) = 0.0;
	gaurdConstraintsMatrix70(1, 2) = 0.0;
	gaurdConstraintsMatrix70(1, 3) = 0.0;
	gaurdConstraintsMatrix70(2, 0) = 0.0;
	gaurdConstraintsMatrix70(2, 1) = -1.0;
	gaurdConstraintsMatrix70(2, 2) = 0.0;
	gaurdConstraintsMatrix70(2, 3) = 0.0;
	gaurdConstraintsMatrix70(3, 0) = 0.0;
	gaurdConstraintsMatrix70(3, 1) = 1.0;
	gaurdConstraintsMatrix70(3, 2) = 0.0;
	gaurdConstraintsMatrix70(3, 3) = 0.0;
	gaurdConstraintsMatrix70(4, 0) = 0.0;
	gaurdConstraintsMatrix70(4, 1) = 0.0;
	gaurdConstraintsMatrix70(4, 2) = -1.0;
	gaurdConstraintsMatrix70(4, 3) = 0.0;
	gaurdConstraintsMatrix70(5, 0) = 0.0;
	gaurdConstraintsMatrix70(5, 1) = 0.0;
	gaurdConstraintsMatrix70(5, 2) = 1.0;
	gaurdConstraintsMatrix70(5, 3) = 0.0;
	gaurdConstraintsMatrix70(6, 0) = 0.0;
	gaurdConstraintsMatrix70(6, 1) = 0.0;
	gaurdConstraintsMatrix70(6, 2) = 0.0;
	gaurdConstraintsMatrix70(6, 3) = -1.0;
	gaurdConstraintsMatrix70(7, 0) = 0.0;
	gaurdConstraintsMatrix70(7, 1) = 0.0;
	gaurdConstraintsMatrix70(7, 2) = 0.0;
	gaurdConstraintsMatrix70(7, 3) = 1.0;

	gaurdBoundValue70.resize(row);
	gaurdBoundValue70[0] = -4.0;
	gaurdBoundValue70[1] = 4.0;
	gaurdBoundValue70[2] = -1.0;
	gaurdBoundValue70[3] = 2.0;
	gaurdBoundValue70[4] = 1000.0;
	gaurdBoundValue70[5] = 1000.0;
	gaurdBoundValue70[6] = 1000.0;
	gaurdBoundValue70[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope70 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix70, gaurdBoundValue70,
					gaurdBoundSign));

	// The transition label ist61

	// Original guard: 4 <= x1 & x1 <= 5 & x2 = 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix71.resize(row, col);
	gaurdConstraintsMatrix71(0, 0) = -1.0;
	gaurdConstraintsMatrix71(0, 1) = 0.0;
	gaurdConstraintsMatrix71(0, 2) = 0.0;
	gaurdConstraintsMatrix71(0, 3) = 0.0;
	gaurdConstraintsMatrix71(1, 0) = 1.0;
	gaurdConstraintsMatrix71(1, 1) = 0.0;
	gaurdConstraintsMatrix71(1, 2) = 0.0;
	gaurdConstraintsMatrix71(1, 3) = 0.0;
	gaurdConstraintsMatrix71(2, 0) = 0.0;
	gaurdConstraintsMatrix71(2, 1) = -1.0;
	gaurdConstraintsMatrix71(2, 2) = 0.0;
	gaurdConstraintsMatrix71(2, 3) = 0.0;
	gaurdConstraintsMatrix71(3, 0) = 0.0;
	gaurdConstraintsMatrix71(3, 1) = 1.0;
	gaurdConstraintsMatrix71(3, 2) = 0.0;
	gaurdConstraintsMatrix71(3, 3) = 0.0;
	gaurdConstraintsMatrix71(4, 0) = 0.0;
	gaurdConstraintsMatrix71(4, 1) = 0.0;
	gaurdConstraintsMatrix71(4, 2) = -1.0;
	gaurdConstraintsMatrix71(4, 3) = 0.0;
	gaurdConstraintsMatrix71(5, 0) = 0.0;
	gaurdConstraintsMatrix71(5, 1) = 0.0;
	gaurdConstraintsMatrix71(5, 2) = 1.0;
	gaurdConstraintsMatrix71(5, 3) = 0.0;
	gaurdConstraintsMatrix71(6, 0) = 0.0;
	gaurdConstraintsMatrix71(6, 1) = 0.0;
	gaurdConstraintsMatrix71(6, 2) = 0.0;
	gaurdConstraintsMatrix71(6, 3) = -1.0;
	gaurdConstraintsMatrix71(7, 0) = 0.0;
	gaurdConstraintsMatrix71(7, 1) = 0.0;
	gaurdConstraintsMatrix71(7, 2) = 0.0;
	gaurdConstraintsMatrix71(7, 3) = 1.0;

	gaurdBoundValue71.resize(row);
	gaurdBoundValue71[0] = -4.0;
	gaurdBoundValue71[1] = 5.0;
	gaurdBoundValue71[2] = -1.0;
	gaurdBoundValue71[3] = 1.0;
	gaurdBoundValue71[4] = 1000.0;
	gaurdBoundValue71[5] = 1000.0;
	gaurdBoundValue71[6] = 1000.0;
	gaurdBoundValue71[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope71 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix71, gaurdBoundValue71,
					gaurdBoundSign));

	// The transition label ist62

	// Original guard: x1 = 4 & 0 <= x2 & x2 <= 1 & -1000 <= v1 & v1 <= 1000 & -1000 <= v2 & v2 <= 1000

	row = 8;
	col = 4;

	gaurdConstraintsMatrix72.resize(row, col);
	gaurdConstraintsMatrix72(0, 0) = -1.0;
	gaurdConstraintsMatrix72(0, 1) = 0.0;
	gaurdConstraintsMatrix72(0, 2) = 0.0;
	gaurdConstraintsMatrix72(0, 3) = 0.0;
	gaurdConstraintsMatrix72(1, 0) = 1.0;
	gaurdConstraintsMatrix72(1, 1) = 0.0;
	gaurdConstraintsMatrix72(1, 2) = 0.0;
	gaurdConstraintsMatrix72(1, 3) = 0.0;
	gaurdConstraintsMatrix72(2, 0) = 0.0;
	gaurdConstraintsMatrix72(2, 1) = -1.0;
	gaurdConstraintsMatrix72(2, 2) = 0.0;
	gaurdConstraintsMatrix72(2, 3) = 0.0;
	gaurdConstraintsMatrix72(3, 0) = 0.0;
	gaurdConstraintsMatrix72(3, 1) = 1.0;
	gaurdConstraintsMatrix72(3, 2) = 0.0;
	gaurdConstraintsMatrix72(3, 3) = 0.0;
	gaurdConstraintsMatrix72(4, 0) = 0.0;
	gaurdConstraintsMatrix72(4, 1) = 0.0;
	gaurdConstraintsMatrix72(4, 2) = -1.0;
	gaurdConstraintsMatrix72(4, 3) = 0.0;
	gaurdConstraintsMatrix72(5, 0) = 0.0;
	gaurdConstraintsMatrix72(5, 1) = 0.0;
	gaurdConstraintsMatrix72(5, 2) = 1.0;
	gaurdConstraintsMatrix72(5, 3) = 0.0;
	gaurdConstraintsMatrix72(6, 0) = 0.0;
	gaurdConstraintsMatrix72(6, 1) = 0.0;
	gaurdConstraintsMatrix72(6, 2) = 0.0;
	gaurdConstraintsMatrix72(6, 3) = -1.0;
	gaurdConstraintsMatrix72(7, 0) = 0.0;
	gaurdConstraintsMatrix72(7, 1) = 0.0;
	gaurdConstraintsMatrix72(7, 2) = 0.0;
	gaurdConstraintsMatrix72(7, 3) = 1.0;

	gaurdBoundValue72.resize(row);
	gaurdBoundValue72[0] = -4.0;
	gaurdBoundValue72[1] = 4.0;
	gaurdBoundValue72[2] = -0.0;
	gaurdBoundValue72[3] = 1.0;
	gaurdBoundValue72[4] = 1000.0;
	gaurdBoundValue72[5] = 1000.0;
	gaurdBoundValue72[6] = 1000.0;
	gaurdBoundValue72[7] = 1000.0;
	gaurdBoundSign = 1;
	gaurd_polytope72 = polytope::ptr(
			new polytope(gaurdConstraintsMatrix72, gaurdBoundValue72,
					gaurdBoundSign));

	// The transition label is   t6

	math::matrix<double> R0;
	row = 4;
	col = 4;
	R0.resize(row, col);
	R0(0, 0) = 1.0;
	R0(0, 1) = 0.0;
	R0(0, 2) = 0.0;
	R0(0, 3) = 0.0;
	R0(1, 0) = 0.0;
	R0(1, 1) = 1.0;
	R0(1, 2) = 0.0;
	R0(1, 3) = 0.0;
	R0(2, 0) = 0.0;
	R0(2, 1) = 0.0;
	R0(2, 2) = 1.0;
	R0(2, 3) = 0.0;
	R0(3, 0) = 0.0;
	R0(3, 1) = 0.0;
	R0(3, 2) = 0.0;
	R0(3, 3) = 1.0;
	std::vector<double> w0(row);
	w0[0] = 0.0;
	w0[1] = 0.0;
	w0[2] = 0.0;
	w0[3] = 0.0;

	Assign assignment0;
	assignment0.Map = R0;
	assignment0.b = w0;

	// The transition label is   t8

	math::matrix<double> R1;
	row = 4;
	col = 4;
	R1.resize(row, col);
	R1(0, 0) = 1.0;
	R1(0, 1) = 0.0;
	R1(0, 2) = 0.0;
	R1(0, 3) = 0.0;
	R1(1, 0) = 0.0;
	R1(1, 1) = 1.0;
	R1(1, 2) = 0.0;
	R1(1, 3) = 0.0;
	R1(2, 0) = 0.0;
	R1(2, 1) = 0.0;
	R1(2, 2) = 1.0;
	R1(2, 3) = 0.0;
	R1(3, 0) = 0.0;
	R1(3, 1) = 0.0;
	R1(3, 2) = 0.0;
	R1(3, 3) = 1.0;
	std::vector<double> w1(row);
	w1[0] = 0.0;
	w1[1] = 0.0;
	w1[2] = 0.0;
	w1[3] = 0.0;

	Assign assignment1;
	assignment1.Map = R1;
	assignment1.b = w1;

	// The transition label is   t7

	math::matrix<double> R2;
	row = 4;
	col = 4;
	R2.resize(row, col);
	R2(0, 0) = 1.0;
	R2(0, 1) = 0.0;
	R2(0, 2) = 0.0;
	R2(0, 3) = 0.0;
	R2(1, 0) = 0.0;
	R2(1, 1) = 1.0;
	R2(1, 2) = 0.0;
	R2(1, 3) = 0.0;
	R2(2, 0) = 0.0;
	R2(2, 1) = 0.0;
	R2(2, 2) = 1.0;
	R2(2, 3) = 0.0;
	R2(3, 0) = 0.0;
	R2(3, 1) = 0.0;
	R2(3, 2) = 0.0;
	R2(3, 3) = 1.0;
	std::vector<double> w2(row);
	w2[0] = 0.0;
	w2[1] = 0.0;
	w2[2] = 0.0;
	w2[3] = 0.0;

	Assign assignment2;
	assignment2.Map = R2;
	assignment2.b = w2;

	// The transition label is   t20

	math::matrix<double> R3;
	row = 4;
	col = 4;
	R3.resize(row, col);
	R3(0, 0) = 1.0;
	R3(0, 1) = 0.0;
	R3(0, 2) = 0.0;
	R3(0, 3) = 0.0;
	R3(1, 0) = 0.0;
	R3(1, 1) = 1.0;
	R3(1, 2) = 0.0;
	R3(1, 3) = 0.0;
	R3(2, 0) = 0.0;
	R3(2, 1) = 0.0;
	R3(2, 2) = 1.0;
	R3(2, 3) = 0.0;
	R3(3, 0) = 0.0;
	R3(3, 1) = 0.0;
	R3(3, 2) = 0.0;
	R3(3, 3) = 1.0;
	std::vector<double> w3(row);
	w3[0] = 0.0;
	w3[1] = 0.0;
	w3[2] = 0.0;
	w3[3] = 0.0;

	Assign assignment3;
	assignment3.Map = R3;
	assignment3.b = w3;

	// The transition label is   t21

	math::matrix<double> R4;
	row = 4;
	col = 4;
	R4.resize(row, col);
	R4(0, 0) = 1.0;
	R4(0, 1) = 0.0;
	R4(0, 2) = 0.0;
	R4(0, 3) = 0.0;
	R4(1, 0) = 0.0;
	R4(1, 1) = 1.0;
	R4(1, 2) = 0.0;
	R4(1, 3) = 0.0;
	R4(2, 0) = 0.0;
	R4(2, 1) = 0.0;
	R4(2, 2) = 1.0;
	R4(2, 3) = 0.0;
	R4(3, 0) = 0.0;
	R4(3, 1) = 0.0;
	R4(3, 2) = 0.0;
	R4(3, 3) = 1.0;
	std::vector<double> w4(row);
	w4[0] = 0.0;
	w4[1] = 0.0;
	w4[2] = 0.0;
	w4[3] = 0.0;

	Assign assignment4;
	assignment4.Map = R4;
	assignment4.b = w4;

	// The transition label is   t19

	math::matrix<double> R5;
	row = 4;
	col = 4;
	R5.resize(row, col);
	R5(0, 0) = 1.0;
	R5(0, 1) = 0.0;
	R5(0, 2) = 0.0;
	R5(0, 3) = 0.0;
	R5(1, 0) = 0.0;
	R5(1, 1) = 1.0;
	R5(1, 2) = 0.0;
	R5(1, 3) = 0.0;
	R5(2, 0) = 0.0;
	R5(2, 1) = 0.0;
	R5(2, 2) = 1.0;
	R5(2, 3) = 0.0;
	R5(3, 0) = 0.0;
	R5(3, 1) = 0.0;
	R5(3, 2) = 0.0;
	R5(3, 3) = 1.0;
	std::vector<double> w5(row);
	w5[0] = 0.0;
	w5[1] = 0.0;
	w5[2] = 0.0;
	w5[3] = 0.0;

	Assign assignment5;
	assignment5.Map = R5;
	assignment5.b = w5;

	// The transition label is   t18

	math::matrix<double> R6;
	row = 4;
	col = 4;
	R6.resize(row, col);
	R6(0, 0) = 1.0;
	R6(0, 1) = 0.0;
	R6(0, 2) = 0.0;
	R6(0, 3) = 0.0;
	R6(1, 0) = 0.0;
	R6(1, 1) = 1.0;
	R6(1, 2) = 0.0;
	R6(1, 3) = 0.0;
	R6(2, 0) = 0.0;
	R6(2, 1) = 0.0;
	R6(2, 2) = 1.0;
	R6(2, 3) = 0.0;
	R6(3, 0) = 0.0;
	R6(3, 1) = 0.0;
	R6(3, 2) = 0.0;
	R6(3, 3) = 1.0;
	std::vector<double> w6(row);
	w6[0] = 0.0;
	w6[1] = 0.0;
	w6[2] = 0.0;
	w6[3] = 0.0;

	Assign assignment6;
	assignment6.Map = R6;
	assignment6.b = w6;

	// The transition label is   t9

	math::matrix<double> R7;
	row = 4;
	col = 4;
	R7.resize(row, col);
	R7(0, 0) = 1.0;
	R7(0, 1) = 0.0;
	R7(0, 2) = 0.0;
	R7(0, 3) = 0.0;
	R7(1, 0) = 0.0;
	R7(1, 1) = 1.0;
	R7(1, 2) = 0.0;
	R7(1, 3) = 0.0;
	R7(2, 0) = 0.0;
	R7(2, 1) = 0.0;
	R7(2, 2) = 1.0;
	R7(2, 3) = 0.0;
	R7(3, 0) = 0.0;
	R7(3, 1) = 0.0;
	R7(3, 2) = 0.0;
	R7(3, 3) = 1.0;
	std::vector<double> w7(row);
	w7[0] = 0.0;
	w7[1] = 0.0;
	w7[2] = 0.0;
	w7[3] = 0.0;

	Assign assignment7;
	assignment7.Map = R7;
	assignment7.b = w7;

	// The transition label is   t11

	math::matrix<double> R8;
	row = 4;
	col = 4;
	R8.resize(row, col);
	R8(0, 0) = 1.0;
	R8(0, 1) = 0.0;
	R8(0, 2) = 0.0;
	R8(0, 3) = 0.0;
	R8(1, 0) = 0.0;
	R8(1, 1) = 1.0;
	R8(1, 2) = 0.0;
	R8(1, 3) = 0.0;
	R8(2, 0) = 0.0;
	R8(2, 1) = 0.0;
	R8(2, 2) = 1.0;
	R8(2, 3) = 0.0;
	R8(3, 0) = 0.0;
	R8(3, 1) = 0.0;
	R8(3, 2) = 0.0;
	R8(3, 3) = 1.0;
	std::vector<double> w8(row);
	w8[0] = 0.0;
	w8[1] = 0.0;
	w8[2] = 0.0;
	w8[3] = 0.0;

	Assign assignment8;
	assignment8.Map = R8;
	assignment8.b = w8;

	// The transition label is   t10

	math::matrix<double> R9;
	row = 4;
	col = 4;
	R9.resize(row, col);
	R9(0, 0) = 1.0;
	R9(0, 1) = 0.0;
	R9(0, 2) = 0.0;
	R9(0, 3) = 0.0;
	R9(1, 0) = 0.0;
	R9(1, 1) = 1.0;
	R9(1, 2) = 0.0;
	R9(1, 3) = 0.0;
	R9(2, 0) = 0.0;
	R9(2, 1) = 0.0;
	R9(2, 2) = 1.0;
	R9(2, 3) = 0.0;
	R9(3, 0) = 0.0;
	R9(3, 1) = 0.0;
	R9(3, 2) = 0.0;
	R9(3, 3) = 1.0;
	std::vector<double> w9(row);
	w9[0] = 0.0;
	w9[1] = 0.0;
	w9[2] = 0.0;
	w9[3] = 0.0;

	Assign assignment9;
	assignment9.Map = R9;
	assignment9.b = w9;

	// The transition label is   t15

	math::matrix<double> R10;
	row = 4;
	col = 4;
	R10.resize(row, col);
	R10(0, 0) = 1.0;
	R10(0, 1) = 0.0;
	R10(0, 2) = 0.0;
	R10(0, 3) = 0.0;
	R10(1, 0) = 0.0;
	R10(1, 1) = 1.0;
	R10(1, 2) = 0.0;
	R10(1, 3) = 0.0;
	R10(2, 0) = 0.0;
	R10(2, 1) = 0.0;
	R10(2, 2) = 1.0;
	R10(2, 3) = 0.0;
	R10(3, 0) = 0.0;
	R10(3, 1) = 0.0;
	R10(3, 2) = 0.0;
	R10(3, 3) = 1.0;
	std::vector<double> w10(row);
	w10[0] = 0.0;
	w10[1] = 0.0;
	w10[2] = 0.0;
	w10[3] = 0.0;

	Assign assignment10;
	assignment10.Map = R10;
	assignment10.b = w10;

	// The transition label is   t14

	math::matrix<double> R11;
	row = 4;
	col = 4;
	R11.resize(row, col);
	R11(0, 0) = 1.0;
	R11(0, 1) = 0.0;
	R11(0, 2) = 0.0;
	R11(0, 3) = 0.0;
	R11(1, 0) = 0.0;
	R11(1, 1) = 1.0;
	R11(1, 2) = 0.0;
	R11(1, 3) = 0.0;
	R11(2, 0) = 0.0;
	R11(2, 1) = 0.0;
	R11(2, 2) = 1.0;
	R11(2, 3) = 0.0;
	R11(3, 0) = 0.0;
	R11(3, 1) = 0.0;
	R11(3, 2) = 0.0;
	R11(3, 3) = 1.0;
	std::vector<double> w11(row);
	w11[0] = 0.0;
	w11[1] = 0.0;
	w11[2] = 0.0;
	w11[3] = 0.0;

	Assign assignment11;
	assignment11.Map = R11;
	assignment11.b = w11;

	// The transition label is   t17

	math::matrix<double> R12;
	row = 4;
	col = 4;
	R12.resize(row, col);
	R12(0, 0) = 1.0;
	R12(0, 1) = 0.0;
	R12(0, 2) = 0.0;
	R12(0, 3) = 0.0;
	R12(1, 0) = 0.0;
	R12(1, 1) = 1.0;
	R12(1, 2) = 0.0;
	R12(1, 3) = 0.0;
	R12(2, 0) = 0.0;
	R12(2, 1) = 0.0;
	R12(2, 2) = 1.0;
	R12(2, 3) = 0.0;
	R12(3, 0) = 0.0;
	R12(3, 1) = 0.0;
	R12(3, 2) = 0.0;
	R12(3, 3) = 1.0;
	std::vector<double> w12(row);
	w12[0] = 0.0;
	w12[1] = 0.0;
	w12[2] = 0.0;
	w12[3] = 0.0;

	Assign assignment12;
	assignment12.Map = R12;
	assignment12.b = w12;

	// The transition label is   t16

	math::matrix<double> R13;
	row = 4;
	col = 4;
	R13.resize(row, col);
	R13(0, 0) = 1.0;
	R13(0, 1) = 0.0;
	R13(0, 2) = 0.0;
	R13(0, 3) = 0.0;
	R13(1, 0) = 0.0;
	R13(1, 1) = 1.0;
	R13(1, 2) = 0.0;
	R13(1, 3) = 0.0;
	R13(2, 0) = 0.0;
	R13(2, 1) = 0.0;
	R13(2, 2) = 1.0;
	R13(2, 3) = 0.0;
	R13(3, 0) = 0.0;
	R13(3, 1) = 0.0;
	R13(3, 2) = 0.0;
	R13(3, 3) = 1.0;
	std::vector<double> w13(row);
	w13[0] = 0.0;
	w13[1] = 0.0;
	w13[2] = 0.0;
	w13[3] = 0.0;

	Assign assignment13;
	assignment13.Map = R13;
	assignment13.b = w13;

	// The transition label is   t32

	math::matrix<double> R14;
	row = 4;
	col = 4;
	R14.resize(row, col);
	R14(0, 0) = 1.0;
	R14(0, 1) = 0.0;
	R14(0, 2) = 0.0;
	R14(0, 3) = 0.0;
	R14(1, 0) = 0.0;
	R14(1, 1) = 1.0;
	R14(1, 2) = 0.0;
	R14(1, 3) = 0.0;
	R14(2, 0) = 0.0;
	R14(2, 1) = 0.0;
	R14(2, 2) = 1.0;
	R14(2, 3) = 0.0;
	R14(3, 0) = 0.0;
	R14(3, 1) = 0.0;
	R14(3, 2) = 0.0;
	R14(3, 3) = 1.0;
	std::vector<double> w14(row);
	w14[0] = 0.0;
	w14[1] = 0.0;
	w14[2] = 0.0;
	w14[3] = 0.0;

	Assign assignment14;
	assignment14.Map = R14;
	assignment14.b = w14;

	// The transition label is   t33

	math::matrix<double> R15;
	row = 4;
	col = 4;
	R15.resize(row, col);
	R15(0, 0) = 1.0;
	R15(0, 1) = 0.0;
	R15(0, 2) = 0.0;
	R15(0, 3) = 0.0;
	R15(1, 0) = 0.0;
	R15(1, 1) = 1.0;
	R15(1, 2) = 0.0;
	R15(1, 3) = 0.0;
	R15(2, 0) = 0.0;
	R15(2, 1) = 0.0;
	R15(2, 2) = 1.0;
	R15(2, 3) = 0.0;
	R15(3, 0) = 0.0;
	R15(3, 1) = 0.0;
	R15(3, 2) = 0.0;
	R15(3, 3) = 1.0;
	std::vector<double> w15(row);
	w15[0] = 0.0;
	w15[1] = 0.0;
	w15[2] = 0.0;
	w15[3] = 0.0;

	Assign assignment15;
	assignment15.Map = R15;
	assignment15.b = w15;

	// The transition label is   t34

	math::matrix<double> R16;
	row = 4;
	col = 4;
	R16.resize(row, col);
	R16(0, 0) = 1.0;
	R16(0, 1) = 0.0;
	R16(0, 2) = 0.0;
	R16(0, 3) = 0.0;
	R16(1, 0) = 0.0;
	R16(1, 1) = 1.0;
	R16(1, 2) = 0.0;
	R16(1, 3) = 0.0;
	R16(2, 0) = 0.0;
	R16(2, 1) = 0.0;
	R16(2, 2) = 1.0;
	R16(2, 3) = 0.0;
	R16(3, 0) = 0.0;
	R16(3, 1) = 0.0;
	R16(3, 2) = 0.0;
	R16(3, 3) = 1.0;
	std::vector<double> w16(row);
	w16[0] = 0.0;
	w16[1] = 0.0;
	w16[2] = 0.0;
	w16[3] = 0.0;

	Assign assignment16;
	assignment16.Map = R16;
	assignment16.b = w16;

	// The transition label is   t35

	math::matrix<double> R17;
	row = 4;
	col = 4;
	R17.resize(row, col);
	R17(0, 0) = 1.0;
	R17(0, 1) = 0.0;
	R17(0, 2) = 0.0;
	R17(0, 3) = 0.0;
	R17(1, 0) = 0.0;
	R17(1, 1) = 1.0;
	R17(1, 2) = 0.0;
	R17(1, 3) = 0.0;
	R17(2, 0) = 0.0;
	R17(2, 1) = 0.0;
	R17(2, 2) = 1.0;
	R17(2, 3) = 0.0;
	R17(3, 0) = 0.0;
	R17(3, 1) = 0.0;
	R17(3, 2) = 0.0;
	R17(3, 3) = 1.0;
	std::vector<double> w17(row);
	w17[0] = 0.0;
	w17[1] = 0.0;
	w17[2] = 0.0;
	w17[3] = 0.0;

	Assign assignment17;
	assignment17.Map = R17;
	assignment17.b = w17;

	// The transition label is   t31

	math::matrix<double> R18;
	row = 4;
	col = 4;
	R18.resize(row, col);
	R18(0, 0) = 1.0;
	R18(0, 1) = 0.0;
	R18(0, 2) = 0.0;
	R18(0, 3) = 0.0;
	R18(1, 0) = 0.0;
	R18(1, 1) = 1.0;
	R18(1, 2) = 0.0;
	R18(1, 3) = 0.0;
	R18(2, 0) = 0.0;
	R18(2, 1) = 0.0;
	R18(2, 2) = 1.0;
	R18(2, 3) = 0.0;
	R18(3, 0) = 0.0;
	R18(3, 1) = 0.0;
	R18(3, 2) = 0.0;
	R18(3, 3) = 1.0;
	std::vector<double> w18(row);
	w18[0] = 0.0;
	w18[1] = 0.0;
	w18[2] = 0.0;
	w18[3] = 0.0;

	Assign assignment18;
	assignment18.Map = R18;
	assignment18.b = w18;

	// The transition label is   t30

	math::matrix<double> R19;
	row = 4;
	col = 4;
	R19.resize(row, col);
	R19(0, 0) = 1.0;
	R19(0, 1) = 0.0;
	R19(0, 2) = 0.0;
	R19(0, 3) = 0.0;
	R19(1, 0) = 0.0;
	R19(1, 1) = 1.0;
	R19(1, 2) = 0.0;
	R19(1, 3) = 0.0;
	R19(2, 0) = 0.0;
	R19(2, 1) = 0.0;
	R19(2, 2) = 1.0;
	R19(2, 3) = 0.0;
	R19(3, 0) = 0.0;
	R19(3, 1) = 0.0;
	R19(3, 2) = 0.0;
	R19(3, 3) = 1.0;
	std::vector<double> w19(row);
	w19[0] = 0.0;
	w19[1] = 0.0;
	w19[2] = 0.0;
	w19[3] = 0.0;

	Assign assignment19;
	assignment19.Map = R19;
	assignment19.b = w19;

	// The transition label is   t29

	math::matrix<double> R20;
	row = 4;
	col = 4;
	R20.resize(row, col);
	R20(0, 0) = 1.0;
	R20(0, 1) = 0.0;
	R20(0, 2) = 0.0;
	R20(0, 3) = 0.0;
	R20(1, 0) = 0.0;
	R20(1, 1) = 1.0;
	R20(1, 2) = 0.0;
	R20(1, 3) = 0.0;
	R20(2, 0) = 0.0;
	R20(2, 1) = 0.0;
	R20(2, 2) = 1.0;
	R20(2, 3) = 0.0;
	R20(3, 0) = 0.0;
	R20(3, 1) = 0.0;
	R20(3, 2) = 0.0;
	R20(3, 3) = 1.0;
	std::vector<double> w20(row);
	w20[0] = 0.0;
	w20[1] = 0.0;
	w20[2] = 0.0;
	w20[3] = 0.0;

	Assign assignment20;
	assignment20.Map = R20;
	assignment20.b = w20;

	// The transition label is   t13

	math::matrix<double> R21;
	row = 4;
	col = 4;
	R21.resize(row, col);
	R21(0, 0) = 1.0;
	R21(0, 1) = 0.0;
	R21(0, 2) = 0.0;
	R21(0, 3) = 0.0;
	R21(1, 0) = 0.0;
	R21(1, 1) = 1.0;
	R21(1, 2) = 0.0;
	R21(1, 3) = 0.0;
	R21(2, 0) = 0.0;
	R21(2, 1) = 0.0;
	R21(2, 2) = 1.0;
	R21(2, 3) = 0.0;
	R21(3, 0) = 0.0;
	R21(3, 1) = 0.0;
	R21(3, 2) = 0.0;
	R21(3, 3) = 1.0;
	std::vector<double> w21(row);
	w21[0] = 0.0;
	w21[1] = 0.0;
	w21[2] = 0.0;
	w21[3] = 0.0;

	Assign assignment21;
	assignment21.Map = R21;
	assignment21.b = w21;

	// The transition label is   t12

	math::matrix<double> R22;
	row = 4;
	col = 4;
	R22.resize(row, col);
	R22(0, 0) = 1.0;
	R22(0, 1) = 0.0;
	R22(0, 2) = 0.0;
	R22(0, 3) = 0.0;
	R22(1, 0) = 0.0;
	R22(1, 1) = 1.0;
	R22(1, 2) = 0.0;
	R22(1, 3) = 0.0;
	R22(2, 0) = 0.0;
	R22(2, 1) = 0.0;
	R22(2, 2) = 1.0;
	R22(2, 3) = 0.0;
	R22(3, 0) = 0.0;
	R22(3, 1) = 0.0;
	R22(3, 2) = 0.0;
	R22(3, 3) = 1.0;
	std::vector<double> w22(row);
	w22[0] = 0.0;
	w22[1] = 0.0;
	w22[2] = 0.0;
	w22[3] = 0.0;

	Assign assignment22;
	assignment22.Map = R22;
	assignment22.b = w22;

	// The transition label is   t4

	math::matrix<double> R23;
	row = 4;
	col = 4;
	R23.resize(row, col);
	R23(0, 0) = 1.0;
	R23(0, 1) = 0.0;
	R23(0, 2) = 0.0;
	R23(0, 3) = 0.0;
	R23(1, 0) = 0.0;
	R23(1, 1) = 1.0;
	R23(1, 2) = 0.0;
	R23(1, 3) = 0.0;
	R23(2, 0) = 0.0;
	R23(2, 1) = 0.0;
	R23(2, 2) = 1.0;
	R23(2, 3) = 0.0;
	R23(3, 0) = 0.0;
	R23(3, 1) = 0.0;
	R23(3, 2) = 0.0;
	R23(3, 3) = 1.0;
	std::vector<double> w23(row);
	w23[0] = 0.0;
	w23[1] = 0.0;
	w23[2] = 0.0;
	w23[3] = 0.0;

	Assign assignment23;
	assignment23.Map = R23;
	assignment23.b = w23;

	// The transition label is   t3

	math::matrix<double> R24;
	row = 4;
	col = 4;
	R24.resize(row, col);
	R24(0, 0) = 1.0;
	R24(0, 1) = 0.0;
	R24(0, 2) = 0.0;
	R24(0, 3) = 0.0;
	R24(1, 0) = 0.0;
	R24(1, 1) = 1.0;
	R24(1, 2) = 0.0;
	R24(1, 3) = 0.0;
	R24(2, 0) = 0.0;
	R24(2, 1) = 0.0;
	R24(2, 2) = 1.0;
	R24(2, 3) = 0.0;
	R24(3, 0) = 0.0;
	R24(3, 1) = 0.0;
	R24(3, 2) = 0.0;
	R24(3, 3) = 1.0;
	std::vector<double> w24(row);
	w24[0] = 0.0;
	w24[1] = 0.0;
	w24[2] = 0.0;
	w24[3] = 0.0;

	Assign assignment24;
	assignment24.Map = R24;
	assignment24.b = w24;

	// The transition label is   t5

	math::matrix<double> R25;
	row = 4;
	col = 4;
	R25.resize(row, col);
	R25(0, 0) = 1.0;
	R25(0, 1) = 0.0;
	R25(0, 2) = 0.0;
	R25(0, 3) = 0.0;
	R25(1, 0) = 0.0;
	R25(1, 1) = 1.0;
	R25(1, 2) = 0.0;
	R25(1, 3) = 0.0;
	R25(2, 0) = 0.0;
	R25(2, 1) = 0.0;
	R25(2, 2) = 1.0;
	R25(2, 3) = 0.0;
	R25(3, 0) = 0.0;
	R25(3, 1) = 0.0;
	R25(3, 2) = 0.0;
	R25(3, 3) = 1.0;
	std::vector<double> w25(row);
	w25[0] = 0.0;
	w25[1] = 0.0;
	w25[2] = 0.0;
	w25[3] = 0.0;

	Assign assignment25;
	assignment25.Map = R25;
	assignment25.b = w25;

	// The transition label is   t1

	math::matrix<double> R26;
	row = 4;
	col = 4;
	R26.resize(row, col);
	R26(0, 0) = 1.0;
	R26(0, 1) = 0.0;
	R26(0, 2) = 0.0;
	R26(0, 3) = 0.0;
	R26(1, 0) = 0.0;
	R26(1, 1) = 1.0;
	R26(1, 2) = 0.0;
	R26(1, 3) = 0.0;
	R26(2, 0) = 0.0;
	R26(2, 1) = 0.0;
	R26(2, 2) = 1.0;
	R26(2, 3) = 0.0;
	R26(3, 0) = 0.0;
	R26(3, 1) = 0.0;
	R26(3, 2) = 0.0;
	R26(3, 3) = 1.0;
	std::vector<double> w26(row);
	w26[0] = 0.0;
	w26[1] = 0.0;
	w26[2] = 0.0;
	w26[3] = 0.0;

	Assign assignment26;
	assignment26.Map = R26;
	assignment26.b = w26;

	// The transition label is   t2

	math::matrix<double> R27;
	row = 4;
	col = 4;
	R27.resize(row, col);
	R27(0, 0) = 1.0;
	R27(0, 1) = 0.0;
	R27(0, 2) = 0.0;
	R27(0, 3) = 0.0;
	R27(1, 0) = 0.0;
	R27(1, 1) = 1.0;
	R27(1, 2) = 0.0;
	R27(1, 3) = 0.0;
	R27(2, 0) = 0.0;
	R27(2, 1) = 0.0;
	R27(2, 2) = 1.0;
	R27(2, 3) = 0.0;
	R27(3, 0) = 0.0;
	R27(3, 1) = 0.0;
	R27(3, 2) = 0.0;
	R27(3, 3) = 1.0;
	std::vector<double> w27(row);
	w27[0] = 0.0;
	w27[1] = 0.0;
	w27[2] = 0.0;
	w27[3] = 0.0;

	Assign assignment27;
	assignment27.Map = R27;
	assignment27.b = w27;

	// The transition label is   t26

	math::matrix<double> R28;
	row = 4;
	col = 4;
	R28.resize(row, col);
	R28(0, 0) = 1.0;
	R28(0, 1) = 0.0;
	R28(0, 2) = 0.0;
	R28(0, 3) = 0.0;
	R28(1, 0) = 0.0;
	R28(1, 1) = 1.0;
	R28(1, 2) = 0.0;
	R28(1, 3) = 0.0;
	R28(2, 0) = 0.0;
	R28(2, 1) = 0.0;
	R28(2, 2) = 1.0;
	R28(2, 3) = 0.0;
	R28(3, 0) = 0.0;
	R28(3, 1) = 0.0;
	R28(3, 2) = 0.0;
	R28(3, 3) = 1.0;
	std::vector<double> w28(row);
	w28[0] = 0.0;
	w28[1] = 0.0;
	w28[2] = 0.0;
	w28[3] = 0.0;

	Assign assignment28;
	assignment28.Map = R28;
	assignment28.b = w28;

	// The transition label is   t28

	math::matrix<double> R29;
	row = 4;
	col = 4;
	R29.resize(row, col);
	R29(0, 0) = 1.0;
	R29(0, 1) = 0.0;
	R29(0, 2) = 0.0;
	R29(0, 3) = 0.0;
	R29(1, 0) = 0.0;
	R29(1, 1) = 1.0;
	R29(1, 2) = 0.0;
	R29(1, 3) = 0.0;
	R29(2, 0) = 0.0;
	R29(2, 1) = 0.0;
	R29(2, 2) = 1.0;
	R29(2, 3) = 0.0;
	R29(3, 0) = 0.0;
	R29(3, 1) = 0.0;
	R29(3, 2) = 0.0;
	R29(3, 3) = 1.0;
	std::vector<double> w29(row);
	w29[0] = 0.0;
	w29[1] = 0.0;
	w29[2] = 0.0;
	w29[3] = 0.0;

	Assign assignment29;
	assignment29.Map = R29;
	assignment29.b = w29;

	// The transition label is   t27

	math::matrix<double> R30;
	row = 4;
	col = 4;
	R30.resize(row, col);
	R30(0, 0) = 1.0;
	R30(0, 1) = 0.0;
	R30(0, 2) = 0.0;
	R30(0, 3) = 0.0;
	R30(1, 0) = 0.0;
	R30(1, 1) = 1.0;
	R30(1, 2) = 0.0;
	R30(1, 3) = 0.0;
	R30(2, 0) = 0.0;
	R30(2, 1) = 0.0;
	R30(2, 2) = 1.0;
	R30(2, 3) = 0.0;
	R30(3, 0) = 0.0;
	R30(3, 1) = 0.0;
	R30(3, 2) = 0.0;
	R30(3, 3) = 1.0;
	std::vector<double> w30(row);
	w30[0] = 0.0;
	w30[1] = 0.0;
	w30[2] = 0.0;
	w30[3] = 0.0;

	Assign assignment30;
	assignment30.Map = R30;
	assignment30.b = w30;

	// The transition label is   t23

	math::matrix<double> R31;
	row = 4;
	col = 4;
	R31.resize(row, col);
	R31(0, 0) = 1.0;
	R31(0, 1) = 0.0;
	R31(0, 2) = 0.0;
	R31(0, 3) = 0.0;
	R31(1, 0) = 0.0;
	R31(1, 1) = 1.0;
	R31(1, 2) = 0.0;
	R31(1, 3) = 0.0;
	R31(2, 0) = 0.0;
	R31(2, 1) = 0.0;
	R31(2, 2) = 1.0;
	R31(2, 3) = 0.0;
	R31(3, 0) = 0.0;
	R31(3, 1) = 0.0;
	R31(3, 2) = 0.0;
	R31(3, 3) = 1.0;
	std::vector<double> w31(row);
	w31[0] = 0.0;
	w31[1] = 0.0;
	w31[2] = 0.0;
	w31[3] = 0.0;

	Assign assignment31;
	assignment31.Map = R31;
	assignment31.b = w31;

	// The transition label is   t22

	math::matrix<double> R32;
	row = 4;
	col = 4;
	R32.resize(row, col);
	R32(0, 0) = 1.0;
	R32(0, 1) = 0.0;
	R32(0, 2) = 0.0;
	R32(0, 3) = 0.0;
	R32(1, 0) = 0.0;
	R32(1, 1) = 1.0;
	R32(1, 2) = 0.0;
	R32(1, 3) = 0.0;
	R32(2, 0) = 0.0;
	R32(2, 1) = 0.0;
	R32(2, 2) = 1.0;
	R32(2, 3) = 0.0;
	R32(3, 0) = 0.0;
	R32(3, 1) = 0.0;
	R32(3, 2) = 0.0;
	R32(3, 3) = 1.0;
	std::vector<double> w32(row);
	w32[0] = 0.0;
	w32[1] = 0.0;
	w32[2] = 0.0;
	w32[3] = 0.0;

	Assign assignment32;
	assignment32.Map = R32;
	assignment32.b = w32;

	// The transition label is   t25

	math::matrix<double> R33;
	row = 4;
	col = 4;
	R33.resize(row, col);
	R33(0, 0) = 1.0;
	R33(0, 1) = 0.0;
	R33(0, 2) = 0.0;
	R33(0, 3) = 0.0;
	R33(1, 0) = 0.0;
	R33(1, 1) = 1.0;
	R33(1, 2) = 0.0;
	R33(1, 3) = 0.0;
	R33(2, 0) = 0.0;
	R33(2, 1) = 0.0;
	R33(2, 2) = 1.0;
	R33(2, 3) = 0.0;
	R33(3, 0) = 0.0;
	R33(3, 1) = 0.0;
	R33(3, 2) = 0.0;
	R33(3, 3) = 1.0;
	std::vector<double> w33(row);
	w33[0] = 0.0;
	w33[1] = 0.0;
	w33[2] = 0.0;
	w33[3] = 0.0;

	Assign assignment33;
	assignment33.Map = R33;
	assignment33.b = w33;

	// The transition label is   t24

	math::matrix<double> R34;
	row = 4;
	col = 4;
	R34.resize(row, col);
	R34(0, 0) = 1.0;
	R34(0, 1) = 0.0;
	R34(0, 2) = 0.0;
	R34(0, 3) = 0.0;
	R34(1, 0) = 0.0;
	R34(1, 1) = 1.0;
	R34(1, 2) = 0.0;
	R34(1, 3) = 0.0;
	R34(2, 0) = 0.0;
	R34(2, 1) = 0.0;
	R34(2, 2) = 1.0;
	R34(2, 3) = 0.0;
	R34(3, 0) = 0.0;
	R34(3, 1) = 0.0;
	R34(3, 2) = 0.0;
	R34(3, 3) = 1.0;
	std::vector<double> w34(row);
	w34[0] = 0.0;
	w34[1] = 0.0;
	w34[2] = 0.0;
	w34[3] = 0.0;

	Assign assignment34;
	assignment34.Map = R34;
	assignment34.b = w34;

	// The transition label is   t40

	math::matrix<double> R35;
	row = 4;
	col = 4;
	R35.resize(row, col);
	R35(0, 0) = 1.0;
	R35(0, 1) = 0.0;
	R35(0, 2) = 0.0;
	R35(0, 3) = 0.0;
	R35(1, 0) = 0.0;
	R35(1, 1) = 1.0;
	R35(1, 2) = 0.0;
	R35(1, 3) = 0.0;
	R35(2, 0) = 0.0;
	R35(2, 1) = 0.0;
	R35(2, 2) = 1.0;
	R35(2, 3) = 0.0;
	R35(3, 0) = 0.0;
	R35(3, 1) = 0.0;
	R35(3, 2) = 0.0;
	R35(3, 3) = 1.0;
	std::vector<double> w35(row);
	w35[0] = 0.0;
	w35[1] = 0.0;
	w35[2] = 0.0;
	w35[3] = 0.0;

	Assign assignment35;
	assignment35.Map = R35;
	assignment35.b = w35;

	// The transition label is   t42

	math::matrix<double> R36;
	row = 4;
	col = 4;
	R36.resize(row, col);
	R36(0, 0) = 1.0;
	R36(0, 1) = 0.0;
	R36(0, 2) = 0.0;
	R36(0, 3) = 0.0;
	R36(1, 0) = 0.0;
	R36(1, 1) = 1.0;
	R36(1, 2) = 0.0;
	R36(1, 3) = 0.0;
	R36(2, 0) = 0.0;
	R36(2, 1) = 0.0;
	R36(2, 2) = 1.0;
	R36(2, 3) = 0.0;
	R36(3, 0) = 0.0;
	R36(3, 1) = 0.0;
	R36(3, 2) = 0.0;
	R36(3, 3) = 1.0;
	std::vector<double> w36(row);
	w36[0] = 0.0;
	w36[1] = 0.0;
	w36[2] = 0.0;
	w36[3] = 0.0;

	Assign assignment36;
	assignment36.Map = R36;
	assignment36.b = w36;

	// The transition label is   t41

	math::matrix<double> R37;
	row = 4;
	col = 4;
	R37.resize(row, col);
	R37(0, 0) = 1.0;
	R37(0, 1) = 0.0;
	R37(0, 2) = 0.0;
	R37(0, 3) = 0.0;
	R37(1, 0) = 0.0;
	R37(1, 1) = 1.0;
	R37(1, 2) = 0.0;
	R37(1, 3) = 0.0;
	R37(2, 0) = 0.0;
	R37(2, 1) = 0.0;
	R37(2, 2) = 1.0;
	R37(2, 3) = 0.0;
	R37(3, 0) = 0.0;
	R37(3, 1) = 0.0;
	R37(3, 2) = 0.0;
	R37(3, 3) = 1.0;
	std::vector<double> w37(row);
	w37[0] = 0.0;
	w37[1] = 0.0;
	w37[2] = 0.0;
	w37[3] = 0.0;

	Assign assignment37;
	assignment37.Map = R37;
	assignment37.b = w37;

	// The transition label is   t36

	math::matrix<double> R38;
	row = 4;
	col = 4;
	R38.resize(row, col);
	R38(0, 0) = 1.0;
	R38(0, 1) = 0.0;
	R38(0, 2) = 0.0;
	R38(0, 3) = 0.0;
	R38(1, 0) = 0.0;
	R38(1, 1) = 1.0;
	R38(1, 2) = 0.0;
	R38(1, 3) = 0.0;
	R38(2, 0) = 0.0;
	R38(2, 1) = 0.0;
	R38(2, 2) = 1.0;
	R38(2, 3) = 0.0;
	R38(3, 0) = 0.0;
	R38(3, 1) = 0.0;
	R38(3, 2) = 0.0;
	R38(3, 3) = 1.0;
	std::vector<double> w38(row);
	w38[0] = 0.0;
	w38[1] = 0.0;
	w38[2] = 0.0;
	w38[3] = 0.0;

	Assign assignment38;
	assignment38.Map = R38;
	assignment38.b = w38;

	// The transition label is   t37

	math::matrix<double> R39;
	row = 4;
	col = 4;
	R39.resize(row, col);
	R39(0, 0) = 1.0;
	R39(0, 1) = 0.0;
	R39(0, 2) = 0.0;
	R39(0, 3) = 0.0;
	R39(1, 0) = 0.0;
	R39(1, 1) = 1.0;
	R39(1, 2) = 0.0;
	R39(1, 3) = 0.0;
	R39(2, 0) = 0.0;
	R39(2, 1) = 0.0;
	R39(2, 2) = 1.0;
	R39(2, 3) = 0.0;
	R39(3, 0) = 0.0;
	R39(3, 1) = 0.0;
	R39(3, 2) = 0.0;
	R39(3, 3) = 1.0;
	std::vector<double> w39(row);
	w39[0] = 0.0;
	w39[1] = 0.0;
	w39[2] = 0.0;
	w39[3] = 0.0;

	Assign assignment39;
	assignment39.Map = R39;
	assignment39.b = w39;

	// The transition label is   t39

	math::matrix<double> R40;
	row = 4;
	col = 4;
	R40.resize(row, col);
	R40(0, 0) = 1.0;
	R40(0, 1) = 0.0;
	R40(0, 2) = 0.0;
	R40(0, 3) = 0.0;
	R40(1, 0) = 0.0;
	R40(1, 1) = 1.0;
	R40(1, 2) = 0.0;
	R40(1, 3) = 0.0;
	R40(2, 0) = 0.0;
	R40(2, 1) = 0.0;
	R40(2, 2) = 1.0;
	R40(2, 3) = 0.0;
	R40(3, 0) = 0.0;
	R40(3, 1) = 0.0;
	R40(3, 2) = 0.0;
	R40(3, 3) = 1.0;
	std::vector<double> w40(row);
	w40[0] = 0.0;
	w40[1] = 0.0;
	w40[2] = 0.0;
	w40[3] = 0.0;

	Assign assignment40;
	assignment40.Map = R40;
	assignment40.b = w40;

	// The transition label is   t38

	math::matrix<double> R41;
	row = 4;
	col = 4;
	R41.resize(row, col);
	R41(0, 0) = 1.0;
	R41(0, 1) = 0.0;
	R41(0, 2) = 0.0;
	R41(0, 3) = 0.0;
	R41(1, 0) = 0.0;
	R41(1, 1) = 1.0;
	R41(1, 2) = 0.0;
	R41(1, 3) = 0.0;
	R41(2, 0) = 0.0;
	R41(2, 1) = 0.0;
	R41(2, 2) = 1.0;
	R41(2, 3) = 0.0;
	R41(3, 0) = 0.0;
	R41(3, 1) = 0.0;
	R41(3, 2) = 0.0;
	R41(3, 3) = 1.0;
	std::vector<double> w41(row);
	w41[0] = 0.0;
	w41[1] = 0.0;
	w41[2] = 0.0;
	w41[3] = 0.0;

	Assign assignment41;
	assignment41.Map = R41;
	assignment41.b = w41;

	// The transition label is   t58

	math::matrix<double> R42;
	row = 4;
	col = 4;
	R42.resize(row, col);
	R42(0, 0) = 1.0;
	R42(0, 1) = 0.0;
	R42(0, 2) = 0.0;
	R42(0, 3) = 0.0;
	R42(1, 0) = 0.0;
	R42(1, 1) = 1.0;
	R42(1, 2) = 0.0;
	R42(1, 3) = 0.0;
	R42(2, 0) = 0.0;
	R42(2, 1) = 0.0;
	R42(2, 2) = 1.0;
	R42(2, 3) = 0.0;
	R42(3, 0) = 0.0;
	R42(3, 1) = 0.0;
	R42(3, 2) = 0.0;
	R42(3, 3) = 1.0;
	std::vector<double> w42(row);
	w42[0] = 0.0;
	w42[1] = 0.0;
	w42[2] = 0.0;
	w42[3] = 0.0;

	Assign assignment42;
	assignment42.Map = R42;
	assignment42.b = w42;

	// The transition label is   t60

	math::matrix<double> R43;
	row = 4;
	col = 4;
	R43.resize(row, col);
	R43(0, 0) = 1.0;
	R43(0, 1) = 0.0;
	R43(0, 2) = 0.0;
	R43(0, 3) = 0.0;
	R43(1, 0) = 0.0;
	R43(1, 1) = 1.0;
	R43(1, 2) = 0.0;
	R43(1, 3) = 0.0;
	R43(2, 0) = 0.0;
	R43(2, 1) = 0.0;
	R43(2, 2) = 1.0;
	R43(2, 3) = 0.0;
	R43(3, 0) = 0.0;
	R43(3, 1) = 0.0;
	R43(3, 2) = 0.0;
	R43(3, 3) = 1.0;
	std::vector<double> w43(row);
	w43[0] = 0.0;
	w43[1] = 0.0;
	w43[2] = 0.0;
	w43[3] = 0.0;

	Assign assignment43;
	assignment43.Map = R43;
	assignment43.b = w43;

	// The transition label is   t59

	math::matrix<double> R44;
	row = 4;
	col = 4;
	R44.resize(row, col);
	R44(0, 0) = 1.0;
	R44(0, 1) = 0.0;
	R44(0, 2) = 0.0;
	R44(0, 3) = 0.0;
	R44(1, 0) = 0.0;
	R44(1, 1) = 1.0;
	R44(1, 2) = 0.0;
	R44(1, 3) = 0.0;
	R44(2, 0) = 0.0;
	R44(2, 1) = 0.0;
	R44(2, 2) = 1.0;
	R44(2, 3) = 0.0;
	R44(3, 0) = 0.0;
	R44(3, 1) = 0.0;
	R44(3, 2) = 0.0;
	R44(3, 3) = 1.0;
	std::vector<double> w44(row);
	w44[0] = 0.0;
	w44[1] = 0.0;
	w44[2] = 0.0;
	w44[3] = 0.0;

	Assign assignment44;
	assignment44.Map = R44;
	assignment44.b = w44;

	// The transition label is   t55

	math::matrix<double> R45;
	row = 4;
	col = 4;
	R45.resize(row, col);
	R45(0, 0) = 1.0;
	R45(0, 1) = 0.0;
	R45(0, 2) = 0.0;
	R45(0, 3) = 0.0;
	R45(1, 0) = 0.0;
	R45(1, 1) = 1.0;
	R45(1, 2) = 0.0;
	R45(1, 3) = 0.0;
	R45(2, 0) = 0.0;
	R45(2, 1) = 0.0;
	R45(2, 2) = 1.0;
	R45(2, 3) = 0.0;
	R45(3, 0) = 0.0;
	R45(3, 1) = 0.0;
	R45(3, 2) = 0.0;
	R45(3, 3) = 1.0;
	std::vector<double> w45(row);
	w45[0] = 0.0;
	w45[1] = 0.0;
	w45[2] = 0.0;
	w45[3] = 0.0;

	Assign assignment45;
	assignment45.Map = R45;
	assignment45.b = w45;

	// The transition label is   t57

	math::matrix<double> R46;
	row = 4;
	col = 4;
	R46.resize(row, col);
	R46(0, 0) = 1.0;
	R46(0, 1) = 0.0;
	R46(0, 2) = 0.0;
	R46(0, 3) = 0.0;
	R46(1, 0) = 0.0;
	R46(1, 1) = 1.0;
	R46(1, 2) = 0.0;
	R46(1, 3) = 0.0;
	R46(2, 0) = 0.0;
	R46(2, 1) = 0.0;
	R46(2, 2) = 1.0;
	R46(2, 3) = 0.0;
	R46(3, 0) = 0.0;
	R46(3, 1) = 0.0;
	R46(3, 2) = 0.0;
	R46(3, 3) = 1.0;
	std::vector<double> w46(row);
	w46[0] = 0.0;
	w46[1] = 0.0;
	w46[2] = 0.0;
	w46[3] = 0.0;

	Assign assignment46;
	assignment46.Map = R46;
	assignment46.b = w46;

	// The transition label is   t54

	math::matrix<double> R47;
	row = 4;
	col = 4;
	R47.resize(row, col);
	R47(0, 0) = 1.0;
	R47(0, 1) = 0.0;
	R47(0, 2) = 0.0;
	R47(0, 3) = 0.0;
	R47(1, 0) = 0.0;
	R47(1, 1) = 1.0;
	R47(1, 2) = 0.0;
	R47(1, 3) = 0.0;
	R47(2, 0) = 0.0;
	R47(2, 1) = 0.0;
	R47(2, 2) = 1.0;
	R47(2, 3) = 0.0;
	R47(3, 0) = 0.0;
	R47(3, 1) = 0.0;
	R47(3, 2) = 0.0;
	R47(3, 3) = 1.0;
	std::vector<double> w47(row);
	w47[0] = 0.0;
	w47[1] = 0.0;
	w47[2] = 0.0;
	w47[3] = 0.0;

	Assign assignment47;
	assignment47.Map = R47;
	assignment47.b = w47;

	// The transition label is   t56

	math::matrix<double> R48;
	row = 4;
	col = 4;
	R48.resize(row, col);
	R48(0, 0) = 1.0;
	R48(0, 1) = 0.0;
	R48(0, 2) = 0.0;
	R48(0, 3) = 0.0;
	R48(1, 0) = 0.0;
	R48(1, 1) = 1.0;
	R48(1, 2) = 0.0;
	R48(1, 3) = 0.0;
	R48(2, 0) = 0.0;
	R48(2, 1) = 0.0;
	R48(2, 2) = 1.0;
	R48(2, 3) = 0.0;
	R48(3, 0) = 0.0;
	R48(3, 1) = 0.0;
	R48(3, 2) = 0.0;
	R48(3, 3) = 1.0;
	std::vector<double> w48(row);
	w48[0] = 0.0;
	w48[1] = 0.0;
	w48[2] = 0.0;
	w48[3] = 0.0;

	Assign assignment48;
	assignment48.Map = R48;
	assignment48.b = w48;

	// The transition label is   t46

	math::matrix<double> R49;
	row = 4;
	col = 4;
	R49.resize(row, col);
	R49(0, 0) = 1.0;
	R49(0, 1) = 0.0;
	R49(0, 2) = 0.0;
	R49(0, 3) = 0.0;
	R49(1, 0) = 0.0;
	R49(1, 1) = 1.0;
	R49(1, 2) = 0.0;
	R49(1, 3) = 0.0;
	R49(2, 0) = 0.0;
	R49(2, 1) = 0.0;
	R49(2, 2) = 1.0;
	R49(2, 3) = 0.0;
	R49(3, 0) = 0.0;
	R49(3, 1) = 0.0;
	R49(3, 2) = 0.0;
	R49(3, 3) = 1.0;
	std::vector<double> w49(row);
	w49[0] = 0.0;
	w49[1] = 0.0;
	w49[2] = 0.0;
	w49[3] = 0.0;

	Assign assignment49;
	assignment49.Map = R49;
	assignment49.b = w49;

	// The transition label is   t49

	math::matrix<double> R50;
	row = 4;
	col = 4;
	R50.resize(row, col);
	R50(0, 0) = 1.0;
	R50(0, 1) = 0.0;
	R50(0, 2) = 0.0;
	R50(0, 3) = 0.0;
	R50(1, 0) = 0.0;
	R50(1, 1) = 1.0;
	R50(1, 2) = 0.0;
	R50(1, 3) = 0.0;
	R50(2, 0) = 0.0;
	R50(2, 1) = 0.0;
	R50(2, 2) = 1.0;
	R50(2, 3) = 0.0;
	R50(3, 0) = 0.0;
	R50(3, 1) = 0.0;
	R50(3, 2) = 0.0;
	R50(3, 3) = 1.0;
	std::vector<double> w50(row);
	w50[0] = 0.0;
	w50[1] = 0.0;
	w50[2] = 0.0;
	w50[3] = 0.0;

	Assign assignment50;
	assignment50.Map = R50;
	assignment50.b = w50;

	// The transition label is   t48

	math::matrix<double> R51;
	row = 4;
	col = 4;
	R51.resize(row, col);
	R51(0, 0) = 1.0;
	R51(0, 1) = 0.0;
	R51(0, 2) = 0.0;
	R51(0, 3) = 0.0;
	R51(1, 0) = 0.0;
	R51(1, 1) = 1.0;
	R51(1, 2) = 0.0;
	R51(1, 3) = 0.0;
	R51(2, 0) = 0.0;
	R51(2, 1) = 0.0;
	R51(2, 2) = 1.0;
	R51(2, 3) = 0.0;
	R51(3, 0) = 0.0;
	R51(3, 1) = 0.0;
	R51(3, 2) = 0.0;
	R51(3, 3) = 1.0;
	std::vector<double> w51(row);
	w51[0] = 0.0;
	w51[1] = 0.0;
	w51[2] = 0.0;
	w51[3] = 0.0;

	Assign assignment51;
	assignment51.Map = R51;
	assignment51.b = w51;

	// The transition label is   t47

	math::matrix<double> R52;
	row = 4;
	col = 4;
	R52.resize(row, col);
	R52(0, 0) = 1.0;
	R52(0, 1) = 0.0;
	R52(0, 2) = 0.0;
	R52(0, 3) = 0.0;
	R52(1, 0) = 0.0;
	R52(1, 1) = 1.0;
	R52(1, 2) = 0.0;
	R52(1, 3) = 0.0;
	R52(2, 0) = 0.0;
	R52(2, 1) = 0.0;
	R52(2, 2) = 1.0;
	R52(2, 3) = 0.0;
	R52(3, 0) = 0.0;
	R52(3, 1) = 0.0;
	R52(3, 2) = 0.0;
	R52(3, 3) = 1.0;
	std::vector<double> w52(row);
	w52[0] = 0.0;
	w52[1] = 0.0;
	w52[2] = 0.0;
	w52[3] = 0.0;

	Assign assignment52;
	assignment52.Map = R52;
	assignment52.b = w52;

	// The transition label is   t50

	math::matrix<double> R53;
	row = 4;
	col = 4;
	R53.resize(row, col);
	R53(0, 0) = 1.0;
	R53(0, 1) = 0.0;
	R53(0, 2) = 0.0;
	R53(0, 3) = 0.0;
	R53(1, 0) = 0.0;
	R53(1, 1) = 1.0;
	R53(1, 2) = 0.0;
	R53(1, 3) = 0.0;
	R53(2, 0) = 0.0;
	R53(2, 1) = 0.0;
	R53(2, 2) = 1.0;
	R53(2, 3) = 0.0;
	R53(3, 0) = 0.0;
	R53(3, 1) = 0.0;
	R53(3, 2) = 0.0;
	R53(3, 3) = 1.0;
	std::vector<double> w53(row);
	w53[0] = 0.0;
	w53[1] = 0.0;
	w53[2] = 0.0;
	w53[3] = 0.0;

	Assign assignment53;
	assignment53.Map = R53;
	assignment53.b = w53;

	// The transition label is   t51

	math::matrix<double> R54;
	row = 4;
	col = 4;
	R54.resize(row, col);
	R54(0, 0) = 1.0;
	R54(0, 1) = 0.0;
	R54(0, 2) = 0.0;
	R54(0, 3) = 0.0;
	R54(1, 0) = 0.0;
	R54(1, 1) = 1.0;
	R54(1, 2) = 0.0;
	R54(1, 3) = 0.0;
	R54(2, 0) = 0.0;
	R54(2, 1) = 0.0;
	R54(2, 2) = 1.0;
	R54(2, 3) = 0.0;
	R54(3, 0) = 0.0;
	R54(3, 1) = 0.0;
	R54(3, 2) = 0.0;
	R54(3, 3) = 1.0;
	std::vector<double> w54(row);
	w54[0] = 0.0;
	w54[1] = 0.0;
	w54[2] = 0.0;
	w54[3] = 0.0;

	Assign assignment54;
	assignment54.Map = R54;
	assignment54.b = w54;

	// The transition label is   t53

	math::matrix<double> R55;
	row = 4;
	col = 4;
	R55.resize(row, col);
	R55(0, 0) = 1.0;
	R55(0, 1) = 0.0;
	R55(0, 2) = 0.0;
	R55(0, 3) = 0.0;
	R55(1, 0) = 0.0;
	R55(1, 1) = 1.0;
	R55(1, 2) = 0.0;
	R55(1, 3) = 0.0;
	R55(2, 0) = 0.0;
	R55(2, 1) = 0.0;
	R55(2, 2) = 1.0;
	R55(2, 3) = 0.0;
	R55(3, 0) = 0.0;
	R55(3, 1) = 0.0;
	R55(3, 2) = 0.0;
	R55(3, 3) = 1.0;
	std::vector<double> w55(row);
	w55[0] = 0.0;
	w55[1] = 0.0;
	w55[2] = 0.0;
	w55[3] = 0.0;

	Assign assignment55;
	assignment55.Map = R55;
	assignment55.b = w55;

	// The transition label is   t52

	math::matrix<double> R56;
	row = 4;
	col = 4;
	R56.resize(row, col);
	R56(0, 0) = 1.0;
	R56(0, 1) = 0.0;
	R56(0, 2) = 0.0;
	R56(0, 3) = 0.0;
	R56(1, 0) = 0.0;
	R56(1, 1) = 1.0;
	R56(1, 2) = 0.0;
	R56(1, 3) = 0.0;
	R56(2, 0) = 0.0;
	R56(2, 1) = 0.0;
	R56(2, 2) = 1.0;
	R56(2, 3) = 0.0;
	R56(3, 0) = 0.0;
	R56(3, 1) = 0.0;
	R56(3, 2) = 0.0;
	R56(3, 3) = 1.0;
	std::vector<double> w56(row);
	w56[0] = 0.0;
	w56[1] = 0.0;
	w56[2] = 0.0;
	w56[3] = 0.0;

	Assign assignment56;
	assignment56.Map = R56;
	assignment56.b = w56;

	// The transition label is   t43

	math::matrix<double> R57;
	row = 4;
	col = 4;
	R57.resize(row, col);
	R57(0, 0) = 1.0;
	R57(0, 1) = 0.0;
	R57(0, 2) = 0.0;
	R57(0, 3) = 0.0;
	R57(1, 0) = 0.0;
	R57(1, 1) = 1.0;
	R57(1, 2) = 0.0;
	R57(1, 3) = 0.0;
	R57(2, 0) = 0.0;
	R57(2, 1) = 0.0;
	R57(2, 2) = 1.0;
	R57(2, 3) = 0.0;
	R57(3, 0) = 0.0;
	R57(3, 1) = 0.0;
	R57(3, 2) = 0.0;
	R57(3, 3) = 1.0;
	std::vector<double> w57(row);
	w57[0] = 0.0;
	w57[1] = 0.0;
	w57[2] = 0.0;
	w57[3] = 0.0;

	Assign assignment57;
	assignment57.Map = R57;
	assignment57.b = w57;

	// The transition label is   t45

	math::matrix<double> R58;
	row = 4;
	col = 4;
	R58.resize(row, col);
	R58(0, 0) = 1.0;
	R58(0, 1) = 0.0;
	R58(0, 2) = 0.0;
	R58(0, 3) = 0.0;
	R58(1, 0) = 0.0;
	R58(1, 1) = 1.0;
	R58(1, 2) = 0.0;
	R58(1, 3) = 0.0;
	R58(2, 0) = 0.0;
	R58(2, 1) = 0.0;
	R58(2, 2) = 1.0;
	R58(2, 3) = 0.0;
	R58(3, 0) = 0.0;
	R58(3, 1) = 0.0;
	R58(3, 2) = 0.0;
	R58(3, 3) = 1.0;
	std::vector<double> w58(row);
	w58[0] = 0.0;
	w58[1] = 0.0;
	w58[2] = 0.0;
	w58[3] = 0.0;

	Assign assignment58;
	assignment58.Map = R58;
	assignment58.b = w58;

	// The transition label is   t44

	math::matrix<double> R59;
	row = 4;
	col = 4;
	R59.resize(row, col);
	R59(0, 0) = 1.0;
	R59(0, 1) = 0.0;
	R59(0, 2) = 0.0;
	R59(0, 3) = 0.0;
	R59(1, 0) = 0.0;
	R59(1, 1) = 1.0;
	R59(1, 2) = 0.0;
	R59(1, 3) = 0.0;
	R59(2, 0) = 0.0;
	R59(2, 1) = 0.0;
	R59(2, 2) = 1.0;
	R59(2, 3) = 0.0;
	R59(3, 0) = 0.0;
	R59(3, 1) = 0.0;
	R59(3, 2) = 0.0;
	R59(3, 3) = 1.0;
	std::vector<double> w59(row);
	w59[0] = 0.0;
	w59[1] = 0.0;
	w59[2] = 0.0;
	w59[3] = 0.0;

	Assign assignment59;
	assignment59.Map = R59;
	assignment59.b = w59;

	// The transition label is   t72

	math::matrix<double> R60;
	row = 4;
	col = 4;
	R60.resize(row, col);
	R60(0, 0) = 1.0;
	R60(0, 1) = 0.0;
	R60(0, 2) = 0.0;
	R60(0, 3) = 0.0;
	R60(1, 0) = 0.0;
	R60(1, 1) = 1.0;
	R60(1, 2) = 0.0;
	R60(1, 3) = 0.0;
	R60(2, 0) = 0.0;
	R60(2, 1) = 0.0;
	R60(2, 2) = 1.0;
	R60(2, 3) = 0.0;
	R60(3, 0) = 0.0;
	R60(3, 1) = 0.0;
	R60(3, 2) = 0.0;
	R60(3, 3) = 1.0;
	std::vector<double> w60(row);
	w60[0] = 0.0;
	w60[1] = 0.0;
	w60[2] = 0.0;
	w60[3] = 0.0;

	Assign assignment60;
	assignment60.Map = R60;
	assignment60.b = w60;

	// The transition label is   t73

	math::matrix<double> R61;
	row = 4;
	col = 4;
	R61.resize(row, col);
	R61(0, 0) = 1.0;
	R61(0, 1) = 0.0;
	R61(0, 2) = 0.0;
	R61(0, 3) = 0.0;
	R61(1, 0) = 0.0;
	R61(1, 1) = 1.0;
	R61(1, 2) = 0.0;
	R61(1, 3) = 0.0;
	R61(2, 0) = 0.0;
	R61(2, 1) = 0.0;
	R61(2, 2) = 1.0;
	R61(2, 3) = 0.0;
	R61(3, 0) = 0.0;
	R61(3, 1) = 0.0;
	R61(3, 2) = 0.0;
	R61(3, 3) = 1.0;
	std::vector<double> w61(row);
	w61[0] = 0.0;
	w61[1] = 0.0;
	w61[2] = 0.0;
	w61[3] = 0.0;

	Assign assignment61;
	assignment61.Map = R61;
	assignment61.b = w61;

	// The transition label is   t70

	math::matrix<double> R62;
	row = 4;
	col = 4;
	R62.resize(row, col);
	R62(0, 0) = 1.0;
	R62(0, 1) = 0.0;
	R62(0, 2) = 0.0;
	R62(0, 3) = 0.0;
	R62(1, 0) = 0.0;
	R62(1, 1) = 1.0;
	R62(1, 2) = 0.0;
	R62(1, 3) = 0.0;
	R62(2, 0) = 0.0;
	R62(2, 1) = 0.0;
	R62(2, 2) = 1.0;
	R62(2, 3) = 0.0;
	R62(3, 0) = 0.0;
	R62(3, 1) = 0.0;
	R62(3, 2) = 0.0;
	R62(3, 3) = 1.0;
	std::vector<double> w62(row);
	w62[0] = 0.0;
	w62[1] = 0.0;
	w62[2] = 0.0;
	w62[3] = 0.0;

	Assign assignment62;
	assignment62.Map = R62;
	assignment62.b = w62;

	// The transition label is   t69

	math::matrix<double> R63;
	row = 4;
	col = 4;
	R63.resize(row, col);
	R63(0, 0) = 1.0;
	R63(0, 1) = 0.0;
	R63(0, 2) = 0.0;
	R63(0, 3) = 0.0;
	R63(1, 0) = 0.0;
	R63(1, 1) = 1.0;
	R63(1, 2) = 0.0;
	R63(1, 3) = 0.0;
	R63(2, 0) = 0.0;
	R63(2, 1) = 0.0;
	R63(2, 2) = 1.0;
	R63(2, 3) = 0.0;
	R63(3, 0) = 0.0;
	R63(3, 1) = 0.0;
	R63(3, 2) = 0.0;
	R63(3, 3) = 1.0;
	std::vector<double> w63(row);
	w63[0] = 0.0;
	w63[1] = 0.0;
	w63[2] = 0.0;
	w63[3] = 0.0;

	Assign assignment63;
	assignment63.Map = R63;
	assignment63.b = w63;

	// The transition label is   t71

	math::matrix<double> R64;
	row = 4;
	col = 4;
	R64.resize(row, col);
	R64(0, 0) = 1.0;
	R64(0, 1) = 0.0;
	R64(0, 2) = 0.0;
	R64(0, 3) = 0.0;
	R64(1, 0) = 0.0;
	R64(1, 1) = 1.0;
	R64(1, 2) = 0.0;
	R64(1, 3) = 0.0;
	R64(2, 0) = 0.0;
	R64(2, 1) = 0.0;
	R64(2, 2) = 1.0;
	R64(2, 3) = 0.0;
	R64(3, 0) = 0.0;
	R64(3, 1) = 0.0;
	R64(3, 2) = 0.0;
	R64(3, 3) = 1.0;
	std::vector<double> w64(row);
	w64[0] = 0.0;
	w64[1] = 0.0;
	w64[2] = 0.0;
	w64[3] = 0.0;

	Assign assignment64;
	assignment64.Map = R64;
	assignment64.b = w64;

	// The transition label is   t66

	math::matrix<double> R65;
	row = 4;
	col = 4;
	R65.resize(row, col);
	R65(0, 0) = 1.0;
	R65(0, 1) = 0.0;
	R65(0, 2) = 0.0;
	R65(0, 3) = 0.0;
	R65(1, 0) = 0.0;
	R65(1, 1) = 1.0;
	R65(1, 2) = 0.0;
	R65(1, 3) = 0.0;
	R65(2, 0) = 0.0;
	R65(2, 1) = 0.0;
	R65(2, 2) = 1.0;
	R65(2, 3) = 0.0;
	R65(3, 0) = 0.0;
	R65(3, 1) = 0.0;
	R65(3, 2) = 0.0;
	R65(3, 3) = 1.0;
	std::vector<double> w65(row);
	w65[0] = 0.0;
	w65[1] = 0.0;
	w65[2] = 0.0;
	w65[3] = 0.0;

	Assign assignment65;
	assignment65.Map = R65;
	assignment65.b = w65;

	// The transition label is   t67

	math::matrix<double> R66;
	row = 4;
	col = 4;
	R66.resize(row, col);
	R66(0, 0) = 1.0;
	R66(0, 1) = 0.0;
	R66(0, 2) = 0.0;
	R66(0, 3) = 0.0;
	R66(1, 0) = 0.0;
	R66(1, 1) = 1.0;
	R66(1, 2) = 0.0;
	R66(1, 3) = 0.0;
	R66(2, 0) = 0.0;
	R66(2, 1) = 0.0;
	R66(2, 2) = 1.0;
	R66(2, 3) = 0.0;
	R66(3, 0) = 0.0;
	R66(3, 1) = 0.0;
	R66(3, 2) = 0.0;
	R66(3, 3) = 1.0;
	std::vector<double> w66(row);
	w66[0] = 0.0;
	w66[1] = 0.0;
	w66[2] = 0.0;
	w66[3] = 0.0;

	Assign assignment66;
	assignment66.Map = R66;
	assignment66.b = w66;

	// The transition label is   t68

	math::matrix<double> R67;
	row = 4;
	col = 4;
	R67.resize(row, col);
	R67(0, 0) = 1.0;
	R67(0, 1) = 0.0;
	R67(0, 2) = 0.0;
	R67(0, 3) = 0.0;
	R67(1, 0) = 0.0;
	R67(1, 1) = 1.0;
	R67(1, 2) = 0.0;
	R67(1, 3) = 0.0;
	R67(2, 0) = 0.0;
	R67(2, 1) = 0.0;
	R67(2, 2) = 1.0;
	R67(2, 3) = 0.0;
	R67(3, 0) = 0.0;
	R67(3, 1) = 0.0;
	R67(3, 2) = 0.0;
	R67(3, 3) = 1.0;
	std::vector<double> w67(row);
	w67[0] = 0.0;
	w67[1] = 0.0;
	w67[2] = 0.0;
	w67[3] = 0.0;

	Assign assignment67;
	assignment67.Map = R67;
	assignment67.b = w67;

	// The transition label is   t63

	math::matrix<double> R68;
	row = 4;
	col = 4;
	R68.resize(row, col);
	R68(0, 0) = 1.0;
	R68(0, 1) = 0.0;
	R68(0, 2) = 0.0;
	R68(0, 3) = 0.0;
	R68(1, 0) = 0.0;
	R68(1, 1) = 1.0;
	R68(1, 2) = 0.0;
	R68(1, 3) = 0.0;
	R68(2, 0) = 0.0;
	R68(2, 1) = 0.0;
	R68(2, 2) = 1.0;
	R68(2, 3) = 0.0;
	R68(3, 0) = 0.0;
	R68(3, 1) = 0.0;
	R68(3, 2) = 0.0;
	R68(3, 3) = 1.0;
	std::vector<double> w68(row);
	w68[0] = 0.0;
	w68[1] = 0.0;
	w68[2] = 0.0;
	w68[3] = 0.0;

	Assign assignment68;
	assignment68.Map = R68;
	assignment68.b = w68;

	// The transition label is   t65

	math::matrix<double> R69;
	row = 4;
	col = 4;
	R69.resize(row, col);
	R69(0, 0) = 1.0;
	R69(0, 1) = 0.0;
	R69(0, 2) = 0.0;
	R69(0, 3) = 0.0;
	R69(1, 0) = 0.0;
	R69(1, 1) = 1.0;
	R69(1, 2) = 0.0;
	R69(1, 3) = 0.0;
	R69(2, 0) = 0.0;
	R69(2, 1) = 0.0;
	R69(2, 2) = 1.0;
	R69(2, 3) = 0.0;
	R69(3, 0) = 0.0;
	R69(3, 1) = 0.0;
	R69(3, 2) = 0.0;
	R69(3, 3) = 1.0;
	std::vector<double> w69(row);
	w69[0] = 0.0;
	w69[1] = 0.0;
	w69[2] = 0.0;
	w69[3] = 0.0;

	Assign assignment69;
	assignment69.Map = R69;
	assignment69.b = w69;

	// The transition label is   t64

	math::matrix<double> R70;
	row = 4;
	col = 4;
	R70.resize(row, col);
	R70(0, 0) = 1.0;
	R70(0, 1) = 0.0;
	R70(0, 2) = 0.0;
	R70(0, 3) = 0.0;
	R70(1, 0) = 0.0;
	R70(1, 1) = 1.0;
	R70(1, 2) = 0.0;
	R70(1, 3) = 0.0;
	R70(2, 0) = 0.0;
	R70(2, 1) = 0.0;
	R70(2, 2) = 1.0;
	R70(2, 3) = 0.0;
	R70(3, 0) = 0.0;
	R70(3, 1) = 0.0;
	R70(3, 2) = 0.0;
	R70(3, 3) = 1.0;
	std::vector<double> w70(row);
	w70[0] = 0.0;
	w70[1] = 0.0;
	w70[2] = 0.0;
	w70[3] = 0.0;

	Assign assignment70;
	assignment70.Map = R70;
	assignment70.b = w70;

	// The transition label is   t61

	math::matrix<double> R71;
	row = 4;
	col = 4;
	R71.resize(row, col);
	R71(0, 0) = 1.0;
	R71(0, 1) = 0.0;
	R71(0, 2) = 0.0;
	R71(0, 3) = 0.0;
	R71(1, 0) = 0.0;
	R71(1, 1) = 1.0;
	R71(1, 2) = 0.0;
	R71(1, 3) = 0.0;
	R71(2, 0) = 0.0;
	R71(2, 1) = 0.0;
	R71(2, 2) = 1.0;
	R71(2, 3) = 0.0;
	R71(3, 0) = 0.0;
	R71(3, 1) = 0.0;
	R71(3, 2) = 0.0;
	R71(3, 3) = 1.0;
	std::vector<double> w71(row);
	w71[0] = 0.0;
	w71[1] = 0.0;
	w71[2] = 0.0;
	w71[3] = 0.0;

	Assign assignment71;
	assignment71.Map = R71;
	assignment71.b = w71;

	// The transition label is   t62

	math::matrix<double> R72;
	row = 4;
	col = 4;
	R72.resize(row, col);
	R72(0, 0) = 1.0;
	R72(0, 1) = 0.0;
	R72(0, 2) = 0.0;
	R72(0, 3) = 0.0;
	R72(1, 0) = 0.0;
	R72(1, 1) = 1.0;
	R72(1, 2) = 0.0;
	R72(1, 3) = 0.0;
	R72(2, 0) = 0.0;
	R72(2, 1) = 0.0;
	R72(2, 2) = 1.0;
	R72(2, 3) = 0.0;
	R72(3, 0) = 0.0;
	R72(3, 1) = 0.0;
	R72(3, 2) = 0.0;
	R72(3, 3) = 1.0;
	std::vector<double> w72(row);
	w72[0] = 0.0;
	w72[1] = 0.0;
	w72[2] = 0.0;
	w72[3] = 0.0;

	Assign assignment72;
	assignment72.Map = R72;
	assignment72.b = w72;

	initial_polytope_I = polytope::ptr(
			new polytope(ConstraintsMatrixI, boundValueI, boundSignI));

	transition::ptr t1 = transition::ptr(
			new transition(1, "t6", 1, 10, gaurd_polytope0, assignment0));
	transition::ptr t2 = transition::ptr(
			new transition(2, "t8", 1, 2, gaurd_polytope1, assignment1));
	transition::ptr t3 = transition::ptr(
			new transition(3, "t7", 1, 4, gaurd_polytope2, assignment2));
	transition::ptr t4 = transition::ptr(
			new transition(4, "t20", 2, 5, gaurd_polytope3, assignment3));
	transition::ptr t5 = transition::ptr(
			new transition(5, "t21", 2, 3, gaurd_polytope4, assignment4));
	transition::ptr t6 = transition::ptr(
			new transition(6, "t19", 2, 1, gaurd_polytope5, assignment5));
	transition::ptr t7 = transition::ptr(
			new transition(7, "t18", 2, 13, gaurd_polytope6, assignment6));
	transition::ptr t8 = transition::ptr(
			new transition(8, "t9", 4, 1, gaurd_polytope7, assignment7));
	transition::ptr t9 = transition::ptr(
			new transition(9, "t11", 4, 5, gaurd_polytope8, assignment8));
	transition::ptr t10 = transition::ptr(
			new transition(10, "t10", 4, 9, gaurd_polytope9, assignment9));
	transition::ptr t11 = transition::ptr(
			new transition(11, "t15", 5, 4, gaurd_polytope10, assignment10));
	transition::ptr t12 = transition::ptr(
			new transition(12, "t14", 5, 2, gaurd_polytope11, assignment11));
	transition::ptr t13 = transition::ptr(
			new transition(13, "t17", 5, 6, gaurd_polytope12, assignment12));
	transition::ptr t14 = transition::ptr(
			new transition(14, "t16", 5, 7, gaurd_polytope13, assignment13));
	transition::ptr t15 = transition::ptr(
			new transition(15, "t32", 6, 3, gaurd_polytope14, assignment14));
	transition::ptr t16 = transition::ptr(
			new transition(16, "t33", 6, 5, gaurd_polytope15, assignment15));
	transition::ptr t17 = transition::ptr(
			new transition(17, "t34", 6, 8, gaurd_polytope16, assignment16));
	transition::ptr t18 = transition::ptr(
			new transition(18, "t35", 6, 18, gaurd_polytope17, assignment17));
	transition::ptr t19 = transition::ptr(
			new transition(19, "t31", 8, 20, gaurd_polytope18, assignment18));
	transition::ptr t20 = transition::ptr(
			new transition(20, "t30", 8, 7, gaurd_polytope19, assignment19));
	transition::ptr t21 = transition::ptr(
			new transition(21, "t29", 8, 6, gaurd_polytope20, assignment20));
	transition::ptr t22 = transition::ptr(
			new transition(22, "t13", 9, 7, gaurd_polytope21, assignment21));
	transition::ptr t23 = transition::ptr(
			new transition(23, "t12", 9, 4, gaurd_polytope22, assignment22));
	transition::ptr t24 = transition::ptr(
			new transition(24, "t4", 10, 1, gaurd_polytope23, assignment23));
	transition::ptr t25 = transition::ptr(
			new transition(25, "t3", 10, 11, gaurd_polytope24, assignment24));
	transition::ptr t26 = transition::ptr(
			new transition(26, "t5", 10, 13, gaurd_polytope25, assignment25));
	transition::ptr t27 = transition::ptr(
			new transition(27, "t1", 11, 10, gaurd_polytope26, assignment26));
	transition::ptr t28 = transition::ptr(
			new transition(28, "t2", 11, 12, gaurd_polytope27, assignment27));
	transition::ptr t29 = transition::ptr(
			new transition(29, "t26", 12, 11, gaurd_polytope28, assignment28));
	transition::ptr t30 = transition::ptr(
			new transition(30, "t28", 12, 14, gaurd_polytope29, assignment29));
	transition::ptr t31 = transition::ptr(
			new transition(31, "t27", 12, 13, gaurd_polytope30, assignment30));
	transition::ptr t32 = transition::ptr(
			new transition(32, "t23", 13, 10, gaurd_polytope31, assignment31));
	transition::ptr t33 = transition::ptr(
			new transition(33, "t22", 13, 12, gaurd_polytope32, assignment32));
	transition::ptr t34 = transition::ptr(
			new transition(34, "t25", 13, 15, gaurd_polytope33, assignment33));
	transition::ptr t35 = transition::ptr(
			new transition(35, "t24", 13, 2, gaurd_polytope34, assignment34));
	transition::ptr t36 = transition::ptr(
			new transition(36, "t40", 14, 12, gaurd_polytope35, assignment35));
	transition::ptr t37 = transition::ptr(
			new transition(37, "t42", 14, 16, gaurd_polytope36, assignment36));
	transition::ptr t38 = transition::ptr(
			new transition(38, "t41", 14, 15, gaurd_polytope37, assignment37));
	transition::ptr t39 = transition::ptr(
			new transition(39, "t36", 15, 14, gaurd_polytope38, assignment38));
	transition::ptr t40 = transition::ptr(
			new transition(40, "t37", 15, 13, gaurd_polytope39, assignment39));
	transition::ptr t41 = transition::ptr(
			new transition(41, "t39", 15, 17, gaurd_polytope40, assignment40));
	transition::ptr t42 = transition::ptr(
			new transition(42, "t38", 15, 3, gaurd_polytope41, assignment41));
	transition::ptr t43 = transition::ptr(
			new transition(43, "t58", 16, 14, gaurd_polytope42, assignment42));
	transition::ptr t44 = transition::ptr(
			new transition(44, "t60", 16, 21, gaurd_polytope43, assignment43));
	transition::ptr t45 = transition::ptr(
			new transition(45, "t59", 16, 17, gaurd_polytope44, assignment44));
	transition::ptr t46 = transition::ptr(
			new transition(46, "t55", 17, 15, gaurd_polytope45, assignment45));
	transition::ptr t47 = transition::ptr(
			new transition(47, "t57", 17, 22, gaurd_polytope46, assignment46));
	transition::ptr t48 = transition::ptr(
			new transition(48, "t54", 17, 16, gaurd_polytope47, assignment47));
	transition::ptr t49 = transition::ptr(
			new transition(49, "t56", 17, 19, gaurd_polytope48, assignment48));
	transition::ptr t50 = transition::ptr(
			new transition(50, "t46", 18, 19, gaurd_polytope49, assignment49));
	transition::ptr t51 = transition::ptr(
			new transition(51, "t49", 18, 24, gaurd_polytope50, assignment50));
	transition::ptr t52 = transition::ptr(
			new transition(52, "t48", 18, 20, gaurd_polytope51, assignment51));
	transition::ptr t53 = transition::ptr(
			new transition(53, "t47", 18, 6, gaurd_polytope52, assignment52));
	transition::ptr t54 = transition::ptr(
			new transition(54, "t50", 19, 17, gaurd_polytope53, assignment53));
	transition::ptr t55 = transition::ptr(
			new transition(55, "t51", 19, 3, gaurd_polytope54, assignment54));
	transition::ptr t56 = transition::ptr(
			new transition(56, "t53", 19, 23, gaurd_polytope55, assignment55));
	transition::ptr t57 = transition::ptr(
			new transition(57, "t52", 19, 18, gaurd_polytope56, assignment56));
	transition::ptr t58 = transition::ptr(
			new transition(58, "t43", 20, 18, gaurd_polytope57, assignment57));
	transition::ptr t59 = transition::ptr(
			new transition(59, "t45", 20, 25, gaurd_polytope58, assignment58));
	transition::ptr t60 = transition::ptr(
			new transition(60, "t44", 20, 8, gaurd_polytope59, assignment59));
	transition::ptr t61 = transition::ptr(
			new transition(61, "t72", 21, 16, gaurd_polytope60, assignment60));
	transition::ptr t62 = transition::ptr(
			new transition(62, "t73", 21, 22, gaurd_polytope61, assignment61));
	transition::ptr t63 = transition::ptr(
			new transition(63, "t70", 22, 17, gaurd_polytope62, assignment62));
	transition::ptr t64 = transition::ptr(
			new transition(64, "t69", 22, 21, gaurd_polytope63, assignment63));
	transition::ptr t65 = transition::ptr(
			new transition(65, "t71", 22, 23, gaurd_polytope64, assignment64));
	transition::ptr t66 = transition::ptr(
			new transition(66, "t66", 23, 22, gaurd_polytope65, assignment65));
	transition::ptr t67 = transition::ptr(
			new transition(67, "t67", 23, 19, gaurd_polytope66, assignment66));
	transition::ptr t68 = transition::ptr(
			new transition(68, "t68", 23, 24, gaurd_polytope67, assignment67));
	transition::ptr t69 = transition::ptr(
			new transition(69, "t63", 24, 23, gaurd_polytope68, assignment68));
	transition::ptr t70 = transition::ptr(
			new transition(70, "t65", 24, 25, gaurd_polytope69, assignment69));
	transition::ptr t71 = transition::ptr(
			new transition(71, "t64", 24, 18, gaurd_polytope70, assignment70));
	transition::ptr t72 = transition::ptr(
			new transition(72, "t61", 25, 24, gaurd_polytope71, assignment71));
	transition::ptr t73 = transition::ptr(
			new transition(73, "t62", 25, 20, gaurd_polytope72, assignment72));

	std::list<transition::ptr> Out_Going_Trans_fromloc3;

	Out_Going_Trans_fromloc3.push_back(t1);
	Out_Going_Trans_fromloc3.push_back(t2);
	Out_Going_Trans_fromloc3.push_back(t3);
	location::ptr l1 = location::ptr(
			new location(1, "loc3", system_dynamics0, invariant0, true,
					Out_Going_Trans_fromloc3));

	std::list<transition::ptr> Out_Going_Trans_fromloc8;

	Out_Going_Trans_fromloc8.push_back(t4);
	Out_Going_Trans_fromloc8.push_back(t5);
	Out_Going_Trans_fromloc8.push_back(t6);
	Out_Going_Trans_fromloc8.push_back(t7);
	location::ptr l2 = location::ptr(
			new location(2, "loc8", system_dynamics1, invariant1, true,
					Out_Going_Trans_fromloc8));

	std::list<transition::ptr> Out_Going_Trans_fromloc13;

	location::ptr l3 = location::ptr(
			new location(3, "BAD", system_dynamics2, invariant2, true,
					Out_Going_Trans_fromloc13));

	std::list<transition::ptr> Out_Going_Trans_fromloc4;

	Out_Going_Trans_fromloc4.push_back(t8);
	Out_Going_Trans_fromloc4.push_back(t9);
	Out_Going_Trans_fromloc4.push_back(t10);
	location::ptr l4 = location::ptr(
			new location(4, "loc4", system_dynamics3, invariant3, true,
					Out_Going_Trans_fromloc4));

	std::list<transition::ptr> Out_Going_Trans_fromloc7;

	Out_Going_Trans_fromloc7.push_back(t11);
	Out_Going_Trans_fromloc7.push_back(t12);
	Out_Going_Trans_fromloc7.push_back(t13);
	Out_Going_Trans_fromloc7.push_back(t14);
	location::ptr l5 = location::ptr(
			new location(5, "loc7", system_dynamics4, invariant4, true,
					Out_Going_Trans_fromloc7));

	std::list<transition::ptr> Out_Going_Trans_fromloc12;

	Out_Going_Trans_fromloc12.push_back(t15);
	Out_Going_Trans_fromloc12.push_back(t16);
	Out_Going_Trans_fromloc12.push_back(t17);
	Out_Going_Trans_fromloc12.push_back(t18);
	location::ptr l6 = location::ptr(
			new location(6, "loc12", system_dynamics5, invariant5, true,
					Out_Going_Trans_fromloc12));

	std::list<transition::ptr> Out_Going_Trans_fromloc6;

	location::ptr l7 = location::ptr(
			new location(7, "GOOD", system_dynamics6, invariant6, true,
					Out_Going_Trans_fromloc6));

	std::list<transition::ptr> Out_Going_Trans_fromloc11;

	Out_Going_Trans_fromloc11.push_back(t19);
	Out_Going_Trans_fromloc11.push_back(t20);
	Out_Going_Trans_fromloc11.push_back(t21);
	location::ptr l8 = location::ptr(
			new location(8, "loc11", system_dynamics7, invariant7, true,
					Out_Going_Trans_fromloc11));

	std::list<transition::ptr> Out_Going_Trans_fromloc5;

	Out_Going_Trans_fromloc5.push_back(t22);
	Out_Going_Trans_fromloc5.push_back(t23);
	location::ptr l9 = location::ptr(
			new location(9, "loc5", system_dynamics8, invariant8, true,
					Out_Going_Trans_fromloc5));

	std::list<transition::ptr> Out_Going_Trans_fromloc2;

	Out_Going_Trans_fromloc2.push_back(t24);
	Out_Going_Trans_fromloc2.push_back(t25);
	Out_Going_Trans_fromloc2.push_back(t26);
	location::ptr l10 = location::ptr(
			new location(10, "loc2", system_dynamics9, invariant9, true,
					Out_Going_Trans_fromloc2));

	std::list<transition::ptr> Out_Going_Trans_fromloc1;

	Out_Going_Trans_fromloc1.push_back(t27);
	Out_Going_Trans_fromloc1.push_back(t28);
	location::ptr l11 = location::ptr(
			new location(11, "loc1", system_dynamics10, invariant10, true,
					Out_Going_Trans_fromloc1));

	std::list<transition::ptr> Out_Going_Trans_fromloc10;

	Out_Going_Trans_fromloc10.push_back(t29);
	Out_Going_Trans_fromloc10.push_back(t30);
	Out_Going_Trans_fromloc10.push_back(t31);
	location::ptr l12 = location::ptr(
			new location(12, "loc10", system_dynamics11, invariant11, true,
					Out_Going_Trans_fromloc10));

	std::list<transition::ptr> Out_Going_Trans_fromloc9;

	Out_Going_Trans_fromloc9.push_back(t32);
	Out_Going_Trans_fromloc9.push_back(t33);
	Out_Going_Trans_fromloc9.push_back(t34);
	Out_Going_Trans_fromloc9.push_back(t35);
	location::ptr l13 = location::ptr(
			new location(13, "loc9", system_dynamics12, invariant12, true,
					Out_Going_Trans_fromloc9));

	std::list<transition::ptr> Out_Going_Trans_fromloc15;

	Out_Going_Trans_fromloc15.push_back(t36);
	Out_Going_Trans_fromloc15.push_back(t37);
	Out_Going_Trans_fromloc15.push_back(t38);
	location::ptr l14 = location::ptr(
			new location(14, "loc15", system_dynamics13, invariant13, true,
					Out_Going_Trans_fromloc15));

	std::list<transition::ptr> Out_Going_Trans_fromloc14;

	Out_Going_Trans_fromloc14.push_back(t39);
	Out_Going_Trans_fromloc14.push_back(t40);
	Out_Going_Trans_fromloc14.push_back(t41);
	Out_Going_Trans_fromloc14.push_back(t42);
	location::ptr l15 = location::ptr(
			new location(15, "loc14", system_dynamics14, invariant14, true,
					Out_Going_Trans_fromloc14));

	std::list<transition::ptr> Out_Going_Trans_fromloc20;

	Out_Going_Trans_fromloc20.push_back(t43);
	Out_Going_Trans_fromloc20.push_back(t44);
	Out_Going_Trans_fromloc20.push_back(t45);
	location::ptr l16 = location::ptr(
			new location(16, "loc20", system_dynamics15, invariant15, true,
					Out_Going_Trans_fromloc20));

	std::list<transition::ptr> Out_Going_Trans_fromloc19;

	Out_Going_Trans_fromloc19.push_back(t46);
	Out_Going_Trans_fromloc19.push_back(t47);
	Out_Going_Trans_fromloc19.push_back(t48);
	Out_Going_Trans_fromloc19.push_back(t49);
	location::ptr l17 = location::ptr(
			new location(17, "loc19", system_dynamics16, invariant16, true,
					Out_Going_Trans_fromloc19));

	std::list<transition::ptr> Out_Going_Trans_fromloc17;

	Out_Going_Trans_fromloc17.push_back(t50);
	Out_Going_Trans_fromloc17.push_back(t51);
	Out_Going_Trans_fromloc17.push_back(t52);
	Out_Going_Trans_fromloc17.push_back(t53);
	location::ptr l18 = location::ptr(
			new location(18, "loc17", system_dynamics17, invariant17, true,
					Out_Going_Trans_fromloc17));

	std::list<transition::ptr> Out_Going_Trans_fromloc18;

	Out_Going_Trans_fromloc18.push_back(t54);
	Out_Going_Trans_fromloc18.push_back(t55);
	Out_Going_Trans_fromloc18.push_back(t56);
	Out_Going_Trans_fromloc18.push_back(t57);
	location::ptr l19 = location::ptr(
			new location(19, "loc18", system_dynamics18, invariant18, true,
					Out_Going_Trans_fromloc18));

	std::list<transition::ptr> Out_Going_Trans_fromloc16;

	Out_Going_Trans_fromloc16.push_back(t58);
	Out_Going_Trans_fromloc16.push_back(t59);
	Out_Going_Trans_fromloc16.push_back(t60);
	location::ptr l20 = location::ptr(
			new location(20, "loc16", system_dynamics19, invariant19, true,
					Out_Going_Trans_fromloc16));

	std::list<transition::ptr> Out_Going_Trans_fromloc25;

	Out_Going_Trans_fromloc25.push_back(t61);
	Out_Going_Trans_fromloc25.push_back(t62);
	location::ptr l21 = location::ptr(
			new location(21, "loc25", system_dynamics20, invariant20, true,
					Out_Going_Trans_fromloc25));

	std::list<transition::ptr> Out_Going_Trans_fromloc24;

	Out_Going_Trans_fromloc24.push_back(t63);
	Out_Going_Trans_fromloc24.push_back(t64);
	Out_Going_Trans_fromloc24.push_back(t65);
	location::ptr l22 = location::ptr(
			new location(22, "loc24", system_dynamics21, invariant21, true,
					Out_Going_Trans_fromloc24));

	std::list<transition::ptr> Out_Going_Trans_fromloc23;

	Out_Going_Trans_fromloc23.push_back(t66);
	Out_Going_Trans_fromloc23.push_back(t67);
	Out_Going_Trans_fromloc23.push_back(t68);
	location::ptr l23 = location::ptr(
			new location(23, "loc23", system_dynamics22, invariant22, true,
					Out_Going_Trans_fromloc23));

	std::list<transition::ptr> Out_Going_Trans_fromloc22;

	Out_Going_Trans_fromloc22.push_back(t69);
	Out_Going_Trans_fromloc22.push_back(t70);
	Out_Going_Trans_fromloc22.push_back(t71);
	location::ptr l24 = location::ptr(
			new location(24, "loc22", system_dynamics23, invariant23, true,
					Out_Going_Trans_fromloc22));

	std::list<transition::ptr> Out_Going_Trans_fromloc21;

	Out_Going_Trans_fromloc21.push_back(t72);
	Out_Going_Trans_fromloc21.push_back(t73);
	location::ptr l25 = location::ptr(
			new location(25, "loc21", system_dynamics24, invariant24, true,
					Out_Going_Trans_fromloc21));

	int dim = initial_polytope_I->getSystemDimension();
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
	Hybrid_Automata.addLocation(l10);
	Hybrid_Automata.addLocation(l11);
	Hybrid_Automata.addLocation(l12);
	Hybrid_Automata.addLocation(l13);
	Hybrid_Automata.addLocation(l14);
	Hybrid_Automata.addLocation(l15);
	Hybrid_Automata.addLocation(l16);
	Hybrid_Automata.addLocation(l17);
	Hybrid_Automata.addLocation(l18);
	Hybrid_Automata.addLocation(l19);
	Hybrid_Automata.addLocation(l20);
	Hybrid_Automata.addLocation(l21);
	Hybrid_Automata.addLocation(l22);
	Hybrid_Automata.addLocation(l23);
	Hybrid_Automata.addLocation(l24);
	Hybrid_Automata.addLocation(l25);
	Hybrid_Automata.setDimension(dim);

	Hybrid_Automata.insert_to_map("x1", 0);
	Hybrid_Automata.insert_to_map("x2", 1);
	Hybrid_Automata.insert_to_map("v1", 2);
	Hybrid_Automata.insert_to_map("v2", 3);

	unsigned int initial_location_id = 17; //the initial Location ID
	symbolic_states::ptr S; //null_pointer as there is no instantiation
	int transition_id = 0; //initial location no transition taken yet
	initial_state::ptr I = initial_state::ptr(
			new initial_state(initial_location_id, initial_polytope_I, S,
					transition_id));

	init_state_list.push_back(I);

}
