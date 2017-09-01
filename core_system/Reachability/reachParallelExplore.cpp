/*
 * reachParallelExplore.cpp
 *
 *  Created on: 16-Nov-2014
 *      Author: amit
 */

#include "core_system/Reachability/reachParallelExplore.h"
#include "core_system/Reachability/reachabilitySequential.h"
#include "core_system/Reachability/NewApproach/Partition_BoundingPolytope.h"
#include "Utilities/invariantBoundaryCheck.h"
#include "application/sf_utility.h"
#include "Utilities/testPolytopePlotting.h"
#include "Utilities/Template_Polyhedra.h"
#include "application/All_PP_Definition.h"


//Optimizing Supp_Func Computation:: Implementation of parallel reachability using OMP approach
const template_polyhedra::ptr reachabilityParallel(unsigned int boundedTotIteration, Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen) {
	//typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;

	//omp_set_num_threads(numVectors);
	//#pragma omp parallel
	//int tid;
	int numVectors = ReachParameters.Directions.size1();
	int dimension = Initial->getSystemDimension();
	unsigned int shm_NewTotalIteration = ReachParameters.Iterations; //Shared Variable for resize iterations number on crossing with invariant
	int Min_Or_Max = 2;

	math::matrix<double> MatrixValue; //Shared Matrix for all child thread
	//size_type row = numVectors, col = shm_NewTotalIteration;
	boost::numeric::ublas::matrix<double>::size_type row = numVectors, col =
			shm_NewTotalIteration;

	if (isInvariantExist == true) { //if invariant exist. Computing
		shm_NewTotalIteration = boundedTotIteration;
		//	cout <<"\nInvariant Exists!!!\n";
	} //End of Invariant Directions
	if (shm_NewTotalIteration == 1) {
		template_polyhedra::ptr poly_emptyp;
		return poly_emptyp;
	}

	col = shm_NewTotalIteration;
	MatrixValue.resize(row, col);

	math::matrix<double> B_trans, phi_tau_Transpose;
	phi_tau_Transpose = ReachParameters.phi_trans;
	B_trans = ReachParameters.B_trans;

	int type = lp_solver_type_choosen;
	//omp_set_dynamic(0);	//handles dynamic adjustment of the number of threads within a team
	//omp_set_nested(3);	//enable nested parallelism
	//omp_set_num_threads(10);
	//omp_set_max_active_levels(3);
#pragma omp parallel for //num_threads(2)
	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {

		/*if (eachDirection==0){
			std::cout<<"\nMax Thread in Inner Level = "<< omp_get_num_threads();
			//std::cout<<"\nMax Active Levels = "<<omp_get_max_active_levels();
		}*/
		//std::cout<<"\n eachDirection = "<<eachDirection<<"\n";
		//std::cout<<"\n Inner threadID = "<<omp_get_thread_num()<<"\n";
		//std::cout<<"\n Inner thread omp_get_nested() = "<<omp_get_nested()<<"\n";
		double res1, result, term2, result1, term1;
		std::vector<double> Btrans_dir, phi_trans_dir, phi_trans_dir1;

		lp_solver s_per_thread_I(type), s_per_thread_U(type), s_per_thread_inv(
				type);
		s_per_thread_I.setMin_Or_Max(2);
		if (!Initial->getIsEmpty()) { //set glpk constraints If not an empty polytope
			s_per_thread_I.setConstraints(ReachParameters.X0->getCoeffMatrix(),
					ReachParameters.X0->getColumnVector(),
					ReachParameters.X0->getInEqualitySign());
		}
		s_per_thread_U.setMin_Or_Max(2);
		if (SystemDynamics.U->getIsEmpty()) { //empty polytope
			//Polytope is empty so no glpk object constraints to be set
		} else {
			s_per_thread_U.setConstraints(SystemDynamics.U->getCoeffMatrix(),
					SystemDynamics.U->getColumnVector(),
					SystemDynamics.U->getInEqualitySign());
		}

		double zIInitial = 0.0, zI = 0.0, zV = 0.0;
		double sVariable, s1Variable;
		std::vector<double> r1Variable(dimension), rVariable(dimension);
		for (int i = 0; i < dimension; i++) {
			rVariable[i] = ReachParameters.Directions(eachDirection, i);
		}
		unsigned int loopIteration = 0;
		double term3, term3a, term3b, res2, term3c = 0.0;
		sVariable = 0.0; //initialize s0
		//  **************    Omega Function   ********************
		res1 = Initial->computeSupportFunction(rVariable, s_per_thread_I);
		if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
			phi_tau_Transpose.mult_vector(rVariable, phi_trans_dir);
			term1 = Initial->computeSupportFunction(phi_trans_dir, s_per_thread_I);
		}

		if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
			B_trans.mult_vector(rVariable, Btrans_dir);

		if (!SystemDynamics.isEmptyMatrixB && !SystemDynamics.U->getIsEmpty())
			term2 = ReachParameters.time_step
					* SystemDynamics.U->computeSupportFunction(Btrans_dir, s_per_thread_U);
		term3a = ReachParameters.result_alfa;
		term3b = support_unitball_infnorm(rVariable);
		if (!SystemDynamics.isEmptyC) {
			term3c = ReachParameters.time_step
					* dot_product(SystemDynamics.C, rVariable); //Added +tau* sf_C(l) 8/11/2015
		}
		term3 = term3a * term3b;
		res2 = term1 + term2 + term3 + term3c; //term3c Added
		if (res1 > res2)
			zIInitial = res1;
		else
			zIInitial = res2;
		//  **************  Omega Function Over  ********************

		MatrixValue(eachDirection, loopIteration) = zIInitial;
		loopIteration++;

		for (; loopIteration < shm_NewTotalIteration;) { //Now stopping condition is only "shm_NewTotalIteration"
			double TempOmega;

			//  **************    W_Support Function   ********************
			//	B_trans.mult_vector(rVariable, Btrans_dir);
			//res1 = ReachParameters.time_step * SystemDynamics.U->computeSupportFunction(Btrans_dir,	s_per_thread_U, s_per_thread_U, Min_Or_Max);
			result1 = term2; //replacement--optimizing
			double beta = ReachParameters.result_beta;
			//double res_beta = beta	* (double) support_unitball_infnorm(rVariable);
			double res_beta = beta * term3b; //replace from Previous step UnitBall
			result = result1 + res_beta + term3c; //Added term3c
			zV = result;
			//  **************  W_Support Function Over  ********************
			s1Variable = sVariable + zV;
			//phi_tau_Transpose.mult_vector(rVariable, r1Variable);
			r1Variable = phi_trans_dir; //replacement--optimizing

			//  **************    Omega Function   ********************
			//res1 = Initial->computeSupportFunction(r1Variable, s_per_thread_I, s_per_thread_U, Min_Or_Max);
			res1 = term1; //replacement--optimizing
			double term3, term3a, res2;


			if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
				phi_tau_Transpose.mult_vector(r1Variable, phi_trans_dir);
				term1 = Initial->computeSupportFunction(phi_trans_dir, s_per_thread_I);
			}

			if (!SystemDynamics.isEmptyMatrixB) { //current_location's SystemDynamics's or ReachParameters
				B_trans.mult_vector(r1Variable, Btrans_dir);
				term2 = ReachParameters.time_step
						* SystemDynamics.U->computeSupportFunction(Btrans_dir, s_per_thread_U);
			}


			term3a = ReachParameters.result_alfa;
			term3b = support_unitball_infnorm(r1Variable);
			if (!SystemDynamics.isEmptyC) {
				term3c = ReachParameters.time_step
						* dot_product(SystemDynamics.C, r1Variable); //Added +tau* sf_C(l) 8/11/2015
			}
			term3 = term3a * term3b;
			res2 = term1 + term2 + term3 + term3c;
			if (res1 > res2)
				zI = res1;
			else
				zI = res2;
			//  **************  Omega Function Over  ********************
			TempOmega = zI + s1Variable; //Y1
			MatrixValue(eachDirection, loopIteration) = TempOmega; //Y1
			rVariable = CopyVector(r1Variable); //source to destination
			sVariable = s1Variable;
			loopIteration++; //for the next Omega-iteration or Time-bound
		} //end of while for each vector
	} //end of pragma omp parallel for


	//todo:: Redundant invariant directional constraints to be removed
	if (isInvariantExist == true) { //if invariant exist. Computing
		math::matrix<double> inv_sfm;
		int num_inv = invariant->getColumnVector().size(); //number of Invariant's constriants
		inv_sfm.resize(num_inv, shm_NewTotalIteration);
		for (int eachInvDirection = 0; eachInvDirection < num_inv;
				eachInvDirection++) {
			for (unsigned int i = 0; i < shm_NewTotalIteration; i++) {
				inv_sfm(eachInvDirection, i) =
						invariant->getColumnVector()[eachInvDirection];
			}
		}
		//cout<<"\nAmit"<<MatrixValue.size2()<<"\n";
		return template_polyhedra::ptr(
				new template_polyhedra(MatrixValue, inv_sfm,
						ReachParameters.Directions,
						invariant->getCoeffMatrix()));
		//return template_polyhedra(MatrixValue, ReachParameters.TotalDirections, ReachParameters.Directions, invariant->getCoeffMatrix());
		//return template_polyhedra(MatrixValue, ReachParameters.TotalDirections);
	} else {
		MatrixValue.resize(numVectors, shm_NewTotalIteration, true); //but writing or resizing only upto the NewTotalIteration
		return template_polyhedra::ptr(
				new template_polyhedra(MatrixValue, ReachParameters.Directions));
	}
}












const template_polyhedra::ptr reachParallelExplore(unsigned int NewTotalIteration, Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant, bool invariantExists, int CORES,
		unsigned int Algorithm_Type, int lp_solver_type_choosen) {
	double T = ReachParameters.TimeBound;
	//double original_time_step = ReachParameters.time_step;

	// ******** New Change after boundaryCheckOutside ***************
	if (invariantExists)	//if invariant exist than NewTotalIteration will be greater than 0
		ReachParameters.Iterations = NewTotalIteration; //Actual number of iteration after boundary check evaluation
	// ******** New Change ***************
	unsigned int original_Iterations = ReachParameters.Iterations;
	double newTimeBound = T / CORES;

	ReachParameters.TimeBound = newTimeBound; //correct one.	//ReachParameters.TimeBound = 1;


	// ******** New Change after boundaryCheckOutside ***************
	//ReachParameters.Iterations = ReachParameters.TimeBound / ReachParameters.time_step; //required in Invarian_Boundary_Check
	ReachParameters.Iterations = ReachParameters.Iterations / CORES; //required in Invarian_Boundary_Check
	//Todo::: please handle when the value of this expression is in Fraction
	// ******** New Change ***************

//	cout << "partition_iterations" << ReachParameters.Iterations << "\n";
//ReachParameters.time_step = newTimeBound/ partition_iterations;	//computed for 1st partition
	template_polyhedra::ptr reachability_region;
	reachability_region = template_polyhedra::ptr(new template_polyhedra());
	//std::list<polytope::ptr> initial_polys_list;
//omp_set_nested(1);		//Enabling the Nested Loop in OpenMP Programming
	//std::cout<<"omp_get_nested() = "<<omp_get_nested()<<std::endl;

//use g++ -fopenmp for compiling #pragma omp



#pragma omp parallel for
	for (int i = 0; i < CORES; i++) { //for (double i = 0; i < T; i += newTimeBound) {

		template_polyhedra::ptr Tpoly;
		math::matrix<double> phi, phi_trans, A_inv_phi, y_matrix, y_trans;
		double START_TIME = i * newTimeBound; //first iteration START_TIME = i = 0 which make beta = 0
		//	std::cout << "\nStart_Time = " << START_TIME << "\n";
		SystemDynamics.MatrixA.matrix_exponentiation(phi, START_TIME); //if MatrixA is empty will not perform this function
		phi.transpose(phi_trans); //phi_trans computed

		math::matrix<double> A_inv; //(SystemDynamics.MatrixA.size1(),SystemDynamics.MatrixA.size2());
		A_inv = ReachParameters.A_inv;
		A_inv.multiply(phi, A_inv_phi);
		A_inv_phi.minus(A_inv, y_matrix);
		y_matrix.transpose(y_trans);
//(phi_trans . X0 + y_trans . U)
//		std::cout << "\ncomputing initial object\n";
		supportFunctionProvider::ptr Initial =transMinkPoly::ptr(new transMinkPoly(ReachParameters.X0, SystemDynamics.U,
						phi_trans, y_trans, 1, 0));
//		cout<<"Done TransMinkPoly \n";
		//Calling Sequential algorithm here and later can mix with parallel for direction
		if (Algorithm_Type == TIME_SLICE) {
//			std::cout << "\nFrom Parallel Only Iterations BEFORE\n";
			Tpoly = reachabilitySequential_For_Parallel_Iterations(ReachParameters.Iterations, SystemDynamics, Initial, ReachParameters, invariant,
					invariantExists, lp_solver_type_choosen);
//			std::cout << "\nFrom Parallel Only Iterations AFTER\n";
		}

//		cout<<"Done creating initial partition \n";
#pragma omp critical
		//	if (i != 0)
		reachability_region = reachability_region->union_TemplatePolytope(Tpoly);
	} //end of pragma for loop
	//this may not be required if Average_number_of_times = 1 otherwise must.
	ReachParameters.Iterations = original_Iterations;
	ReachParameters.TimeBound = T; //restore the original timebound for next transition or location
	//GenerateInitialPolytopePlotter(initial_polys_list);
	return reachability_region;
}











