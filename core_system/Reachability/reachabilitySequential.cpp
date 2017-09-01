/*
 * reachabilitySequential.cpp
 *
 *  Created on: 24-Nov-2014
 *      Author: amit
 */


#include "core_system/Reachability/reachabilitySequential.h"


//Reachability Algorithm after optimization of the duplicate support function computation
template_polyhedra::ptr reachabilitySequential(unsigned int boundedTotIteration, Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant, bool isInvariantExist, int lp_solver_type_choosen) {

	int numVectors = ReachParameters.Directions.size1();
	int dimension = Initial->getSystemDimension();
	unsigned int shm_NewTotalIteration = ReachParameters.Iterations; //Shared Variable for resize iterations number on crossing with invariant

	int Min_Or_Max = 2;

	math::matrix<double> MatrixValue; //Shared Matrix for all child thread
	size_type row = numVectors, col = shm_NewTotalIteration;
//	cout << "\nBefore calling InvariantBoundaryCheck"<< "\n";
	if (isInvariantExist == true) { //if invariant exist. Computing
		shm_NewTotalIteration = boundedTotIteration;

	} //End of Invariant Directions
	//cout << "\nNew shm_NewTotalIteration = " << shm_NewTotalIteration << "\n";
	if (shm_NewTotalIteration < 1) {
		template_polyhedra::ptr poly_emptyp;
		return poly_emptyp;
	}

	col = shm_NewTotalIteration; //if invariant exist col will be resized
	MatrixValue.resize(row, col);
	int solver_type = lp_solver_type_choosen;
	lp_solver s_per_thread_I(solver_type), s_per_thread_U(solver_type), s_per_thread_inv(solver_type);
	s_per_thread_I.setMin_Or_Max(2);
	if (!ReachParameters.X0->getIsEmpty()) //set glpk constraints If not an empty polytope
		s_per_thread_I.setConstraints(ReachParameters.X0->getCoeffMatrix(),
				ReachParameters.X0->getColumnVector(), ReachParameters.X0->getInEqualitySign());

	s_per_thread_U.setMin_Or_Max(2);
	if (!SystemDynamics.U->getIsEmpty()) { //empty polytope
		s_per_thread_U.setConstraints(SystemDynamics.U->getCoeffMatrix(),
				SystemDynamics.U->getColumnVector(),
				SystemDynamics.U->getInEqualitySign());
	}
	double res1, result, term2 = 0.0, result1, term1 = 0.0;
	std::vector<double> Btrans_dir, phi_trans_dir, phi_trans_dir1;
	math::matrix<double> B_trans, phi_tau_Transpose;
	if (!SystemDynamics.isEmptyMatrixA) //current_location's SystemDynamics's or ReachParameters
		phi_tau_Transpose = ReachParameters.phi_trans;
	if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
		B_trans = ReachParameters.B_trans;

	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {
		double zIInitial = 0.0, zI = 0.0, zV = 0.0;
		double sVariable, s1Variable;
		std::vector<double> r1Variable(dimension), rVariable(dimension);
		for (int i = 0; i < dimension; i++) {
			rVariable[i] = ReachParameters.Directions(eachDirection, i);
		}
		unsigned int loopIteration = 0;
		double term3, term3a = 0.0, term3b = 0.0, res2, term3c = 0.0;
		sVariable = 0.0; //initialize s0
		//  **************    Omega Function   ********************
		res1 = Initial->computeSupportFunction(rVariable, s_per_thread_I);
	//	cout<<"res1 = "<<res1 <<"\n";

		if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
			phi_tau_Transpose.mult_vector(rVariable, phi_trans_dir);
			term1 = Initial->computeSupportFunction(phi_trans_dir, s_per_thread_I);
		}else if (SystemDynamics.isEmptyMatrixA) { //if A is empty :: {tau.A}' reduces to zero so, e^{tau.A}' reduces to 1
													// so, 1 * rVariable give only rVariable
			term1 = Initial->computeSupportFunction(rVariable, s_per_thread_I);
		}//handling constant dynamics

		if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
			B_trans.mult_vector(rVariable, Btrans_dir);

		if (!SystemDynamics.isEmptyMatrixB && !SystemDynamics.U->getIsEmpty())
			term2 = ReachParameters.time_step * SystemDynamics.U->computeSupportFunction(Btrans_dir,s_per_thread_U);
		term3a = ReachParameters.result_alfa;
		term3b = (double) support_unitball_infnorm(rVariable);
	//	cout<<"term3b = "<<term3b<<"\n";

		if (!SystemDynamics.isEmptyC) {
			term3c = ReachParameters.time_step * dot_product(SystemDynamics.C, rVariable); //Added +tau* sf_C(l) 8/11/2015
		//	cout<<"term3c = "<<term3c<<"\n";
		//	cout<<"dot_product(SystemDynamics.C, rVariable) = "<<dot_product(SystemDynamics.C, rVariable)<<"\n";
		}
		term3 = term3a * term3b;
		res2 = term1 + term2 + term3 + term3c; //term3c Added
		//cout<<"res2 = "<<res2<<"\n";
		/*if (res1<0 && res2 < 0){	//if both negative
			if (res1 < res2)
				zIInitial = res1;
			else
				zIInitial = res2;
		}else{	//otherwise normal convention
			if (res1 > res2)
				zIInitial = res1;
			else
				zIInitial = res2;
		}*/

		if (res1 > res2)
			zIInitial = res1;
		else
			zIInitial = res2;


		//  **************  Omega Function Over  ********************
		MatrixValue(eachDirection, loopIteration) = zIInitial;
		//cout<<"zIInitial = "<< zIInitial<<std::endl;
		loopIteration++;
		for (; loopIteration < shm_NewTotalIteration;) { //Now stopping condition is only "shm_NewTotalIteration"
			double TempOmega;
			//  **************    W_Support Function   ********************
			//	std::vector<double> trans_dir;
			//	B_trans.mult_vector(rVariable, Btrans_dir);
			//res1 = ReachParameters.time_step * SystemDynamics.U->computeSupportFunction(Btrans_dir,	s_per_thread_U, s_per_thread_U, Min_Or_Max);
			result1 = term2;
			double beta = ReachParameters.result_beta;
			//double res_beta = beta * (double) support_unitball_infnorm(rVariable);
			double res_beta = beta * term3b; //Replacing term3b from previous step
			result = result1 + res_beta + term3c; //Added term3c
			zV = result;
			//  **************  W_Support Function Over  ********************
			s1Variable = sVariable + zV;
			//phi_tau_Transpose.mult_vector(rVariable, r1Variable);
			//r1Variable = phi_trans_dir;

			if (SystemDynamics.isEmptyMatrixA) { //Matrix A is empty for constant dynamics
				r1Variable = rVariable;
			} else {
				r1Variable = phi_trans_dir;
			}


			//  **************    Omega Function   ********************
			//res1 = Initial->computeSupportFunction(r1Variable, s_per_thread_I, s_per_thread_U, Min_Or_Max);

			//res1 = term1;
			if (SystemDynamics.isEmptyMatrixA) { //Matrix A is empty for constant dynamics
				//res1 = res1; //A is empty than r1Variable is NOT computable and so is term1. Hence res1 is previous  res1
			} else {
				res1 = term1; //A is not empty than r1Variable is computable and so is term1
			}

			double term3, term3a, res2;
			if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
				phi_tau_Transpose.mult_vector(r1Variable, phi_trans_dir);
				term1 = Initial->computeSupportFunction(phi_trans_dir, s_per_thread_I);
			}

			if (!SystemDynamics.isEmptyMatrixB) { //current_location's SystemDynamics's or ReachParameters
				B_trans.mult_vector(r1Variable, Btrans_dir);
				term2 = ReachParameters.time_step
						* SystemDynamics.U->computeSupportFunction(Btrans_dir,s_per_thread_U);
			}

			term3a = ReachParameters.result_alfa;
			term3b = support_unitball_infnorm(r1Variable);
			//cout<<"term3b = "<<term3b<<"\n";
			if (!SystemDynamics.isEmptyC) {
				term3c = ReachParameters.time_step
						* dot_product(SystemDynamics.C, r1Variable); //Added +tau* sf_C(l) 8/11/2015
			//	cout<<"dot_product(SystemDynamics.C, r1Variable) = "<<dot_product(SystemDynamics.C, r1Variable)<<"\n";
			}
			term3 = term3a * term3b;
			res2 = term1 + term2 + term3 + term3c;



			/*if (res1<0 && res2 < 0){	//if both negative
				if (res1 < res2)
					zI = res1;
				else
					zI = res2;
			}else{	//otherwise normal convention
				if (res1 > res2)
					zI = res1;
				else
					zI = res2;
			}*/


			if (res1 > res2)
				zI = res1;
			else
				zI = res2;





			//  **************  Omega Function Over  ********************
			TempOmega = zI + s1Variable; //Y1
		//	std::cout<<"TempOmega = "<< TempOmega<<std::endl;
			MatrixValue(eachDirection, loopIteration) = TempOmega; //Y1
			rVariable = CopyVector(r1Variable); //source to destination
			sVariable = s1Variable;
			loopIteration++; //for the next Omega-iteration or Time-bound
		} //end of while for each vector
	}

	//todo:: Redundant invariant directional constraints to be removed

	if (isInvariantExist == true) { //if invariant exist. Computing
		math::matrix<double> inv_sfm;
		int num_inv = invariant->getColumnVector().size(); //number of Invariant's constraints
		inv_sfm.resize(num_inv, shm_NewTotalIteration);
		for (int eachInvDirection = 0; eachInvDirection < num_inv;
				eachInvDirection++) {
			for (unsigned int i = 0; i < shm_NewTotalIteration; i++) {
				inv_sfm(eachInvDirection, i) =
						invariant->getColumnVector()[eachInvDirection];
			}
		}
		return template_polyhedra::ptr( new template_polyhedra(MatrixValue, inv_sfm,
				ReachParameters.Directions, invariant->getCoeffMatrix()));
	} else {
		return template_polyhedra::ptr( new template_polyhedra(MatrixValue, ReachParameters.Directions));
	}
}


/*
 * Code AFTER optimising the support function computation
 * Same Code as above but with critical region in InvariantExist block
 */
template_polyhedra::ptr reachabilitySequential_For_Parallel_Iterations(unsigned int boundedTotIteration,
		Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial, ReachabilityParameters& ReachParameters,
		polytope::ptr invariant, bool isInvariantExist, int lp_solver_type_choosen) {

	int numVectors = ReachParameters.Directions.size1();
	int dimension = Initial->getSystemDimension();
	unsigned int shm_NewTotalIteration = ReachParameters.Iterations; //Shared Variable for resize iterations number on crossing with invariant
	int Min_Or_Max = 2;

	math::matrix<double> MatrixValue; //Shared Matrix for all child thread
	size_type row = numVectors, col = shm_NewTotalIteration;
	MatrixValue.resize(row, col);

	if (isInvariantExist == true) { //if invariant exist. Computing
		//shm_NewTotalIteration = InvariantBoundaryCheck(SystemDynamics, Initial, ReachParameters, invariant, lp_solver_type_choosen);
		shm_NewTotalIteration = boundedTotIteration;
	} //End of Invariant Directions
	if (shm_NewTotalIteration == 1) {
		template_polyhedra::ptr poly_emptyp;
		return poly_emptyp;
	}
//cout<<"shm_NewTotalIteration = " <<shm_NewTotalIteration<<std::endl;
	//	cout<<"OK 1 \n";
	int solver_type = lp_solver_type_choosen;

	lp_solver s_per_thread_I(solver_type), s_per_thread_U(solver_type),
			s_per_thread_inv(solver_type);
	s_per_thread_I.setMin_Or_Max(2);
	if (!ReachParameters.X0->getIsEmpty()) //set glpk constraints If not an empty polytope
		s_per_thread_I.setConstraints(ReachParameters.X0->getCoeffMatrix(),
				ReachParameters.X0->getColumnVector(), ReachParameters.X0->getInEqualitySign());

	s_per_thread_U.setMin_Or_Max(2);
	if (SystemDynamics.U->getIsEmpty()) { //empty polytope
		//Polytope is empty so no glpk object constraints to be set
	} else {
		s_per_thread_U.setConstraints(SystemDynamics.U->getCoeffMatrix(),
				SystemDynamics.U->getColumnVector(), SystemDynamics.U->getInEqualitySign());
	}
//	cout<<"OK 2 \n";
	double res1, result, term2, result1, term1;
	std::vector<double> Btrans_dir, phi_trans_dir, phi_trans_dir1;
	math::matrix<double> B_trans, phi_tau_Transpose;
	phi_tau_Transpose = ReachParameters.phi_trans;
	B_trans = ReachParameters.B_trans;

	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {

		double zIInitial = 0.0, zI = 0.0, zV = 0.0;
		double sVariable, s1Variable;
		std::vector<double> r1Variable(dimension), rVariable(dimension);
		for (int i = 0; i < dimension; i++) {
			rVariable[i] = ReachParameters.Directions(eachDirection, i);
		}
		unsigned int loopIteration = 0;
		double term3, term3a, term3b, res2, term3c = 0.0;
		sVariable = 0.0; //initialize s0
//		cout<<"OK 3 \n";
		//  **************    Omega Function   ********************
		res1 = Initial->computeSupportFunction(rVariable, s_per_thread_I);
//		cout<<"OK 4 \n";
		if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
			phi_tau_Transpose.mult_vector(rVariable, phi_trans_dir);
			term1 = Initial->computeSupportFunction(phi_trans_dir, s_per_thread_I);
		}
		if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
			B_trans.mult_vector(rVariable, Btrans_dir);

		if (!SystemDynamics.isEmptyMatrixB && !SystemDynamics.U->getIsEmpty())
			term2 = ReachParameters.time_step * SystemDynamics.U->computeSupportFunction(Btrans_dir,s_per_thread_U);

		term3a = ReachParameters.result_alfa;
		term3b = support_unitball_infnorm(rVariable);
		if (!SystemDynamics.isEmptyC) {
			term3c = ReachParameters.time_step * dot_product(SystemDynamics.C, rVariable); //Added +tau* sf_C(l) 8/11/2015
		}
		term3 = term3a * term3b;
		res2 = term1 + term2 + term3 + term3c;
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
			//	std::vector<double> trans_dir;
			//	B_trans.mult_vector(rVariable, Btrans_dir);
			//res1 = ReachParameters.time_step * SystemDynamics.U->computeSupportFunction(Btrans_dir,	s_per_thread_U, s_per_thread_U, Min_Or_Max);
			result1 = term2;
			double beta = ReachParameters.result_beta;
			//double res_beta = beta	* (double) support_unitball_infnorm(rVariable);
			double res_beta = beta * term3b; //replace from previous steps UnitBall
			result = result1 + res_beta + term3c; //Added term3c
			zV = result;
			//  **************  W_Support Function Over  ********************
			s1Variable = sVariable + zV;

			//phi_tau_Transpose.mult_vector(rVariable, r1Variable);
			r1Variable = phi_trans_dir; //replacement

			//  **************    Omega Function   ********************
			//res1 = Initial->computeSupportFunction(r1Variable, s_per_thread_I, s_per_thread_U, Min_Or_Max);
			res1 = term1; //replacement

			double term3, term3a, res2;

			if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
				phi_tau_Transpose.mult_vector(r1Variable, phi_trans_dir);
				term1 = Initial->computeSupportFunction(phi_trans_dir, s_per_thread_I);
			}

			if (!SystemDynamics.isEmptyMatrixB) { //current_location's SystemDynamics's or ReachParameters
				B_trans.mult_vector(r1Variable, Btrans_dir);
				term2 = ReachParameters.time_step * SystemDynamics.U->computeSupportFunction(Btrans_dir,s_per_thread_U);
			}

			term3a = ReachParameters.result_alfa;
			term3b = support_unitball_infnorm(r1Variable);
			if (!SystemDynamics.isEmptyC) {
				term3c = ReachParameters.time_step * dot_product(SystemDynamics.C, r1Variable); //Added +tau* sf_C(l) 8/11/2015
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
		} //end of for each iteration
	} //end of for all directions
	//cout<<"OK 3 \n";
	template_polyhedra::ptr tpolys;
	tpolys = template_polyhedra::ptr(new template_polyhedra());
#pragma omp critical
	{	/* This critical is used when we call the module PAR_ITER as we are updating the variable ReachParameters
	 and MatrixValue if this function is called from the module SEQ this #pragma will be ignored */
		//Todo:: Redundant invariant directional constraints to be removed
		if (isInvariantExist == true) { //if invariant exist. Computing
			math::matrix<double> inv_sfm;
			int num_inv = invariant->getColumnVector().size(); //number of Invariant's constriants
			inv_sfm.resize(num_inv, shm_NewTotalIteration);
			for (int eachInvDirection = 0; eachInvDirection < num_inv; eachInvDirection++) {
				for (unsigned int i = 0; i < shm_NewTotalIteration; i++)
					inv_sfm(eachInvDirection, i) = invariant->getColumnVector()[eachInvDirection];
			}
			tpolys->setTemplateDirections(ReachParameters.Directions);
			tpolys->setMatrixSupportFunction(MatrixValue);
			tpolys->setInvariantDirections(invariant->getCoeffMatrix());
			tpolys->setMatrix_InvariantBound(inv_sfm);
			//return template_polyhedra(MatrixValue, inv_sfm,ReachParameters.Directions, invariant->getCoeffMatrix());
			//return template_polyhedra(MatrixValue, ReachParameters.TotalDirections);
		} else {
		//	cout<<"OK 4 \n";
			tpolys->setTemplateDirections(ReachParameters.Directions);
		//	cout<<"OK 5 \n";
			tpolys->setMatrixSupportFunction(MatrixValue);
		//	cout<<"OK 6 \n";
			//return template_polyhedra(MatrixValue, ReachParameters.Directions);
		}
	} //end of the Critical Section
	return tpolys;
}
