/*
 * AGJHGPU.cpp
 *
 *  Created on: Nov 30, 2016
 *      Author: tomas
 */

#include "AGJHGPU.cuh"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort)
			exit(code);
	}
}

const int blockSize = 1;
const int threadSize = 1;
const int threadCount = blockSize * threadSize;
const int threadQueueSize = 16;

__constant__ CudaReachabilityParameters d_reach_params;

__device__ double support_unitball_infnorm(CudaVector* dir) {
	double sum = 0.0;
	for (unsigned int i = 0; i < dir->size; i++)
		sum += abs(dir->get(i));
	return sum;
}
__device__ void gpuReachabilitySequential(CudaState * s, unsigned int &newTotIters, double result_alfa, double result_beta, CudaPolyhedra * res) {
	/*template_polyhedra::ptr reachabilitySequential(unsigned int boundedTotIteration, Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
	 ReachabilityParameters& ReachParameters, polytope::ptr invariant, bool isInvariantExist, int lp_solver_type_choosen) {*/
//	int numVectors = ReachParameters.Directions.size1();
	int numVectors = d_reach_params.Directions_rows;
	int dimension = s->initialSet.systemDimension;
	unsigned int shm_NewTotalIteration = d_reach_params.Iterations; //Shared Variable for resize iterations number on crossing with invariant

	int Min_Or_Max = 2;

//	CudaMatrix MatrixValue; //Shared Matrix for all child thread moved dows
	size_type row = numVectors, col = shm_NewTotalIteration;
//	cout << "\nBefore calling InvariantBoundaryCheck"<< "\n";
//	if (isInvariantExist == true) { //if invariant exist. Computing
	if (s->location->invariantExists == true) { //if invariant exist. Computing
		shm_NewTotalIteration = d_reach_params.Iterations;

	} //End of Invariant Directions
	  //cout << "\nNew shm_NewTotalIteration = " << shm_NewTotalIteration << "\n";
	if (shm_NewTotalIteration < 1) {
		CudaPolyhedra poly_emptyp;
//		return poly_emptyp;
		res->assign(&poly_emptyp);
		return;
	}

	col = shm_NewTotalIteration; //if invariant exist col will be resized
//	MatrixValue.resize(row, col);
	CudaMatrix MatrixValue(row, col); //Shared Matrix for all child thread
	//int solver_type = lp_solver_type_choosen;
//	lp_solver s_per_thread_I(solver_type), s_per_thread_U(solver_type), s_per_thread_inv(solver_type);
	CudaLpSolver s_per_thread_I, s_per_thread_U, s_per_thread_inv;

	s_per_thread_I.setMin_Or_Max(2);
//	if (!d_reach_params.X0->getIsEmpty()) //set glpk constraints If not an empty polytope
//		s_per_thread_I.setConstraints(ReachParameters.X0->getCoeffMatrix(),
//				ReachParameters.X0->getColumnVector(), ReachParameters.X0->getInEqualitySign());

	if (!s->initialSet.isEmpty) { //set glpk constraints If not an empty polytope
		CudaMatrix *coeffMatrix = &s->initialSet.coeffMatrix;
		s_per_thread_I.setConstraints(coeffMatrix, &s->initialSet.columnVector, s->initialSet.InEqualitySign);
	}
	s_per_thread_U.setMin_Or_Max(2);
//	if (!SystemDynamics->U->getIsEmpty()) { //empty polytope
//		s_per_thread_U.setConstraints(SystemDynamics->U->getCoeffMatrix(),
//				SystemDynamics->U->getColumnVector(),
//				SystemDynamics->U->getInEqualitySign());
//	}
	CudaDynamics* SystemDynamics = &s->location->system_dynamics;
	if (!s->location->system_dynamics.U.isEmpty) { //empty polytope
		s_per_thread_U.setConstraints(&s->location->system_dynamics.U.coeffMatrix, &s->location->system_dynamics.U.columnVector, s->location->system_dynamics.U.InEqualitySign);
	}
	double res1, result, term2 = 0.0, result1, term1 = 0.0;
	CudaVector Btrans_dir, phi_trans_dir, phi_trans_dir1;
//	CudaMatrix B_trans, phi_tau_Transpose;
//	if (!SystemDynamics->isEmptyMatrixA) //current_location's SystemDynamics's or ReachParameters
//		phi_tau_Transpose = ReachParameters.phi_trans;
//	if (!SystemDynamics->isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
//		B_trans = ReachParameters.B_trans;
	CudaMatrix phi_tau_Transpose;
	CudaMatrix B_trans;

	if (!SystemDynamics->isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
//		phi_tau_Transpose.assign(d_reach_params.phi_trans_DATA_rows, d_reach_params.phi_trans_DATA_cols, d_reach_params.phi_trans_DATA); //= d_reach_params.phi_trans;

		if (!SystemDynamics->calculatedphiTrans) {
			SystemDynamics->matrixA.matrix_exponentiation(d_reach_params.time_step, &SystemDynamics->phi_trans);

		}
		phi_tau_Transpose.assign(&SystemDynamics->phi_trans);
	}
	if (!SystemDynamics->isEmptyMatrixB) { //current_location's SystemDynamics's or ReachParameters
		B_trans.assign(&SystemDynamics->calculateBtrans()); // = ReachParameters.B_trans;
	}
	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {
		double zIInitial = 0.0, zI = 0.0, zV = 0.0;
		double sVariable, s1Variable;
		CudaVector r1Variable(dimension), rVariable(dimension);
		for (int i = 0; i < dimension; i++) {
//			rVariable[i] = ReachParameters.Directions(eachDirection, i);
			rVariable.set(i, d_reach_params.Directions_DATA[IDX2F(eachDirection, i, d_reach_params.Directions_rows)]);
		}
		unsigned int loopIteration = 0;
		double term3, term3a = 0.0, term3b = 0.0, res2, term3c = 0.0;
		sVariable = 0.0; //initialize s0
		//  **************    Omega Function   ********************
//		res1 = Initial->computeSupportFunction(rVariable, s_per_thread_I);
		res1 = s->initialSet.computeSupportFunction(&rVariable, &s_per_thread_I);
		//	cout<<"res1 = "<<res1 <<"\n";

		if (!SystemDynamics->isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
			phi_tau_Transpose.mult_vector(&rVariable, &phi_trans_dir);
//			term1 = Initial->computeSupportFunction(phi_trans_dir, s_per_thread_I);
			term1 = s->initialSet.computeSupportFunction(&phi_trans_dir, &s_per_thread_I);
		} else if (SystemDynamics->isEmptyMatrixA) { //if A is empty :: {tau.A}' reduces to zero so, e^{tau.A}' reduces to 1
													 // so, 1 * rVariable give only rVariable
//			term1 = Initial->computeSupportFunction(rVariable, s_per_thread_I);
			term1 = s->initialSet.computeSupportFunction(&rVariable, &s_per_thread_I);
		}													//handling constant dynamics

		if (!SystemDynamics->isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
			B_trans.mult_vector(&rVariable, &Btrans_dir);

		if (!SystemDynamics->isEmptyMatrixB && !SystemDynamics->U.isEmpty)
			term2 = d_reach_params.time_step * SystemDynamics->U.computeSupportFunction(&Btrans_dir, &s_per_thread_U);
		term3a = result_alfa;
		term3b = (double) support_unitball_infnorm(&rVariable);
		//	cout<<"term3b = "<<term3b<<"\n";

		if (!SystemDynamics->isEmptyC) {
//			term3c = d_reach_params.time_step * dot_product(SystemDynamics->C, rVariable); //Added +tau* sf_C(l) 8/11/2015
			term3c = d_reach_params.time_step * CudaVector::dot_product(&SystemDynamics->C, &rVariable); //Added +tau* sf_C(l) 8/11/2015
			//	cout<<"term3c = "<<term3c<<"\n";
			//	cout<<"dot_product(SystemDynamics->C, rVariable) = "<<dot_product(SystemDynamics->C, rVariable)<<"\n";
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
//		MatrixValue(eachDirection, loopIteration) = zIInitial;
		MatrixValue.set(eachDirection, loopIteration, zIInitial);
		//cout<<"zIInitial = "<< zIInitial<<std::endl;
		loopIteration++;

		for (; loopIteration < shm_NewTotalIteration;) { //Now stopping condition is only "shm_NewTotalIteration"

			double TempOmega;
			//  **************    W_Support Function   ********************
			//	std::vector<double> trans_dir;
			//	B_trans.mult_vector(rVariable, Btrans_dir);
			//res1 = ReachParameters.time_step * SystemDynamics->U->computeSupportFunction(Btrans_dir,	s_per_thread_U, s_per_thread_U, Min_Or_Max);
			result1 = term2;
			double beta = result_beta;
			//double res_beta = beta * (double) support_unitball_infnorm(rVariable);
			double res_beta = beta * term3b; //Replacing term3b from previous step
			result = result1 + res_beta + term3c; //Added term3c
			zV = result;
			//  **************  W_Support Function Over  ********************
			s1Variable = sVariable + zV;
			//phi_tau_Transpose.mult_vector(rVariable, r1Variable);
			//r1Variable = phi_trans_dir;

			if (SystemDynamics->isEmptyMatrixA) { //Matrix A is empty for constant dynamics
				r1Variable.assign(&rVariable);
			} else {
				r1Variable.assign(&phi_trans_dir);
			}

			//  **************    Omega Function   ********************
			//res1 = Initial->computeSupportFunction(r1Variable, s_per_thread_I, s_per_thread_U, Min_Or_Max);

			//res1 = term1;
			if (SystemDynamics->isEmptyMatrixA) { //Matrix A is empty for constant dynamics
				//res1 = res1; //A is empty than r1Variable is NOT computable and so is term1. Hence res1 is previous  res1
			} else {
				res1 = term1; //A is not empty than r1Variable is computable and so is term1
			}

			double term3, term3a, res2;
			if (!SystemDynamics->isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
				phi_tau_Transpose.mult_vector(&r1Variable, &phi_trans_dir);
				term1 = s->initialSet.computeSupportFunction(&phi_trans_dir, &s_per_thread_I);
//				term1 = Initial->computeSupportFunction(phi_trans_dir, s_per_thread_I);
			}

			if (!SystemDynamics->isEmptyMatrixB) { //current_location's SystemDynamics's or ReachParameters
				B_trans.mult_vector(&r1Variable, &Btrans_dir);
				term2 = d_reach_params.time_step * SystemDynamics->U.computeSupportFunction(&Btrans_dir, &s_per_thread_U);
			}

			term3a = result_alfa;
			term3b = support_unitball_infnorm(&r1Variable);
			//cout<<"term3b = "<<term3b<<"\n";
			if (!SystemDynamics->isEmptyC) {
				term3c = d_reach_params.time_step * CudaVector::dot_product(&SystemDynamics->C, &r1Variable); //Added +tau* sf_C(l) 8/11/2015
				//	cout<<"dot_product(SystemDynamics->C, r1Variable) = "<<dot_product(SystemDynamics->C, r1Variable)<<"\n";
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
			MatrixValue.set(eachDirection, loopIteration, TempOmega);
//			rVariable = CopyVector(r1Variable); //source to destination
			CudaVector::copyVector(&r1Variable, &rVariable); //source to destination
			sVariable = s1Variable;
			loopIteration++; //for the next Omega-iteration or Time-bound
		} //end of while for each vector
	}

	CudaMatrix Directions;
	Directions.assign(d_reach_params.Directions_rows, d_reach_params.Directions_cols, d_reach_params.Directions_DATA);
	if (s->location->invariantExists) { //if invariant exist. Computing
//		CudaMatrix inv_sfm;
		int num_inv = s->location->invariant.columnVector.size; //number of Invariant's constraints
//		inv_sfm.resize(num_inv, shm_NewTotalIteration);
		CudaMatrix inv_sfm(num_inv, shm_NewTotalIteration);
		for (int eachInvDirection = 0; eachInvDirection < num_inv; eachInvDirection++) {
			for (unsigned int i = 0; i < shm_NewTotalIteration; i++) {
//				inv_sfm(eachInvDirection, i) =
//						invariant->getColumnVector()[eachInvDirection];
				inv_sfm.set(eachInvDirection, i, s->location->invariant.columnVector.get(eachInvDirection));
			}
		}
//		return template_polyhedra::ptr( new template_polyhedra(MatrixValue, inv_sfm,
//				ReachParameters.Directions, invariant->getCoeffMatrix()));
		CudaPolyhedra ret(&MatrixValue, &inv_sfm, &Directions, &s->location->invariant.coeffMatrix);
		res->assign(&ret);
		return;
	} else {
//		return template_polyhedra::ptr( new template_polyhedra(MatrixValue, ReachParameters.Directions));
		CudaPolyhedra ret(&MatrixValue, &Directions);
		res->assign(&ret);
		return;
	}
}

__device__ void cudaInvariantBoundaryCheck(CudaState * s, unsigned int &newTotIters, double result_alfa, double result_beta) {
	/*Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial, ReachabilityParameters& ReachParameters,
	 polytope::ptr invariant, int lp_solver_type_choosen, unsigned int &newTotIters) {*/
	CudaDynamics *SystemDynamics = &s->location->system_dynamics;
	CudaPolytope *Initial = &s->initialSet;
	CudaPolytope *invariant = &s->location->invariant;
	unsigned int shm_NewTotalIteration = d_reach_params.Iterations; //Shared Variable for resize iterations number on crossing with invariant
	int dimension = Initial->systemDimension;
	int Min_Or_Max = 2;
	int numberOfInvariants = invariant->columnVector.size; //total number of Invariant's constraints
	// ******************* Probable race condition variables:: can be  placed inside for-loop *******************
	int foundStart = 0, intersection_start, intersection_end;
	//bool invariantCrossed = false;
	// *************************** For Negative ************************************
	double res1_minus, term2_minus = 0.0, result1_minus, term1_minus = 0.0, result_minus;
	CudaVector Btrans_dir_minus, phi_trans_dir_minus, phi_trans_dir1_minus;
	CudaMatrix phi_tau_Transpose;
	CudaMatrix B_trans;
	if (!SystemDynamics->isEmptyMatrixA) {		//current_location's SystemDynamics's or ReachParameters
//			phi_tau_Transpose .assign(d_reach_params.phi_trans_DATA_rows,d_reach_params.phi_trans_DATA_cols,d_reach_params.phi_trans_DATA);//= d_reach_params.phi_trans;
		if (!SystemDynamics->calculatedphiTrans) {
			SystemDynamics->matrixA.matrix_exponentiation(d_reach_params.time_step, &SystemDynamics->phi_trans);

		}
		phi_tau_Transpose.assign(&SystemDynamics->phi_trans);
	}
	if (!SystemDynamics->isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
		B_trans.assign(&SystemDynamics->calculateBtrans()); // = ReachParameters.B_trans;
	// *******************************************************************************************************************
	//std::vector<int> boundaryIterations(numberOfInvariants, shm_NewTotalIteration); // size(dimension_size,initial_value)
	CudaVector boundaryIterations(numberOfInvariants, shm_NewTotalIteration);

	//cout<<"Test 1.1 \n";
	//int type = lp_solver_type_choosen;
	//#pragma omp parallel for //num_threads(numberOfInvariants)
	for (int eachInvariantDirection = 0; eachInvariantDirection < numberOfInvariants; eachInvariantDirection++) {
		double TempOmega_min;
		double zI_min = 0.0, zV_min = 0.0;
		double sVariable_min, s1Variable_min; //For Minimization of First Vector only

		CudaVector r1Variable_minus(dimension);
//			r1Variable_minus.resize(dimension);

		CudaVector rVariable_minus(dimension);
//			rVariable_minus.resize(dimension);

		for (int i = 0; i < dimension; i++) {
//				rVariable_minus[i] = -1 * invariant->getCoeffMatrix()(eachInvariantDirection, i); //Second vector negative of First vector
			rVariable_minus.set(i, -1 * invariant->coeffMatrix.get(eachInvariantDirection, i));
		}
		//	cout<<"Test 1.1 \n";
		// ******** Lamda Computation for each invariant's directions/constraints **********
		double invariant_SupportFunction;
//			invariant_SupportFunction = invariant->getColumnVector()[eachInvariantDirection];
		invariant_SupportFunction = invariant->columnVector.get(eachInvariantDirection);
		// **********************XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX******************
		CudaLpSolver s_per_thread_I_minus;
		CudaLpSolver s_per_thread_U_minus;
		//	lp_solver s_per_thread_I_minus(type), s_per_thread_U_minus(type);
		// ******************************************* For Negative Direction Starts *******************************************
		s_per_thread_I_minus.setMin_Or_Max(2);
		if (!Initial->isEmpty) //set glpk constraints If not an empty polytope
			//s_per_thread_I_minus.setConstraints(Initial.X0->getCoeffMatrix(), ReachParameters.X0->getColumnVector(), ReachParameters.X0->getInEqualitySign());
			s_per_thread_I_minus.setConstraints(&Initial->coeffMatrix, &Initial->columnVector, Initial->InEqualitySign);

		//cout<<"Test 18\n";
		s_per_thread_U_minus.setMin_Or_Max(2);
		if (!SystemDynamics->U.isEmpty) { //empty polytope
			s_per_thread_U_minus.setConstraints(&SystemDynamics->U.coeffMatrix, &SystemDynamics->U.columnVector, SystemDynamics->U.InEqualitySign);
		}
		// ******************************************* Negative Direction Ends *******************************************
		unsigned int loopIteration = 0;
		sVariable_min = 0.0;

		// ******************************************* For Negative Direction Starts *******************************************
		//zIInitial = Omega_Support(ReachParameters, rVariable_minus, Initial,SystemDynamics, s_per_thread_I_minus, s_per_thread_U_minus, Min_Or_Max);
		double term3_minus, term3a_minus, term3b_minus, res2_minus, term3c_minus = 0.0;
		//  **************    Omega Function   ********************
		res1_minus = Initial->computeSupportFunction(&rVariable_minus, &s_per_thread_I_minus);

		if (!SystemDynamics->isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
//				phi_tau_Transpose.mult_vector(rVariable_minus, phi_trans_dir_minus);
			phi_tau_Transpose.mult_vector(&rVariable_minus, &phi_trans_dir_minus);
			term1_minus = Initial->computeSupportFunction(&phi_trans_dir_minus, &s_per_thread_I_minus);
		} else if (SystemDynamics->isEmptyMatrixA) { //if A is empty :: {tau.A}' reduces to zero so, e^{tau.A}' reduces to 1
			// so, 1 * rVariable give only rVariable

			//some code here

			term1_minus = Initial->computeSupportFunction(&rVariable_minus, &s_per_thread_I_minus);

		}				//handling constant dynamics

		if (!SystemDynamics->isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
//				B_trans.mult_vector(rVariable_minus, Btrans_dir_minus);
			B_trans.mult_vector(&rVariable_minus, &Btrans_dir_minus);
//			if (!SystemDynamics->isEmptyMatrixB && !SystemDynamics->U->getIsEmpty())
		if (!SystemDynamics->isEmptyMatrixB && !SystemDynamics->U.isEmpty)
			term2_minus = d_reach_params.time_step * (SystemDynamics->U.computeSupportFunction(&Btrans_dir_minus, &s_per_thread_U_minus));
		term3a_minus = result_alfa;
		term3b_minus = (double) support_unitball_infnorm(&rVariable_minus);
		//cout<<"Test Before dot_product \n";
		if (!SystemDynamics->isEmptyC) {
			term3c_minus = d_reach_params.time_step * CudaVector::dot_product(&SystemDynamics->C, &rVariable_minus); //Added +tau* sf_C(l) 8/11/2015
			//std::cout<<"No Error here "<<std::endl;
		}
		term3_minus = term3a_minus * term3b_minus;
		res2_minus = term1_minus + term2_minus + term3_minus + term3c_minus;
		if (res1_minus > res2_minus)
			TempOmega_min = res1_minus;	//zIInitial_minus = res1_minus;
		else
			TempOmega_min = res2_minus;	//zIInitial_minus = res2_minus;

		//  **************  Omega Function Over  ********************
		// ******************************************* Negative Direction Ends *******************************************
		loopIteration++;	//needed for returning 1 for Omega_0
		for (; loopIteration < shm_NewTotalIteration;) { //Now stopping condition is only "shm_NewTotalIteration"

			if ((-1 * TempOmega_min) > invariant_SupportFunction) { // Should have been correct
				intersection_end = loopIteration;
//					boundaryIterations[eachInvariantDirection] = loopIteration; //Made Changes here due to circle
				boundaryIterations.set(eachInvariantDirection, loopIteration); //Made Changes here due to circle
				break; //no need to compute reachable set anymore due to out of invariant boundary
			}
			// ******************************************* Intersection Detection Section Ends *******************************************
			// ************************************************  For Negative Direction Starts *******************************************
			//double TempOmega_min;
			double term3_minus, term3a_minus, res2_minus, beta_minus, res_beta_minus;

			//	zV = W_Support(ReachParameters, SystemDynamics, rVariable,s_per_thread_U, Min_Or_Max);
			//  **************    W_Support Function   ********************
			result1_minus = term2_minus;
			beta_minus = result_beta;
			res_beta_minus = beta_minus * term3b_minus; //Replacing term3b from previous step
			result_minus = result1_minus + res_beta_minus + term3c_minus;

			zV_min = result_minus;
			//  **************  W_Support Function Over  ********************
			s1Variable_min = sVariable_min + zV_min;
			if (SystemDynamics->isEmptyMatrixA) { //Matrix A is empty for constant dynamics
				r1Variable_minus.assign(&rVariable_minus);
			} else {
				r1Variable_minus.assign(&phi_trans_dir_minus);
			}

			//zI = Omega_Support(ReachParameters, r1Variable, Initial,SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
			//  **************    Omega Function   ********************

			res1_minus = term1_minus;

			if (!SystemDynamics->isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
//					phi_tau_Transpose.mult_vector(r1Variable_minus, phi_trans_dir_minus);
				phi_tau_Transpose.mult_vector(&r1Variable_minus, &phi_trans_dir_minus);
//					term1_minus = Initial->computeSupportFunction(phi_trans_dir_minus, s_per_thread_I_minus);
				term1_minus = Initial->computeSupportFunction(&phi_trans_dir_minus, &s_per_thread_I_minus);
			}
			if (!SystemDynamics->isEmptyMatrixB) { //current_location's SystemDynamics's or ReachParameters
				B_trans.mult_vector(&r1Variable_minus, &Btrans_dir_minus);
				term2_minus = d_reach_params.time_step * SystemDynamics->U.computeSupportFunction(&Btrans_dir_minus, &s_per_thread_U_minus);
			}
			term3a_minus = result_alfa;
			term3b_minus = support_unitball_infnorm(&r1Variable_minus);
			//std::cout<<"No Error here "<<loopIteration<<std::endl;

			if (!SystemDynamics->isEmptyC) {
				//std::cout<<"SystemDynamics->C.size() = "<<SystemDynamics->C.size()<<std::endl;
				//std::cout<<"rVariable_minus.size() = "<<rVariable_minus.size()<<std::endl;
				term3c_minus = d_reach_params.time_step * CudaVector::dot_product(&SystemDynamics->C, &rVariable_minus); //Added +tau* sf_C(l) 8/11/2015
			}
			term3_minus = term3a_minus * term3b_minus;
			res2_minus = term1_minus + term2_minus + term3_minus + term3c_minus;
			if (res1_minus > res2_minus)
				zI_min = res1_minus;
			else
				zI_min = res2_minus;
			//  **************  Omega Function Over  ********************
			TempOmega_min = zI_min + s1Variable_min; //Y1
			//		cout<<"TempOmega_min = "<<TempOmega_min<<"\t";
			// ************************************************   Negative Direction Ends *******************************************
			loopIteration++; // Placed here as Omega_0 and Omega_1 is computed so loopIteration value == 2

//				rVariable_minus = CopyVector(r1Variable_minus);
			rVariable_minus.assign(&r1Variable_minus);
			sVariable_min = s1Variable_min;

		} //end of iterations
		  //	if (invariantCrossed)	//We need to check all invariants as we are taking the min of all the directions
		  //		break;//boundary crossed so no need to check other invariant's directions
	} //end of parallel for each Iterations or Time-Bound
	  // At the end of the For-Loop or all invariant_Directions we have boundaryIterations vector with the different limit to stop iterations
	  // The Minimum(boundaryIterations) will be the final Boundary Limit for the number of iterations
//		unsigned int min_Total_Iteration;
//		min_Total_Iteration = *min_element(boundaryIterations.begin(), boundaryIterations.end()); // - 1 ;//excluding the last outside Omega
//		newTotIters = min_Total_Iteration;

	newTotIters = boundaryIterations.minimalElement();
}
__device__ void getInitialSet(double START_TIME, CudaDynamics * systemDynamics, CudaPolytope * x0, CudaTransMinkPoly*res) {

	CudaMatrix phi, phi_trans;

	if ((systemDynamics->isEmptyMatrixB || systemDynamics->U.isEmpty) && systemDynamics->isEmptyC) {	//both B and C is empty, so we have x'(t) = Ax(t)
		//cout <<"Matrix B and Vector C is Empty!!, Dynamics is x'(t) = Ax(t)\n";
		systemDynamics->matrixA.matrix_exponentiation(START_TIME, &phi); //if MatrixA is empty will not perform this function
		phi.transpose(&phi_trans); //phi_trans computed
		//	std::cout << "\ncomputing initial object\n";
		CudaTransMinkPoly Initial = CudaTransMinkPoly(x0, &systemDynamics->U, &phi_trans, &phi_trans, 0, 0);		//when  Bu(t) and C are empty
//		return Initial;
		res->assign(&Initial);
		return;
	}
	printf("Dynamics is NOT of type --:  x'(t) = Ax(t)\n");
	CudaMatrix A_inv_phi, y_matrix, y_trans;
	systemDynamics->matrixA.matrix_exponentiation(START_TIME, &phi); //if MatrixA is empty will not perform this function
	phi.transpose(&phi_trans); //phi_trans computed	//cout << "phi_trans" << phi_trans << std::endl;
//	CudaMatrix A_inv(d_reach_params.A_inv_DATA_rows,d_reach_params.A_inv_DATA_cols,d_reach_params.A_inv_DATA); //(SystemDynamics->MatrixA.size1(),SystemDynamics->MatrixA.size2());
	CudaMatrix A_inv;
	A_inv.assign(d_reach_params.A_inv_DATA_rows, d_reach_params.A_inv_DATA_cols, d_reach_params.A_inv_DATA); //(SystemDynamics->MatrixA.size1(),SystemDynamics->MatrixA.size2());;

	CudaMatrix::times(&A_inv, &phi, &A_inv_phi);	//cout<<"A_inv_phi = "<<A_inv_phi<<std::endl;
	CudaMatrix::minus(&A_inv_phi, &A_inv, &y_matrix);	//cout << "y_matrix = " << y_matrix << std::endl;
	y_matrix.transpose(&y_trans);
	//	std::cout << "\ncomputing initial object\n";

	if (systemDynamics->isEmptyC) {
		printf("C is Empty\n");
		CudaTransMinkPoly Initial = CudaTransMinkPoly(x0, &systemDynamics->U, &phi_trans, &y_trans, 1, 0);	//when only C is empty
//		return Initial;
		res->assign(&Initial);
		return;
	} else if (!systemDynamics->isEmptyC) {
		printf("C is NOT Empty\n");
		CudaTransMinkPoly Initial = CudaTransMinkPoly(x0, &systemDynamics->U, &systemDynamics->C, &phi_trans, &y_trans, 1, 0);
//		return Initial;
		res->assign(&Initial);
		return;
	}
}

__device__ double invariantCrossingCheck1(double START_TIME, double time_step, double time_horizon, CudaPolytope * Initial, CudaPolytope *invariant, CudaDynamics *SystemDynamics) {

	// ************************************************ object creations ************************************************
	CudaLpSolver neg_lp;

	if (!Initial->isEmpty) //set glpk constraints If not an empty polytope
//		X0 is the initial set  based on the code in postc
	//neg_lp.setConstraints(ReachParameters.X0->getCoeffMatrix(), ReachParameters.X0->getColumnVector(), ReachParameters.X0->getInEqualitySign());

		neg_lp.setConstraints(&Initial->coeffMatrix, &Initial->columnVector, Initial->InEqualitySign);
	//cout<<"problem here\n";
	// ************************************************ object creations ************************************************
	int dimension = Initial->systemDimension;
	int numberOfInvariants = invariant->columnVector.size;	//total number of Invariant's constraints
	double neg_val;
	bool flag;

	for (double t = START_TIME; t <= time_horizon; t = t + time_step) {		//cout<<"test a\n";
		flag = false;

		CudaTransMinkPoly p;
		getInitialSet(t, SystemDynamics, Initial, &p);	//initial set at time = t
		double invBound;
		CudaVector neg_dir(dimension);
		for (int eachInvariant = 0; eachInvariant < numberOfInvariants; eachInvariant++) {
			for (int j = 0; j < dimension; j++)
				neg_dir.set(j, -1 * invariant->coeffMatrix.get(eachInvariant, j));
			neg_val = -1 * p.computeSupportFunction(&neg_dir, &neg_lp);
			invBound = invariant->columnVector.get(eachInvariant);	//if this does not works use LP solver
			if (neg_val > invBound) {
				flag = true; //crossed the invariant bound
				break;	//return t;
			}
		}
		/*
		 //Not required at runtime only for Printing and Debugging
		 polytope::ptr bound_poly;
		 bound_poly = getBoundedConvexSet(p,ReachParameters);	//creating boundedPolytope over template directions using convex set

		 math::matrix<double> vertices_list = bound_poly->get_2dVertices(0, 1); // 1 and 2 are the output variables
		 // ------------- Printing the vertices on the Output File -------------
		 for (int i = 0; i < vertices_list.size1(); i++) {
		 for (int j = 0; j < vertices_list.size2(); j++) {
		 outFile << vertices_list(i, j) << " ";
		 }
		 outFile << std::endl;
		 }
		 outFile << std::endl; // 1 gap after each polytope plotted
		 // ------------- Printing the vertices on the Output File -------------
		 */
		if (flag) {	//this condition is created only to be able to print the crossing convex set
			//outFile.close();//Not required at runtime only for Printing and Debugging
			return t;
		}
	}
	//outFile.close();//Not required at runtime only for Printing and Debugging
	return time_horizon; //crossed time-horizon without crossing the invariant bound or didnot find intersecting with invariant face
}

__device__ void cudaJumpInvariantBoundaryCheck(CudaState * s, unsigned int &newTotIters) {
	/*Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
	 ReachabilityParameters& ReachParameters, polytope::ptr invariant, int lp_solver_type_choosen, unsigned int &newTotIters*/
	CudaDynamics *SystemDynamics = &s->location->system_dynamics;
	CudaPolytope *Initial = &s->initialSet;
	CudaPolytope *invariant = &s->location->invariant;

	//-- get the time-horizon and time-step and discritization_factor=10 times
	// compute the iteration using this new time-step
	double widening_Factor = 5;	//10 times
	double time_step = d_reach_params.time_step * widening_Factor;
	double time_horizon = d_reach_params.TimeBound;
	double start_time = 0, time_crossed, next_search_time, actual_time_bound;	//last time when it is still inside //cout<<"test 1\n";

	//time_crossed = invariantCrossingCheck(start_time, time_step, time_horizon, Initial, ReachParameters, invariant, SystemDynamics, lp_solver_type_choosen);
	time_crossed = invariantCrossingCheck1(start_time, time_step, time_horizon, Initial, invariant, SystemDynamics);
	//cout<<"First coarse time-step completed at time = "<< time_crossed<<"\n";
	if (time_crossed == time_horizon) {	//Search failed to find a crossing. Flowpipe completely inside Invariant
		newTotIters = d_reach_params.Iterations;	//orginal iteration
		//cout<<"Invariant Condition completed at coarse time-step\n";
	} else {	//if (time_crossed < time_horizon) meaning it has just crossed with widening_Factor but now to get precise crossing-time we use fine-time-step
		//cout<<"Fine time-step testing invoked!!!\n";
		next_search_time = (time_crossed - time_step) + d_reach_params.time_step;
		time_step = d_reach_params.time_step;	//setting the original time-step as fine-time-step

		actual_time_bound = invariantCrossingCheck1(next_search_time, time_step, time_crossed, Initial, invariant, SystemDynamics);
		//Debugging the correct time
		//Since I am following '0-indexing' throughout. So it must return the value/index N as its indices are (0...(N-1))
		//since the function invariantCrossingCheck1(...) return crossing time bound 'N' starting from time at '0'
		//so, we should at one more time-step to return the size as 'N+1' for array indexing
		//For eg., return 1 if only Omega_0 is inside, 2 if Omega_0 and Omega_1 are inside
		//::Debugging CONCLUDED that math::ceil() function rounds up when values are >=5 to its ceiling. which creates problem
		//But if we do not use math::ceil() then normal division on double data type gives in-correct result in g++/gcc.
		//So we use math::ceil() as an extra iterations or over-approximation returned is acceptable in our case.
		//---
		//cout << "actual_time_bound = " << actual_time_bound << " ReachParameters.time_step = " << ReachParameters.time_step;
		//cout << "  actual_time_bound = " << (actual_time_bound + time_step);
		newTotIters = ceil((actual_time_bound + time_step) / d_reach_params.time_step);
		//cout << "  Iterations = " << newTotIters << "\n";
	}
}
__device__ double CudaComputeBeta(CudaState * s) {
	double tau = d_reach_params.time_step;
	CudaDynamics* SysD = &s->location->system_dynamics;
	double norm_A = 0.0, result;
	unsigned int dim_for_Max_norm;
	if (!SysD->isEmptyMatrixA) { //if Not Empty
		norm_A = SysD->matrixA.norm_inf();
	}
	CudaMatrix Btrans;
	double V_max_norm = 0.0;

	if (!SysD->isEmptyMatrixB) { //if NOT Empty
		SysD->matrixB.transpose(&Btrans);
		dim_for_Max_norm = SysD->matrixB.rows;	//dimension for computing Max_Norm(V): V=(B)29x6 . (u)6x1 = (dim of V)29x1
		CudaTransMinkPoly *Vptr = new CudaTransMinkPoly(&SysD->U, &Btrans);
		V_max_norm = Vptr->max_norm(dim_for_Max_norm);
	}

	if (SysD->isEmptyMatrixA) { //if A is Empty
		result = 0;	//norm_A will be zero and which is common term
	} else {
		result = (exp(tau * norm_A) - 1 - tau * norm_A) * (V_max_norm / norm_A);
	}
	//result = (exp(tau * norm_A) - 1 - tau * norm_A) * (V_max_norm / norm_A);
	//	cout<<"\nBeta = "<<(double)result<<endl;
	return result;

}
__device__ double CudaComputeAlfa(CudaState * s) {
	double tau = d_reach_params.time_step;
	double norm_A = 0.0;
	double result;
	CudaPolytope *continuous_initial_polytope = &s->initialSet;
	CudaDynamics *system_dynamics = &s->location->system_dynamics;

	unsigned int dim_for_Max_norm = 0;
	double V_max_norm = 0.0, I_max_norm = 0.0;
	if (!system_dynamics->isEmptyMatrixA) { //if Not Empty
		norm_A = system_dynamics->matrixA.norm_inf();
	}

	dim_for_Max_norm = continuous_initial_polytope->systemDimension;
	I_max_norm = continuous_initial_polytope->max_norm(dim_for_Max_norm); //R_X_o ie max_norm of the Initial polytope

	CudaMatrix Btrans;
	if (!system_dynamics->isEmptyMatrixB) { //if NOT Empty
		system_dynamics->matrixB.transpose(&Btrans);

		CudaTransMinkPoly Vptr(&system_dynamics->U, &Btrans);

		dim_for_Max_norm = system_dynamics->matrixB.rows;
		V_max_norm = Vptr.max_norm(dim_for_Max_norm);

	}

	if (system_dynamics->isEmptyMatrixA) {
		result = 0;
	} else {

		result = (exp(tau * norm_A) - 1 - tau * norm_A) * (I_max_norm + (V_max_norm / norm_A));
	}
	return result;

}

//std::list<initial_state::ptr> agjh::postD(symbolic_states::ptr symb)
__device__ void cudaPostD(CudaState * symb, CudaLocation * d_location, CudaVectorTemplate<CudaState*> *ret) {
//	template_polyhedra::ptr reach_region= symb->getContinuousSetptr();
	CudaPolyhedra * reach_region = &symb->continuousSet;
//	int locId = *(symb->getDiscreteSet().getDiscreteElements().begin());
	int locId = symb->discreteSet.discrete_elements[0];
//	location::ptr current_location = H.getLocation(locId);
	CudaLocation* current_location = &d_location[locId];
//	CudaVectorTemplate<CudaState>res(50);
	ret->size = 50;
	ret->deAlocData();
	ret->alocData();

	int retCount = 0;

	if (reach_region->getTotalIterations() != 0) { //computed reach_region is empty && optimize transition BreadthLevel-wise
		//	cout<<"1\n";
//		for (std::list<transition::ptr>::iterator t = current_location->getOut_Going_Transitions().begin();
//				t != current_location->getOut_Going_Transitions().end(); t++) { // get each destination_location_id and push into the pwl.waiting_list
		for (int i = 0; i < current_location->Out_Going_Transitions_size; i++) {
//			int transition_id = (*t)->getTransitionId();

			CudaTransition *t = &current_location->Out_Going_Transitions[i];

//			CudaLocation current_destination;
			CudaAssign current_assignment;
			CudaPolytope gaurd_polytope;
			gaurd_polytope.assign(current_location->Out_Going_Transitions[i].getGaurd());
			//std::list < template_polyhedra::ptr > intersected_polyhedra;
			CudaPolytope intersectedRegion;	//created two objects here
			//CudaDiscreteSet ds;
//			current_destination = H.getLocation((*t)->getDestination_Location_Id());
			CudaLocation * current_destination = t->destination_location;

			int locName = current_destination->name;

			//cout<<"2\n";
			CudaVectorTemplate<CudaPolytope*> polys;
			//intersected_polyhedra = reach_region->polys_intersectionSequential(gaurd_polytope, lp_solver_type_choosen); //, intersection_start_point);
			reach_region->flowpipe_intersectionSequential(&gaurd_polytope, &polys);

			//cout<<"3\n";
			// to make is even procedure with Sequential procedure.... so intersection is done first and then decide to skip this loc
			if ((locName == CudaLocation::BAD) || (locName == CudaLocation::GOOD) || (locName == CudaLocation::FINAL) || (locName == CudaLocation::UNSAVE)) {
				continue;			//do not push into the waitingList
				printf("skiping %d\n", locName);
			}

			current_assignment.assign(&t->assignT);
			// *** interesected_polyhedra included with invariant_directions also ******
			//	cout<<"size = "<< intersected_polyhedra.size();

//			int destination_locID = (*t)->getDestination_Location_Id();
			int destination_locID = t->destination_location->id;

			//int destination_locID =t.destination_location_id;
			//	ds.insert_element(destination_locID);

//			for (std::list<polytope::ptr>::iterator i = polys.begin(); i != polys.end(); i++) {

			for (int i = 0; i < polys.size; i++) {


				intersectedRegion.assign(polys.get(i));
				//intersectedRegion = (*i)->getTemplate_approx(lp_solver_type_choosen);
				//Returns a single over-approximated polytope from the list of intersected polytopes
				//	GeneratePolytopePlotter(intersectedRegion);
				CudaPolytope newShiftedPolytope;				//created an object here
				CudaPolytope newPolytope;
				intersectedRegion.GetPolytope_Intersection(&gaurd_polytope, &newPolytope);				//Retuns the intersected region as a single newpolytope. **** with added directions
				//newShiftedPolytope = post_assign_exact(newPolytope, current_assignment.Map, current_assignment.b);//initial_polytope_I = post_assign_exact(newPolytope, R, w);

				CudaMatrix test(current_assignment.Map.rows, current_assignment.Map.cols);
				if (current_assignment.Map.inverse(&test))	//invertible?
						{

					//std::cout << "Exact Post Assignment\n";
					CudaPolytope::post_assign_exact(&newPolytope, &current_assignment.Map, &current_assignment.b, &newShiftedPolytope);

						} else {
					//std::cout << "Approximate Post Assignment\n";
//					newShiftedPolytope = CudaPolytope::post_assign_approx_deterministic(
//							newPolytope, current_assignment.Map,
//							current_assignment.b, reach_parameters.Directions,
//							lp_solver_type_choosen);
					CudaMatrix Directions;
					Directions.assign(d_reach_params.Directions_rows, d_reach_params.Directions_cols, d_reach_params.Directions_DATA);
					CudaPolytope::post_assign_approx_deterministic(&newPolytope, &current_assignment.Map, &current_assignment.b, &Directions, &newShiftedPolytope);
				}

				//	newShiftedPolytope->print2file(newInitSet,0,1); //printing the New Initial Set

				CudaState* newState = new CudaState(destination_locID, &newShiftedPolytope);

				newState->transition = t;	// keeps track of the transition_ID
				if (retCount == ret->size) {
					ret->resize(ret->size + 10);
				}
				ret->set(retCount, newState);
				retCount++;

//				res.push_back(newState);
			} //end of multiple intersected region with guard

		} //end of multiple transaction
	} //end-if of flowpipe computation not empty
	ret->trim(retCount);
	return;

} //end of while loop checking waiting_list != empty

__device__ void cudaPostC(CudaState * iS, CudaPolyhedra * ret) {

	CudaDiscreteSet d;
	d.discrete_elements[0] = iS->location->id;

	//d_reach_params;
//	CudaState S;
	CudaLocation* current_location = iS->location;
	double result_alfa = CudaComputeAlfa(iS);
	double result_beta = CudaComputeBeta(iS);
//
//	CudaMatrix  phi_trans;
//
//	if (!current_location->system_dynamics.isEmptyMatrixA) { //if A not Empty
//		CudaMatrix phi_matrix =current_location->system_dynamics.matrixA.matrix_exponentiation( d_reach_params.time_step);
//		phi_trans.assign(phi_matrix.transpose());
//
//	}
//	moved to  initreachabilityparameters
//	CudaMatrix B_trans;
//	// transpose to be done once and stored in the structure of parameters
//	if (!current_location.system_dynamics.isEmptyMatrixB) { //if B not Empty
//		current_location.system_dynamics.matrixB.transpose(B_trans);
//		//d_reach_parameters.B_trans = B_trans
//	}

	// ******************* Computing Parameters *******************************
	unsigned int NewTotalIteration = d_reach_params.Iterations;
	if (current_location->invariantExists) {

		if (current_location->system_dynamics.isEmptyMatrixB == true && current_location->system_dynamics.isEmptyC == true) {
			//Approach of Coarse-time-step and Fine-time-step
			cudaJumpInvariantBoundaryCheck(iS, NewTotalIteration);
		} else {
			//Approach of Sequential invariant check will work for all case
			cudaInvariantBoundaryCheck(iS, NewTotalIteration, result_alfa, result_beta);
		}
	}
	gpuReachabilitySequential(iS, NewTotalIteration, result_alfa, result_beta, ret);

//	return reach_region;

}
__global__ void testOnGpu(CudaMatrix * testMatrix) {

	int id = blockIdx.x * blockDim.x + threadIdx.x;
//	CudaState *tQueue = d_queue + (threadQueueSize * sizeof(CudaState) * id);
//	//dequeue
//	cudaError_t cudaStat;

	//test convertMatrix
	printf("\nconvertMatrix test");
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {

			if (!testMatrix->get(i, j) == i * j) {
				printf("Matrix test error convertMatrix has a bug \n");
				break;
			}
		}
	}

	printf("OK\n");
	//test aloc dealoc
	CudaMatrix mx = CudaMatrix(10, 10);
	mx.deallocateData();
	mx.alocData();
	mx.deallocateData();
	mx.deallocateData();
	mx.alocData();
	mx.set(0, 0, 1);

	mx.deallocateData();
	mx.alocData();

	//test asign , set ,get
	CudaMatrix m1 = CudaMatrix(10, 10);
	CudaMatrix m2 = CudaMatrix(10, 10);
	CudaMatrix m3 = CudaMatrix(10, 10);

	for (int i = 0; i < m1.rows; i++) {
		for (int j = 0; j < m1.cols; j++) {
			m1.set(i, j, 1);
			m2.set(i, j, 2);
			m3.set(i, j, 10);

		}

	}
	printf("get/set test");
	for (int i = 0; i < m1.rows; i++) {
		for (int j = 0; j < m1.cols; j++) {
			if (m1.get(i, j) != 1) {
				printf("Matrix test error either set or get have a bug\n");
				break;
			}
			if (m2.get(i, j) != 2) {
				printf("Matrix test error either set or get have a bug\n");
				break;
			}
		}

	}

	printf("OK\n");
	printf("assign test");
	m1.assign(&m3);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("\n");
	m2.assign(&m3);
	for (int i = 0; i < m1.rows; i++) {
		for (int j = 0; j < m1.cols; j++) {
			if (m1.get(i, j) != 10) {
				printf("Matrix test error assign has a bug\n");
				break;
			}
			if (m2.get(i, j) != 10) {
				printf("Matrix test error  assign has  a bug\n");
				break;
			}
		}

	}

	printf("OK\n");
	printf("max Norm test ");
	CudaMatrix normTest = CudaMatrix(10, 10);

	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {

			normTest.set(i, j, j * i);

		}

	}

	double norm = normTest.norm_inf();
	if (norm != 81) {
		printf("Matrix test error norm_inf has a bug norm was %f\n", norm);
	}
	printf("OK\n");
	printf("multiplyByScalar test ");

	CudaMatrix ms1;
	CudaMatrix::multiplyByScalar(&m1, 10, &ms1);
	CudaMatrix ms2;
	CudaMatrix::multiplyByScalar(&m2, 2, &ms2);
	for (int i = 0; i < ms1.rows; i++) {
		for (int j = 0; j < ms1.cols; j++) {
			if (ms1.get(i, j) != 100) {
				printf("Matrix test error multiplyByScalar has a bug\n");
				break;
			}
			if (ms2.get(i, j) != 20) {
				printf("Matrix test error  multiplyByScalar has  a bug\n");
				break;
			}
		}

	}

	printf("OK\n");

	printf("add scalar test ");

	CudaMatrix ma1 = CudaMatrix(10, 10);
	CudaMatrix ma2 = CudaMatrix(10, 10);
	for (int i = 0; i < ma1.rows; i++) {
		for (int j = 0; j < ma2.cols; j++) {

			ma1.set(i, j, 3);
			ma2.set(i, j, 4);
		}

	}
	CudaMatrix mAdd1;
	CudaMatrix::addScalar(&ma1, 10, &mAdd1);
	CudaMatrix mAdd2;
	CudaMatrix::addScalar(&ma2, 2, &mAdd2);
	for (int i = 0; i < mAdd1.rows; i++) {
		for (int j = 0; j < mAdd1.cols; j++) {
			if (mAdd1.get(i, j) != 13) {
				printf("mAddtrix test error addScalar has a bug 13!= %f\n", mAdd1.get(i, j));
				break;
			}
			if (mAdd2.get(i, j) != 6) {
				printf("mAddtrix test error addScalar has a bug 6!= %f\n", mAdd2.get(i, j));
				break;
			}
		}

	}

	printf("OK\n");
	printf("times test ");

	CudaMatrix m9x2 = CudaMatrix(9, 2);
	CudaMatrix m2x2 = CudaMatrix(2, 2);
	for (int i = 0; i < m9x2.rows; i++) {
		for (int j = 0; j < m9x2.cols; j++) {
//			printf(" %d",i);
			m9x2.set(i, j, i);
		}
//		printf("\n 9x2 ");

	}
	for (int i = 0; i < m2x2.rows; i++) {
		for (int j = 0; j < m2x2.cols; j++) {

			m2x2.set(i, j, i);
		}

	}

	CudaMatrix timesRes;
	CudaMatrix::times(&m9x2, &m2x2, &timesRes);
	for (int i = 0; i < timesRes.rows; i++) {
			for (int j = 0; j < timesRes.cols; j++) {

				m9x2.set(i, j, i);
			}


		}
	if (timesRes.rows != 9) {
		printf("TIMES ERROR EXPECTED ROWS 9 BUT WAS %d ", timesRes.rows);
	}
	if (timesRes.cols != 2) {
		printf("TIMES ERROR EXPECTED cols 2 BUT WAS %d ", timesRes.cols);
	}
	for (int i = 0; i < m2x2.rows; i++) {
		for (int j = 0; j < m2x2.cols; j++) {

			if (timesRes.get(i, j) != i) {
				printf("matrix multipl. error should be %d, is %f", i, timesRes.get(i, j));
			}
		}

	}

	CudaMatrix mt3 = CudaMatrix(10, 10);
	CudaMatrix mt4 = CudaMatrix(10, 10);
	for (int i = 0; i < mt3.rows; i++) {
		for (int j = 0; j < mt3.cols; j++) {

			mt3.set(i, j, 3);
			mt4.set(i, j, 4);
		}

	}
	CudaMatrix mt;
	CudaMatrix::times(&mt3, &mt4, &mt);

	for (int i = 0; i < mt.rows; i++) {
		for (int j = 0; j < mt.cols; j++) {
			if (mt.get(i, j) != 120) {
				printf("Matrix test error times has a bug result was %f\n", mt.get(i, j));
				break;
			}
		}

	}

	printf("OK\n");
	printf("minus test");

	CudaMatrix mm3 = CudaMatrix(10, 10);
	CudaMatrix mm4 = CudaMatrix(10, 10);
	for (int i = 0; i < mm3.rows; i++) {
		for (int j = 0; j < mm3.cols; j++) {

			mm3.set(i, j, 3 * i * j);
			mm4.set(i, j, 4 * i * j);
		}

	}

	CudaVectorTemplate<double> vec1(10);
	if (vec1.size != 10) {
		printf("VECTOR SIZE MISSMATCH\n");
	}
	for (int i = 0; i < vec1.size; i++) {
		vec1.set(i, i);
	}
	vec1.resize(15);

	if (vec1.size != 15) {
		printf("VECTOR SIZE MISSMATCH\n");
	}
	for (int i = 0; i < vec1.size; i++) {
		vec1.set(i, i);
	}
//	CudaMatrix mm = CudaMatrix::minus(&mm3, &mm4);
//
//	for (int i = 0; i < mm.rows; i++) {
//		for (int j = 0; j < mm.cols; j++) {
//			if (mm.get(i, j) != -1 * i * j) {
//				printf("Matrix test error minus has a bug result was %f\n", mm.get(i, j));
//				break;
//			}
//		}
//
//	}
//	printf("OK\n");
//	printf("minusEquals test");
//
//	CudaMatrix mmeq3 = CudaMatrix(10, 10);
//	CudaMatrix mmeq4 = CudaMatrix(10, 10);
//	for (int i = 0; i < mmeq3.rows; i++) {
//		for (int j = 0; j < mmeq3.cols; j++) {
//
//			mmeq3.set(i, j, 3 * i * j);
//			mmeq4.set(i, j, 4 * i * j);
//		}
//
//	}
//	mmeq3.minusEquals(&mmeq4);
//
//	for (int i = 0; i < mmeq3.rows; i++) {
//		for (int j = 0; j < mmeq3.cols; j++) {
//			if (mmeq3.get(i, j) != -1 * i * j) {
//				printf("Matrix test error minus has a bug result was %f\n", mmeq3.get(i, j));
//				break;
//			}
//		}
//
//	}
//	printf("OK\n");
//	printf("plus test");
//
//	CudaMatrix mp3 = CudaMatrix(10, 10);
//	CudaMatrix mp4 = CudaMatrix(10, 10);
//	for (int i = 0; i < mp3.rows; i++) {
//		for (int j = 0; j < mp3.cols; j++) {
//
//			mp3.set(i, j, 3 * i * j);
//			mp4.set(i, j, 4 * i * j);
//		}
//
//	}
//	CudaMatrix mp = CudaMatrix::plus(&mp3, &mp4);
//
//	for (int i = 0; i < mp.rows; i++) {
//		for (int j = 0; j < mp.cols; j++) {
//			if (mp.get(i, j) != 7 * i * j) {
//				printf("Matrix test error plus has a bug result was %f\n", mp.get(i, j));
//				break;
//			}
//		}
//
//	}
//	printf("OK\n");
//	printf("identity_matrix test");
//
//	CudaMatrix im = CudaMatrix::identity_matrix(3);
//	double data[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
//	CudaMatrix test = CudaMatrix(3, 3, data);
//	for (int i = 0; i < im.rows; i++) {
//		for (int j = 0; j < im.cols; j++) {
//			if (im.get(i, j) != test.get(i, j)) {
//				printf("Matrix test error identity_matrix has a bug %d,%d was %f and not %f\n", i, j, im.get(i, j), test.get(i, j));
//				break;
//			}
//		}
//
//	}
//	printf("OK\n");
//	printf("result is %d <\n",res);
//	printf("test 3 <\n");

}
__global__ void agjhGPu(CudaState *d_queue, CudaLocation * d_location, int size) {

	int id = blockIdx.x * blockDim.x + threadIdx.x;
	CudaState *tQueue = d_queue + (threadQueueSize * sizeof(CudaState) * id);
	//dequeue
	cudaError_t cudaStat;

	if (id == 0) {
		//for (int i = 0; i < threadQueueSize; i++) {

		int i = 0;
		//search
		CudaPolyhedra R;
		cudaPostC(&tQueue[i], &R);

//				 symbolic_states::ptr R1 = symbolic_states::ptr(new symbolic_states());
//				 R1->setContinuousSetptr(R);
		CudaState R1;
		R1.continuousSet.assign(&R);

		CudaDiscreteSet d;
//				 discrete_set d;
//				 d.insert_element(s->getLocationId());
		d.discrete_elements[0] = tQueue[i].location_id;
//				 R1->setDiscreteSet(d);
	R1.discreteSet = d;
		//PASSED.push_back(R1); TODO r tree and stuff
//
		CudaVectorTemplate<CudaState*> R2;
		cudaPostD(&R1, d_location, &R2);
		printf("\nfound %d states", R2.size);
		//	}
	}
}
CudaMatrix convertMatrix(math::matrix<double> m) {
	CudaMatrix ret(m.size1(), m.size2());
	for (int i = 0; i < ret.rows; i++) {
		for (int j = 0; j < ret.cols; j++) {
			ret.set(i, j, m(i, j));
		}

	}
	double *d;
	gpuErrchk(cudaMalloc(&d, (ret.rows * ret.cols) * sizeof(double*)));
	gpuErrchk(cudaMemcpy(d, ret.data, (ret.rows * ret.cols) * sizeof(double*), cudaMemcpyHostToDevice));
	ret.data = d;
	return ret;

}

CudaVector convertVector(math::vector<double> m) {

	CudaVector ret(m.size());
	for (int i = 0; i < ret.size; i++) {

		ret.set(i, m.at(i));
	}

	double *d;
	gpuErrchk(cudaMalloc(&d, (ret.size) * sizeof(double*)));
	gpuErrchk(cudaMemcpy(d, ret.data, (ret.size) * sizeof(double*), cudaMemcpyHostToDevice));
	ret.data = d;

	return ret;

}
void loadPolytope(CudaPolytope* cpol, polytope* pol) {

	cpol->InEqualitySign = pol->getInEqualitySign();
	cpol->coeffMatrix = convertMatrix(pol->getCoeffMatrix());
	cpol->columnVector = convertVector(pol->getColumnVector());

	cpol->isEmpty = pol->getIsEmpty();
	cpol->isUniverse = pol->getIsUniverse();
	cpol->number_facets = pol->getNumberFacets();
	cpol->systemDimension = pol->getSystemDimension();

}
CudaState AGJHGPU::convertStateToCuda(initial_state::ptr state) {
	CudaState ret;
	ret.transition_id = state->getTransitionId();
	ret.location_id = state->getLocationId();
//	ret.location pointer is set in the copy to gpu method
	//loadLocation(&(ret.location),H.getLocation(state.get()->getLocationId()).get());

	loadPolytope(&ret.initialSet, state.get()->getInitialSet().get());
//TODO CUDA VECTOR AS LIST : SIZE DOES NOT EQUAL CAPACity
	return ret;

}
void loadSystemDynamics(CudaDynamics* cDyn, Dynamics * dyn) {

	cDyn->isEmptyMatrixA = dyn->isEmptyMatrixA;	//True if empty otherwise False
	cDyn->matrixA = convertMatrix(dyn->MatrixA);
	cDyn->isEmptyMatrixB = dyn->isEmptyMatrixB;	//True if empty otherwise False
	cDyn->matrixB = convertMatrix(dyn->MatrixB);
	loadPolytope(&(cDyn->U), dyn->U.get());
	cDyn->C = convertVector(dyn->C);
	cDyn->isEmptyC = dyn->isEmptyC;

}

void loadTransition(CudaTransition *cTrans, transition* trans) {

//		cTrans->assignT.Map.assign(convertMatrix(trans->getAssignT().Map));
	cTrans->assignT.Map = convertMatrix(trans->getAssignT().Map);
	cTrans->assignT.b = convertVector(trans->getAssignT().b);

	loadPolytope(&cTrans->Gaurd, trans->getGaurd().get());
	cTrans->trans_id = trans->getTransitionId();

}
void AGJHGPU::loadLocation(CudaLocation *cLocation, location* loc) {
	string locName = loc->getName();
	if ((locName.compare("BAD") == 0)) {
		cLocation->name = CudaLocation::BAD;
	}
	if ((locName.compare("GOOD") == 0)) {
		cLocation->name = CudaLocation::GOOD;
	}
	if ((locName.compare("FINAL") == 0)) {
		cLocation->name = CudaLocation::FINAL;
	}
	if ((locName.compare("UNSAVE") == 0)) {
		cLocation->name = CudaLocation::UNSAVE;
	}

	cLocation->id = loc->getLocId();
//	cState->location.invariant
	loadPolytope(&(cLocation->invariant), loc->getInvariant().get());
//	cState->location.invariantExists
	cLocation->invariantExists = loc->isInvariantExists();
//	cState->location.system_dynamics
	loadSystemDynamics(&(cLocation->system_dynamics), &(loc->getSystem_Dynamics()));

}
void AGJHGPU::copyDataToGpu() {

	//copy locations
	int transitionCounter = 0;
	int locationCounter = 0;
	std::map<int, int> transitionMap;
	//std::map<int, int> locationMap;

	int maxLocId=0;

	//calculate pointers
	for (std::map<int, location::ptr>::iterator it = H.list_locations.begin(); it != H.list_locations.end(); ++it) {
		std::list<transition::ptr> tList = it->second.get()->getOut_Going_Transitions();
	//	locationMap[it->second.get()->getLocId()] = locationCounter;
		if(it->second.get()->getLocId()<0){
			cout<<"ERROR GPU IMPLEMENTATION DOES NOT SUPORT NEGATIVE VALUET TRANSITION IDS";
		}
		maxLocId=max(maxLocId,it->second.get()->getLocId());
		locationCounter++;
		for (std::list<transition::ptr>::iterator tr = tList.begin(); tr != tList.end(); ++tr) {
			transitionMap[(*tr)->getTransitionId()] = transitionCounter;
			transitionCounter++;
		}
	}
	int j = 0;
	gpuErrchk(cudaMalloc((void ** ) &d_location,  (maxLocId+1) * sizeof(CudaLocation)));
	for (std::map<int, location::ptr>::iterator it = H.list_locations.begin(); it != H.list_locations.end(); ++it) {
		CudaLocation l;
		loadLocation(&l, it->second.get());
		std::list<transition::ptr> tList = it->second.get()->getOut_Going_Transitions();
		CudaTransition * d_loc_transition;
		gpuErrchk(cudaMalloc((void ** ) &d_loc_transition, transitionCounter* sizeof(CudaTransition)));
		int k = 0;
		for (std::list<transition::ptr>::iterator tr = tList.begin(); tr != tList.end(); ++tr) {
//			CudaTransition*  target;
//			target=d_transitions;
//			target+=(transitionMap[(*tr).get()->getTransitionId()]);
//			gpuErrchk(cudaMemcpy(d_loc_transition + k , &target, sizeof(CudaTransition*), cudaMemcpyHostToDevice));

			CudaTransition t;//TODO TRANSITION ID A POSITION IN TRANSITION ARRAy`
			loadTransition(&t, (*tr).get());
			t.destination_location = d_location +(*tr).get()->getDestination_Location_Id();// locationMap[(*tr).get()->getDestination_Location_Id()];
			t.source_location = d_location + (*tr).get()->getDestination_Location_Id();//locationMap[(*tr).get()->getSource_Location_Id()];
			//int pos=transitionMap[t.trans_id];
			gpuErrchk(cudaMemcpy(d_loc_transition + k, &t, sizeof(CudaTransition), cudaMemcpyHostToDevice));
			int pos = transitionMap[t.trans_id];
		//	gpuErrchk(cudaMemcpy(d_location + pos, &t, sizeof(CudaLocation), cudaMemcpyHostToDevice));
			k++;
		}
		l.Out_Going_Transitions = d_loc_transition;
		l.Out_Going_Transitions_size = tList.size();
		int pos =l.id;// locationMap[l.id];
		gpuErrchk(cudaMemcpy(d_location + pos, &l, sizeof(CudaLocation), cudaMemcpyHostToDevice));
		j++;

	}

	j = 0;
	for (std::list<initial_state::ptr>::iterator i = I.begin(); i != I.end(); i++) {
		CudaState s = convertStateToCuda((*i));

		//update transition pointer
		int pos = transitionMap[s.transition_id];
		s.transition = d_transitions + pos;
//		update location pointer

		s.location = d_location + s.location_id	;

		gpuErrchk(cudaMemcpy(d_queue + j, &s, sizeof(CudaState), cudaMemcpyHostToDevice));
		j++;
	}

	initReachabilityParameters();
}

void runUnitTest() {
	CudaMatrix *d_m;

	cudaMalloc((void **) &d_m, sizeof(CudaMatrix));
	matrix<double> m(10, 10);
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			m(i, j) = i * j;

		}
	}
	/*	matrix<double> m9x2(9, 2);
	matrix<double> m2x2(2, 2);
	for (int i = 0; i < m9x2.size1(); i++) {
		for (int j = 0; j < m9x2.size2(); j++) {

			m9x2(i, j) = i;
		}

	}
	for (int i = 0; i < m2x2.size1(); i++) {
		for (int j = 0; j < m2x2.size2(); j++) {

			m2x2(i, j) = i;
		}

	}
	matrix<double> res;
	m9x2.multiply(m2x2, res);
	for (int i = 0; i < res.size1(); i++) {
		for (int j = 0; j < res.size2(); j++) {

			printf(" %f", res(i, j));
		}
		printf("\n");

	}*/
	CudaMatrix h_m = convertMatrix(m);
	gpuErrchk(cudaMemcpy(d_m, &h_m, sizeof(CudaMatrix), cudaMemcpyHostToDevice));
	h_m.data = NULL;
	testOnGpu<<<blockSize,threadSize>>>(d_m);
}
void AGJHGPU::initReachabilityParameters() {

	CudaReachabilityParameters params;
	params.Iterations = this->reach_parameters.Iterations;
	params.TimeBound = this->reach_parameters.TimeBound;
	params.time_step = this->reach_parameters.time_step;
	//DIRETIONS
	matrix<double> m = this->reach_parameters.Directions;
	for (int i = 0; i < m.size1(); i++) {
		for (int j = 0; j < m.size2(); j++) {
			int p = IDX2F(i, j, m.size1());

			params.Directions_DATA[p] = m(i, j);
		}

	}
	params.Directions_rows = m.size1();
	params.Directions_cols = m.size2();
	//TOTAL DIRECTIONS
	m = this->reach_parameters.TotalDirections;
	params.TotalDirections_rows = m.size1();
	params.TotalDirections_cols = m.size2();
	for (int i = 0; i < m.size1(); i++) {
		for (int j = 0; j < m.size2(); j++) {
			int p = IDX2F(i, j, m.size1());

			params.TotalDirections_DATA[p] = m(i, j);
		}

	}
	//A_INV
	m = this->reach_parameters.A_inv;
	params.A_inv_DATA_rows = m.size1();
	params.A_inv_DATA_cols = m.size2();
	for (int i = 0; i < m.size1(); i++) {
		for (int j = 0; j < m.size2(); j++) {
			int p = IDX2F(i, j, m.size1());

			params.A_inv_DATA[p] = m(i, j);
		}

	}

	gpuErrchk(cudaMemcpyToSymbol(d_reach_params, &params, sizeof(CudaReachabilityParameters)));

}

std::list<symbolic_states::ptr> AGJHGPU::gpuReachability(std::list<abstractCE::ptr> &ce_candidates) {
	std::list<symbolic_states::ptr> ret; //	template_polyhedra::ptr reach_region;

	int device = 0;
	cudaDeviceProp props;
	cudaGetDevice(&device);
	cudaGetDeviceProperties(&props, device);

	size_t threadQueueMemSize = threadQueueSize * threadCount * sizeof(CudaState);
	size_t locationArraySize = H.list_locations.size() * sizeof(CudaLocation);
	int transitionsCount = 0;
	for (std::map<int, location::ptr>::iterator it = H.list_locations.begin(); it != H.list_locations.end(); ++it) {
		transitionsCount += it->second.get()->getOut_Going_Transitions().size();

	}

	size_t transitionsSize = sizeof(CudaTransition) * transitionsCount;

	size_t hashTableMemSize = props.totalGlobalMem - transitionsSize - threadQueueMemSize - sizeof(CudaReachabilityParameters) - (1024 * 1024 * 1024); //<TODO FIND MISSING  MEMORY

	//alocate and fill queue

	cout << "using " << transitionsSize << " B  for transition array. Size of one transition:" << sizeof(CudaTransition) << " transition count:" << transitionsCount << "\n";
	gpuErrchk(cudaMalloc((void ** ) &d_transitions, transitionsSize * sizeof(CudaTransition)));

	cout << "using " << locationArraySize << " B  for location array. Size of one location:" << sizeof(CudaLocation) << " location count:" << H.list_locations.size() << "\n";

	cout << "using " << threadQueueMemSize << " B  for queue. Size of one state:" << sizeof(CudaState) << "\n";
	gpuErrchk(cudaMalloc((void ** ) &d_queue, threadQueueSize * threadCount * sizeof(CudaState)));

	//alocate  hashTable

	cout << "using " << hashTableMemSize / (1024 * 1024) << " mb for hashTable\n";
	//gpuErrchk(cudaMalloc((void ** ) &d_hashTable, hashTableMemSize));

	int queueItemCount = threadQueueSize * threadCount;
	runUnitTest();
	copyDataToGpu();

	//BFS
	cout << "RUNING CUDA KERNEL WITH BLOCK SIZE: " << blockSize << " THREAD SIZE: " << threadSize;

	agjhGPu<<<blockSize,threadSize>>>(d_queue,d_location,queueItemCount);
	gpuErrchk(cudaPeekAtLastError());

	cudaFree(d_queue);
//	cudaFree(d_hashTable);
	/*	symbolic_states c;
	 c.transition_id=555;
	 // create class storage on device and copy top level class
	 symbolic_states *d_c;
	 cudaMalloc((void **)&d_c, sizeof(symbolic_states));

	 cudaMemcpy(d_c, &c, sizeof(symbolic_states), cudaMemcpyHostToDevice);

	 test<<<1,32>>>(c);*/

	/*

	 int a, b, c;
	 int *d_a, *d_b, *d_c;
	 int size = sizeof(int);
	 // host copies of a, b, c
	 // device copies of a, b, c
	 // Allocate space for device copies of a, b, c
	 cudaMalloc((void **)&d_a, size);
	 cudaMalloc((void **)&d_b, size);
	 cudaMalloc((void **)&d_c, size);
	 // Setup input values
	 a = 2;
	 b = 7;

	 // Copy inputs to device
	 cudaMemcpy(d_a, &a, size, cudaMemcpyHostToDevice);
	 cudaMemcpy(d_b, &b, size, cudaMemcpyHostToDevice);
	 // Launch add() kernel on GPU
	 //	add<<<blockSize,threadSize>>>(d_a, d_b, d_c);
	 // Copy result back to host
	 cudaMemcpy(&c, d_c, size, cudaMemcpyDeviceToHost);

	 cout<<"result\n";
	 cout<<c;
	 cout<<"\n";
	 cout<<a+b;*/
	// Cleanup
	cudaDeviceSynchronize();

	return ret;
}

AGJHGPU::AGJHGPU() {

}

AGJHGPU::~AGJHGPU() {

}
