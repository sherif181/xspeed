/*
 * invariantCheck.cpp
 *
 *  Created on: 16-Nov-2014
 *      Author: amit
 */

#include "invariantBoundaryCheck.h"

void InvariantBoundaryCheck(Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial, ReachabilityParameters& ReachParameters,
		polytope::ptr invariant, int lp_solver_type_choosen, unsigned int &newTotIters) {

	unsigned int shm_NewTotalIteration = ReachParameters.Iterations; //Shared Variable for resize iterations number on crossing with invariant
	int dimension = Initial->getSystemDimension();
	int Min_Or_Max = 2;
	int numberOfInvariants = invariant->getColumnVector().size(); //total number of Invariant's constraints
	// ******************* Probable race condition variables:: can be  placed inside for-loop *******************
	int foundStart = 0, intersection_start, intersection_end;
	//bool invariantCrossed = false;
	// *************************** For Negative ************************************
	double res1_minus, term2_minus = 0.0, result1_minus, term1_minus = 0.0, result_minus;
	std::vector<double> Btrans_dir_minus, phi_trans_dir_minus, phi_trans_dir1_minus;
	math::matrix<double> B_trans, phi_tau_Transpose;
	if (!SystemDynamics.isEmptyMatrixA) //current_location's SystemDynamics's or ReachParameters
		phi_tau_Transpose = ReachParameters.phi_trans;
	if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
		B_trans = ReachParameters.B_trans;
	// *******************************************************************************************************************
	std::vector<int> boundaryIterations(numberOfInvariants, shm_NewTotalIteration); // size(dimension_size,initial_value)
	//cout<<"Test 1.1 \n";
int type = lp_solver_type_choosen;
//#pragma omp parallel for //num_threads(numberOfInvariants)
	for (int eachInvariantDirection = 0; eachInvariantDirection < numberOfInvariants; eachInvariantDirection++) {
		double TempOmega_min;
		double zI_min = 0.0, zV_min = 0.0;
		double sVariable_min, s1Variable_min; //For Minimization of First Vector only

		std::vector<double> r1Variable_minus;
		r1Variable_minus.resize(dimension);

		std::vector<double> rVariable_minus;
		rVariable_minus.resize(dimension);

		for (int i = 0; i < dimension; i++) {
			rVariable_minus[i] = -1 * invariant->getCoeffMatrix()(eachInvariantDirection, i); //Second vector negative of First vector
		}
		//	cout<<"Test 1.1 \n";
		// ******** Lamda Computation for each invariant's directions/constraints **********
		double invariant_SupportFunction;
		invariant_SupportFunction = invariant->getColumnVector()[eachInvariantDirection];
		// **********************XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX******************

		lp_solver s_per_thread_I_minus(type), s_per_thread_U_minus(type);
// ******************************************* For Negative Direction Starts *******************************************
		s_per_thread_I_minus.setMin_Or_Max(2);
		if (!Initial->getIsEmpty()) //set glpk constraints If not an empty polytope
			s_per_thread_I_minus.setConstraints(ReachParameters.X0->getCoeffMatrix(), ReachParameters.X0->getColumnVector(), ReachParameters.X0->getInEqualitySign());

		//cout<<"Test 18\n";
		s_per_thread_U_minus.setMin_Or_Max(2);
		if (!SystemDynamics.U->getIsEmpty()) { //empty polytope
			s_per_thread_U_minus.setConstraints(SystemDynamics.U->getCoeffMatrix(), SystemDynamics.U->getColumnVector(), SystemDynamics.U->getInEqualitySign());
		}
// ******************************************* Negative Direction Ends *******************************************
		unsigned int loopIteration = 0;
		sVariable_min = 0.0;

		// ******************************************* For Negative Direction Starts *******************************************
		//zIInitial = Omega_Support(ReachParameters, rVariable_minus, Initial,SystemDynamics, s_per_thread_I_minus, s_per_thread_U_minus, Min_Or_Max);
		double term3_minus, term3a_minus, term3b_minus, res2_minus, term3c_minus = 0.0;
		//  **************    Omega Function   ********************
		res1_minus = Initial->computeSupportFunction(rVariable_minus, s_per_thread_I_minus);
		if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
			phi_tau_Transpose.mult_vector(rVariable_minus, phi_trans_dir_minus);
			term1_minus = Initial->computeSupportFunction(phi_trans_dir_minus, s_per_thread_I_minus);
		}else if (SystemDynamics.isEmptyMatrixA) { //if A is empty :: {tau.A}' reduces to zero so, e^{tau.A}' reduces to 1
			// so, 1 * rVariable give only rVariable
			term1_minus = Initial->computeSupportFunction(rVariable_minus, s_per_thread_I_minus);
		}//handling constant dynamics

		if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
			B_trans.mult_vector(rVariable_minus, Btrans_dir_minus);
		if (!SystemDynamics.isEmptyMatrixB && !SystemDynamics.U->getIsEmpty())
			term2_minus = ReachParameters.time_step * SystemDynamics.U->computeSupportFunction(Btrans_dir_minus, s_per_thread_U_minus);
		term3a_minus = ReachParameters.result_alfa;
		term3b_minus = (double) support_unitball_infnorm(rVariable_minus);
		//cout<<"Test Before dot_product \n";
		if (!SystemDynamics.isEmptyC) {
			term3c_minus = ReachParameters.time_step * dot_product(SystemDynamics.C, rVariable_minus); //Added +tau* sf_C(l) 8/11/2015
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
				boundaryIterations[eachInvariantDirection] = loopIteration; //Made Changes here due to circle
				break;//no need to compute reachable set anymore due to out of invariant boundary
			}
			// ******************************************* Intersection Detection Section Ends *******************************************
// ************************************************  For Negative Direction Starts *******************************************
			//double TempOmega_min;
			double term3_minus, term3a_minus, res2_minus, beta_minus, res_beta_minus;

			//	zV = W_Support(ReachParameters, SystemDynamics, rVariable,s_per_thread_U, Min_Or_Max);
			//  **************    W_Support Function   ********************
			result1_minus = term2_minus;
			beta_minus = ReachParameters.result_beta;
			res_beta_minus = beta_minus * term3b_minus; //Replacing term3b from previous step
			result_minus = result1_minus + res_beta_minus + term3c_minus;	// **********************************TODO :: VERIFY res_beta_minus from res_beta
			zV_min = result_minus;
			//  **************  W_Support Function Over  ********************
			s1Variable_min = sVariable_min + zV_min;
			if (SystemDynamics.isEmptyMatrixA) { //Matrix A is empty for constant dynamics
				r1Variable_minus = rVariable_minus;
			}else{
				r1Variable_minus = phi_trans_dir_minus;
			}

			//zI = Omega_Support(ReachParameters, r1Variable, Initial,SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
			//  **************    Omega Function   ********************
			res1_minus = term1_minus;
			if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
				phi_tau_Transpose.mult_vector(r1Variable_minus, phi_trans_dir_minus);
				term1_minus = Initial->computeSupportFunction(phi_trans_dir_minus, s_per_thread_I_minus);
			}
			if (!SystemDynamics.isEmptyMatrixB) { //current_location's SystemDynamics's or ReachParameters
				B_trans.mult_vector(r1Variable_minus, Btrans_dir_minus);
				term2_minus = ReachParameters.time_step * SystemDynamics.U->computeSupportFunction(Btrans_dir_minus, s_per_thread_U_minus);
			}
			term3a_minus = ReachParameters.result_alfa;
			term3b_minus = support_unitball_infnorm(r1Variable_minus);
			//std::cout<<"No Error here "<<loopIteration<<std::endl;
			if (!SystemDynamics.isEmptyC) {
				//std::cout<<"SystemDynamics.C.size() = "<<SystemDynamics.C.size()<<std::endl;
				//std::cout<<"rVariable_minus.size() = "<<rVariable_minus.size()<<std::endl;
				term3c_minus = ReachParameters.time_step * dot_product(SystemDynamics.C, rVariable_minus); //Added +tau* sf_C(l) 8/11/2015
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

			rVariable_minus = CopyVector(r1Variable_minus);
			sVariable_min = s1Variable_min;
		} //end of iterations
		//	if (invariantCrossed)	//We need to check all invariants as we are taking the min of all the directions
		//		break;//boundary crossed so no need to check other invariant's directions
	} //end of parallel for each Iterations or Time-Bound
// At the end of the For-Loop or all invariant_Directions we have boundaryIterations vector with the different limit to stop iterations
// The Minimum(boundaryIterations) will be the final Boundary Limit for the number of iterations
	unsigned int min_Total_Iteration;
	min_Total_Iteration = *min_element(boundaryIterations.begin(), boundaryIterations.end()); // - 1 ;//excluding the last outside Omega
	newTotIters = min_Total_Iteration;
}

//In-efficient method of checking each invariants face one after another
void InvariantBoundaryCheck1(Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial, ReachabilityParameters& ReachParameters,
		polytope::ptr invariant, int lp_solver_type_choosen, unsigned int &newTotIters) {

	unsigned int shm_NewTotalIteration = ReachParameters.Iterations; //Shared Variable for resize iterations number on crossing with invariant
//	cout<<"shm_NewTotalIteration = "<<shm_NewTotalIteration<<"OK\n";
	int dimension = Initial->getSystemDimension();
	int Min_Or_Max = 2;
	int numberOfInvariants = invariant->getColumnVector().size(); //total number of Invariant's constraints
	/*cout <<"invariant->getCoeffMatrix() = "<<invariant->getCoeffMatrix()<<"\n";
	 for (int i=0;i<invariant->getColumnVector().size();i++)
	 cout<<invariant->getColumnVector()[i]<<"\t";*/

	// ******************* Probable race condition variables:: can be  placed inside for-loop *******************
	int foundStart = 0, intersection_start, intersection_end;
	//bool invariantCrossed = false;
	// *************************** For Positive ************************************
	double res1, term2 = 0.0, result1, term1 = 0.0, result;
	std::vector<double> Btrans_dir, phi_trans_dir, phi_trans_dir1;
	// *************************** For Negative ************************************
	double res1_minus, term2_minus = 0.0, result1_minus, term1_minus = 0.0,
			result_minus;
	std::vector<double> Btrans_dir_minus, phi_trans_dir_minus,
			phi_trans_dir1_minus;
	math::matrix<double> B_trans, phi_tau_Transpose;
	if (!SystemDynamics.isEmptyMatrixA) //current_location's SystemDynamics's or ReachParameters
		phi_tau_Transpose = ReachParameters.phi_trans;
	if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
		B_trans = ReachParameters.B_trans;
	// *******************************************************************************************************************
	std::vector<int> boundaryIterations(numberOfInvariants,
			shm_NewTotalIteration); // size(dimension_size,initial_value)

//#pragma omp parallel for //num_threads(numberOfInvariants)
	for (int eachInvariantDirection = 0;
			eachInvariantDirection < numberOfInvariants;
			eachInvariantDirection++) {

		double TempOmega, TempOmega_min;
		//double zIInitial, zIInitial_minus;

		double zI = 0.0, zV = 0.0;
		double zI_min = 0.0, zV_min = 0.0;
		double sVariable, s1Variable;
		double sVariable_min, s1Variable_min; //For Minimization of First Vector only

		std::vector<double> r1Variable, r1Variable_minus;
		r1Variable.resize(dimension);
		r1Variable_minus.resize(dimension);

		std::vector<double> rVariable, rVariable_minus;
		rVariable.resize(dimension);
		rVariable_minus.resize(dimension);
		//cout<<"Invariant Number = "<<eachInvariantDirection<<endl;
		//	cout<<"Test dim = "<<dimension<<"\n";
		for (int i = 0; i < dimension; i++) {
			rVariable[i] = invariant->getCoeffMatrix()(eachInvariantDirection,
					i); //First vector positive
			rVariable_minus[i] = -1
					* invariant->getCoeffMatrix()(eachInvariantDirection, i); //Second vector negative of First vector
		}
		//	cout<<"Test 1.1 \n";
		// ******** Lamda Computation for each invariant's directions/constraints **********
		double invariant_SupportFunction;
		int type = lp_solver_type_choosen;
		//	cout<<"Test 2.a \n";
		lp_solver lpSolver(type), lp_U_dummy(type);
		//	cout<<"Test 2 \n";
		lpSolver.setMin_Or_Max(2);
		if (!invariant->getIsEmpty() && !invariant->getIsUniverse()) { //set glpk constraints If not an empty polytope
			lpSolver.setConstraints(invariant->getCoeffMatrix(),
					invariant->getColumnVector(),
					invariant->getInEqualitySign());
			invariant_SupportFunction = invariant->computeSupportFunction(
					rVariable, lpSolver);
			//	std::cout << "\neachInvariantDirection = "<<eachInvariantDirection<<" Invariant SupportFunction = " << invariant_SupportFunction << endl;
		}
		// **********************XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX******************
		lp_solver s_per_thread_I(type), s_per_thread_U(type);
// ************************************************ For Positive Direction Starts ************************************************
		s_per_thread_I.setMin_Or_Max(2);
		if (!Initial->getIsEmpty()) { //set glpk constraints If not an empty polytope
			s_per_thread_I.setConstraints(ReachParameters.X0->getCoeffMatrix(),
					ReachParameters.X0->getColumnVector(),
					ReachParameters.X0->getInEqualitySign());
			//s_per_thread_I.setConstraints(Initial->getCoeffMatrix(),ReachParameters.X0->getColumnVector(),ReachParameters.X0->getInEqualitySign());
		}
		s_per_thread_U.setMin_Or_Max(2);
		if (SystemDynamics.U->getIsEmpty()) { //empty polytope
			//Polytope is empty so no glpk object constraints to be set
			//std::cout << "InvariantBoundCheck: Input set U is empty\n";
		} else {
			s_per_thread_U.setConstraints(SystemDynamics.U->getCoeffMatrix(),
					SystemDynamics.U->getColumnVector(),
					SystemDynamics.U->getInEqualitySign());
		}
		//cout<<"Test 16\n";
// ************************************************   Positive Direction Ends *******************************************
		lp_solver s_per_thread_I_minus(type), s_per_thread_U_minus(type);
// ******************************************* For Negative Direction Starts *******************************************
		s_per_thread_I_minus.setMin_Or_Max(2);
		if (!Initial->getIsEmpty()) //set glpk constraints If not an empty polytope
			s_per_thread_I_minus.setConstraints(
					ReachParameters.X0->getCoeffMatrix(),
					ReachParameters.X0->getColumnVector(),
					ReachParameters.X0->getInEqualitySign());

		//cout<<"Test 18\n";
		s_per_thread_U_minus.setMin_Or_Max(2);
		if (SystemDynamics.U->getIsEmpty()) { //empty polytope
			//Polytope is empty so no glpk object constraints to be set
		} else {
			s_per_thread_U_minus.setConstraints(
					SystemDynamics.U->getCoeffMatrix(),
					SystemDynamics.U->getColumnVector(),
					SystemDynamics.U->getInEqualitySign());
		}
// ******************************************* Negative Direction Ends *******************************************
		unsigned int loopIteration = 0;
		sVariable = 0.0; //initialize s0
		sVariable_min = 0.0;
		// ************************************************  For Positive Direction Starts *******************************************
		//zIInitial = Omega_Support(ReachParameters, rVariable, Initial,SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
		double term3, term3a, term3b, res2, term3c = 0.0;
		//  **************    Omega Function   ********************
		res1 = Initial->computeSupportFunction(rVariable, s_per_thread_I);
		if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
			phi_tau_Transpose.mult_vector(rVariable, phi_trans_dir);
			term1 = Initial->computeSupportFunction(phi_trans_dir,
					s_per_thread_I);
		}

		if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
			B_trans.mult_vector(rVariable, Btrans_dir);
		if (!SystemDynamics.isEmptyMatrixB && !SystemDynamics.U->getIsEmpty())
			term2 = ReachParameters.time_step
					* SystemDynamics.U->computeSupportFunction(Btrans_dir,
							s_per_thread_U);
		term3a = ReachParameters.result_alfa;
		term3b = (double) support_unitball_infnorm(rVariable);
		if (!SystemDynamics.isEmptyC) {
			term3c = ReachParameters.time_step
					* dot_product(SystemDynamics.C, rVariable); //Added +tau* sf_C(l) 8/11/2015
		}
		term3 = term3a * term3b;
		res2 = term1 + term2 + term3 + term3c; //term3c Added
		if (res1 > res2)
			TempOmega = res1;	//zIInitial = res1;
		else
			TempOmega = res2;	//zIInitial = res2;
		//	cout<<"res2 = "<<res2<<" and ";
		//  **************  Omega Function Over  ********************
		// ************************************************   Positive Direction Ends *******************************************

		// ******************************************* For Negative Direction Starts *******************************************
		//zIInitial = Omega_Support(ReachParameters, rVariable_minus, Initial,SystemDynamics, s_per_thread_I_minus, s_per_thread_U_minus, Min_Or_Max);
		double term3_minus, term3a_minus, term3b_minus, res2_minus,
				term3c_minus = 0.0;
		//  **************    Omega Function   ********************
		res1_minus = Initial->computeSupportFunction(rVariable_minus,
				s_per_thread_I_minus);
		if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
			phi_tau_Transpose.mult_vector(rVariable_minus, phi_trans_dir_minus);
			term1_minus = Initial->computeSupportFunction(phi_trans_dir_minus,
					s_per_thread_I_minus);
		}
		if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
			B_trans.mult_vector(rVariable_minus, Btrans_dir_minus);
		if (!SystemDynamics.isEmptyMatrixB && !SystemDynamics.U->getIsEmpty())
			term2_minus = ReachParameters.time_step
					* SystemDynamics.U->computeSupportFunction(Btrans_dir_minus,
							s_per_thread_U_minus);
		term3a_minus = ReachParameters.result_alfa;
		term3b_minus = (double) support_unitball_infnorm(rVariable_minus);
		if (!SystemDynamics.isEmptyC) {
			term3c_minus = ReachParameters.time_step
					* dot_product(SystemDynamics.C, rVariable_minus); //Added +tau* sf_C(l) 8/11/2015
		}
		term3_minus = term3a_minus * term3b_minus;
		res2_minus = term1_minus + term2_minus + term3_minus + term3c_minus;
		if (res1_minus > res2_minus)
			TempOmega_min = res1_minus;	//zIInitial_minus = res1_minus;
		else
			TempOmega_min = res2_minus;	//zIInitial_minus = res2_minus;
		//	cout<<"res2_minus = "<<res2_minus<<"\n";
//  **************  Omega Function Over  ********************
// ******************************************* Negative Direction Ends *******************************************
		loopIteration++;
		for (; loopIteration < shm_NewTotalIteration;) { //Now stopping condition is only "shm_NewTotalIteration"

			//Debugging ---inserted here invariant checking
			//		cout<<"-1 * TempOmega_min = " << -1 * TempOmega_min <<" invariant_SupportFunction = "<<invariant_SupportFunction;
			//		cout<<"TempOmega = " << TempOmega<<std::endl;
			// ******************************************* Intersection Detection Section Starts *****************************************
			if (((-1 * TempOmega_min) < invariant_SupportFunction)
					&& (invariant_SupportFunction < TempOmega)) { //Should have been correct
				intersection_start = loopIteration; //Partially inside the region
				//	cout<<"\nINTERSECTION_start = "<<intersection_start <<endl;
				//}
			}
			if ((-1 * TempOmega_min) > invariant_SupportFunction) { // Should have been correct
				//if ( ((-1 * TempOmega_min) < invariant_SupportFunction) && (TempOmega < invariant_SupportFunction)) {
				intersection_end = loopIteration;
				//	cout << "\nIntersection_End = " << intersection_end << endl;
				//	cout<<"TempOmega_min = "<<TempOmega_min<< " and TempOmega = "<<TempOmega<<"\t";
				boundaryIterations[eachInvariantDirection] = loopIteration; //Made Changes here due to circle
				//	invariantCrossed = true;//Omega is out of the invariant's boundary
				break;//no need to compute reachable set anymore due to out of invariant boundary
			}
			// ******************************************* Intersection Detection Section Ends *******************************************
			//---inserted here
// ************************************************  For Positive Direction Starts *******************************************
			//double TempOmega;
			double term3, term3a, res2, beta, res_beta;
			//	zV = W_Support(ReachParameters, SystemDynamics, rVariable,s_per_thread_U, Min_Or_Max);
			//  **************    W_Support Function   ********************
			result1 = term2;
			beta = ReachParameters.result_beta;
			res_beta = beta * term3b; //Replacing term3b from previous step
			result = result1 + res_beta + term3c; //Added term3c
			zV = result;
			//  **************  W_Support Function Over  ********************
			s1Variable = sVariable + zV;
			r1Variable = phi_trans_dir;
			//zI = Omega_Support(ReachParameters, r1Variable, Initial,SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
			//  **************    Omega Function   ********************
			res1 = term1;
			if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
				phi_tau_Transpose.mult_vector(r1Variable, phi_trans_dir);
				term1 = Initial->computeSupportFunction(phi_trans_dir,
						s_per_thread_I);
			}
			if (!SystemDynamics.isEmptyMatrixB) { //current_location's SystemDynamics's or ReachParameters
				B_trans.mult_vector(r1Variable, Btrans_dir);
				term2 = ReachParameters.time_step
						* SystemDynamics.U->computeSupportFunction(Btrans_dir,
								s_per_thread_U);
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
			//	cout<<"TempOmega = "<<TempOmega<<"\t";
// ************************************************   Positive Direction Ends *******************************************
// ************************************************  For Negative Direction Starts *******************************************
			//double TempOmega_min;
			double term3_minus, term3a_minus, res2_minus, beta_minus,
					res_beta_minus;

			//	zV = W_Support(ReachParameters, SystemDynamics, rVariable,s_per_thread_U, Min_Or_Max);
			//  **************    W_Support Function   ********************
			result1_minus = term2_minus;
			beta_minus = ReachParameters.result_beta;
			res_beta_minus = beta_minus * term3b_minus; //Replacing term3b from previous step
			result_minus = result1_minus + res_beta + term3c_minus;
			zV_min = result_minus;
			//  **************  W_Support Function Over  ********************
			s1Variable_min = sVariable_min + zV_min;
			r1Variable_minus = phi_trans_dir_minus;
			//zI = Omega_Support(ReachParameters, r1Variable, Initial,SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
			//  **************    Omega Function   ********************
			res1_minus = term1_minus;
			if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
				phi_tau_Transpose.mult_vector(r1Variable_minus,
						phi_trans_dir_minus);
				term1_minus = Initial->computeSupportFunction(
						phi_trans_dir_minus, s_per_thread_I_minus);
			}
			if (!SystemDynamics.isEmptyMatrixB) { //current_location's SystemDynamics's or ReachParameters
				B_trans.mult_vector(r1Variable_minus, Btrans_dir_minus);
				term2_minus = ReachParameters.time_step
						* SystemDynamics.U->computeSupportFunction(
								Btrans_dir_minus, s_per_thread_U_minus);
			}
			term3a_minus = ReachParameters.result_alfa;
			term3b_minus = support_unitball_infnorm(r1Variable_minus);
			if (!SystemDynamics.isEmptyC) {
				term3c_minus = ReachParameters.time_step
						* dot_product(SystemDynamics.C, rVariable_minus); //Added +tau* sf_C(l) 8/11/2015
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
			//Debug --- Invariant check
			//***** Intersection Detection Section *****
			//---taken from here

			rVariable = CopyVector(r1Variable); //source to destination		//Also works   rVariable=r1Variable
			sVariable = s1Variable;
			rVariable_minus = CopyVector(r1Variable_minus);
			sVariable_min = s1Variable_min;
		} //end of iterations
		//	if (invariantCrossed)	//We need to check all invariants as we are taking the min of all the directions
		//		break;//boundary crossed so no need to check other invariant's directions
	} //end of parallel for each Iterations or Time-Bound
// At the end of the For-Loop or all invariant_Directions we have boundaryIterations vector with the different limit to stop iterations
// The Minimum(boundaryIterations) will be the final Boundary Limit for the number of iterations
	unsigned int min_Total_Iteration;
	min_Total_Iteration = *min_element(boundaryIterations.begin(),
			boundaryIterations.end()); // - 1 ;//excluding the last outside Omega
//#pragma omp critical
//	{
	newTotIters = min_Total_Iteration;
//	}
//	cout<<"\nmin_Total_Iteration = "<<min_Total_Iteration<<endl;
//	return min_Total_Iteration;
}

void quickInvariantBoundaryCheck(Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		int lp_solver_type_choosen, unsigned int &newTotIters) {

	int dimension = ReachParameters.X0->getSystemDimension();
	int numberOfInvariants = invariant->getColumnVector().size(); //total number of Invariant's constraints
	std::vector<int> boundaryIters(numberOfInvariants,
			ReachParameters.Iterations); // size(dimension_size,initial_value)
	for (int eachInvariant = 0; eachInvariant < numberOfInvariants;
			eachInvariant++) {
		std::vector<double> pos_dir(dimension), neg_dir(dimension);
		//pos_dir.resize(dimension);	neg_dir.resize(dimension);
		//cout<<"Invariant Number = "<<eachInvariantDirection<<endl;
		//	cout<<"Invariant Directions ";
		for (int i = 0; i < dimension; i++) {
			pos_dir[i] = invariant->getCoeffMatrix()(eachInvariant, i); //First vector positive
			cout << pos_dir[i] << "\t";
			neg_dir[i] = -1 * invariant->getCoeffMatrix()(eachInvariant, i); //Second vector negative of First vector
		}
		// ******** Lamda Computation for each invariant's directions/constraints **********
		//Todo:: just get the boundvalue for the invariant direction without using LP_Solver

		double invariant_SupportFunction;
		lp_solver lpSolver(lp_solver_type_choosen);
		//	cout<<"Test 2 \n";
		lpSolver.setMin_Or_Max(2);
		if (!invariant->getIsEmpty()) { //set glpk constraints If not an empty polytope
			lpSolver.setConstraints(invariant->getCoeffMatrix(),
					invariant->getColumnVector(),
					invariant->getInEqualitySign());
			invariant_SupportFunction = invariant->computeSupportFunction(
					pos_dir, lpSolver);
			//std::cout << "\neachInvariantDirection = "<<eachInvariant<<" Invariant SupportFunction = " << invariant_SupportFunction << endl;
		}
		// **********************XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX******************
		boundaryIters[eachInvariant] = invariantCheck(pos_dir, neg_dir,
				invariant_SupportFunction, ReachParameters, SystemDynamics,
				Initial, lp_solver_type_choosen);
		cout << "boundaryIters[eachInvariant] = "
				<< boundaryIters[eachInvariant] << std::endl;
	}	//end of each Invariant directions

	//check if critical section required
	newTotIters = *min_element(boundaryIters.begin(), boundaryIters.end());
}

unsigned int invariantCheck(std::vector<double>& pos_dir, std::vector<double>& neg_dir, double SearchKey, ReachabilityParameters& ReachParameters, Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial, int lp_solver_type_choosen) {
	supportFunctionProvider::ptr p;
	unsigned int start_iter = 0, end_iter = ReachParameters.Iterations;
// ************************************************ object creations ************************************************
	lp_solver pos_lp(lp_solver_type_choosen);
	pos_lp.setMin_Or_Max(2);
	if (!Initial->getIsEmpty()) { //set glpk constraints If not an empty polytope
		pos_lp.setConstraints(ReachParameters.X0->getCoeffMatrix(), ReachParameters.X0->getColumnVector(),
				ReachParameters.X0->getInEqualitySign());
	}
	lp_solver neg_lp(lp_solver_type_choosen);
	neg_lp.setMin_Or_Max(2);

	if (!Initial->getIsEmpty()) //set glpk constraints If not an empty polytope
		neg_lp.setConstraints(ReachParameters.X0->getCoeffMatrix(), ReachParameters.X0->getColumnVector(), ReachParameters.X0->getInEqualitySign());
	// ************************************************ object creations ************************************************

	// *******************************************   Quick Check  *******************************************
	p = getInitialSet(ReachParameters.TimeBound, ReachParameters,
			SystemDynamics);	//initial set at time = timeBound
	double pos_val = p->computeSupportFunction(pos_dir, pos_lp);
	double neg_val = p->computeSupportFunction(neg_dir, neg_lp);
	cout << "pos_val = " << pos_val << "\t";
	cout << "neg_val = " << -1 * neg_val << "\t";
	cout << "SearchKey = " << SearchKey << std::endl;
	if (pos_val <= SearchKey && (-1 * neg_val) <= SearchKey) { //Last Omega completely inside the invariant
		cout << "Flowpipe completely inside Invariant!!" << std::endl;
		return ReachParameters.Iterations;//total iterations supplied by the user
	}
	//cout <<"Not coming here"<<std::endl;
	//other wise perform binary Search to find out the actual iterations that just crosses the invariant
	// *******************************************   Quick Check  *******************************************
	int searchSteps = 0;
	while (start_iter <= end_iter) {
		searchSteps++;
		unsigned int mid = (start_iter + end_iter - 1) / 2;	//iterations to be returned
		cout << "mid = " << mid << "\t";
		double timeBound = mid * ReachParameters.time_step;	//computing the timebound at that point
		p = getInitialSet(timeBound, ReachParameters, SystemDynamics);//initial set at time = timeBound
		double pos_val = p->computeSupportFunction(pos_dir, pos_lp);

		double neg_val = -1 * p->computeSupportFunction(neg_dir, neg_lp);
		cout << "pos_val = " << pos_val << "\t";
		cout << "neg_val = " << neg_val << "\t";
		if (pos_val > SearchKey && neg_val == SearchKey) {//fortunate Condition. May not occur always
			cout << "Search stopped at = " << searchSteps << " steps" << std::endl;
			//return (mid - 1);	//iterations found where neg_dir is just out on the line of invariant
			return mid;	//iterations found where neg_dir is just out on the line of invariant
			//so (mid - 1) will be the iteration inside or touching the invariant
		}//if only (pos_val > SearchKey) is true then need the next condition check
		if (neg_val > SearchKey && pos_val > SearchKey) {
			end_iter = mid - 1;
		}
		if (neg_val < SearchKey) {		//  && pos_val > SearchKey){
			start_iter = mid + 1;
		}
		int diff = end_iter - start_iter;//when they both converge close to eachother near the invariant bound
		if (abs(diff) <= 1) {//sufficiently closed to the boundValue and log(N) iterations over
			cout << "Search crossed at = " << searchSteps << " steps" << std::endl;
			return end_iter;	//or start_iter
		}
	}	//endwhile
	cout << "Complete Search Steps = " << searchSteps << std::endl;
}

supportFunctionProvider::ptr getInitialSet(double START_TIME, ReachabilityParameters& ReachParameters, Dynamics& SystemDynamics) {

	math::matrix<double> phi, phi_trans;

	if ((SystemDynamics.isEmptyMatrixB || SystemDynamics.U->getIsEmpty()) && SystemDynamics.isEmptyC){	//both B and C is empty, so we have x'(t) = Ax(t)
		//cout <<"Matrix B and Vector C is Empty!!, Dynamics is x'(t) = Ax(t)\n";
		SystemDynamics.MatrixA.matrix_exponentiation(phi, START_TIME); //if MatrixA is empty will not perform this function
		phi.transpose(phi_trans); //phi_trans computed
		//	std::cout << "\ncomputing initial object\n";
		supportFunctionProvider::ptr Initial = transMinkPoly::ptr( new transMinkPoly(ReachParameters.X0, SystemDynamics.U,
						phi_trans, phi_trans, 0, 0));//when  Bu(t) and C are empty
		return Initial;
	}
	cout <<"Dynamics is NOT of type --:  x'(t) = Ax(t)\n";
	math::matrix<double> A_inv_phi, y_matrix, y_trans;
	SystemDynamics.MatrixA.matrix_exponentiation(phi, START_TIME); //if MatrixA is empty will not perform this function
	phi.transpose(phi_trans); //phi_trans computed	//cout << "phi_trans" << phi_trans << std::endl;
	math::matrix<double> A_inv; //(SystemDynamics.MatrixA.size1(),SystemDynamics.MatrixA.size2());
	A_inv = ReachParameters.A_inv;
	A_inv.multiply(phi, A_inv_phi);	//cout<<"A_inv_phi = "<<A_inv_phi<<std::endl;
	A_inv_phi.minus(A_inv, y_matrix);	//cout << "y_matrix = " << y_matrix << std::endl;
	y_matrix.transpose(y_trans);
	//	std::cout << "\ncomputing initial object\n";

	if (SystemDynamics.isEmptyC) {
		cout << "C is Empty\n";
		supportFunctionProvider::ptr Initial = transMinkPoly::ptr(new transMinkPoly(ReachParameters.X0, SystemDynamics.U, phi_trans, y_trans, 1, 0));	//when only C is empty
		return Initial;
	} else if (!SystemDynamics.isEmptyC) {
		cout << "C is NOT Empty\n";
		supportFunctionProvider::ptr Initial = transMinkPoly::ptr(new transMinkPoly(ReachParameters.X0, SystemDynamics.U, SystemDynamics.C, phi_trans, y_trans, 1, 0));
		return Initial;
	}
}

//NOTE :: this approach of increasing the time-step does not work in support-function algorithm
void SlowStartInvariantBoundaryCheck(Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParametersOld, polytope::ptr invariant,
		int lp_solver_type_choosen, unsigned int &newTotIters) {

	ReachabilityParameters ReachParameters;	//local variable to work with
	ReachParameters = ReachParametersOld;	//copy of Original

	bool slowStarted = true;	//slowStart started
	unsigned int shm_NewTotalIteration = ReachParameters.Iterations; //Shared Variable for resize iterations number on crossing with invariant
	//	cout<<"shm_NewTotalIteration = "<<shm_NewTotalIteration<<"OK\n";
	int dimension = Initial->getSystemDimension();
	int Min_Or_Max = 2;
	int numberOfInvariants = invariant->getColumnVector().size(); //total number of Invariant's constraints
	/*cout <<"invariant->getCoeffMatrix() = "<<invariant->getCoeffMatrix()<<"\n";
	 for (int i=0;i<invariant->getColumnVector().size();i++)
	 cout<<invariant->getColumnVector()[i]<<"\t";*/

	// ******************* Probable race condition variables:: can be  placed inside for-loop *******************
	int foundStart = 0, intersection_start, intersection_end;
	//bool invariantCrossed = false;
	// *************************** For Positive ************************************
	double res1, term2 = 0.0, result1, term1 = 0.0, result;
	std::vector<double> Btrans_dir, phi_trans_dir, phi_trans_dir1;
	// *************************** For Negative ************************************
	double res1_minus, term2_minus = 0.0, result1_minus, term1_minus = 0.0,
			result_minus;
	std::vector<double> Btrans_dir_minus, phi_trans_dir_minus,
			phi_trans_dir1_minus;
	math::matrix<double> B_trans, phi_tau_Transpose;
	if (!SystemDynamics.isEmptyMatrixA) //current_location's SystemDynamics's or ReachParameters
		phi_tau_Transpose = ReachParameters.phi_trans;
	if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
		B_trans = ReachParameters.B_trans;

	// *******************************************************************************************************************

	std::vector<int> boundaryIterations(numberOfInvariants,
			shm_NewTotalIteration); // size(dimension_size,initial_value)

//#pragma omp parallel for //num_threads(numberOfInvariants)
	for (int eachInvariantDirection = 0;
			eachInvariantDirection < numberOfInvariants;
			eachInvariantDirection++) {

		bool TimeBoundFound = false;//final result or actual flowpipe_cost found

		double zIInitial, zIInitial_minus, zI = 0.0, zV = 0.0;
		double zI_min = 0.0, zV_min = 0.0;
		double sVariable, s1Variable;
		double sVariable_min, s1Variable_min; //For Minimization of First Vector only

		std::vector<double> r1Variable, r1Variable_minus;
		r1Variable.resize(dimension);
		r1Variable_minus.resize(dimension);

		std::vector<double> rVariable, rVariable_minus;
		rVariable.resize(dimension);
		rVariable_minus.resize(dimension);
		//	cout<<"Invariant Number = "<<eachInvariantDirection<<endl;
		for (int i = 0; i < dimension; i++) {
			rVariable[i] = invariant->getCoeffMatrix()(eachInvariantDirection,
					i); //First vector positive
			rVariable_minus[i] = -1
					* invariant->getCoeffMatrix()(eachInvariantDirection, i); //Second vector negative of First vector
		}
		// ******** Lamda Computation for each invariant's directions/constraints **********
		double invariant_SupportFunction;
		int type = lp_solver_type_choosen;
		lp_solver lpSolver(type), lp_U_dummy(type);
		//	cout<<"Test 2 \n";
		lpSolver.setMin_Or_Max(2);
		if (!invariant->getIsEmpty() && !invariant->getIsUniverse()) { //set glpk constraints If not an empty polytope
			lpSolver.setConstraints(invariant->getCoeffMatrix(),
					invariant->getColumnVector(),
					invariant->getInEqualitySign());
			invariant_SupportFunction = invariant->computeSupportFunction(
					rVariable, lpSolver);
			std::cout << "\neachInvariantDirection = " << eachInvariantDirection
					<< " Invariant SupportFunction = "
					<< invariant_SupportFunction << endl;
		}
		// **********************XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX******************

		lp_solver s_per_thread_I(type), s_per_thread_U(type);
// ************************************************ For Positive Direction Starts ************************************************
		s_per_thread_I.setMin_Or_Max(2);
		if (!Initial->getIsEmpty()) { //set glpk constraints If not an empty polytope
			s_per_thread_I.setConstraints(ReachParameters.X0->getCoeffMatrix(),
					ReachParameters.X0->getColumnVector(),
					ReachParameters.X0->getInEqualitySign());
			//s_per_thread_I.setConstraints(Initial->getCoeffMatrix(),ReachParameters.X0->getColumnVector(),ReachParameters.X0->getInEqualitySign());
		}
		s_per_thread_U.setMin_Or_Max(2);
		if (SystemDynamics.U->getIsEmpty()) { //empty polytope
			//Polytope is empty so no glpk object constraints to be set
			//std::cout << "InvariantBoundCheck: Input set U is empty\n";
		} else {
			s_per_thread_U.setConstraints(SystemDynamics.U->getCoeffMatrix(),
					SystemDynamics.U->getColumnVector(),
					SystemDynamics.U->getInEqualitySign());
		}
// ************************************************   Positive Direction Ends *******************************************
		lp_solver s_per_thread_I_minus(type), s_per_thread_U_minus(type);
// ******************************************* For Negative Direction Starts *******************************************
		s_per_thread_I_minus.setMin_Or_Max(2);
		if (!Initial->getIsEmpty()) //set glpk constraints If not an empty polytope
			s_per_thread_I_minus.setConstraints(
					ReachParameters.X0->getCoeffMatrix(),
					ReachParameters.X0->getColumnVector(),
					ReachParameters.X0->getInEqualitySign());

		s_per_thread_U_minus.setMin_Or_Max(2);
		if (SystemDynamics.U->getIsEmpty()) { //empty polytope
			//Polytope is empty so no glpk object constraints to be set
		} else {
			s_per_thread_U_minus.setConstraints(
					SystemDynamics.U->getCoeffMatrix(),
					SystemDynamics.U->getColumnVector(),
					SystemDynamics.U->getInEqualitySign());
		}
// ******************************************* Negative Direction Ends *******************************************
		unsigned int loopIteration = 0;
		sVariable = 0.0; //initialize s0
		sVariable_min = 0.0;
		// ************************************************  For Positive Direction Starts *******************************************
		//zIInitial = Omega_Support(ReachParameters, rVariable, Initial,SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
		double term3, term3a, term3b, res2, term3c = 0.0;
		//  **************    Omega Function   ********************
		res1 = Initial->computeSupportFunction(rVariable, s_per_thread_I);
		if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
			phi_tau_Transpose.mult_vector(rVariable, phi_trans_dir);
			term1 = Initial->computeSupportFunction(phi_trans_dir,
					s_per_thread_I);
		}
		if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
			B_trans.mult_vector(rVariable, Btrans_dir);
		if (!SystemDynamics.isEmptyMatrixB && !SystemDynamics.U->getIsEmpty())
			term2 = ReachParameters.time_step
					* SystemDynamics.U->computeSupportFunction(Btrans_dir,
							s_per_thread_U);
		term3a = ReachParameters.result_alfa;
		term3b = (double) support_unitball_infnorm(rVariable);
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
		//	cout<<"res2 = "<<res2<<" and ";
		//  **************  Omega Function Over  ********************
		// ************************************************   Positive Direction Ends *******************************************

		// ******************************************* For Negative Direction Starts *******************************************
		//zIInitial = Omega_Support(ReachParameters, rVariable_minus, Initial,SystemDynamics, s_per_thread_I_minus, s_per_thread_U_minus, Min_Or_Max);
		double term3_minus, term3a_minus, term3b_minus, res2_minus,
				term3c_minus = 0.0;
		//  **************    Omega Function   ********************
		res1_minus = Initial->computeSupportFunction(rVariable_minus,
				s_per_thread_I_minus);
		if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
			phi_tau_Transpose.mult_vector(rVariable_minus, phi_trans_dir_minus);
			term1_minus = Initial->computeSupportFunction(phi_trans_dir_minus,
					s_per_thread_I_minus);
		}
		if (!SystemDynamics.isEmptyMatrixB) //current_location's SystemDynamics's or ReachParameters
			B_trans.mult_vector(rVariable_minus, Btrans_dir_minus);
		if (!SystemDynamics.isEmptyMatrixB && !SystemDynamics.U->getIsEmpty())
			term2_minus = ReachParameters.time_step
					* SystemDynamics.U->computeSupportFunction(Btrans_dir_minus,
							s_per_thread_U_minus);
		term3a_minus = ReachParameters.result_alfa;
		term3b_minus = (double) support_unitball_infnorm(rVariable_minus);
		if (!SystemDynamics.isEmptyC) {
			term3c_minus = ReachParameters.time_step
					* dot_product(SystemDynamics.C, rVariable_minus); //Added +tau* sf_C(l) 8/11/2015
		}
		term3_minus = term3a_minus * term3b_minus;
		res2_minus = term1_minus + term2_minus + term3_minus + term3c_minus;
		if (res1_minus > res2_minus)
			zIInitial_minus = res1_minus;
		else
			zIInitial_minus = res2_minus;
		//	cout<<"res2_minus = "<<res2_minus<<"\n";
//  **************  Omega Function Over  ********************
// ******************************************* Negative Direction Ends *******************************************
		loopIteration++;
		//for (; loopIteration < shm_NewTotalIteration ;) { //Now stopping condition is NOT only "shm_NewTotalIteration"

		double time_stepBeforeCrossing, time_stepAfterCrossing;	//time_step is used to determine the crossing
		int foundIntersectionEnd_After_jump = 0;

		for (;
				ReachParameters.time_step <= ReachParametersOld.TimeBound
						&& TimeBoundFound != true;) {
			//Now stopping condition is the omega that actually cross.  Trick is dealing with time-step
// ************************************************  For Positive Direction Starts *******************************************
			double TempOmega, term3, term3a, res2, beta, res_beta;
			//	zV = W_Support(ReachParameters, SystemDynamics, rVariable,s_per_thread_U, Min_Or_Max);
			//  **************    W_Support Function   ********************
			result1 = term2;
			beta = ReachParameters.result_beta;
			res_beta = beta * term3b; //Replacing term3b from previous step
			result = result1 + res_beta + term3c; //Added term3c
			zV = result;
			//  **************  W_Support Function Over  ********************
			s1Variable = sVariable + zV;
			r1Variable = phi_trans_dir;
			//zI = Omega_Support(ReachParameters, r1Variable, Initial,SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
			//  **************    Omega Function   ********************
			res1 = term1;
			if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
				phi_tau_Transpose.mult_vector(r1Variable, phi_trans_dir);
				term1 = Initial->computeSupportFunction(phi_trans_dir,
						s_per_thread_I);
			}
			if (!SystemDynamics.isEmptyMatrixB) { //current_location's SystemDynamics's or ReachParameters
				B_trans.mult_vector(r1Variable, Btrans_dir);
				term2 = ReachParameters.time_step
						* SystemDynamics.U->computeSupportFunction(Btrans_dir,
								s_per_thread_U);
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
			cout << "TempOmega = " << TempOmega << "\t";
// ************************************************   Positive Direction Ends *******************************************
// ************************************************  For Negative Direction Starts *******************************************
			double TempOmega_min, term3_minus, term3a_minus, res2_minus,
					beta_minus, res_beta_minus;
			//	zV = W_Support(ReachParameters, SystemDynamics, rVariable,s_per_thread_U, Min_Or_Max);
			//  **************    W_Support Function   ********************
			result1_minus = term2_minus;
			beta_minus = ReachParameters.result_beta;
			res_beta_minus = beta_minus * term3b_minus; //Replacing term3b from previous step
			result_minus = result1_minus + res_beta + term3c_minus;
			zV_min = result_minus;
			//  **************  W_Support Function Over  ********************
			s1Variable_min = sVariable_min + zV_min;
			r1Variable_minus = phi_trans_dir_minus;
			//zI = Omega_Support(ReachParameters, r1Variable, Initial,SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
			//  **************    Omega Function   ********************
			res1_minus = term1_minus;
			if (!SystemDynamics.isEmptyMatrixA) { //current_location's SystemDynamics's or ReachParameters
				phi_tau_Transpose.mult_vector(r1Variable_minus,
						phi_trans_dir_minus);
				term1_minus = Initial->computeSupportFunction(
						phi_trans_dir_minus, s_per_thread_I_minus);
			}
			if (!SystemDynamics.isEmptyMatrixB) { //current_location's SystemDynamics's or ReachParameters
				B_trans.mult_vector(r1Variable_minus, Btrans_dir_minus);
				term2_minus = ReachParameters.time_step
						* SystemDynamics.U->computeSupportFunction(
								Btrans_dir_minus, s_per_thread_U_minus);
			}
			term3a_minus = ReachParameters.result_alfa;
			term3b_minus = support_unitball_infnorm(r1Variable_minus);
			if (!SystemDynamics.isEmptyC) {
				term3c_minus = ReachParameters.time_step
						* dot_product(SystemDynamics.C, rVariable_minus); //Added +tau* sf_C(l) 8/11/2015
			}
			term3_minus = term3a_minus * term3b_minus;
			res2_minus = term1_minus + term2_minus + term3_minus + term3c_minus;
			if (res1_minus > res2_minus)
				zI_min = res1_minus;
			else
				zI_min = res2_minus;
			//  **************  Omega Function Over  ********************
			TempOmega_min = zI_min + s1Variable_min; //Y1
			cout << "TempOmega_min = " << TempOmega_min << "\t";
			// ************************************************   Negative Direction Ends *******************************************
			loopIteration++; // Placed here as Omega_0 and Omega_1 is computed so loopIteration value == 2
// ******************************************* Intersection Detection Section Starts *******************************************
			if (((-1 * TempOmega_min) < invariant_SupportFunction)
					&& (invariant_SupportFunction < TempOmega)) { //Should have been correct
				intersection_start = loopIteration; //Partially inside the region
				//	cout<<"\nINTERSECTION_start = "<<intersection_start <<endl;
				//}
			}
			bool foundIntersectionEnd_jump = false;
			cout << "\time_stepBeforeCrossing = " << time_stepBeforeCrossing
					<< " ReachParameters.time_step = "
					<< ReachParameters.time_step << "\n";
			if ((-1 * TempOmega_min) > invariant_SupportFunction) { // Should have been correct
				//if ( ((-1 * TempOmega_min) < invariant_SupportFunction) && (TempOmega < invariant_SupportFunction)) {
				intersection_end = loopIteration;
				foundIntersectionEnd_jump = true;
				time_stepAfterCrossing = ReachParameters.time_step;	//the time_step at which it crossed
				foundIntersectionEnd_After_jump++;//1st found == 1, 2nd found is sequential found
				cout << "\nIntersection_End = " << intersection_end
						<< " foundIntersectionEnd_After_jump = "
						<< foundIntersectionEnd_After_jump << endl;

				//	cout<<"TempOmega_min = "<<TempOmega_min<< " and TempOmega = "<<TempOmega<<"\t";
				//????	boundaryIterations[eachInvariantDirection] = loopIteration; //Made Changes here due to circle
				//	invariantCrossed = true;//Omega is out of the invariant's boundary
				//????	break; //no need to compute reachable set anymore due to out of invariant boundary
			}

			if (foundIntersectionEnd_jump) {
				/*
				 * when Slow Start finds the crossing point
				 * Now call a separate/may be same function to either perform sequential increment of time-step
				 * or half on reverse directions
				 */
				//	break; //no need to compute reachable set anymore due to out of invariant boundary
				slowStarted = false;//STOP slow start as it jumped the invariant
			} else {	//if not crossed
				time_stepBeforeCrossing = ReachParameters.time_step; //the time_step at which it did NOT cross
			}
			if (foundIntersectionEnd_After_jump == 2) {	//found sequential search
				boundaryIterations[eachInvariantDirection] =
						ReachParameters.time_step
								/ ReachParametersOld.time_step;
				TimeBoundFound = true;
			}

			if (slowStarted) {	// slow started procedure
				slowStarted = true;	//slow start started
				ReachParameters.time_step = ReachParameters.time_step * 2;//doubled every time. initially 2 steps done sequentially
				if (ReachParameters.time_step >= ReachParametersOld.TimeBound)//increment cross the given timeBound
					ReachParameters.time_step = ReachParametersOld.TimeBound;//so handle it by assigning the given timeBound
			} else { //slow start has stopped so time-step to be set as original
				cout << "NOOOOOOOTTTT Entered!!! === ";
				ReachParameters.time_step = time_stepBeforeCrossing
						+ ReachParametersOld.time_step;
				cout << " ReachParameters.time_step = "
						<< ReachParameters.time_step << "\n";

				//Method 1) it is increasing sequentially
			}

// ******************************************* Intersection Detection Section Ends *******************************************
			rVariable = CopyVector(r1Variable); //source to destination		//Also works   rVariable=r1Variable
			sVariable = s1Variable;
			rVariable_minus = CopyVector(r1Variable_minus);
			sVariable_min = s1Variable_min;
		} //end of iterations
		//	if (invariantCrossed)	//We need to check all invariants as we are taking the min of all the directions
		//		break;//boundary crossed so no need to check other invariant's directions
	} //end of parallel for each Iterations or Time-Bound
// At the end of the For-Loop or all invariant_Directions we have boundaryIterations vector with the different limit to stop iterations
// The Minimum(boundaryIterations) will be the final Boundary Limit for the number of iterations
	unsigned int min_Total_Iteration;
	min_Total_Iteration = *min_element(boundaryIterations.begin(),
			boundaryIterations.end()); // - 1 ;//excluding the last outside Omega
//#pragma omp critical
//	{
	newTotIters = min_Total_Iteration;
//	}
//	cout<<"\nmin_Total_Iteration = "<<min_Total_Iteration<<endl;
//	return min_Total_Iteration;
}

void jumpInvariantBoundaryCheck(Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant, int lp_solver_type_choosen, unsigned int &newTotIters) {

	//-- get the time-horizon and time-step and discritization_factor=10 times
	// compute the iteration using this new time-step
	double widening_Factor = 5;	//10 times
	double time_step = ReachParameters.time_step * widening_Factor;
	double time_horizon = ReachParameters.TimeBound;
	double start_time = 0, time_crossed, next_search_time, actual_time_bound;//last time when it is still inside //cout<<"test 1\n";
std::string filename="./coarse.out";

	//time_crossed = invariantCrossingCheck(start_time, time_step, time_horizon, Initial, ReachParameters, invariant, SystemDynamics, lp_solver_type_choosen);
	time_crossed = invariantCrossingCheck1(start_time, time_step, time_horizon, Initial, ReachParameters, invariant,
			SystemDynamics, lp_solver_type_choosen, filename);
	//cout<<"First coarse time-step completed at time = "<< time_crossed<<"\n";
	if (time_crossed == time_horizon){	//Search failed to find a crossing. Flowpipe completely inside Invariant
		newTotIters = ReachParameters.Iterations;	//orginal iteration
		//cout<<"Invariant Condition completed at coarse time-step\n";
	}else {	//if (time_crossed < time_horizon) meaning it has just crossed with widening_Factor but now to get precise crossing-time we use fine-time-step
		//cout<<"Fine time-step testing invoked!!!\n";
		next_search_time =(time_crossed - time_step) + ReachParameters.time_step;
		time_step = ReachParameters.time_step;	//setting the original time-step as fine-time-step
		filename="./fine.out";
		actual_time_bound = invariantCrossingCheck1(next_search_time, time_step, time_crossed, Initial, ReachParameters,
				invariant, SystemDynamics, lp_solver_type_choosen, filename);
		//Debugging the correct time
			//Since I am following '0-indexing' throughout. So it must return the value/index N as its indices are (0...(N-1))
			//since the function invariantCrossingCheck1(...) return crossing time bound 'N' starting from time at '0'
			//so, we should at one more time-step to return the size as 'N+1' for array indexing
			//For eg., return 1 if only Omega_0 is inside, 2 if Omega_0 and Omega_1 are inside
			//TODO::Debugging CONCLUDED that math::ceil() function rounds up when values are >=5 to its ceiling. which creates problem
					//But if we do not use math::ceil() then normal division on double data type gives in-correct result in g++/gcc.
				//So we use math::ceil() as an extra iterations or over-approximation returned is acceptable in our case.
		//---
		//cout << "actual_time_bound = " << actual_time_bound << " ReachParameters.time_step = " << ReachParameters.time_step;
		//cout << "  actual_time_bound = " << (actual_time_bound + time_step);
		newTotIters = math::ceil((actual_time_bound + time_step)/ ReachParameters.time_step);
		//cout << "  Iterations = " << newTotIters << "\n";
	}
}

//Returns THE time when a convex set just crosses the invariant boundary one after another of a invariant polyhedra.
//Note: This is an in-efficient method
double invariantCrossingCheck(double START_TIME, double time_step, double time_horizon, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant, Dynamics& SystemDynamics, int lp_solver_type_choosen) {
	int dimension = ReachParameters.X0->getSystemDimension();
	int numberOfInvariants = invariant->getColumnVector().size(); //total number of Invariant's constraints
	//unsigned int tot_iters = math::ceil(time_horizon / time_step);
	std::vector<double> boundaryIters(numberOfInvariants, time_horizon); // size(dimension_size,initial_value)
	for (int eachInvariant = 0; eachInvariant < numberOfInvariants; eachInvariant++) {
		std::vector<double> pos_dir(dimension), neg_dir(dimension);		//	neg_dir.resize(dimension);		//cout<<"Invariant Number = "<<eachInvariantDirection<<endl;		//	cout<<"Invariant Directions ";
	//	cout<<"eachInvariant = "<<eachInvariant<<"\n";
		for (int i = 0; i < dimension; i++) {
			pos_dir[i] = invariant->getCoeffMatrix()(eachInvariant, i); //First vector positive			cout << pos_dir[i] << "\t";
			neg_dir[i] = -1 * invariant->getCoeffMatrix()(eachInvariant, i); //Second vector negative of First vector
		}
		double invariant_SupportFunction;
		//invariant_SupportFunction = invariant->getColumnVector(eachInvariant);	//this is same as solving LP but better use support_function concept
		lp_solver lpSolver(lp_solver_type_choosen);
		lpSolver.setMin_Or_Max(2);
		if (!invariant->getIsEmpty()) { //set glpk constraints If not an empty polytope
			lpSolver.setConstraints(invariant->getCoeffMatrix(), invariant->getColumnVector(), invariant->getInEqualitySign());
			invariant_SupportFunction = invariant->computeSupportFunction(pos_dir, lpSolver);	//std::cout << "\neachInvariantDirection = "<<eachInvariant<<" Invariant SupportFunction = " << invariant_SupportFunction << endl;
		}
		//cout<<"Test 2:: invariant_SupportFunction = "<<invariant_SupportFunction<<"\n";
		boundaryIters[eachInvariant] = invariantFaceCrossingCheck(neg_dir, invariant_SupportFunction, START_TIME, time_step, time_horizon, Initial, ReachParameters,
				SystemDynamics, lp_solver_type_choosen);//cout << "boundaryIters[eachInvariant] = " << boundaryIters[eachInvariant] << std::endl;
	}	//end of each Invariant directions
	return *min_element(boundaryIters.begin(), boundaryIters.end());
}

double invariantFaceCrossingCheck(std::vector<double>& neg_dir, double invBound, double START_TIME, double time_step, double time_horizon,
		supportFunctionProvider::ptr Initial, ReachabilityParameters& ReachParameters, Dynamics& SystemDynamics, int lp_solver_type_choosen){

	supportFunctionProvider::ptr p;
	// ************************************************ object creations ************************************************
	lp_solver neg_lp(lp_solver_type_choosen);
	neg_lp.setMin_Or_Max(2);
	if (!Initial->getIsEmpty()) //set glpk constraints If not an empty polytope
		neg_lp.setConstraints(ReachParameters.X0->getCoeffMatrix(), ReachParameters.X0->getColumnVector(), ReachParameters.X0->getInEqualitySign());
	//cout<<"problem here\n";
	// ************************************************ object creations ************************************************
	double neg_val;
	for(double t=START_TIME; t<=time_horizon; t = t + time_step){
		//cout<<"test a\n";
		p = getInitialSet(t, ReachParameters, SystemDynamics);	//initial set at time = t
		neg_val = -1 * p->computeSupportFunction(neg_dir, neg_lp);
		if (neg_val > invBound) {
			return t;	//crossed the invariant bound
		}
	}
	return time_horizon; //crossed time-horizon without crossing the invariant bound or didnot find intersecting with invariant face
}

double invariantCrossingCheck1(double START_TIME, double time_step, double time_horizon, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant, Dynamics& SystemDynamics, int lp_solver_type_choosen, std::string fileName){
	supportFunctionProvider::ptr p;
	// ************************************************ object creations ************************************************
	lp_solver neg_lp(lp_solver_type_choosen);
	neg_lp.setMin_Or_Max(2);
	if (!Initial->getIsEmpty()) //set glpk constraints If not an empty polytope
		neg_lp.setConstraints(ReachParameters.X0->getCoeffMatrix(), ReachParameters.X0->getColumnVector(), ReachParameters.X0->getInEqualitySign());
	//cout<<"problem here\n";
	// ************************************************ object creations ************************************************
	int dimension = ReachParameters.X0->getSystemDimension();
	int numberOfInvariants = invariant->getColumnVector().size(); //total number of Invariant's constraints
	double neg_val;
	bool flag;
/*
//Not required at runtime only for Printing and Debugging
	std::ofstream outFile;
	outFile.open(fileName.c_str());
*/
	for(double t=START_TIME; t<=time_horizon; t = t + time_step){		//cout<<"test a\n";
		flag=false;
		p = getInitialSet(t, ReachParameters, SystemDynamics);	//initial set at time = t
		double invBound;
		std::vector<double> neg_dir(dimension);
		for (int eachInvariant = 0; eachInvariant < numberOfInvariants; eachInvariant++) {
			for (int j = 0; j < dimension; j++)
				neg_dir[j] = -1 * invariant->getCoeffMatrix()(eachInvariant, j);
			neg_val = -1 * p->computeSupportFunction(neg_dir, neg_lp);
			invBound = invariant->getColumnVector()[eachInvariant];	//if this does not works use LP solver
			if (neg_val > invBound) {
				flag=true; //crossed the invariant bound
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
		if (flag){	//this condition is created only to be able to print the crossing convex set
			//outFile.close();//Not required at runtime only for Printing and Debugging
			return t;
		}
	}
	//outFile.close();//Not required at runtime only for Printing and Debugging
	return time_horizon; //crossed time-horizon without crossing the invariant bound or didnot find intersecting with invariant face
}

polytope::ptr getBoundedConvexSet(supportFunctionProvider::ptr& p, ReachabilityParameters& ReachParameters){
	math::matrix<double> tempDirs;
	std::vector<double> bounds(ReachParameters.Directions.size1());
	tempDirs=ReachParameters.Directions;
//	cout<<"	ReachParameters.Directions = "<< ReachParameters.Directions<<std::endl;
	lp_solver lp(GLPK_SOLVER);
	lp.setMin_Or_Max(2);//Maximizing

	if (!ReachParameters.X0->getIsEmpty()){	//checking if (!p->getIsEmpty()) may return false as only X0 is not empty but U is empty
		lp.setConstraints(ReachParameters.X0->getCoeffMatrix(),ReachParameters.X0->getColumnVector(),ReachParameters.X0->getInEqualitySign());
	}
	std::vector<double> dir(ReachParameters.Directions.size2());
	for (int d = 0; d < tempDirs.size1(); d++) {
		for (int j = 0; j < tempDirs.size2(); j++){
			dir[j] = ReachParameters.Directions(d, j);
			//cout<<"  "<<dir[j]<<"\t";
		}
		bounds[d]=p->computeSupportFunction(dir,lp);
		//cout<<"		bounds[d] = "<<	bounds[d]<<std::endl;
	}
	return (polytope::ptr(new polytope(tempDirs,bounds,1)));	//1 is <= inEquality sign
}

