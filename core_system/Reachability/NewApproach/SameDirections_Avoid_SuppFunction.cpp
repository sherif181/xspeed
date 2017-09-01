/*
 * SameDirections_Avoid_SuppFunction.cpp
 *
 *  Created on: 12-Mar-2015
 *      Author: amit
 */

#include "core_system/Reachability/NewApproach/SameDirections_Avoid_SuppFunction.h"

template_polyhedra::ptr reachabilitySameDirection(Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen) {

	int numVectors = ReachParameters.Directions.size1();
	int dimension = Initial->getSystemDimension();
	unsigned int shm_NewTotalIteration = ReachParameters.Iterations; //Shared Variable for resize iterations number on crossing with invariant
	int Min_Or_Max = 2;

	math::matrix<double> MatrixValue;	//Shared Matrix for all child thread
	size_type row = numVectors, col = shm_NewTotalIteration;
	if (isInvariantExist == true) {	//if invariant exist. Computing

		InvariantBoundaryCheck(SystemDynamics, Initial,
				ReachParameters, invariant, lp_solver_type_choosen, shm_NewTotalIteration);
	}	//End of Invariant Directions

	if (shm_NewTotalIteration == 1) {
		template_polyhedra::ptr poly_emptyp;
		return poly_emptyp;
	}

	col = shm_NewTotalIteration;	//if invariant exist col will be resized
	MatrixValue.resize(row, col);
	int solver_type = lp_solver_type_choosen;
	lp_solver s_per_thread_I(solver_type), s_per_thread_U(solver_type),
			s_per_thread_inv(solver_type);

	s_per_thread_I.setMin_Or_Max(2);
	if (!ReachParameters.X0->getIsEmpty())//set glpk constraints If not an empty polytope
		s_per_thread_I.setConstraints(ReachParameters.X0->getCoeffMatrix(),
				ReachParameters.X0->getColumnVector(),
				ReachParameters.X0->getInEqualitySign());

	s_per_thread_U.setMin_Or_Max(2);
	if (SystemDynamics.U->getIsEmpty()) {	//empty polytope
		//Polytope is empty so no glpk object constraints to be set
	} else {
		s_per_thread_U.setConstraints(SystemDynamics.U->getCoeffMatrix(),
				SystemDynamics.U->getColumnVector(),
				SystemDynamics.U->getInEqualitySign());
	}
unsigned int temp=0;
//cout<<"\nTESTing Parallel"<<endl;
//#pragma omp parallel for
	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {
		// *************************
		bool r_same_r1 = false;
		Optimize_Omega_Support Optimize_Omega;
		Optimize_Omega.input_epsilon = 0.01;	//difference of angle between r and r1;	//can take inputs from USER
		// *************************
		double zIInitial = 0.0, zI = 0.0, zV = 0.0;
		double sVariable, s1Variable;
		//std::cout<<"\nCheck 1"<<endl;
		std::vector<double> r1Variable;	//now single dimension
		r1Variable.resize(dimension);
		std::vector<double> rVariable;
		rVariable.resize(dimension);
		//	cout<<"(";
		for (int i = 0; i < dimension; i++) {
			rVariable[i] = ReachParameters.Directions(eachDirection, i);
			//	cout<<ReachParameters.Directions(eachDirection, i)<<" , ";
		}
		//	cout<<")"<<endl;
		unsigned int loopIteration = 0;
		sVariable = 0.0; 		//initialize s0
		//cout<<"\nOmega_Support(EachDirection) = "<<eachDirection<<endl;
		r_same_r1 = false;		//before entering inside the LOOP
		zIInitial = Omega_Support_Similar_Direction(ReachParameters, rVariable, Initial,
				SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max, r_same_r1, Optimize_Omega);
		MatrixValue(eachDirection, loopIteration) = zIInitial;
		loopIteration++;
		for (; loopIteration < shm_NewTotalIteration;) { //Now stopping condition is only "shm_NewTotalIteration"
			double TempOmega;
		//	std::cout<<"\tHello = "<<loopIteration;
			ReachParameters.phi_trans.mult_vector(rVariable, r1Variable);	//r1 computed
			/*
			 * now compare if r1 is similar to r  for computing Omega_Support (later can think of W_Support)
			 */
			double compute_epsilon = compute_theta(Optimize_Omega.dir1, r1Variable);
			if (compute_epsilon <= Optimize_Omega.input_epsilon){
				Optimize_Omega.computed_epsilon = compute_epsilon;
				r_same_r1 = true;
			//	temp++;
			}else{
				r_same_r1 = false;
			}

			/** Precompute beta and send it as parameter */
			zV = W_Support(ReachParameters, SystemDynamics, rVariable,
					s_per_thread_U, Min_Or_Max);
			//	cout<<"zV= "<<zV<<"\t";
			s1Variable = sVariable + zV;
			zI = Omega_Support_Similar_Direction(ReachParameters, r1Variable, Initial,
					SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max, r_same_r1, Optimize_Omega);
			//		cout<<"zi= "<<zI<<"\t";
			TempOmega = zI + s1Variable; 		//Y1
			MatrixValue(eachDirection, loopIteration) = TempOmega; 		//Y1
			rVariable = CopyVector(r1Variable);		//source to destination
			sVariable = s1Variable;
			loopIteration++;	//for the next Omega-iteration or Time-bound
		}	//end of while for each vector
	}
	//cout<<"temp = "<<temp<<endl;
//cout<<"Sequential Algorithm :: Time Step :="<<ReachParameters.time_step<<endl;
	//cout<<"Outside Parallel\n";
	if (isInvariantExist == true) {		//if invariant exist. Computing
		math::matrix<double> inv_sfm;
		int num_inv = invariant->getColumnVector().size();//number of Invariant's constriants
		inv_sfm.resize(num_inv, shm_NewTotalIteration);
		for (int eachInvDirection = 0; eachInvDirection < num_inv;
				eachInvDirection++) {
			for (unsigned int i = 0; i < shm_NewTotalIteration; i++) {
				inv_sfm(eachInvDirection, i) = invariant->getColumnVector()[eachInvDirection];
			}
		}
		return template_polyhedra::ptr (new template_polyhedra(MatrixValue, inv_sfm,
				ReachParameters.Directions, invariant->getCoeffMatrix()));
	} else {
		return template_polyhedra::ptr (new template_polyhedra(MatrixValue, ReachParameters.Directions));
	}
}
