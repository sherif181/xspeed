/*
 * reachability_Sequential_AllDirections.cpp
 *
 *  Created on: 23-Mar-2015
 *      Author: amit
 */

#include "core_system/Reachability/reachability_Sequential_AllDirections.h"
/*

template_polyhedra reachabilitySequentialAllDirections(Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen) {

	std::vector<D> AllDirections;	//if invariant exist then need not compute all transposed directions

	int numVectors = ReachParameters.Directions.size1();
	int dimension = Initial->getSystemDimension();
	unsigned int shm_NewTotalIteration = ReachParameters.Iterations;//Shared Variable for resize iterations number on crossing with invariant
	int Min_Or_Max = 2;

	math::matrix<double> MatrixValue;	//Shared Matrix for all child thread
	size_type row = numVectors, col = shm_NewTotalIteration;
	if (isInvariantExist == true) {	//if invariant exist. Computing
		shm_NewTotalIteration = InvariantBoundaryCheck(SystemDynamics, Initial,
				ReachParameters, invariant, lp_solver_type_choosen);
	}	//End of Invariant Directions

if (shm_NewTotalIteration==1){
	template_polyhedra poly_emptyp;
	return poly_emptyp;
}

	col = shm_NewTotalIteration;	//if invariant exist Iterations value will change
	AllDirections = get_directions(ReachParameters, shm_NewTotalIteration);	//compute all directions



//	MatrixValue.resize(row, col);

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

	math::matrix<double> M1_for_Omega(numVectors, col);	//declaration of the two matrix M1 and M2
	math::matrix<double> M2_for_WSupp(numVectors, col);
	std::vector<double> rVariable(dimension);
	double res_Omega = 0.0, res_WSupp = 0.0;
	unsigned int row_id=0,col_id=0;
	for(unsigned int eachDir= 0; eachDir < AllDirections.size(); eachDir++){

		rVariable = AllDirections[eachDir].v;		//get the direction
		row_id = AllDirections[eachDir].R;
		col_id = AllDirections[eachDir].C;
		res_Omega = Omega_Support(ReachParameters, rVariable, Initial,
				SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
		res_WSupp =


		M1_for_Omega[row_id][col_id] = res_Omega;
		M2_for_WSupp[row_id][col_id] = res_WSupp;




	}

	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {
		double zIInitial = 0.0, zI = 0.0, zV = 0.0;
		double sVariable, s1Variable;
		//	std::cout<<"\nCheck 1"<<endl;
		std::vector<double> r1Variable;	//now single dimension
		r1Variable.resize(dimension);
		std::vector<double> rVariable;
		rVariable.resize(dimension);
		for (int i = 0; i < dimension; i++) {
			rVariable[i] = ReachParameters.Directions(eachDirection, i);
		}
		unsigned int loopIteration = 0;
		sVariable = 0.0; 		//initialize s0
		zIInitial = Omega_Support(ReachParameters, rVariable, Initial,
				SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
		MatrixValue(eachDirection, loopIteration) = zIInitial;
		loopIteration++;
		for (; loopIteration < shm_NewTotalIteration;) { //Now stopping condition is only "shm_NewTotalIteration"
			double TempOmega;
			ReachParameters.phi_trans.mult_vector(rVariable, r1Variable);
			zV = W_Support(ReachParameters, SystemDynamics, rVariable,
					s_per_thread_U, Min_Or_Max);
			s1Variable = sVariable + zV;
			zI = Omega_Support(ReachParameters, r1Variable, Initial,
					SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
			TempOmega = zI + s1Variable; 		//Y1
			MatrixValue(eachDirection, loopIteration) = TempOmega; 		//Y1
			rVariable = CopyVector(r1Variable);		//source to destination
			sVariable = s1Variable;
			loopIteration++;	//for the next Omega-iteration or Time-bound
		}	//end of while for each vector
	}









	if (isInvariantExist == true) {		//if invariant exist. Computing
		math::matrix<double> inv_sfm;
		int num_inv = invariant->getColumnVector().size();//number of Invariant's constriants
		inv_sfm.resize(num_inv, shm_NewTotalIteration);
		for (int eachInvDirection = 0; eachInvDirection < num_inv;
				eachInvDirection++) {
			for (unsigned int i = 0; i < shm_NewTotalIteration; i++) {
				inv_sfm(eachInvDirection, i) =
						invariant->getColumnVector()[eachInvDirection];
			}
		}
		return template_polyhedra(MatrixValue, inv_sfm,
				ReachParameters.Directions, invariant->getCoeffMatrix());
	} else {
		return template_polyhedra(MatrixValue, ReachParameters.Directions);
	}
}

*/
