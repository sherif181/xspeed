/*
 * computeInitialPolytope.cpp
 *
 *  Created on: 24-Nov-2014
 *      Author: amit
 */

#include "Utilities/computeInitialPolytope.h"

/*

polytope::ptr create_polytope_set(supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, Dynamics& SystemDynamics) {
	//polytope p;
	std::vector<double> columnvector(ReachParameters.Directions.size1());
	glpk_lp_solver lp, lp_U;
	lp.setMin_Or_Max(2);
	lp_U.setMin_Or_Max(2);

	int dim = ReachParameters.Directions.size2();
	std::vector<double> direction(dim);
	std::cout << "Dim =  " << dim;
	std::cout << "Direction = " << ReachParameters.Directions.size1()
			<< std::endl;
	for (unsigned int i = 0; i < ReachParameters.Directions.size1(); i++) {
		std::cout << "loop 1 \n";
		for (int j = 0; j < dim; j++) {
			direction[j] = ReachParameters.Directions(i, j);
			std::cout << "loop 2 \n";
		}
		if (!Initial->getIsEmpty()) //set glpk constraints If not an empty polytope
			lp.setConstraints(ReachParameters.X0->getCoeffMatrix(),
					ReachParameters.X0->getColumnVector(),
					ReachParameters.X0->getInEqualitySign());

		if (!SystemDynamics.U->getIsEmpty()) {	//empty polytope
			lp_U.setConstraints(SystemDynamics.U->getCoeffMatrix(),
					SystemDynamics.U->getColumnVector(),
					SystemDynamics.U->getInEqualitySign());

			columnvector[i] = Initial->computeSupportFunction(direction, lp,
					lp_U, 2);
		}
	}
	/ *
	 p.setCoeffMatrix(ReachParameters.Directions);
	 p.setColumnVector(columnvector);
	 p.setInEqualitySign(1);
	 * /
	polytope::ptr p = polytope::ptr(
			new polytope(ReachParameters.Directions, columnvector, 1));

	return p;
}


template_polyhedra Reach_Parallel_Iteration(Dynamics& SystemDynamics,
		supportFunctionProvider::ptr initial,
		ReachabilityParameters& reach_parameters, polytope::ptr invariant,
		bool isInvariantExist) {

	template_polyhedra reach_region, Reach;
	int NCores = 4;
	std::list<polytope::ptr> polys;

	polys = compute_initial_polytopes_transMink(SystemDynamics, initial,
			reach_parameters, invariant, isInvariantExist, NCores);
	std::list<polytope::ptr>::iterator it;
	//double START_TIME = ReachParameters.time_step * col;//same as newTimeBound
	unsigned int newTimeBound = reach_parameters.TimeBound / NCores;
	reach_parameters.TimeBound = newTimeBound;

	int totIters;
	if (isInvariantExist == true) {	//if invariant exist. Computing
		totIters = InvariantBoundaryCheck(SystemDynamics, initial,
				reach_parameters, invariant);
	}	//End of Invariant Directions

	//int partition_iterations = reach_parameters.Iterations / NCores;
	int partition_iterations = totIters / NCores;
	reach_parameters.Iterations = partition_iterations;	//required in Invariant_Boundary_Check
	//cout << "partition_iterations" << partition_iterations << "\n";

	int count = 1;
	for (it = polys.begin(); it != polys.end(); it++) {
		std::cout << "Number of Initial Polytopes = " << count << "\n";
		if (count == 2)
			GeneratePolytopePlotter((*it));
		reach_parameters.X0 = (*it);		//changing the initial polytope
		reach_region = reachabilityParallel(SystemDynamics, (*it),
				reach_parameters, invariant, isInvariantExist);
		reach_parameters.TimeBound += newTimeBound;
		if (count == 1)
			Reach = reach_region;
		else
			Reach = Reach.union_TemplatePolytope(reach_region);
		count++;
	}

	return Reach;
}


/ *
 * Trying Idea :: Method(i) described in my note book
 * /
std::list<polytope::ptr> compute_initial_polytopes_transMink(
		Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int NCores) {

	std::list<polytope::ptr> initial_polytopes;
	std::vector<D> AllDirections;
	AllDirections = ReachParameters.AdvanceTransposeDirection;
	polytope::ptr original_initial_poly;
	original_initial_poly = ReachParameters.X0;
	int numVectors = ReachParameters.Directions.size1();
	int dimension = Initial->getSystemDimension();
	unsigned int shm_NewTotalIteration = ReachParameters.Iterations; //Shared Variable for resize iterations number on crossing with invariant
	unsigned int Actual_Iteration_size = ReachParameters.Iterations; //needed to know the direction_no in Transposed_Direction_Matrix
	int Min_Or_Max = 2;

	//std::cout << "\nNumber of Iterations to be executed = " << shm_NewTotalIteration << std::endl;
	int PartSize = shm_NewTotalIteration / NCores;//partition size to be handled in each cores
	std::cout << "\nPartSize  = " << PartSize << std::endl;
	// Becarefull about the fractional partsize

	int iterations_cross_boundary=0;

	for (int part_no = 0; part_no < NCores; part_no++) {
		std::vector<double> columnVector(ReachParameters.Directions.size1());
		int col = part_no * PartSize;
		ReachParameters.TimeBound = ReachParameters.time_step * col;	//timebound changes for each CORES

		for (int eachDirection = 0; eachDirection < numVectors;
				eachDirection++) {
		/ *	if (part_no > 0)
				iterations_cross_boundary = iterations_cross_boundary -1;
			int position = eachDirection * Actual_Iteration_size + iterations_cross_boundary;
		* /
			int position = eachDirection * Actual_Iteration_size + col;

			// IF I CAN KNOW WHEN THE PREVIOUS ITERATION STOPS
			// OR WHEN THE PREVIOUS DYNAMICS CROSSED THE INVARIANT_BOUNDARY THEN THAT WILL BE THE POSITION OF THE TRANSPOSED DIRECTION
			// ie if i can knonw the col in the above case
			std::vector<double> transdirection(dimension);
			//std::cout<<"(";
//std::cout<<"Test 1\n";
			for (unsigned int i = 0; i < dimension; i++) {
				transdirection[i] = AllDirections[position].v[i];
				//std::cout<<transdirection[i]<<",";
			}
			//std::cout<<")";
//std::cout<<"Test 2\n";
			glpk_lp_solver s_per_thread_I, s_per_thread_U, s_per_thread_inv;
			s_per_thread_I.setMin_Or_Max(2);
			//if (!Initial->getIsEmpty())	//set glpk constraints If not an empty polytope
			if (!ReachParameters.X0->getIsEmpty())
				s_per_thread_I.setConstraints(
						ReachParameters.X0->getCoeffMatrix(),
						ReachParameters.X0->getColumnVector(),
						ReachParameters.X0->getInEqualitySign());

			s_per_thread_U.setMin_Or_Max(2);
			if (!SystemDynamics.U->getIsEmpty()) {	//empty polytope
				s_per_thread_U.setConstraints(
						SystemDynamics.U->getCoeffMatrix(),
						SystemDynamics.U->getColumnVector(),
						SystemDynamics.U->getInEqualitySign());
			}
			double START_TIME = ReachParameters.time_step * col;//same as newTimeBound
			math::matrix<double> phi, phi_trans;
			SystemDynamics.MatrixA.matrix_exponentiation(phi, START_TIME);
			phi.transpose(phi_trans);
			/ *
			supportFunctionProvider::ptr Initial = transMinkPoly::ptr(
					new transMinkPoly(ReachParameters.X0, SystemDynamics.U,
							phi_trans, ReachParameters.B_trans, START_TIME));
			* /

			columnVector[eachDirection] = Initial->computeSupportFunction(
					transdirection, s_per_thread_I, s_per_thread_U, 2);
			//	std::cout<<"columnVector[eachDirection] = "<<columnVector[eachDirection] <<"\t";
		}
//		std::cout << "\n";
		polytope::ptr p = polytope::ptr(
				new polytope(ReachParameters.Directions, columnVector, 1));
		initial_polytopes.push_back(p);
//		std::cout<<"Test 3\n";
		/ *
		 * Before computing the next polytope find out where the previous polytope "p" will meet the invariant_boundary
		 * /
/ *

		if (isInvariantExist == true) {	//if invariant exist. Computing
			ReachParameters.X0 = p;
			col = (part_no+1) * PartSize;	//for next time bound
			//ReachParameters.TimeBound = ReachParameters.time_step * col;	//timebound changes for each CORES
			//ReachParameters.Iterations = Actual_Iteration_size / NCores;
			iterations_cross_boundary = InvariantBoundaryCheck(SystemDynamics, p,
						ReachParameters, invariant);
			std::cout<<"Stoped after iterations = " <<iterations_cross_boundary <<std::endl;
			}	//End of Invariant Directions
		ReachParameters.X0 = original_initial_poly;	//restoring the original initial_polytope
* /

	}
	return initial_polytopes;
}


std::list<polytope::ptr> compute_initial_polytopes(Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int NCores) {
	std::list<polytope::ptr> initial_polytopes;
	std::vector<D> AllDirections;
	AllDirections = ReachParameters.AdvanceTransposeDirection;

	polytope::ptr initial_poly;

	int numVectors = ReachParameters.Directions.size1();

	int dimension = Initial->getSystemDimension();
	unsigned int shm_NewTotalIteration = ReachParameters.Iterations; //Shared Variable for resize iterations number on crossing with invariant
	unsigned int Actual_Iteration_size = ReachParameters.Iterations; //needed to know the direction_no in Transposed_Direction_Matrix
	int Min_Or_Max = 2;

	if (isInvariantExist == true) {	//if invariant exist. Computing
		shm_NewTotalIteration = InvariantBoundaryCheck(SystemDynamics, Initial,
				ReachParameters, invariant);
	}	//End of Invariant Directions

	int flag = 0;
	//std::cout << "\nNumber of Iterations to be executed = " << shm_NewTotalIteration << std::endl;

	int PartSize = shm_NewTotalIteration / NCores;//partition size to be handled in each cores
	//std::cout << "\nPartSize  = " << PartSize << std::endl;
	//Becarefull about the fractional partsize
	for (int part_no = 0; part_no < NCores; part_no++) {
		std::vector<double> columnVector(ReachParameters.Directions.size1());

		for (int eachDirection = 0; eachDirection < numVectors;
				eachDirection++) {
			glpk_lp_solver s_per_thread_I, s_per_thread_U, s_per_thread_inv;

			s_per_thread_I.setMin_Or_Max(2);
			if (!Initial->getIsEmpty())	//set glpk constraints If not an empty polytope
				s_per_thread_I.setConstraints(
						ReachParameters.X0->getCoeffMatrix(),
						ReachParameters.X0->getColumnVector(),
						ReachParameters.X0->getInEqualitySign());

			s_per_thread_U.setMin_Or_Max(2);
			if (SystemDynamics.U->getIsEmpty()) {	//empty polytope
				//Polytope is empty so no glpk object constraints to be set
			} else {
				s_per_thread_U.setConstraints(
						SystemDynamics.U->getCoeffMatrix(),
						SystemDynamics.U->getColumnVector(),
						SystemDynamics.U->getInEqualitySign());
			}
			double zIInitial = 0.0, zI = 0.0, zV = 0.0;
			double sVariable, s1Variable;

			std::vector<double> r1Variable;	//now single dimension
			r1Variable.resize(dimension);
			std::vector<double> rVariable;
			rVariable.resize(dimension);
			//	std::cout << "\ntest 1 "<< std::endl;
			int li = part_no * PartSize;//computes the Column no of the Transposed_Direction_Matrix
			//	std::cout << "\n li = "<<li<< std::endl;
			int position = eachDirection * Actual_Iteration_size + li;//	Row Major representation of Transposed_Direction_Matrix
			//	std::cout << "\nposition = "<<position<< std::endl;
			//int position = eachDirection * shm_NewTotalIteration + 0;
			/ *	if (AllDirections[position].R == eachDirection
			 && AllDirections[position].C == li)
			 std::cout << "\nPosition of the Data Structure of All Directions Correctly mapped!!!\n";
			 else
			 std::cout
			 << "\nPosition of the Data Structure of All Directions NOT mapped Correctly ??????XXXXX\n";
			 * /
			 for (int i = 0; i < dimension; i++) {
				rVariable[i] = AllDirections[position].v[i];
			}

			unsigned int loopIteration = 0;
			sVariable = 0.0; 		//initialize s0

			zIInitial = Omega_Support(ReachParameters, rVariable, Initial,
					SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);

			loopIteration++;

			double TempOmega;
			for (; loopIteration < 2;) { //for Second Omega (or just the next Omega)

				double original_tau = ReachParameters.time_step;
				double new_delta_tau = ReachParameters.time_step * li; //shifted dt
				ReachParameters.time_step = new_delta_tau;

				//int same_position = eachDirection * Actual_Iteration_size + li;
				for (int i = 0; i < dimension; i++) {
					//r1Variable[i] = AllDirections[same_position].v[i];		//transposed direction
					r1Variable[i] = AllDirections[position].v[i]; //transposed direction
				}
				zV = W_Support(ReachParameters, SystemDynamics, rVariable,
						s_per_thread_U, Min_Or_Max);

				s1Variable = sVariable + zV;
				zI = Omega_Support(ReachParameters, r1Variable, Initial,
						SystemDynamics, s_per_thread_I, s_per_thread_U,
						Min_Or_Max);
				TempOmega = zI + s1Variable; 		//Y1
				loopIteration++;	//for the next Omega-iteration or Time-bound

				ReachParameters.time_step = original_tau;//restore original for next direction

			}	//end of while for each vector
			//cout<<endl;
			if (part_no == 0) {			//First Omega
				columnVector[eachDirection] = zIInitial;
			} else {		//Second Omega onwards
				columnVector[eachDirection] = TempOmega;
			}	//end of pragma omp parallel for
		}
		initial_poly = polytope::ptr(
				new polytope(ReachParameters.Directions, columnVector, 1));
		initial_polytopes.push_back(initial_poly);
	}
	return initial_polytopes;
}

*/

