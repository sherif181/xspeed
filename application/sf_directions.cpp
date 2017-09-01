/*
 * sf_directions.cpp
 *
 *  Created on: 14-Dec-2014
 *      Author: amit
 */

#include "application/sf_directions.h"

std::vector<std::vector<double> > generate_axis_directions(unsigned int N) {
	std::vector<std::vector<double> > directions;
	for (unsigned int i = 0; i < N; i++) {
		std::vector<double> d1(N, 0.0);
		d1[i] = 1;
		directions.push_back(d1);
		std::vector<double> d2(N, 0.0);
		d2[i] = -1;
		directions.push_back(d2);
	}
	return directions;
}

std::vector<std::vector<double> > get_octagonal_directions(unsigned int dim) {
	std::vector<std::vector<double> > mydirs;

	for (unsigned int i = 0; i < dim; i++) {
		for (unsigned int j = i; j < dim; j++) {
			if (i == j) {
				std::vector<double> v1(dim, 0);
				v1[i] = 1;
				mydirs.push_back(v1);
				std::vector<double> v2(dim, 0);
				v2[i] = -1;
				mydirs.push_back(v2);
			} else {

				std::vector<double> v1(dim, 0);
				v1[i] = 1;
				v1[j] = 1;
				mydirs.push_back(v1);

				std::vector<double> v2(dim, 0);
				v2[i] = 1;
				v2[j] = -1;
				mydirs.push_back(v2);

				std::vector<double> v3(dim, 0);
				v3[i] = -1;
				v3[j] = 1;
				mydirs.push_back(v3);

				std::vector<double> v4(dim, 0);
				v4[i] = -1;
				v4[j] = -1;
				mydirs.push_back(v4);
			}

		}
	}
	return mydirs;
}

/*
 * Directly generating separate Direction List for polytope -- X0 and U
 * Optimising the Computation within the Loop to avoid duplicate directions at the same time
 *
 * Now making the index suitable for parallelizing using #pragma OMP parallel for each numVectors
 */
void getDirectionList_X0_and_U(int numCoresAvail, ReachabilityParameters &ReachParameters,
		unsigned int newiters, math::matrix<float> &list_X0, math::matrix<float> &list_U, bool U_empty, Dynamics& SystemDynamics) {

	int numVectors = ReachParameters.Directions.size1();
	int dimension = ReachParameters.Directions.size2();

	math::matrix<double> B_trans, phi_tau_Transpose;
	if (!SystemDynamics.isEmptyMatrixB) //if Matrix B is not Empty
		B_trans = ReachParameters.B_trans;
	if (!SystemDynamics.isEmptyMatrixA) //if Matrix A is not Empty
		phi_tau_Transpose = ReachParameters.phi_trans;

	unsigned int total_list_X0 = list_X0.size1(); //total number of directions for X0 is [ numDirs * (iters + 1) ]
	unsigned int total_list_U = list_U.size1(); //total number of directions for U is [ numDirs * iters ]

	int cores;
	if (numVectors >= numCoresAvail)
		cores = numVectors;
	else
		cores = numCoresAvail;

//	omp_set_dynamic(0);	//handles dynamic adjustment of the number of threads within a team
//#pragma omp parallel for num_threads(cores)
	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {
		unsigned int index_X, indexU; //making the index suitable for parallelizing
		if (!U_empty) {
			indexU = eachDirection * newiters;
		}
		if (eachDirection == 0) { //only starting loop begins with 0
			index_X = eachDirection * newiters;
		} else {
			index_X = eachDirection * newiters + eachDirection; //only X0(list_X0) has 2 directions for first-iter
		}
		std::vector<double> rVariable(dimension), r1Variable(dimension); //now single dimension
		for (int i = 0; i < dimension; i++) {
			rVariable[i] = ReachParameters.Directions(eachDirection, i);
			list_X0(index_X, i) = (float) rVariable[i];
		}
		index_X++; //1st
		unsigned int loopIteration = 0;
		// ********** Omega's Directions  **********************
		std::vector<double> phi_trans_dir, B_trans_dir;
		if (!SystemDynamics.isEmptyMatrixA){ //if Matrix A is not Empty
			phi_tau_Transpose.mult_vector(rVariable, phi_trans_dir);
		}

		if (!SystemDynamics.isEmptyMatrixB) { //if not only than will be required to multiply
			B_trans.mult_vector(rVariable, B_trans_dir);
		}
		//std::cout << index_X0 << " ";
		for (unsigned int x = 0; x < list_X0.size2(); x++) { //dimension of direction
			if (!SystemDynamics.isEmptyMatrixA){ //if Matrix A is not Empty
				list_X0(index_X, x) = (float) phi_trans_dir[x]; //X0 and U both has the same dimension
			}else {
				list_X0(index_X, x) = (float) rVariable[x]; //Since A is empty :: {tau.A}' reduces to zero so, e^{tau.A}' reduces to 1
				// so, 1 * rVariable give only rVariable
			}//handling Constant Dynamics
			//list_X0(index_X, x) = (float) phi_trans_dir[x]; //X0 and U both has the same dimension

			if (!SystemDynamics.isEmptyMatrixB && !SystemDynamics.U->getIsEmpty()) { //if not only than will be required to multiply
				list_U(indexU, x) = (float) B_trans_dir[x]; //optimizing in a single loop
			}
		}
		index_X++; //2nd
		if (!U_empty) { //if not only than will be required to multiply
			indexU++; //for next entry
		}
		// ********** Omega's Directions End **********************
		loopIteration++;
		for (; loopIteration < newiters;) { //Now stopping condition is only "shm_NewTotalIteration"
			//		ReachParameters.phi_trans.mult_vector(rVariable, r1Variable);	//replacement from previous step
			//r1Variable = phi_trans_dir; //direct replacement from previous computation
			if (!SystemDynamics.isEmptyMatrixA){ //if Matrix A is not Empty
				r1Variable = phi_trans_dir; //direct replacement from previous computation
			} else { //if Matrix A is Empty in case of constant dynamics
				r1Variable = rVariable;
			}//handling Constant Dynamics

			// ********** W_Support's Directions  **********************
			//		std::vector<double> Btrans_dir;
			//		B_trans.mult_vector(rVariable, Btrans_dir);	//replacement from previous step
			// ********** W_Support's Directions End **********************
			// ********** Omega's Directions  **********************
			std::vector<double> B_trans_dir1;
			if (!SystemDynamics.isEmptyMatrixA) //if Matrix A is not Empty
				phi_tau_Transpose.mult_vector(r1Variable, phi_trans_dir);

			if (!U_empty) { //if not only than will be required to multiply
				B_trans.mult_vector(r1Variable, B_trans_dir1);
			}
			for (unsigned int x = 0; x < list_X0.size2(); x++) { //dimension of direction
				//list_X0(index_X, x) = (float) phi_trans_dir[x]; //X0 and U both has the same dimension
				if (!SystemDynamics.isEmptyMatrixA){ //if Matrix A is not Empty
					list_X0(index_X, x) = (float) phi_trans_dir[x]; //X0 and U both has the same dimension
				}else {
					list_X0(index_X, x) = (float) r1Variable[x]; //Since A is empty :: {tau.A}' reduces to zero so, e^{tau.A}' reduces to 1
					// so, 1 * r1Variable give only r1Variable
				}//handling Constant Dynamics

				if (!U_empty) { //if not only than will be required to multiply
					list_U(indexU, x) = (float) B_trans_dir1[x]; //optimizing in a single loop
				}
			}
			index_X++; //for next entry
			if (!U_empty) { //if not only than will be required to multiply
				indexU++; //for next entry
			}
			// ********** Omega's Directions End **********************
			rVariable = CopyVector(r1Variable); //source to destination
			loopIteration++; //for the next iteration
		} //end of iteration for each direction
		  //std::cout<<"Over eachDirection = "<< eachDirection<<std::endl;
	} //end of all directions
	  //std::cout<<"what is the Problem?"<<std::endl;
	  //std::cout<<std::endl;
	  //return Direction_List;
}

//Only for profiling GLPK solver Time for comparison with boundary value implementation
void getDirectionList_X0_and_U_OnlyForGLPK(
		ReachabilityParameters &ReachParameters, unsigned int newiters,
		std::list<std::vector<double> > &list_X0,
		std::list<std::vector<double> > &list_U, bool U_empty,
		Dynamics& SystemDynamics) {
	int numVectors = ReachParameters.Directions.size1();
	int dimension = ReachParameters.Directions.size2();
//std::cout<<"\nCalling the GLPK directions Vector List\n";

	math::matrix<double> B_trans, phi_tau_Transpose;
	if (!SystemDynamics.isEmptyMatrixB) //if Matrix B is not Empty
		B_trans = ReachParameters.B_trans;
	if (!SystemDynamics.isEmptyMatrixA) //if Matrix A is not Empty
		phi_tau_Transpose = ReachParameters.phi_trans;

//	unsigned int total_list_X0 = list_X0.size(); //total number of directions for X0 is [ numDirs * (iters+1) ]
//	unsigned int total_list_U = list_U.size(); //total number of directions for U is [ numDirs * iters ]
	unsigned int total_list_X0 = numVectors * (newiters + 1); //1 extra for loop1
	unsigned int total_list_U = numVectors * newiters; //'n' dirs for each 'n' loops

	std::list<std::vector<double> >::iterator index_X;
	std::list<std::vector<double> >::iterator indexU;
	index_X = list_X0.begin(); //starts with 0
	indexU = list_U.begin(); //starts with 0
	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {
		std::vector<double> rVariable(dimension), r1Variable(dimension); //now single dimension
		for (int i = 0; i < dimension; i++) {
			rVariable[i] = ReachParameters.Directions(eachDirection, i);
			//list_X0(index_X, i) = (float) rVariable[i];
		}
		list_X0.insert(index_X, rVariable);
		unsigned int loopIteration = 0;
		// ********** Omega's Directions  **********************
		std::vector<double> phi_trans_dir, B_trans_dir;
		if (!SystemDynamics.isEmptyMatrixA) //if Matrix A is not Empty
			phi_tau_Transpose.mult_vector(rVariable, phi_trans_dir);
		if (!U_empty) { //if not only than will be required to multiply
			B_trans.mult_vector(rVariable, B_trans_dir);
		}
		list_X0.insert(index_X, phi_trans_dir); //X0 and U both has the same dimension
		if (!U_empty) { //if not only than will be required to multiply
			list_U.insert(indexU, B_trans_dir); //optimizing in a single loop
		}
		// ********** Omega's Directions End **********************
		loopIteration++;
		for (; loopIteration < newiters;) { //Now stopping condition is only "shm_NewTotalIteration"
			//		ReachParameters.phi_trans.mult_vector(rVariable, r1Variable);	//replacement from previous step
			r1Variable = phi_trans_dir; //direct replacement from previous computation
			// ********** W_Support's Directions  **********************
			//		std::vector<double> Btrans_dir;
			//		B_trans.mult_vector(rVariable, Btrans_dir);	//replacement from previous step
			// ********** W_Support's Directions End **********************
			// ********** Omega's Directions  **********************
			std::vector<double> B_trans_dir1;
			if (!SystemDynamics.isEmptyMatrixA) //if Matrix A is not Empty
				phi_tau_Transpose.mult_vector(r1Variable, phi_trans_dir);
			if (!U_empty) { //if not only than will be required to multiply
				B_trans.mult_vector(r1Variable, B_trans_dir1);
			}

			list_X0.insert(index_X, phi_trans_dir); //X0 and U both has the same dimension
			if (!U_empty) { //if not only than will be required to multiply
				list_U.insert(indexU, B_trans_dir1); //optimizing in a single loop
			}
			// ********** Omega's Directions End **********************
			rVariable = CopyVector(r1Variable); //source to destination
			loopIteration++; //for the next iteration
		} //end of iteration for each direction
	} //end of all directions
}
