/*
 * simplex.cuh
 *
 *  Created on: 11-Apr-2015
 */

#ifndef SIMPLEX_CUH_
#define SIMPLEX_CUH_
#include <iostream>
#include <stdio.h>
#include <cstdlib>
//#include <cuda_runtime.h>
#include <ctime>
#include <malloc.h>
#include <vector>
#include <list>
#include "core_system/math/matrix.h"
#include "application/DataStructureDirections.h"
#include <climits>

//interface is   Simplex(A<B<C,N_S,No_O,No_C)

struct block_lp {
	math::matrix<float> block_obj_coeff;
};
struct block_lp_result {
	std::vector<float> results;
	std::vector<float> results_UnitBall;
};

class Simplex {
private:
	unsigned int number_of_LPs; //total number of LPs to be solved per instance
	unsigned int number_of_LPs_for_X; //total number of LPs to be solved for X
	unsigned int number_of_LPs_for_U; //total number of LPs to be solved for U

	math::matrix<double> orig_CoefficientMatrix;
	std::vector<double> BoundValue;

	std::vector<double> BoundValue_for_X;
	std::vector<double> BoundValue_for_U;
	//unsigned int number_constraints; //can be obtained from orig_CoefficientMatrix.size1();
	math::matrix<float> C;
	unsigned int number_of_Constraint;
	bool single_bound_flag_result;
	float *R_X, *R_U, *R_UnitBall, *R_dotProduct; //return value from Device/Kernel
	float *G_R_X, *G_R_U, *G_R_UnitBall, *G_R_dotProduct;

	float *bound_result_for_X, *d_bound_result_for_X, *d_obj_coeff_for_X,
			*obj_coeff_for_X,*U_C, *d_U_C;

public:
	float *bound_result, *d_bound_result, *d_obj_coeff, *obj_coeff;
	//math::matrix<float> A, C;
	//std::vector<float> B;
	int *Sel, *G_Sel;
	float a;
	int M, N, i, j;
	float *MAT, *G_MAT, *G_R, *R, *N_MAT;
	int NB, f, c, No_c;

	__host__ Simplex(unsigned int N_S); //Old implementation

//get status of particular simplex
	__host__ int getStatus(int n); //get status of particular simplex

//get the No of simplex the object is ruuning on GPU
	__host__ int getNo_OF_Simplx(); //get the No of simplex the object is ruuning on GPU

//get the result of all simplex
	__host__ std::vector<float> getResultAll(); //OLD implementation

//get the result of all simplex

	__host__ float getResult(int n); // get result of particular simplex

	__host__ std::vector<int> getStatusAll(); //get the status of all simplex

	__host__ void setConstratint(math::matrix<double> A1,
			std::vector<double> B1); //OLD implementation setting constraints of simplex

	__host__ void ComputeLP(math::matrix<float> &C1, unsigned int streams); //Old implementation

	// ************ All New Implementation ********
	__host__ Simplex(int UnitBall ,unsigned int N_S_for_X); //New implementation

	__host__ void setConstratint(math::matrix<double> coeffMatrix_for_X,
			std::vector<double> columVector_for_X,
			int UnitBall);//New implementation

	__host__ void getResult_X(std::vector<float> &result); //New implementation
	__host__ void getResult_U(std::vector<float> &result); //New implementation
	__host__ void getResult_UnitBall(std::vector<float> &result); //New implementation
	__host__ void getResult_dotProduct(std::vector<float> &result); //New implementation


	//here obj_funs_for_X/U/UnitBall is the List of all the directions on which support function is to be computed
	__host__ void ComputeLP(math::matrix<float> &obj_funs_for_X,
			int UnitBall, unsigned int streams); //New implementation

	//here obj_funs_for_X/dotProduct(U.C * directions)/UnitBall is the List of all the directions on which support function is to be computed
	__host__ void ComputeLP(math::matrix<float> &obj_funs_for_X,
			int UnitBall, unsigned int streams, std::vector<double> U_C); //New implementation


	~Simplex() {
	}

// Compute the LP. argument is Objective function(S)
};

#endif /* SIMPLEX_CUH_ */
