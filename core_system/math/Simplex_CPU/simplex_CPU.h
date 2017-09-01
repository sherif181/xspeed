/*
 * simplex.cuh
 *
 *  Created on: 11-Apr-2015
 *      Author: bilzcuda
 */

#ifndef SIMPLEX_CPUH_
#define SIMPLEX_CPUH_

#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <malloc.h>
#include <vector>
#include <list>
#include "core_system/math/matrix.h"
#include <climits>
#include <boost/shared_ptr.hpp>

//interface is   Simplex(A<B<C,N_S,No_O,No_C)

class Simplex_CPU {
private:
	unsigned int number_of_LPs;	//total number of LPs to be solved per instance
	math::matrix<double> orig_CoefficientMatrix;
	std::vector<double> BoundValue;
	//unsigned int number_constraints; //can be obtained from orig_CoefficientMatrix.size1();
	//math::matrix<float> C;
	unsigned int number_of_Constraint;
public:
	typedef boost::shared_ptr<Simplex_CPU> simplex_ptr;
	//math::matrix<float> A, C;
	//std::vector<float> B;
	int *Sel, *G_Sel;
	float a;
	int M, N, i, j;
	float *MAT, *G_MAT, *G_R, *R, *N_MAT;
	float *orig_MAT;

	int NB, f, c, No_c;
	Simplex_CPU();

	void setConstratint_CPU(math::matrix<double> &A1, std::vector<double> &B1);	//setting constraints of simplex
	float ComputeLP_CPU(std::vector<double> &C1);	//GPU computations
	~Simplex_CPU() {
	}

// Compute the LP. argument is Objective function(S)
};

#endif /* SIMPLEX_CPUH_ */
