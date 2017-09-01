/*
 * matrix_vector.cuh
 *
 *  Created on: 13-Apr-2015
 *      Author: amit
 */

#ifndef MATRIX_VECTOR_CUH_
#define MATRIX_VECTOR_CUH_

#include "core_system/math/matrix.h"
//#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>
#include <cstdlib>
#include <cuda_runtime.h>
#include <ctime>
#include <cublas_v2.h>
#include <curand.h>

/*void gpu_blas_mmul(const float *A, const float *B, float *C, const int m,
		const int k, const int n);*/

void gpu_blas_mmul(cublasHandle_t handle, const float *A, const float *B, float *C, const int m,
		const int k, const int n);

/*thrust::device_vector<float> GPU_Multiply_Matrix_Vector(
		thrust::device_vector<float> A, int rowA, int colA,
		thrust::device_vector<float> B, int rowB, int colB);*/

void GPU_Multiply_Matrix_Vector(cublasHandle_t handle,
		thrust::device_vector<float> &A, int rowA, int colA,
		thrust::device_vector<float> &B, int rowB, int colB, thrust::device_vector<float> &C);

//thrust::device_vector<float> convert_thrust_matrix(math::matrix<double> A);
void convert_thrust_matrix(math::matrix<double> A, thrust::device_vector<float> &d_A);

/*
int main() {
	// Allocate 3 arrays on CPU
	int nr_rows_A, nr_cols_A, nr_rows_B, nr_cols_B, nr_rows_C, nr_cols_C;
	// for simplicity we are going to use square arrays
	nr_rows_A = nr_cols_A = nr_rows_B = nr_cols_B = nr_rows_C = nr_cols_C = 3;
	thrust::device_vector<float> d_A(nr_rows_A * nr_cols_A), d_B(
			nr_rows_B * nr_cols_B), d_C(nr_rows_C * nr_cols_C);
	thrust::host_vector<float> A(nr_rows_A * nr_cols_A), B(
			nr_rows_A * nr_cols_A), C(nr_rows_A * nr_cols_A);
	// Fill the arrays A and B on GPU with random numbers
//	GPU_fill_rand(thrust::raw_pointer_cast(&d_A[0]), nr_rows_A, nr_cols_A);
//	GPU_fill_rand(thrust::raw_pointer_cast(&d_B[0]), nr_rows_B, nr_cols_B);
	A = d_A;
	B = d_B;
	// Multiply A and B on GPU
	gpu_blas_mmul(thrust::raw_pointer_cast(&d_A[0]),
			thrust::raw_pointer_cast(&d_B[0]),
			thrust::raw_pointer_cast(&d_C[0]), nr_rows_A, nr_cols_A, nr_cols_B);
	C = d_C;

	return 0;
}
*/

#endif /* MATRIX_VECTOR_CUH_ */
