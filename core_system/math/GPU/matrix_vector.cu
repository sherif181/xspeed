#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>
#include <cstdlib>
#include <cuda_runtime.h>
#include <ctime>
#include <cublas_v2.h>
#include <curand.h>
#include "core_system/math/matrix.h"

//void gpu_blas_mmul(cublasHandle_t handle, const float *A, const float *B, float *C, const int m,
//		const int k, const int n) {

void gpu_blas_mmul(cublasHandle_t handle, const float *A, const float *B, float *C, const int m,
		const int k, const int n) {
	int lda = m, ldb = k, ldc = m;
	const float alf = 1;
	const float bet = 0;
	const float *alpha = &alf;
	const float *beta = &bet;

// Create a handle for CUBLAS
	/*cublasHandle_t handle;
	cublasCreate(&handle);*/
// Do the actual multiplication
	cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, A, lda, B,
			ldb, beta, C, ldc);
//	std::cout<<"Yes called";
// Destroy the handle
	//cublasDestroy(handle);
}

void GPU_Multiply_Matrix_Vector(cublasHandle_t handle,
		thrust::device_vector<float> &A, int rowA, int colA,
		thrust::device_vector<float> &B, int rowB, int colB, thrust::device_vector<float> &C) {

/*thrust::device_vector<float> GPU_Multiply_Matrix_Vector(thrust::device_vector<float> A, int rowA, int colA,
		thrust::device_vector<float> B, int rowB, int colB) {*/
	//*A and *B -- are Matrices in column-major format
	//*A will be a Matrix and *B a vector for Matrix to vector multiplication
	//rowA,colA,rowB,colB -- are the transformed values of original Matrices

	gpu_blas_mmul(handle, thrust::raw_pointer_cast(&A[0]),
			thrust::raw_pointer_cast(&B[0]), thrust::raw_pointer_cast(&C[0]),
			colA, colB, rowB);
//	return C;	//rows=rowB and col=colA all outputs are in column-major format including row,col
}

void convert_thrust_matrix(math::matrix<double> A, thrust::device_vector<float> &d_A) {
//converting matrix A from row-major to thrust::device_vector in column-major format
		for (int i = 0; i < A.size1(); i++) {
		for(int j=0;j<A.size2();j++){
			d_A[j * A.size2() + i] = (float) A(i,j);
		}
	}
	//return d_A;
}
