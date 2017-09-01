#include "simplex.cuh"
#include <omp.h>
#include "iostream"
//Trying Pivot row as "First positive coefficient Rule"
//Tested Working for Both Helicopter and Five Dimensional system
//Streaming working for General GPU Simplex but only for Helicopter Model but not for Five Dimensional System
__global__ void bound_simplex(float *bound_value, float *C, float *result,
		int dimension, unsigned int N_S, int offset) {

	/*
	 //Each thread solves one LP
	 unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	 if (i < N_S) {
	 float sum = 0.0;
	 for (int j = 0; j < dimension; j++) { //every column is obj_coeff
	 float obj_value = C[i * dimension + j];
	 if (obj_value < 0) {
	 sum = sum + -1 * obj_value * bound_value[j * 2 + 1];
	 } else {
	 sum = sum + obj_value * bound_value[j * 2];
	 }
	 }
	 result[i] = sum;
	 }
	 */

	//Each block solves one LP with only one thread working in each block
	int index = blockIdx.x + offset; //LP to be solved

	if (index < N_S) {

		//	int i = index * blockDim.x; //address of each LP(object_function) is (index x size of one LP)
		if (threadIdx.x == 0) {
			float sum = 0.0;
			for (int j = 0; j < blockDim.x; j++) { //every column is obj_coeff
				float obj_value = C[index * blockDim.x + j];
				if (obj_value < 0) {
					sum = sum + -1 * obj_value * bound_value[j * 2 + 1];
				} else {
					sum = sum + obj_value * bound_value[j * 2];
				}
			}
			result[index] = sum;
		} //end of working thread=0
	}
}

//New Implementation
/*
 * N_S_for_X :: number of simplex for polytope X
 * Kernel Launch::
 * 	GridSize:: nos. of LPs to be solved: N_S_for_X
 * 	BlockSize/threads:: is the dimension of the System/Model/direction.size()
 */__global__ void bounded_Solver_for_X_UnitBall(float *bound_value_for_X,
		float *C, float *result_for_X, unsigned int N_S_for_X,
		float *result_for_UnitBall, int offset) {
	//Each block solves one LP with only one thread working in each block
	int index = blockIdx.x + offset; //LP to be solved
	if (threadIdx.x == 0) {
		if (index < N_S_for_X) { //if index is greater(for UnitBall case) then it will not be evaluated
			int i = index * blockDim.x; //address of each LP(object_function) is (index x size of one LP)
			float sum_for_X = 0.0;
			float sum_for_UnitBall = 0.0;
			float obj_value_for_X;
			for (int j = 0; j < blockDim.x; j++) { //every column is obj_coeff:: Here blockDim.x is the dimension
				obj_value_for_X = C[index * blockDim.x + j];

				//	printf("obj_value_for_X = %f  ",obj_value_for_X);
				if (obj_value_for_X < 0) {
					sum_for_X = sum_for_X
							+ -1 * obj_value_for_X
									* bound_value_for_X[j * 2 + 1];
					sum_for_UnitBall += (-1 * obj_value_for_X);	//making abs(UnitBall)
				} else {
					sum_for_X = sum_for_X
							+ obj_value_for_X * bound_value_for_X[j * 2];
					sum_for_UnitBall += obj_value_for_X;//making abs(UnitBall)
				}
			}
			result_for_X[index] = sum_for_X;
			result_for_UnitBall[index] = sum_for_UnitBall;
			//printf("result_for_X = %f  result_for_UnitBall = %f\n",result_for_X[index], sum_for_UnitBall);
		} //end of thread/block check for N_S_for_X
	} //end of threadIdx.x==0
}

// New Implementation
//Again modified to include the + C of the equation Ax + Bu + C
/*
 * N_S_for_X :: number of simplex for polytope X
 * Kernel Launch::
 * 	GridSize:: nos. of LPs to be solved: N_S_for_X
 * 	BlockSize/threads:: is the dimension of the System/Model/direction.size()
 *
 *
 */__global__ void bounded_Solver_for_X_UnitBall_dotProduct(
		float *bound_value_for_X, float *C, float *U_C, float *result_for_X,
		unsigned int N_S_for_X, float *result_for_UnitBall,
		float *result_for_dotProduct, int offset) {
	//Each block solves one LP with only one thread working in each block
	int index = blockIdx.x + offset; //LP to be solved
	if (threadIdx.x == 0) {
		if (index < N_S_for_X) { //if index is greater(for UnitBall case) then it will not be evaluated
			int i = index * blockDim.x; //address of each LP(object_function) is (index x size of one LP)
			float sum_for_X = 0.0;
			float sum_for_UnitBall = 0.0;
			float prod_for_dotProd = 0.0;
			float obj_value_for_X;
			for (int j = 0; j < blockDim.x; j++) { //every column is obj_coeff:: Here blockDim.x is the dimension
				obj_value_for_X = C[index * blockDim.x + j];//retrieves each element of the directions/vector

				prod_for_dotProd = prod_for_dotProd + obj_value_for_X * U_C[j];

				//	printf("obj_value_for_X = %f  ",obj_value_for_X);
				if (obj_value_for_X < 0) {
					sum_for_X = sum_for_X
							+ -1 * obj_value_for_X
									* bound_value_for_X[j * 2 + 1];
					sum_for_UnitBall += (-1 * obj_value_for_X);	//making abs(UnitBall)
				} else {
					sum_for_X = sum_for_X
							+ obj_value_for_X * bound_value_for_X[j * 2];
					sum_for_UnitBall += obj_value_for_X;//making abs(UnitBall)
				}
			}
			result_for_X[index] = sum_for_X;
			result_for_UnitBall[index] = sum_for_UnitBall;
			result_for_dotProduct[index] = prod_for_dotProd;
			//printf("result_for_X = %f  result_for_UnitBall = %f\n",result_for_X[index], sum_for_UnitBall);
		} //end of thread/block check for N_S_for_X
	} //end of threadIdx.x==0
}

//__global__ void mykernel(float *S_MAT, int S_row, int S_col, float *Result, int S_N, int *S_Sel, float *R_data, int *R_index) {
__global__ void mykernel(float *S_MAT, int S_row, int S_col, float *Result,
		int S_N, float *R_data, int *R_index, int offset_res) {
//offset_res :: is the pointer to the LP number(start of the stream_block) to be solved
//int index = threadIdx.x + (blockIdx.x * blockDim.x);
//int index = blockIdx.x;
	int index = offset_res + blockIdx.x;
	if (index < (offset_res + S_N)) { //offset_size of lps plus block_size(ie., S_N)
//	if (index < S_N) {
		int tid;
		int i; // used for for index
		unsigned int temp_index;
		unsigned int temp_index1;
		int base = index * S_row * S_col;
		int R_base = index * blockDim.x; // blockDim.x = 96
		__shared__ int c;
		__shared__ int rm;
		__shared__ int row; //pivotRow
		__shared__ int pivotCol; //pivotCol this can remove global variable S_Sel

		int col = 1;
		__shared__ int remember[96]; //Found a column which is negative but theta/Min has no positive value
		__shared__ float col1[60]; //pivotColumn

		if (threadIdx.x == 0) {
			c = 0;
			rm = 0;
			row = -1; //pivotRow
			pivotCol = -1;
		}
		__syncthreads();
		while (!c) {
			__syncthreads();
			int minValue = 0;
			__shared__ int newpivotcol;
			int Last_row = S_row - 1; // Last row
			if (threadIdx.x == 0) {
				newpivotcol = -1;
				c = 1;
			}
			__syncthreads(); //make sure that c and newpivotcol is initialized
			//   ***************** Get_Pivot function begins  *****************
			//for (int j = 2; j < S_col - 1; j++) {//only last row but all column
			if (threadIdx.x >= 2 && threadIdx.x < (S_col - 1)) {
				int j = threadIdx.x;
				unsigned int temp_index1 = Last_row + j * S_row + base; //avoiding re-computation
				if (S_MAT[temp_index1] < minValue) {
					//minValue = S_MAT[temp_index1];
					newpivotcol = j; //"First positive coefficient rule"
					int local_NewPivotCol;
					local_NewPivotCol = *(volatile int*) &newpivotcol;
					//atomicCAS(&newpivotcol, local_NewPivotCol, j);
					/*
					 http://stackoverflow.com/questions/27616417/cuda-is-there-any-way-to-prevent-other-threads-from-changing-a-shared-or-global
					 if (atomicCAS(&newpivotcol, local_NewPivotCol, j)==local_NewPivotCol){
					 //this thread won the write
					 printf("Thread ID = %d ",threadIdx.x);
					 }*/
					//break;
				}
			}
			__syncthreads(); //here we have min and newpivotcol

			if (threadIdx.x == 0) {
				if (newpivotcol == -1) {
					//return -2;
					row = -2;
				} else {
					float row_min = INT_MAX;
					float row_num = -1;
					//TODO:: this temp_res can be an array of value computed in parallel
					//TODO:: row_min and row_num can then be computed using reduction
					for (i = 0; i < S_row - 1; i++) {
						temp_index = newpivotcol * S_row + i + base; //avoiding re-computation
						temp_index1 = i + (S_col - 1) * S_row + base; //avoiding re-computation
						if ((S_MAT[temp_index] > 0)
								&& (S_MAT[temp_index1] > 0)) {
							float temp_res = S_MAT[temp_index1]
									/ S_MAT[temp_index]; //avoiding re-computation
							if (temp_res <= row_min) {
								row_min = temp_res;
								row_num = i;
							}
						}
					}
					// *S_Sel = newpivotcol;
					pivotCol = newpivotcol;
					//S_Sel[index] = newpivotcol;
					if (row_min == INT_MAX) {
						//return -1;
						row = -1;
					}
					if (row_num != -1) {
						//return row_num;
						row = row_num;
					}
				}
			} //end of one thread
			__syncthreads();

			//   ***************** Get_Pivot function ends  *****************

			//col = S_Sel[index];
			//col = *S_Sel;
			col = pivotCol;
			if (row > -1) {
				tid = threadIdx.x;
				if (threadIdx.x >= 2 && threadIdx.x < S_col) {
					//for (int i1 = 2; i1 < S_col; i1++) {		//Data Parallel section 1
					if (tid == remember[tid - 2]) {
						temp_index = (S_row - 1) + (tid * S_row) + base; //avoiding re-computation
						S_MAT[temp_index] = -1 * S_MAT[temp_index]; //replacing back to original
					}
				} //Data Parallel section 1 done
				__syncthreads();
				tid = threadIdx.x;
				if (threadIdx.x >= 0 && threadIdx.x < S_row) {
					//for (int i = 0; i < S_row; i++) {	//Data Parallel section 2
					col1[tid] = S_MAT[(tid + col * S_row) + base]; //keeping the old pivotcol coeff
				} //Data Parallel section 2 done
				__syncthreads();

				unsigned int temp_row_base = row + base; //avoiding re-computation
				S_MAT[temp_row_base + S_row] =
						S_MAT[temp_row_base + col * S_row];
				//S_MAT[temp_row_base] = col - 1;
				S_MAT[row + base] = col - 1; //now temp_row_base is not required
				tid = threadIdx.x;
				if (threadIdx.x >= 2 && threadIdx.x < S_col) {
					//for (int j = 2; j < S_col; j++){		//Data Parallel section 3
					unsigned int row_base = row + base; //avoiding re-computation
					temp_index = row_base + (tid * S_row); //avoiding re-computation
					S_MAT[temp_index] = S_MAT[temp_index] / col1[row]; //S_MAT[row_base + S_row];
					//S_MAT[temp_index] = S_MAT[temp_index] / S_MAT[row_base + S_row];
				} //Data Parallel section 3 done
				__syncthreads();
				//printf("Row here = %d",row);
				tid = threadIdx.x;
				if (threadIdx.x >= 0 && threadIdx.x < S_row) {
					//for (int i = 0; i < S_row; i++) {	//Data parallel section 4
					for (i = 2; i < S_col; i++) {
						if (tid != row) {
							temp_index1 = i * S_row + base;
							temp_index = tid + temp_index1;
							S_MAT[temp_index] = S_MAT[temp_index]
									- (col1[tid] * S_MAT[row + temp_index1]);
						} else {
							break;
						}
					}
				} //Data Parallel section 4 done
				__syncthreads();

				if (threadIdx.x >= 2 && threadIdx.x < (S_col - 1)) {
					i = threadIdx.x;
					//if (threadIdx.x == 0) {
					//for (i = 2; i < (S_col - 1); i++) {	//Data Parallel section 5
					if (S_MAT[((S_row - 1) + i * S_row) + base] < 0) {

						c = false; // check needed for race condition here.
						int local_c;
						local_c = *(volatile int*) &c;
						//atomicCAS(&c, local_c, 0);
						//break;	not needed for parallel
					}
					//}//Data Parallel section 5 done here
				}
				__syncthreads();

			} else if (row == -1) {
				if (threadIdx.x == 0) {
					remember[rm] = col;
					c = 1;
					rm++;
				}
				__syncthreads();

				temp_index = (S_row - 1) + (col * S_row) + base; //if col==-1 than problem for base==0 i.e. temp_index==-1
				S_MAT[temp_index] = -1 * S_MAT[temp_index]; //remembering by making positive

				if (threadIdx.x >= 2 && threadIdx.x < (S_col - 1)) {
					i = threadIdx.x;
					//if (threadIdx.x == 0) {
					//for (i = 2; i < (S_col - 1); i++) {		//Data parallel section 6
					if ((S_MAT[((S_row - 1) + i * S_row) + base] < 0)) {
						c = false; // check needed for race condition here.
						int local_c;
						local_c = *(volatile int*) &c;
						//atomicCAS(&c, local_c, 0);
						//break;	not needed for parallel
					}
				} //Data parallel section 6 done
				__syncthreads();
			}
		} //end of while
		__syncthreads();
		if (threadIdx.x == 0) {
			//printf("Value = %f ",S_MAT[(S_row - 1 + (S_col - 1) * S_row) + base]);
			Result[index] = S_MAT[(S_row - 1 + (S_col - 1) * S_row) + base];
		}
	} //endif
}

__host__ Simplex::Simplex(unsigned int N_S) {
	number_of_LPs = N_S;
	M = 0;
	N = 0;
	c = 0;
	No_c = 0;
	/*unsigned int memSize = N_S * sizeof(float);
	 R = (float*) malloc(memSize);*/
}
// ************ All New Implementation ********
__host__ Simplex::Simplex(int UnitBall, unsigned int N_S_for_X) {
	//New implementation
	number_of_LPs_for_X = N_S_for_X;
}

//get status of particular simplex
__host__ int Simplex::getStatus(int n) {
	int s;
	for (int i = 0; i < C.size1(); i++) {
		if (i == (n - 1)) {
			if (R[i] == -1) {
				s = 6; // 6 = Simplex Is Unbounded
			} else if (R[i] > 0) {
				s = 2; // 2= Simplex has feasible and Optimal solution
			}
		}
	}
	return s;
} //get status of particular simplex

//get the No of simplex the object is ruuning on GPU
__host__ int Simplex::getNo_OF_Simplx() {
	return C.size1();
} //get the No of simplex the object is ruuning on GPU

//get the result of all simplex
__host__ std::vector<float> Simplex::getResultAll() {

	std::vector<float> Res(C.size1());
	for (int i = 0; i < C.size1(); i++) {
		Res[i] = R[i];
	}
	return Res;
}
__host__ void Simplex::getResult_X(std::vector<float> &result) { //New implementation
	result.resize(number_of_LPs_for_X);
	for (int i = 0; i < number_of_LPs_for_X; i++) {
		result[i] = R_X[i];
	}

}
__host__ void Simplex::getResult_U(std::vector<float> &result) { //New implementation
	result.resize(number_of_LPs_for_U);
	for (int i = 0; i < number_of_LPs_for_U; i++) {
		result[i] = R_U[i];
	}
}
__host__ void Simplex::getResult_UnitBall(std::vector<float> &result) { //New implementation
	result.resize(number_of_LPs_for_X);
	for (int i = 0; i < number_of_LPs_for_X; i++) {
		result[i] = R_UnitBall[i];
	}
}
__host__ void Simplex::getResult_dotProduct(std::vector<float> &result) {
	result.resize(number_of_LPs_for_X);
	for (int i = 0; i < number_of_LPs_for_X; i++) {
		result[i] = R_dotProduct[i];
	}
}

//get the result of all simplex

__host__ float Simplex::getResult(int n) {
// get result of particular simplex
	float r;
	for (int i = 0; i < C.size1(); i++) {
		if (i == (n - 1)) {
			r = R[i];
		}
	}
	return r;
} // get result of particular simplex

__host__ std::vector<int> Simplex::getStatusAll() {

	std::vector<int> Status(C.size1());
	for (int i = 0; i < C.size1(); i++) {
		if (R[i] == -1)
			Status[i] = 6;
		else
			Status[i] = 2;
	}
	return Status;
} //get the status of all simplex

__host__ void Simplex::setConstratint(math::matrix<double> A,
		std::vector<double> B) { //Old implementation
	unsigned int N_S = number_of_LPs;
	orig_CoefficientMatrix = A;
	BoundValue = B;

	bool single_bound_flag = false;
	int count = 0;

	for (int i = 0; i < A.size1(); i++) {
		single_bound_flag = false;
		for (int j = 0; j < A.size2(); j++) {
			if (A(i, j) != 0) {
				count++; //keeping track of the single variable
			}
		}
		if (count == 1) {
			single_bound_flag = true;
			count = 0;
		}
	}
//single_bound_flag = false;
	single_bound_flag_result = single_bound_flag;

	if (!single_bound_flag) { //If not single bound value
		//std::cout << "\nNot a Boundary Check Model!!!!\n";

		int No_O = A.size2();
		int No_C = A.size1();
		M = No_C + 1;
		N = No_O + 3 + No_C;
		c = 1 + No_O;
		//MAT = (float *) calloc(N_S * M * N, sizeof(float));
		unsigned int memSize = N_S * M * N * sizeof(float);
		cudaError_t err;
		err = cudaMallocHost(&MAT, memSize); //Pinned memory Syntax:: cudaMallocHost(&h_ptr,bytes);
		printf("CUDA cudaMallocHost-- MAT: %s\n", cudaGetErrorString(err));
		cudaMemset(MAT, 0, memSize); //initializing all elements to zero
		//std::cout << "Before OMP for MAT!!!" << std::endl;
#pragma omp parallel for
		for (int s = 0; s < N_S; s++) {
			for (int i = 0; i < M - 1; i++) {
				for (int j = 0; j < N; j++) {
					if (j == 0) {
						MAT[(int) ((i + j * M) + (M * N * s))] = c + i;
					} else if (j > 1) {
						if (j < (No_O + 2)) {
							//Coefficient of A
							MAT[(int) ((i + j * M) + (M * N * s))] = (float) A(
									i, j - 2);
						} else if (j == N - 1) {
							MAT[(int) ((i + j * M) + (M * N * s))] =
									(float) B[i];
						} else if (j < N - 1) {
							MAT[(int) ((i + (No_O + 2 + i) * M) + (M * N * s))] =
									1;

						}
					}
				}
			}
		} //pragma for-loop
		  //	std::cout << "Going out from setConstraints()!!!" << std::endl;
	} //endif
} //setting constraints of simplex

__host__ void Simplex::setConstratint(math::matrix<double> coeffMatrix_for_X,
		std::vector<double> columVector_for_X, int UnitBall) { //New implementation
	BoundValue_for_X = columVector_for_X;

	bool single_bound_flag_for_X = false;
	int count = 0;

	for (int i = 0; i < coeffMatrix_for_X.size1(); i++) {
		single_bound_flag_for_X = false;
		for (int j = 0; j < coeffMatrix_for_X.size2(); j++) {
			if (coeffMatrix_for_X(i, j) != 0) {
				count++; //keeping track of the single variable
			}
		}
		if (count == 1) {
			single_bound_flag_for_X = true;
			count = 0;
		}
	}

//single_bound_flag = false;
	if (single_bound_flag_for_X) {
		single_bound_flag_result = true;
	} else {
		single_bound_flag_result = false;
		std::cout << "\nModel does not has Bounded Variable!!!\n"; //Todo::Take care of this condition
	}
}

__host__ void Simplex::ComputeLP(math::matrix<float> &C1,
		unsigned int number_of_streams) { //OLD implementation
	cudaError_t err;
//If the size of threads exceeds 32(warp-size) be careful of (unknown error/deadlock/abort) situation
//For helicopter the variable=28+slack(56 constraints) + 0 (artificial) + 3(extra) = 87 NEEDS 32 * INT(87/32) + 32
//For Five Dim System the variable=5+slack(10 constraints) + 1 (artificial) + 3(extra) = 19 NEEDS 32 * INT(19/32) + 32
	unsigned int threads_per_block; //Maximum threads depends on CC 1.x =512 2.x and > = 1024

	unsigned int number_of_blocks; //depends on our requirements (better to be much more than the number of SMs)
	int N_S = C1.size1();
	int No_C = orig_CoefficientMatrix.size1();
	C = math::matrix<float>(C1);
//single_bound_flag_result = true;

	unsigned int memSize = N_S * sizeof(float);
//R = (float*) malloc(memSize);
	cudaMallocHost(&R, memSize); //PINNED Memory		//Doing here First Time

	if (single_bound_flag_result) {
		//	printf("\nOur New Optimizations Result %d\n", single_bound_flag_result);
		unsigned int mysize = N_S * sizeof(float);
		//cudaMallocHost(&R, mysize);	//PINNED Memory		//Doing here First Time
		err = cudaMalloc((void **) &G_R, mysize); //Doing here First Time

		mysize = BoundValue.size() * sizeof(float);
		err = cudaMalloc((void **) &d_bound_result, mysize); //C1.size1() * 96 being the maximum threads
		//	printf("CUDA malloc d_bound_result: %s\n", cudaGetErrorString(err));
		mysize = C1.size1() * C1.size2() * sizeof(float);
		err = cudaMalloc((void **) &d_obj_coeff, mysize); //C1.size1() * 96 being the maximum threads
		//	printf("CUDA malloc d_obj_coeff: %s\n", cudaGetErrorString(err));

		mysize = BoundValue.size() * sizeof(float);
		//only single variable per constraints
		//bound_result = (float *) malloc(BoundValue.size() * sizeof(float));	//creating bound_values for kernel compute
		cudaMallocHost(&bound_result, mysize); //Pinned memory Syntax:: cudaMallocHost(&h_ptr,bytes);
		//printf("CUDA cudaMallocHost-- bound_result: %s\n", cudaGetErrorString(err));
		for (int i = 0; i < BoundValue.size(); i++)
			bound_result[i] = BoundValue[i];

		mysize = C1.size1() * C1.size2() * sizeof(float);
		//obj_coeff = (float *) malloc(C1.size1() * C1.size2() * sizeof(float));//creating obj_coeff_values for kernel compute
		cudaMallocHost(&obj_coeff, mysize); //Pinned memory Syntax:: cudaMallocHost(&h_ptr,bytes);
		//printf("CUDA cudaMallocHost-- obj_coeff: %s\n", cudaGetErrorString(err));

		threads_per_block = C1.size2(); //dimension size is the number of threads per block
		//threads_per_block = 512;	//1024;	//Maximum on CC 3.0

		// **** Begin of Stream Processing *******	//Using Asynchronous Memory copy:: needs //MAT to be a PINNED memory
		//100% Occupancy with 1024 threads using 4 registers shows 2 active thread blocks per Multiprocessor
		// we have 7 SM = 7 x 2 = 14 blocks
		//14 blocks x 1024 threads(or LPs) = 14336 LPs per Stream
		int num_streams = number_of_streams; //number of streams desired to create ::Note check for odd numbers
		int num_LPs_perStream;

		if (N_S % num_streams == 0) {
			num_LPs_perStream = (N_S / num_streams);
		} else {
			num_LPs_perStream = (N_S / num_streams); //last will not be the same will be less
			num_streams = num_streams + 1; //one extra streams. though nos of LPs to be solved will be less;
		}
		cudaStream_t stream[num_streams];
		cudaError_t result;
		//int Each_LP_size = M * N;	// * sizeof(float);

		//Creation of Streams
		for (int i = 0; i < num_streams; i++) {
			result = cudaStreamCreate(&stream[i]);
		}

		cudaMemcpy(d_bound_result, bound_result,
				(BoundValue.size() * sizeof(float)), cudaMemcpyHostToDevice);
		//cudaMemcpy(d_obj_coeff, obj_coeff,(C1.size1() * C1.size2() * sizeof(float)),cudaMemcpyHostToDevice);
#pragma omp parallel for
		for (unsigned int s = 0; s < C1.size1(); s++) {
			for (int i = 0; i < C1.size2(); i++) {
				obj_coeff[s * C1.size2() + i] = C1(s, i); //Row major Copy
			}
		}

		//cudaMemcpy(G_MAT, MAT, (N_S * M * N * sizeof(float)),	cudaMemcpyHostToDevice);
		//Stream -- memcopy Host to Device
		for (int i = 0; i < num_streams; i++) {
			int offset = i * num_LPs_perStream * C1.size2();
			/*cudaMemcpyAsync(&d_bound_result, &bound_result,
			 (BoundValue.size() * sizeof(float)), cudaMemcpyHostToDevice,
			 stream[i]);*/

			cudaMemcpyAsync(&d_obj_coeff[offset], &obj_coeff[offset],
					(num_LPs_perStream * C1.size2() * sizeof(float)),
					cudaMemcpyHostToDevice, stream[i]); //TODO:: What happens to the Last stream
		}
		/*for (int i = 0; i < num_streams; i++) {
		 cudaMemcpyAsync(&d_bound_result, &bound_result,
		 (BoundValue.size() * sizeof(float)), cudaMemcpyHostToDevice,
		 stream[i]);
		 }*/

		int dimension = C1.size2();
		//std::cout << "\nnum_LPs_perStream  = " << num_LPs_perStream << "\n";

		//Stream -- Kernel		//streaming is easy if each block solves one LP
		for (int i = 0; i < num_streams; i++) {
			int offset = i * num_LPs_perStream; //for result here offset_res is a pointer to the LP number
			//bound_simplex<<<grid_size, threads_per_block>>>(d_bound_result, d_obj_coeff, G_R, dimension, C1.size1());
			bound_simplex<<<num_LPs_perStream, threads_per_block, 0, stream[i]>>>(
					d_bound_result, d_obj_coeff, G_R, dimension, N_S, offset);
		}

		//Stream -- memcopy Device to Host
		for (int i = 0; i < num_streams; i++) {
			int offset = i * num_LPs_perStream; //for result here offset_res is a pointer to the LP number
			cudaMemcpyAsync(&R[offset], &G_R[offset],
					num_LPs_perStream * sizeof(float), cudaMemcpyDeviceToHost,
					stream[i]); //TODO:: What happens to the Last stream
		}
		// **** End of Stream Processing *******
		//cudaMemcpy(R, G_R, C1.size1() * sizeof(float), cudaMemcpyDeviceToHost);
		cudaFree(G_R);
		cudaFree(d_bound_result);
		cudaFree(d_obj_coeff);

		cudaFreeHost(obj_coeff); //Newly Added but not tested
		cudaFreeHost(bound_result); //Newly Added but not tested
	} else { //otherwise normal simplex algorithm

		unsigned int threads_per_block; //Maximum threads depends on CC 1.x =512 2.x and > = 1024
		unsigned int number_of_blocks; //depends on our requirements (better to be much more than the number of SMs)
		int device;
		cudaDeviceProp props;
		cudaGetDevice(&device);
		cudaGetDeviceProperties(&props, device);

		int No_O = C.size2();
		M = No_C + 1, N = No_O + 3 + No_C;
		int N_C = No_C;
		c = 1 + No_O;
#pragma omp parallel for
		for (int s = 0; s < N_S; s++) {
			for (int i = M - 1; i < M; i++) { //int i=M-1; //should be enough
				for (int j = 2; j < N; j++) {
					if (j < 2 + No_O) {
						MAT[(int) ((i + j * M) + (M * N * s))] = -C(s, j - 2);
					}
				}
			}
		}
		std::vector<int> rem;
		for (int i = 0; i < N_C; i++) {
			if (BoundValue[i] < 0) {
				rem.push_back(i);
			}
		}
		//std::cout << "C= " << rem.size() << "\n";
		int nc = N + rem.size();
		threads_per_block = 32 * (nc / 32) + 32; //if count equal 0 than nc=N so works for for Model
		if (threads_per_block > props.maxThreadsPerBlock) //Assuming maximum threads supported by CC is 1024
			threads_per_block = props.maxThreadsPerBlock;

		int offset; //temporary fixing for this if- part
		//	int num_LPs_perStream;	//temporary fixing for this if- part

		int *R_index; //reduction data
		float *R_data; //reduction index
		err = cudaMalloc((void **) &R_data,
				C1.size1() * threads_per_block * sizeof(float)); //C1.size1() * 96 being the maximum threads
		err = cudaMalloc((void **) &R_index,
				C1.size1() * threads_per_block * sizeof(int)); //C1.size1() being the number of LPs

		err = cudaMalloc((void **) &G_R, N_S * sizeof(float)); //Doing it here for the First Time

		if (rem.size() > 0) { //Helicopter model has no negative bound so count=0
			N_MAT = (float *) calloc(N_S * M * nc, sizeof(float));

#pragma omp parallel for
			for (int i = 0; i < N_S; i++) {
				for (int j = 0; j < M; j++) {
					int base = i * M * N;
					int basen = i * M * nc;
					N_MAT[j + ((nc - 1) * M) + basen] = MAT[j + ((N - 1) * M)
							+ base];
				}
			}

//Creating Artificial Variables
#pragma omp parallel for
			for (int k = 0; k < N_S; k++) {
				int base = k * M * N;
				int basen = k * M * nc;
				for (int i = 0; i < rem.size(); i++) {
					//ch = 0;
					for (int j = 0; j < nc; j++) {
						if (MAT[rem[i] + ((N - 1) * M) + base] < 0) {
							if ((j >= (N - 1)) && (j < (nc - 1))) {
								N_MAT[(rem[i] + j * M) + basen] = 1;
							} else if (j == (nc - 1)) {
								N_MAT[(rem[i] + j * M) + basen] = -1
										* N_MAT[(rem[i] + j * M) + basen];
							} else if (j == 1) {
								N_MAT[(rem[i] + j * M) + basen] = -1;
							} else if (j == 0) {
								N_MAT[((rem[i] + j * M)) + basen] = ((N - M)
										+ rem[i] + 3 + No_O);
							} else if (j > 1) {
								N_MAT[(rem[i] + j * M) + basen] = -1
										* (MAT[(rem[i] + j * M) + base]);
							}
						}
					}
				}
			}

#pragma omp parallel for
			for (int k = 0; k < N_S; k++) {
				int base = k * M * N;
				int basen = k * M * nc;
				for (int j = (N - 1); j < (nc - 1); j++) {
					N_MAT[((M - 1) + j * M) + basen] = -1;
				}
			}

//Creation of Last Row or Z-Value(Z-C)
#pragma omp parallel for
			for (int k = 0; k < N_S; k++) {
				int basen = k * M * nc;
				for (int k1 = 2; k1 < nc; k1++) {
					float sum = 0;
					for (int j = 0; j < (M - 1); j++) {
						sum = sum
								+ (N_MAT[(j + k1 * M) + basen]
										* N_MAT[(j + 1 * M) + basen]);
					}
					N_MAT[((M - 1) + k1 * M) + basen] = sum
							- N_MAT[((M - 1) + k1 * M) + basen];
				}
			}

			cudaMalloc((void **) &G_MAT, (N_S * M * nc * sizeof(float)));
			//printf("CUDA malloc G_MAT: %s\n", cudaGetErrorString(err));
			cudaMemcpy(G_MAT, N_MAT, (N_S * M * nc * sizeof(float)),
					cudaMemcpyHostToDevice);
			//printf("CUDA malloc N_MAT: %s\n", cudaGetErrorString(err));
			//	cudaMemcpy(G_Sel, Sel, sizeof(int), cudaMemcpyHostToDevice);
			//printf("CUDA malloc G_Sel: %s\n", cudaGetErrorString(err));
			//mykernel<<<N_S, threads_per_block>>>(G_MAT, M, nc, G_R, N_S, G_Sel, R_data, R_index);
			//mykernel<<<N_S, threads_per_block>>>(G_MAT, M, nc, G_R, N_S, R_data, R_index);
			mykernel<<<N_S, threads_per_block>>>(G_MAT, M, nc, G_R, N_S, R_data,
					R_index, offset);

			cudaDeviceSynchronize();
			cudaMemcpy(R, G_R, N_S * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(N_MAT, G_MAT, (N_S * M * nc * sizeof(float)),
					cudaMemcpyDeviceToHost);

			//	std::cout << "Result for Artificial\n";
#pragma omp parallel for
			for (int i = 0; i < N_S; i++) {
				int base = i * M * N;
				int basen = i * M * nc;
				for (int j = 0; j < M; j++) {
					for (int k = 0; k < N; k++) {
						if (N_MAT[j + basen] == k + 1) {
							N_MAT[j + M + basen] = MAT[(M - 1) + (2 + k) * M
									+ base];
						}
					}
				}
			}

#pragma omp parallel for
			for (int s = 0; s < N_S; s++) {
				if (R[s] == 0) {
					int base = s * M * N;
					int basen = s * M * nc;
					for (int i = 0; i < N; i++) {
						float sum = 0;
						for (int j = 0; j < M; j++) {
							if ((j < M - 1)) {
								if (i != (N - 1)) {
									sum = sum
											+ (N_MAT[(j + (i * M)) + basen]
													* N_MAT[(j + (1 * M))
															+ basen]);
									MAT[(j + (i * M)) + base] = N_MAT[(j
											+ (i * M)) + basen];
								} else if (i == N - 1) {
									sum = sum
											+ (N_MAT[(j + (nc - 1) * M) + basen]
													* N_MAT[(j + (1 * M))
															+ basen]);
									MAT[(j + (i * M)) + base] = N_MAT[(j
											+ (nc - 1) * M) + basen];
								}
							}
							if (j == (M - 1)) {
								if (i > 1) {
									MAT[(j + (i * M)) + base] =
											MAT[(j + (i * M)) + base]
													+ (-1) * sum;
								}
							}
						}
					}
				}
			}

			cudaFree(G_MAT);
			//		cudaFree(G_R);
			//cudaFree(G_Sel);
			//cudaDeviceSynchronize();
			//		cudaMalloc((void **) &G_R, N_S * sizeof(float));
			cudaMalloc((void **) &G_MAT, (N_S * M * N * sizeof(float)));
			// printf("CUDA malloc G_MAT: %s\n", cudaGetErrorString(err));
			cudaMemcpy(G_MAT, MAT, (N_S * M * N * sizeof(float)),
					cudaMemcpyHostToDevice);
			//printf("CUDA malloc N_MAT: %s\n", cudaGetErrorString(err));
			//	cudaMemcpy(G_Sel, Sel, sizeof(int), cudaMemcpyHostToDevice);
			//printf("CUDA malloc G_Sel: %s\n", cudaGetErrorString(err));

			//mykernel<<<N_S, threads_per_block>>>(G_MAT, M, N, G_R, N_S, G_Sel, R_data, R_index);
			//mykernel<<<N_S, threads_per_block>>>(G_MAT, M, N, G_R, N_S, R_data, R_index);
			mykernel<<<N_S, threads_per_block>>>(G_MAT, M, N, G_R, N_S, R_data,
					R_index, offset);
			cudaDeviceSynchronize();
			cudaMemcpy(R, G_R, N_S * sizeof(float), cudaMemcpyDeviceToHost);

			/*cudaMemcpy(MAT, G_MAT, (N_S * M * N * sizeof(float)),
			 cudaMemcpyDeviceToHost);*/
			//std::cout
			// << "***********Final SIMPLEX from GPU*************\n Time took:\n";*/
		} else {

			err = cudaMalloc((void **) &G_MAT, (N_S * M * N * sizeof(float)));
			//printf("CUDA malloc G_MAT : %s\n", cudaGetErrorString(err));
			// **** Begin of Stream Processing *******
			//Using Asynchronous Memory copy:: needs //MAT to be a PINNED memory
			int num_streams = number_of_streams; //number of streams desired to create ::Note check for odd numbers
			int Each_LP_size = M * N; // * sizeof(float);
			int num_LPs_perStream;
			bool equal_stream = true;
			if (N_S % num_streams == 0) {
				num_LPs_perStream = (N_S / num_streams);
				equal_stream = true;
			} else {
				num_LPs_perStream = (N_S / num_streams); //last stream will not be of the same size
				num_streams = num_streams + 1; //one extra streams.where nos of LPs to be solved will be less;
				equal_stream = false;
			}
			cudaStream_t stream[num_streams];
			cudaError_t result;

			//Creation of Streams
			for (int i = 0; i < num_streams; i++) {
				result = cudaStreamCreate(&stream[i]);
			}

			//err = cudaMemcpy(G_MAT, MAT, (N_S * M * N * sizeof(float)),cudaMemcpyHostToDevice);
			//Stream -- memcopy Host to Device
			std::cout << "\nNumber of LPs_perStream = " << num_LPs_perStream
					<< std::endl;
			unsigned int lastBlock_size;
			if (equal_stream == false) {
				lastBlock_size = N_S
						- (N_S / (num_streams - 1)) * (num_streams - 1); //LAST Stream Size
				std::cout << "\nAmit Last Block size (LPs is )= "
						<< lastBlock_size << std::endl;
			}

			for (int i = 0; i < num_streams; i++) {
				if (equal_stream == false && i == (num_streams - 1)) { //last stream
					int offset = i * Each_LP_size * lastBlock_size; //for memory copy
					cudaMemcpyAsync(&G_MAT[offset], &MAT[offset],
							(lastBlock_size * M * N * sizeof(float)),
							cudaMemcpyHostToDevice, stream[i]);
				} else {
					int offset = i * Each_LP_size * num_LPs_perStream; //for memory copy
					cudaMemcpyAsync(&G_MAT[offset], &MAT[offset],
							(num_LPs_perStream * M * N * sizeof(float)),
							cudaMemcpyHostToDevice, stream[i]);
				}
			}
			//mykernel<<<N_S, threads_per_block>>>(G_MAT, M, N, G_R, N_S, G_Sel, R_data, R_index);
			//	mykernel<<<N_S, threads_per_block>>>(G_MAT, M, N, G_R, N_S, R_data, R_index);
			//std::cout << "Before Kernel Call!!!" << std::endl;
			//Stream -- Kernel
			for (int i = 0; i < num_streams; i++) {
				if (equal_stream == false && i == (num_streams - 1)) { //last stream
					int offset_res = i * lastBlock_size; //for result here offset_res is a pointer to the LP number
					//mykernel<<<num_LPs_perStream, 256, 0, stream[i]>>>(G_MAT, M, N, G_R, G_Sel, num_LPs_perStream, offset_res);
					mykernel<<<lastBlock_size, threads_per_block, 0, stream[i]>>>(
							G_MAT, M, N, G_R, lastBlock_size, R_data, R_index,
							offset_res);
				} else {
					int offset_res = i * num_LPs_perStream; //for result here offset_res is a pointer to the LP number
					//mykernel<<<num_LPs_perStream, 256, 0, stream[i]>>>(G_MAT, M, N, G_R, G_Sel, num_LPs_perStream, offset_res);
					mykernel<<<num_LPs_perStream, threads_per_block, 0,
							stream[i]>>>(G_MAT, M, N, G_R, num_LPs_perStream,
							R_data, R_index, offset_res);
				}

			}
			//	std::cout << "After Kernel Call!!!" << std::endl;
			//cudaDeviceSynchronize();//removed as hopping that cudaFree will handle it
			//err = cudaMemcpy(R, G_R, N_S * sizeof(float),cudaMemcpyDeviceToHost);

			//Stream -- memcopy Device to Host
			for (int i = 0; i < num_streams; i++) {
				if (equal_stream == false && i == (num_streams - 1)) { //last stream
					int offset_res = i * lastBlock_size; //for result here offset_res is a pointer to the LP number
					cudaMemcpyAsync(&R[offset_res], &G_R[offset_res],
							(lastBlock_size * sizeof(float)),
							cudaMemcpyDeviceToHost, stream[i]);
				} else {
					int offset_res = i * num_LPs_perStream; //for result here offset_res is a pointer to the LP number
					cudaMemcpyAsync(&R[offset_res], &G_R[offset_res],
							(num_LPs_perStream * sizeof(float)),
							cudaMemcpyDeviceToHost, stream[i]);
				}
			}
			// **** End of Stream Processing *******

			//printf("CUDA memcpy G_R: %s\n", cudaGetErrorString(err));
			//	std::cout << "N_S = " << N_S << std::endl;
		}
		cudaFree(G_MAT);
		cudaFree(G_R);
		cudaFree(R_index); //Only to synchronise with the cudamemcpy
		cudaFree(R_data); //Only to synchronise with the cudamemcpy
		cudaFreeHost(MAT); //This is needed to avoid Segmentation fault error
		//cudaDeviceReset();
	}
}

/*	Used here due to lack of proper place
 * numVectors:: number of Template Directions/Vectors
 * NewTotalIteration :: Iterations
 */

__host__ void Simplex::ComputeLP(math::matrix<float> &obj_funs_for_X,
		int UnitBall, unsigned int number_of_streams) { //NEW implementation
	cudaError_t err;
	unsigned int threads_per_block; //Maximum threads depends on CC 1.x =512 2.x and > = 1024
	unsigned int number_of_blocks; //depends on our requirements (better to be much more than the number of SMs)
	int N_S = this->number_of_LPs_for_X;
//std::cout<<"N_S = "<<N_S<<"\n";
//single_bound_flag_result = true;

	unsigned int memSize = N_S * sizeof(float);
//R = (float*) malloc(memSize);
	cudaMallocHost(&R_X, memSize); //PINNED Memory		//Doing here First Time
	cudaMallocHost(&R_UnitBall, memSize); //PINNED Memory		//Doing here First Time

	if (single_bound_flag_result) {
		//printf("\nOur New Optimizations Result %d\n", single_bound_flag_result);
		unsigned int mysize = N_S * sizeof(float);
		//cudaMallocHost(&R, mysize);	//PINNED Memory		//Doing here First Time
		err = cudaMalloc((void **) &G_R_X, mysize); //Doing here First Time
		err = cudaMalloc((void **) &G_R_UnitBall, mysize); //Doing here First Time

		mysize = BoundValue_for_X.size() * sizeof(float);
		err = cudaMalloc((void **) &d_bound_result_for_X, mysize); //C1.size1() * 96 being the maximum threads
		//	printf("CUDA malloc d_bound_result: %s\n", cudaGetErrorString(err));
		mysize = obj_funs_for_X.size1() * obj_funs_for_X.size2()
				* sizeof(float);
		err = cudaMalloc((void **) &d_obj_coeff_for_X, mysize); //C1.size1() * 96 being the maximum threads
		//	printf("CUDA malloc d_obj_coeff: %s\n", cudaGetErrorString(err));
		//mysize = obj_funs_for_X.size1() * obj_funs_for_X.size2()* sizeof(float);
		//obj_coeff = (float *) malloc(C1.size1() * C1.size2() * sizeof(float));//creating obj_coeff_values for kernel compute
		cudaMallocHost(&obj_coeff_for_X, mysize); //Pinned memory Syntax:: cudaMallocHost(&h_ptr,bytes);
		//printf("CUDA cudaMallocHost-- obj_coeff: %s\n", cudaGetErrorString(err));
		int dimension = obj_funs_for_X.size2();

		mysize = BoundValue_for_X.size() * sizeof(float);
		//only single variable per constraints
		//bound_result = (float *) malloc(BoundValue.size() * sizeof(float));	//creating bound_values for kernel compute
		cudaMallocHost(&bound_result_for_X, mysize); //Pinned memory Syntax:: cudaMallocHost(&h_ptr,bytes);
		//printf("CUDA cudaMallocHost-- bound_result: %s\n", cudaGetErrorString(err));
		for (int i = 0; i < BoundValue_for_X.size(); i++) {
			bound_result_for_X[i] = BoundValue_for_X[i];
		}
		threads_per_block = obj_funs_for_X.size2(); //dimension size is the number of threads per block
		//threads_per_block = 512;	//1024;	//Maximum on CC 3.0

		// **** Begin of Stream Processing *******	//Using Asynchronous Memory copy:: needs //MAT to be a PINNED memory
		//100% Occupancy with 1024 threads using 4 registers shows 2 active thread blocks per Multiprocessor

		int num_streams = number_of_streams; //number of streams desired to create ::Note check for odd numbers
		int num_LPs_perStream;

		if (N_S % num_streams == 0) { //Maximum nos. of LPs is of UnitBall_infinity_norm
			num_LPs_perStream = (N_S / num_streams);
		} else {
			num_LPs_perStream = (N_S / num_streams); //last will not be the same will be less
			num_streams = num_streams + 1; //one extra streams. though nos of LPs to be solved will be less in the extra stream;
		}
		//std::cout<<"num_LPs_perStream = " <<num_LPs_perStream<<"\n";
		cudaStream_t stream[num_streams];
		cudaError_t result;
		//int Each_LP_size = M * N;	// * sizeof(float);

		//Creation of Streams
		for (int i = 0; i < num_streams; i++) {
			result = cudaStreamCreate(&stream[i]);
		}

		cudaMemcpy(d_bound_result_for_X, bound_result_for_X,
				(BoundValue_for_X.size() * sizeof(float)),
				cudaMemcpyHostToDevice);
		//cudaMemcpy(d_obj_coeff, obj_coeff,(C1.size1() * C1.size2() * sizeof(float)),cudaMemcpyHostToDevice);
#pragma omp parallel for
		for (unsigned int s = 0; s < obj_funs_for_X.size1(); s++) {
			for (int i = 0; i < obj_funs_for_X.size2(); i++) {
				obj_coeff_for_X[s * obj_funs_for_X.size2() + i] =
						obj_funs_for_X(s, i); //Row major Copy
				//std::cout<<" obj_coeff_for_X = "<<obj_coeff_for_X[s * obj_funs_for_X.size2() + i]<<"\t";
			}
		}
		//cudaMemcpy(G_MAT, MAT, (N_S * M * N * sizeof(float)),	cudaMemcpyHostToDevice);
		//Stream -- memcopy Host to Device
		for (int i = 0; i < num_streams; i++) {
			int offset = i * num_LPs_perStream * obj_funs_for_X.size2();
			cudaMemcpyAsync(&d_obj_coeff_for_X[offset],
					&obj_coeff_for_X[offset],
					(num_LPs_perStream * obj_funs_for_X.size2() * sizeof(float)),
					cudaMemcpyHostToDevice, stream[i]);
		}
		//std::cout << "\nnum_LPs_perStream  = " << num_LPs_perStream << "\n";
		//Stream -- Kernel		//streaming is easy if each block solves one LP
		for (int i = 0; i < num_streams; i++) {
			int offset = i * num_LPs_perStream; //for result here offset_res is a pointer to the LP number
			//bound_simplex<<<grid_size, threads_per_block>>>(d_bound_result, d_obj_coeff, G_R, dimension, C1.size1());
			bounded_Solver_for_X_UnitBall<<<num_LPs_perStream,
					threads_per_block, 0, stream[i]>>>(d_bound_result_for_X,
					d_obj_coeff_for_X, G_R_X, N_S, G_R_UnitBall, offset);
		}

		//Stream -- memcopy Device to Host
		for (int i = 0; i < num_streams; i++) {
			int offset = i * num_LPs_perStream; //for result here offset_res is a pointer to the LP number
			cudaMemcpyAsync(&R_X[offset], &G_R_X[offset],
					num_LPs_perStream * sizeof(float), cudaMemcpyDeviceToHost,
					stream[i]);
		}
		mysize = N_S * sizeof(float); //dimension is same for UnitBall's direction also
		cudaMemcpy(R_UnitBall, G_R_UnitBall, mysize, cudaMemcpyDeviceToHost); //TODO:: make this as Async memcopy this is Big Chunk
		printf("CUDA cudaMemcpy-- R_UnitBall: %s\n", cudaGetErrorString(err));

		// **** End of Stream Processing *******
		//cudaMemcpy(R, G_R, C1.size1() * sizeof(float), cudaMemcpyDeviceToHost);

		cudaFree(G_R_X);
		cudaFree(d_bound_result_for_X);
		cudaFree(d_obj_coeff_for_X);
		cudaFree(G_R_UnitBall);

		cudaFreeHost(obj_coeff_for_X); //Newly Added but not tested
		cudaFreeHost(bound_result_for_X); //Newly Added but not tested
	} else { //otherwise normal simplex algorithm
		std::cout << "\nModel does not has Bounded Variable!!!\n";
		//TODO:: take care of this condition
	}
}

//New implementation to replace the above implementation
//here obj_funs_for_X/dotProduct(U.C * directions)/UnitBall is the List of all the directions on which support function is to be computed
__host__ void Simplex::ComputeLP(math::matrix<float> &obj_funs_for_X,
		int UnitBall, unsigned int number_of_streams,
		std::vector<double> sysU_C) {
	/*for (int i=0;i<sysU_C.size();i++){
		cout<<sysU_C[i]<<"\t";
	}
	cout<<endl;*/

	cudaError_t err;
	unsigned int threads_per_block; //Maximum threads depends on CC 1.x =512 2.x and > = 1024
	unsigned int number_of_blocks; //depends on our requirements (better to be much more than the number of SMs)
	int N_S = this->number_of_LPs_for_X;
	//std::cout<<"N_S = "<<N_S<<"\n";
	//single_bound_flag_result = true;

	unsigned int memSize = N_S * sizeof(float);
	//R = (float*) malloc(memSize);
	cudaMallocHost(&R_X, memSize); //PINNED Memory		//Doing here First Time
	cudaMallocHost(&R_UnitBall, memSize); //PINNED Memory		//Doing here First Time
	cudaMallocHost(&R_dotProduct, memSize); //PINNED Memory		//Doing here First Time

	if (single_bound_flag_result) {
		//printf("\nOur New Optimizations Result %d\n", single_bound_flag_result);
		unsigned int mysize = N_S * sizeof(float);
		//cudaMallocHost(&R, mysize);	//PINNED Memory		//Doing here First Time
		err = cudaMalloc((void **) &G_R_X, mysize); //Doing here First Time
		err = cudaMalloc((void **) &G_R_UnitBall, mysize); //Doing here First Time

		err = cudaMalloc((void **) &G_R_dotProduct, mysize); //Doing here First Time
		printf("CUDA malloc G_R_dotProduct: %s\n", cudaGetErrorString(err));
		mysize = BoundValue_for_X.size() * sizeof(float);
		err = cudaMalloc((void **) &d_bound_result_for_X, mysize); //C1.size1() * 96 being the maximum threads
		//	printf("CUDA malloc d_bound_result: %s\n", cudaGetErrorString(err));
		mysize = obj_funs_for_X.size1() * obj_funs_for_X.size2()
				* sizeof(float);
		err = cudaMalloc((void **) &d_obj_coeff_for_X, mysize); //C1.size1() * 96 being the maximum threads

		err = cudaMalloc((void **) &d_U_C, sysU_C.size() * sizeof(float)); //C1.size1() * 96 being the maximum threads
		printf("CUDA malloc d_U_C: %s\n", cudaGetErrorString(err));
		//mysize = obj_funs_for_X.size1() * obj_funs_for_X.size2()* sizeof(float);
		//obj_coeff = (float *) malloc(C1.size1() * C1.size2() * sizeof(float));//creating obj_coeff_values for kernel compute
		cudaMallocHost(&obj_coeff_for_X, mysize); //Pinned memory Syntax:: cudaMallocHost(&h_ptr,bytes);
		//printf("CUDA cudaMallocHost-- obj_coeff: %s\n", cudaGetErrorString(err));
		int dimension = obj_funs_for_X.size2();

		mysize = BoundValue_for_X.size() * sizeof(float);
		//only single variable per constraints
		//bound_result = (float *) malloc(BoundValue.size() * sizeof(float));	//creating bound_values for kernel compute
		cudaMallocHost(&bound_result_for_X, mysize); //Pinned memory Syntax:: cudaMallocHost(&h_ptr,bytes);
		//printf("CUDA cudaMallocHost-- bound_result: %s\n", cudaGetErrorString(err));
		for (int i = 0; i < BoundValue_for_X.size(); i++) {
			bound_result_for_X[i] = BoundValue_for_X[i];
		}
		cudaMallocHost(&U_C, sysU_C.size() * sizeof(float)); //Pinned memory Syntax:: cudaMallocHost(&h_ptr,bytes);
		printf("CUDA cudaMallocHost-- U_C: %s\n", cudaGetErrorString(err));
		std::cout<<"U_C[i]\n";
		for (int i = 0; i < sysU_C.size(); i++) {
			U_C[i] = sysU_C[i];
			std::cout<<U_C[i]<<"\t";
		}

		threads_per_block = obj_funs_for_X.size2(); //dimension size is the number of threads per block
		//threads_per_block = 512;	//1024;	//Maximum on CC 3.0

		// **** Begin of Stream Processing *******	//Using Asynchronous Memory copy:: needs //MAT to be a PINNED memory
		//100% Occupancy with 1024 threads using 4 registers shows 2 active thread blocks per Multiprocessor

		int num_streams = number_of_streams; //number of streams desired to create ::Note check for odd numbers
		int num_LPs_perStream;

		if (N_S % num_streams == 0) { //Maximum nos. of LPs is of UnitBall_infinity_norm
			num_LPs_perStream = (N_S / num_streams);
		} else {
			num_LPs_perStream = (N_S / num_streams); //last will not be the same will be less
			num_streams = num_streams + 1; //one extra streams. though nos of LPs to be solved will be less in the extra stream;
		}
		//std::cout<<"num_LPs_perStream = " <<num_LPs_perStream<<"\n";
		cudaStream_t stream[num_streams];
		cudaError_t result;
		//int Each_LP_size = M * N;	// * sizeof(float);

		//Creation of Streams
		for (int i = 0; i < num_streams; i++) {
			result = cudaStreamCreate(&stream[i]);
		}

		cudaMemcpy(d_bound_result_for_X, bound_result_for_X,
				(BoundValue_for_X.size() * sizeof(float)),
				cudaMemcpyHostToDevice);	//this is in the default Stream

		cudaMemcpy(d_U_C, U_C, (sysU_C.size() * sizeof(float)),
				cudaMemcpyHostToDevice);	//this is in the default Stream

		//cudaMemcpy(d_obj_coeff, obj_coeff,(C1.size1() * C1.size2() * sizeof(float)),cudaMemcpyHostToDevice);
#pragma omp parallel for
		for (unsigned int s = 0; s < obj_funs_for_X.size1(); s++) {
			for (int i = 0; i < obj_funs_for_X.size2(); i++) {
				obj_coeff_for_X[s * obj_funs_for_X.size2() + i] =
						obj_funs_for_X(s, i); //Row major Copy
				//std::cout<<" obj_coeff_for_X = "<<obj_coeff_for_X[s * obj_funs_for_X.size2() + i]<<"\t";
			}
		}
		//cudaMemcpy(G_MAT, MAT, (N_S * M * N * sizeof(float)),	cudaMemcpyHostToDevice);
		//Stream -- memcopy Host to Device
		for (int i = 0; i < num_streams; i++) {
			int offset = i * num_LPs_perStream * obj_funs_for_X.size2();
			cudaMemcpyAsync(&d_obj_coeff_for_X[offset],
					&obj_coeff_for_X[offset],
					(num_LPs_perStream * obj_funs_for_X.size2() * sizeof(float)),
					cudaMemcpyHostToDevice, stream[i]);
		}
		//std::cout << "\nnum_LPs_perStream  = " << num_LPs_perStream << "\n";
		//Stream -- Kernel		//streaming is easy if each block solves one LP
		for (int i = 0; i < num_streams; i++) {
			int offset = i * num_LPs_perStream; //for result here offset_res is a pointer to the LP number
			//bound_simplex<<<grid_size, threads_per_block>>>(d_bound_result, d_obj_coeff, G_R, dimension, C1.size1());
			bounded_Solver_for_X_UnitBall_dotProduct<<<num_LPs_perStream,
					threads_per_block, 0, stream[i]>>>(d_bound_result_for_X,
					d_obj_coeff_for_X, d_U_C, G_R_X, N_S, G_R_UnitBall,
					G_R_dotProduct, offset);
		}

		//Stream -- memcopy Device to Host
		for (int i = 0; i < num_streams; i++) {
			int offset = i * num_LPs_perStream; //for result here offset_res is a pointer to the LP number
			cudaMemcpyAsync(&R_X[offset], &G_R_X[offset],
					num_LPs_perStream * sizeof(float), cudaMemcpyDeviceToHost,
					stream[i]);
		}
		mysize = N_S * sizeof(float); //dimension is same for UnitBall's direction also
		cudaMemcpy(R_UnitBall, G_R_UnitBall, mysize, cudaMemcpyDeviceToHost); //TODO:: make this as Async memcopy this is Big Chunk
		printf("CUDA cudaMemcpy-- R_UnitBall: %s\n", cudaGetErrorString(err));

		cudaMemcpy(R_dotProduct, G_R_dotProduct, mysize,
				cudaMemcpyDeviceToHost); //TODO:: make this as Async memcopy this is Big Chunk
		printf("CUDA cudaMemcpy-- R_dotProduct: %s\n", cudaGetErrorString(err));

		// **** End of Stream Processing *******
		//cudaMemcpy(R, G_R, C1.size1() * sizeof(float), cudaMemcpyDeviceToHost);

		cudaFree(G_R_X);
		cudaFree(d_bound_result_for_X);
		cudaFree(d_obj_coeff_for_X);
		cudaFree(G_R_UnitBall);
		cudaFree(G_R_dotProduct);

		cudaFreeHost(obj_coeff_for_X); //Newly Added but not tested
		cudaFreeHost(bound_result_for_X); //Newly Added but not tested
	} else { //otherwise normal simplex algorithm
		std::cout << "\nModel does not has Bounded Variable!!!\n";
		//TODO:: take care of this condition
	}
}
