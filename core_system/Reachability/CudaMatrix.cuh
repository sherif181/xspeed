#define IDX2F(i,j,ld) ((((j))*(ld))+((i)))
#include "cublas_v2.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
template <class type = double>  class CudaPair {
public:
	type data[2];
};
template <class type = double>  class CudaVectorTemplate {
public:

	__device__ static void   vector_join(CudaVectorTemplate* v1,
			CudaVectorTemplate* v2,CudaVectorTemplate* res) {

		CudaVectorTemplate result();
		res->deAlocData();
		res->size=v1->size + v2->size;
		res->alocData();

//		result = v1;
//		unsigned int tot_size;
//		tot_size = result.size + v2.size;
//		result.resize(tot_size);

//		for (unsigned int i = v1.size, j = 0; j < v2.size; i++, j++)
//			result[i] = v2[j];
		for (unsigned int j = 0; j < v1->size; j++){
			res->set(j,v1->get(j));
		}
		for (unsigned int j = v1->size, i=0; j < v2->size; i++,j++){
				res->set(j,v2->get(i));
		}


	}
	__host__ __device__ void assign(CudaVectorTemplate* v) {
			deAlocData();
			size=v->size;
			alocData();
			for (int i = 0; i < size; i++) {

					set(i, v->get(i));


			}

		}
	/// \brief Resize the vector
		/// Resize the vector to a new size. If \c preserve is true, data are copied otherwise data are lost. If the new size is bigger, the remaining values are filled in with the initial value (0 by default) in the case of \c unbounded_array, which is the container by default. If the new size is smaller, last values are lost. This behaviour can be different if you explicitely specify another type of container.
		/// \param size new size of the vector
		/// \param preserve if true, keep values
	__host__ __device__ void resize (int nsize) {
		type* nd = new type[nsize];
		 for (unsigned int i = 0; i < size; i++){
			 nd[i]=get(i);
		 }
		 size=nsize;
		 deAlocData();
		 data=nd;

	}
	__device__ void trim(int nsize){
		 type* nd = new type[nsize];
		 for (unsigned int i = 0; i < min(size,nsize); i++){
			 nd[i]=get(i);
		 }
		 deAlocData();
		 size=nsize;
		 data=nd;
	}
	__device__  static void vector_add(CudaVectorTemplate * v1,
			CudaVectorTemplate*  v2,CudaVectorTemplate* res) {

		res->deAlocData();
		res->size=v2->size;
		res->alocData();

		for (unsigned int i = 0; i < v2->size; i++)
//			result[i] = v1[i] + v2[i];
		    res->set(i,v1->get(i)+v2->get(i));

	}

	__device__ static void copyVector(CudaVectorTemplate* v,CudaVectorTemplate* res) {
		res->deAlocData();
		res->size=v->size;
		res->alocData();

		for (int i = 0; i < v->size; ++i) {
			res->set(i,v->get(i));
		}

	}
	__device__ static type dot_product(CudaVectorTemplate * vector1, CudaVectorTemplate* vector2) {
		type res = 0;
		assert(vector1->size == vector2->size);

		for (int i = 0; i < vector1->size; i++) {
				res = res + vector1->get(i) * vector2->get(i);
		}
		return res;
	}
	type * data=NULL;
	int size=0;
	__host__ __device__ type get(int i) {

		return data[i];

	}
	__host__ __device__ void set(int i, type val) {

		data[i] = val;

	}
	__device__ __host__ void alocData() {
		data = new type[size]();

	}


	__host__ __device__ CudaVectorTemplate() {
		this->size = 0;

	}
	__host__  __device__  CudaVectorTemplate(int size) {
		this->size = size;
		alocData();
	}

	 __device__ CudaVectorTemplate(int size, type initValue) {
		this->size = size; //TODO TEST
		alocData();
		for (int i = 0; i < size; i++) {
			data[i] = initValue;
		}


	}
	__host__ __device__ type minimalElement() {
		type ret = DBL_MAX; //TODO TEST
		for (int i = 0; i < this->size; i++) {
			ret = min(data[i], ret);
		}

		return ret;
	}
//	__host__ __device__ 	void deAlocData();
	__host__ __device__ 	void deAlocData() {
			if (data != NULL) {
						delete[] data;
						data = NULL;
			}

		}
	//TODO DESTRUCTOR

	//TODO THIS CLASS DOES NOT HAVE A DESTRUCTOR !!
};

//template <class type>
//__host__ __device__ void CudaVectorTemplate<type>::deAlocData() {
//		for (int i = 0; i < this->size; i++) {
//			data[i]->deAlocData();
//		}
//		if (data != NULL) {
//					delete[] data;
//					data = NULL;
//		}
//TODO DEALOCATE IF POINTER TYPE
//	}

//template <>
//	__host__ __device__ 	void CudaVectorTemplate<double>::deAlocData() {
//			if (data != NULL) {
//						delete[] data;
//						data = NULL;
//			}
//
//		}





typedef CudaVectorTemplate<> CudaVector;
class CudaMatrix {
public:
    // Resizing
	  /** Resize a matrix to new dimensions
	   * If data are preserved, then if the size if bigger at least on one dimension, extra values are filled with zeros.
	   * If data are not preserved, then nothing has to be assumed regarding the content of the matrix after resizing.
	   * \param size1 the new number of rows
	   * \param size2 the new number of colums
	   * \param preserve a boolean to say if one wants the data to be preserved during the resizing. Default is true.
	   */


	__device__ void resize(int nr,int nc,bool preserve){
			//TODO test
			if(preserve){
				resize(nr,nc);
				return;

			}
			//do not preserve
			deallocateData();
			rows=nr;
			cols=nc;
			alocData();
	}
	__device__ void resize(int nr,int nc){
				//TODO test
				CudaMatrix temp(nr,nc);
				for (int i = 0; i < rows; i++) {
							for (int j = 0; j < cols; j++) {
								temp.set(i, j, get(i, j));
						}

				}
				assign(&temp);
		}

	__device__ bool  inverse(CudaMatrix*  inverse) {
		//bool NonSingular=false;
		using namespace boost::numeric::ublas;
		// create a working copy of the input
//		ublas_matrix_impl A(this->size1(), this->size2(), this->data());
		CudaMatrix A;
		A.assign(this);
		// make A same as the current matrix.
		/*
		 cout << "\nAmit here\n";
		 for (int i = 0; i < A.size1(); i++) {
		 for (int j = 0; j < A.size2(); j++)
		 cout << A(i, j) << "\t";
		 cout << endl;
		 }
		 cout << "\nAmit there\n";
		 */
		// create a permutation matrix for the LU-factorization
//		pmatrix pm(A.size1());
		CudaVector pm(A.rows);

		// perform LU-factorization
		int res = lu_factorize(&A,&pm);
		if (res != 0) {
			//NonSingular=false;
			return false;
		}
		// create identity matrix of "inverse"
//		inverse.assign(identity_matrix<scalar_type>(A.size1()));
		CudaMatrix inRes;
		CudaMatrix::identity_matrix(A.rows,&inRes);
		inverse->assign(&inRes);
		// backsubstitute to get the inverse
		lu_substitute(&A, &pm, inverse);
		return true;	//	NonSingular=true;
	}
	__device__ void     matrix_join(CudaMatrix* mat2,CudaMatrix* res) {
		int row, col;
		row = rows;
		col =cols;
	//	std::cout <<"sfm-directions = "<<col<<" and invariant/mat2.size2() = "<<mat2.size2()<<"\n";
		if (mat2->cols==0){	//second matrix is empty
//			joined_matrix = matrix(row,col,this->data());

			res->assign(this);
			return;
		} else if (col==0){
//			joined_matrix = mat2;
			res->assign(mat2);
			return;
		} else if (col == mat2->cols) {
//			ublas_matrix_impl m(this->size1(), this->size2(), this->data());
			CudaMatrix m;
			m.assign(this);
//			CudaMatrix mat1(m.rows, m.cols, m.data);
			CudaMatrix mat1;
			mat1.assign(&m);
			row = row + mat2->rows; //only row will increase as the col is the system_dimension, so will not change
			CudaMatrix mat2_temp(row, col); //cout << "This is the new Rows = " << mat2_temp.size1();
//			mat1.matrix_copy(mat2_temp);
			mat2_temp.assign(&mat1);
			mat2_temp.resize(row, col, true); //cout << "This is the new Rows after matrix_copy and resize = " << mat2_temp.size1();
			//joined_matrix.resize(row, col);		//mat1.matrix_copy(joined_matrix);		//not working here
			for (int i = mat1.rows, index_i = 0; i < mat2_temp.rows;i++, index_i++) {
				for (int j = 0; j < mat2_temp.cols; j++) {
					mat2_temp.set(i, j, mat2->get(index_i, j));
				}
			}
//			joined_matrix = math::matrix<scalar_type>(mat2_temp.size1(),
//					mat2_temp.size2(), mat2_temp.data());
			res->assign(&mat2_temp);
			return;
		} else {
			printf ("Matices are not of the same dimension:: number of columns do not match!!!");
		}
		//TODO TEST
	}

	/**
	 * Matrix vector multiplication. The result is stored in the passed argument res
	 */
	 __device__ void mult_vector(CudaVector* v,CudaVector * res) const {

		res->deAlocData();
		res->size=rows;
		res->alocData();
//TODO UNIT TEST
		cublasStatus_t stat;
		cublasHandle_t handle;
		stat = cublasCreate(&handle);
		double alpha = 1;

		double beta = 1;

		cublasDgemv(handle, CUBLAS_OP_N, rows, cols, &alpha, data, rows, v->data, 1, &beta, res->data, 1);
		if (CUBLAS_STATUS_SUCCESS != stat) {
			printf("cublas error %d\n", stat);
		}


		cublasDestroy(handle);
	}

	__host__ __device__ static int lu_factorize(CudaMatrix * m,CudaVector* P) {
		/* INPUT: A - array of pointers to rows of a square matrix having dimension N
		 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
		 * OUTPUT: Matrix A is changed, it contains both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
		 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
		 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
		 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
		 */

		int N = m->rows;
		P->resize(N+1);
		if (m->rows != m->cols) {
			printf("matrix not square");
			return 1;
		}


		double **A = new double*[m->rows];
		for (int i = 0; i < N; i++) {
			A[i]= new double[m->cols];
			for (int j = 0; j < N; j++) {
				A[i][j]=m->get(i,j);
			}
		}

		double Tol = 0.2;
		int k, imax;
		double maxA, *ptr, absA;

		for (int i = 0; i <= N; i++)
			P->set(i,i); //Unit permutation matrix, P[N] initialized with N
		int piv;
		for (int i = 0; i < N; i++) {
			maxA = 0.0;
			imax = i;

			for (k = i; k < N; k++)
				if ((absA = fabs(A[k][i])) > maxA) {
					maxA = absA;
					imax = k;
				}

			if (maxA < Tol) {
				for (int i = 0; i < N; i++) {
							delete[] A[i];
				}
				return 1; //failure, matrix is degenerate
			}

			if (imax != i) {
				//pivoting P
				piv = P->get(i);
				P[i] = P[imax];
				P[imax] = piv;

				//pivoting rows of A
				ptr = A[i];
				A[i] = A[imax];
				A[imax] = ptr;

				//counting pivots starting from N (for determinant)
//				P[N]++;
				P->set(N,P->get(N)+1);
			}

			for (int j = i + 1; j < N; j++) {
				A[j][i] /= A[i][i];

				for (k = i + 1; k < N; k++)
					A[j][k] -= A[j][i] * A[i][k];
			}
		}
		//copy to m.data
		for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						m->set(i,j,A[i][j]);
					}
		}
		//free
		for (int i = 0; i < N; i++) {
			delete[] A[i];
		}
		delete[] A;
		return 0;//TODO IS THIS THE CORRECT RETUN VALUE ASK?  //decomposition done

	}
//	__host__ __device__ static int lu_factorize(CudaMatrix m, CudaMatrix pm) {
//
//	}
	__host__ __device__ static void lu_substitute(CudaMatrix * m, CudaVector * P, CudaMatrix * ret) {  //,CudaMatrix*res) {
//TODO ASK IF CORRECT IMPL

		int N=m->rows;
		ret->reinit(N,N);
		for (int j = 0; j < N; j++) {
	        for (int i = 0; i < N; i++) {
	            if (P->get(i) == j)
	                ret->set(i,j,1.0);
	            else
//	                IA[i][j] = 0.0;
	            	ret->set(i,j,0.0);

	            for (int k = 0; k < i; k++)
//	                IA[i][j] -= A[i][k] * IA[k][j];
	            	ret->set(i,j,ret->get(i,j)-(m->get(i,k)*ret->get(k,j)));
	        }

	        for (int i = N - 1; i >= 0; i--) {
	            for (int k = i + 1; k < N; k++)
//	                IA[i][j] -= A[i][k] * IA[k][j];
	            	ret->set(i,j,ret->get(i,j)-(m->get(i,k)*ret->get(k,j)));

//	            IA[i][j] = IA[i][j] / A[i][i];
	            ret->set(i,j,ret->get(i,j)/m->get(i,i));
	        }
	    }

		}
	__host__ __device__ static void prod(CudaMatrix *m1, CudaMatrix* m2,CudaMatrix *res) {

//TODO ist this the same ?
		return times(m1,m2,res);

	}

	__host__ __device__ static void identity_matrix(int n,CudaMatrix* res) {
		CudaMatrix im = CudaMatrix(n, n);
		res->reinit(n,n);

		for (int row = 0; row < res->rows; row++) {
			for (int col = 0; col < res->cols; col++) {
				if (col == row) {
					res->set(row, col, 1);
				}
			}

		}


	}
	__host__ __device__ void reinit(int r,int c){
		deallocateData();
		cols=c;
		rows=r;
		alocData();
	 }
	__host__ __device__ static void plus(CudaMatrix * m1, CudaMatrix *m2,CudaMatrix *res) {
		CudaMatrix c = CudaMatrix(m1->rows, m1->cols);

		res->assign(m1);
		cublasStatus_t stat;
		cublasHandle_t handle;
		stat = cublasCreate(&handle);
		double alpha = 1;
		cublasDaxpy(handle, m2->cols * m2->rows, &alpha, m2->data, 1, res->data, 1);
		if (CUBLAS_STATUS_SUCCESS != stat) {
			printf("cublas error %d\n", stat);
		}

		cublasDestroy(handle);

	}
	__host__ __device__ static void addScalar(CudaMatrix *m, double scalar,CudaMatrix * res) {
		res->reinit(m->rows, m->cols);
		for (int j = 0; j < m->rows * m->cols; j++) {

			res->data[j] = scalar;
		}

		cublasStatus_t stat;
		cublasHandle_t handle;
		stat = cublasCreate(&handle);
		double alpha = 1;
		cublasDaxpy(handle, m->cols * m->rows, &alpha, m->data, 1, res->data, 1);
		if (CUBLAS_STATUS_SUCCESS != stat) {
			printf("cublas error %d\n", stat);
		}

		cublasDestroy(handle);

	}

	__host__ __device__ static void minus(CudaMatrix * m1, CudaMatrix* m2,CudaMatrix*res) {

		res->assign(m1);
		cublasStatus_t stat;
		cublasHandle_t handle;
		stat = cublasCreate(&handle);
		double alpha = -1;
		cublasDaxpy(handle, m2->cols * m2->rows, &alpha, m2->data, 1, res->data, 1);
		if (CUBLAS_STATUS_SUCCESS != stat) {
			printf("cublas error %d\n", stat);
		}

		cublasDestroy(handle);

	}

	__host__ __device__ void deallocateData() {

		if (data != NULL) {
			delete[] data;
			data = NULL;
		}
	}
	__host__ __device__ static void multiplyByScalar(CudaMatrix* m1, double scalar,CudaMatrix *res) {
		res->reinit(m1->rows, m1->cols);
		cublasStatus_t stat;
		cublasHandle_t handle;
		stat = cublasCreate(&handle);
		cublasDaxpy(handle, m1->cols * m1->rows, &scalar, m1->data, 1, res->data, 1);
		if (CUBLAS_STATUS_SUCCESS != stat) {
			printf("cublas error %d\n", stat);
		}

		cublasDestroy(handle);


	}
	__host__ __device__ static void  times(CudaMatrix * m1, CudaMatrix *m2,CudaMatrix * res) {
		res->reinit(m1->rows,m1->cols);//TODO DIFERENT SIZES

		cublasStatus_t stat;
		cublasHandle_t handle;
		stat = cublasCreate(&handle);
		double alpha = 1.0;
		double beta = 0.0;
		//asserts m1->rows==m2->rows
		// res(m,n) = m1(m m1->rows,m1->cols k) * m2(k m2->rows,m2->cols n)
// int lda=m,ldb=k,ldc=m;
		//// C(m,n) = A(m,k) * B(k,n)
//		cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m1->rows, m2->cols, m1->cols, alpha, A,         m1->rows,        B, m2->rows, beta, C, ldc);

		cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m1->rows, m2->cols, m1->cols, &alpha, m1->data, m1->rows, m2->data, m2->rows, &beta, res->data, res->rows);

		if (CUBLAS_STATUS_SUCCESS != stat) {
			printf("cublas error %d\n", stat);
		}


		cublasDestroy(handle);


	}

	__host__ __device__ void minusEquals(CudaMatrix* m) {
		CudaMatrix r;
		CudaMatrix::minus(this, m,&r);
		this->assign(&r);
	}
	/**
	 * THis method internally allocates a new array  and assigns the values according to the matrix m
	 */
	__host__ __device__ void assign(CudaMatrix* m) {
		deallocateData();
		rows = m->rows;
		cols = m->cols;
		alocData();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {


				set(i, j, m->get(i, j));
			}

		}

	}
//TODO TEST
	__host__ __device__ void assign(int prows, int pcols, double* data) {
		deallocateData();
		rows = prows;
		cols = pcols;

		alocData();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				set(i, j, data[IDX2F(i, j, rows)]);
			}

		}

	}
	__host__ __device__ void expm_pad(CudaMatrix *H, double t,CudaMatrix *res) {
		const int p = 6;
		const int n = H->rows;
		CudaMatrix I ;
		CudaMatrix::identity_matrix(n,&I);
		CudaMatrix U(n, n), H2(n, n), P(n, n), Q(n, n);
		double norm = 0.0;
		// Calcuate Pade coefficients
		double c[p + 1];
		c[0] = 1;
		for (int i = 0; i < p; ++i)
			c[i + 1] = c[i] * ((p - i) / ((i + 1.0) * (2.0 * p - i)));
		// Calcuate the infinty norm of H, which is defined as the largest row sum of a matrix
		/*
		 for(size_type i=0; i<n; ++i) {
		 real_value_type temp = 0.0;
		 for(size_type j = 0; j < n; j++)
		 temp += std::abs(H(i, j));
		 norm = t * std::max<real_value_type>(norm, temp);
		 }
		 */
		norm = H->norm_inf();
		// If norm = 0, and all H elements are not NaN or infinity but zero,
		// then U should be identity.
		if (norm == 0.0) {
			bool all_H_are_zero = true;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					if (H->get(i, j) != 0.0)
						all_H_are_zero = false;
			if (all_H_are_zero == true){
//				return I;
				res->assign(&I);
				return;
			}
			// Some error happens, H has elements which are NaN or infinity.
			printf("TODO !!! Null input error in the template expm_pad.\n");
			printf("Null INPUT : "); // << H <<"\n";

		}
		// Scaling, seek s such that || H*2^(-s) || < 1/2, and set scale = 2^(-s)
		int s = 0;
		double scale = 1.0;
		if (norm > 0.5) {
			using std::log;
			s = max(0, static_cast<int>((log(norm) / log(2.0) + 2.0)));
			scale /= pow(2.0, s);
			CudaMatrix r;
			CudaMatrix::multiplyByScalar(H, (scale * t),&r);
			U.assign(&r); // Here U is used as temp value due to that H is const
		} else
			U.assign(H);

		// Horner evaluation of the irreducible fraction, see the following ref above.
		// Initialise P (numerator) and Q (denominator)
		prod(&U, &U,&H2);
//		H2.assign(&h2Prod);
		CudaMatrix::multiplyByScalar(&I, c[p],&Q);
//		Q.assign(&Qres);
		CudaMatrix::multiplyByScalar(&I, c[p - 1],&P);
//		P.assign(&Pres);
		size_type odd = 1;
		for (size_type k = p - 1; k > 0; --k) {
			if(odd == 1) {
					CudaMatrix prodRes;
					CudaMatrix::prod(&Q, &H2,&prodRes);
					CudaMatrix mulScalRes;
					CudaMatrix::multiplyByScalar(&I, c[k - 1],&mulScalRes);
				CudaMatrix::plus(&prodRes, &mulScalRes,&Q);
//					Q .assign(&plusRes);

			}

			else{
					 CudaMatrix prodRes;
					 CudaMatrix::prod(&P, &H2,&prodRes);
					 CudaMatrix mulScalRes;
					 CudaMatrix::multiplyByScalar(&I, c[k - 1],&mulScalRes);
					 CudaMatrix::plus(&prodRes, &mulScalRes,&  P);

			}

			odd = 1 - odd;
		}
		if(odd == 1) {
			(CudaMatrix::prod(&Q, &U,&Q));
		}else {
			CudaMatrix::prod(&P, &U,&P);

//			(P = CudaMatrix::prod(&P, &U));
		}
		Q.minusEquals(&P);
		// In origine expokit package, they use lapack ZGESV to obtain inverse matrix,
		// and in that ZGESV routine, it uses LU decomposition for obtaing inverse matrix.
		// Since in ublas, there is no matrix inversion template, I simply use the build-in
		// LU decompostion package in ublas, and back substitute by myself.

		// Implement Matrix Inversion
		CudaVector pm(n) ;

//		int res = CudaMatrix::lu_factorize(Q, pm);
		int resFac=CudaMatrix::lu_factorize(&Q,&pm);
		if (resFac != 0) {
			printf("Matrix inversion error in the template expm_pad.\n");

		}
		// H2 is not needed anymore, so it is temporary used as identity matrix for substituting.
		H2.assign(&I);
		CudaMatrix::lu_substitute(&Q, &pm, &H2);
		if(odd == 1) {

					CudaMatrix addScal;
					CudaMatrix::addScalar(&I, 2.0,&addScal);
					CudaMatrix prod;
					CudaMatrix::prod(&H2, &P,&prod);
					CudaMatrix times;
					CudaMatrix::times(&addScal, &prod,&times);
					CudaMatrix::multiplyByScalar(&times, -1,&U);
//					U.assign(&res);
				}else{
					CudaMatrix mulByScal;
					CudaMatrix::multiplyByScalar(&I, 2.0,&mulByScal);
					CudaMatrix prod;
					CudaMatrix::prod(&H2, &P,&prod);

					CudaMatrix::times(&mulByScal,&prod,&U);
//					U.assign(&times)
				}
		// Squaring
		for (size_type i = 0; i < (size_type) s; ++i){
			prod(&U, &U,&U);
//			U.assign(&p);
		}

//		return U;
		res->assign(&U);
	}

	__host__ __device__ void alocData() {

		data = new double[rows * cols]();

	}
//destructor
	__host__ __device__ ~CudaMatrix() {
		#if defined(__CUDA_ARCH__)
		 // Device code here
		 deallocateData();
		#else
		 // Host code here
		 // no dealocation on host as the data has been  alocated on the device and used for execution
		#endif

	}

	__host__ __device__ CudaMatrix() {
		rows = 0;
		cols = 0;
	}
	__host__ __device__ CudaMatrix(int r, int c) {
		rows = r;
		cols = c;
		alocData();
	}



	__host__ __device__ void set(int i, int j, double val) {
		int p = IDX2F(i, j, this->rows);
		this->data[p] = val;

	}
	__host__ __device__ double get(int i, int j) {
		return data[IDX2F(i, j, rows)];

	}
	__host__ __device__ void matrix_exponentiation( double time_tau,CudaMatrix *res) {

		CudaMatrix m;
		m.assign(this);
		 expm_pad(&m, time_tau,res);

	}
	__host__ __device__  void transpose(CudaMatrix*res) {
		//TODO JUNIT ?
		int r = this->cols;
		int c = this->rows;
		res->reinit(r, c);
		for (size_type i = 0; i < r; i++) {
			for (size_type j = 0; j < c; j++) {
				res->set(i, j, this->get(j, i));
			}
		}

	}
	/**
	 * returns the infinity norm of the matrix. inf norm of a matrix m is defined as
	 * inf_norm = max(a_{i,j}), 0<=i<r, 0<=j<c
	 */

	int rows; //size 1

	int cols; //size 2
	__host__ __device__ double abs(double x) {
		if (x < 0)
			return -x;
		else
			return x;
	}

	double *data = NULL;
	__host__ __device__ double norm_inf() {
		double norm = 0;
		for (int i = 0; i < this->rows; i++) {
			for (int j = 0; j < this->cols; j++) {
				if (abs(get(i, j)) > norm) {
					norm = abs(this->get(i, j));
				}


			}
		}
		return norm;

	}

};
