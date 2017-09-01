/*
 * AGJHGPU.h
 *
 *  Created on: Nov 30, 2016
 *      Author: tomas
 */

#ifndef AGJHGPU_H_
#define AGJHGPU_H_

#include "reachability.h"
#include "core_system/symbolic_states/symbolic_states.h"
#include "cublas_v2.h"
#include "CudaMatrix.cuh"

__constant__ struct CudaReachabilityParameters {
public:
	/*math::matrix<double> Directions;		//List of Template Directions
	 math::matrix<double> TotalDirections;		//Directions Template + Invariant Directions(excluding negative of invairants)
	 math::matrix<double> phi_trans;
	 math::matrix<double> B_trans;// transpose of the matrix B of Dynamics's polytope : placed here for above reason
	 math::matrix<double> A_inv;		//inverse of the matrix A
	 std::vector<int> Stop_locID;	//List of Locations for which Flowpipe need not be computed
	 double result_alfa;		// placed here for the same above reason
	 double result_beta;		// placed here for the same above reason

	 */
	double TimeBound;
	unsigned int Iterations;
	double time_step; // TimeBound/Iterations

	double A_inv_DATA[128]; //MAXIMUM IS 32KB
	int A_inv_DATA_rows;
	int A_inv_DATA_cols;

//	double phi_trans_DATA[128]; //MAXIMUM IS 32KB
//	int phi_trans_DATA_rows;
//	int phi_trans_DATA_cols;

//	math::matrix<double> Directions;		//List of Template Directions
	double Directions_DATA[128];
	int Directions_rows;
	int Directions_cols;

//	math::matrix<double> TotalDirections;		//Directions Template + Invariant Directions(excluding negative of invairants)
	double TotalDirections_DATA[128];
	int TotalDirections_rows;
	int TotalDirections_cols;

};

class CudaDiscreteSet {
public:
	int discrete_elements[1];

};

class CudaLpSolver {

public:
	int *Sel, *G_Sel;
	float a;

	double *MAT=NULL, *G_MAT, *G_R, *R, *N_MAT=NULL;
	 double *orig_MAT=NULL;
	 double * BoundValue;
	 double * BoundValue_for_X;
	bool single_bound_flag_result;
	CudaMatrix orig_CoefficientMatrix;
	int NB, f, c, No_c;__device__ unsigned int TestConstraints() {
		//TODO use a different solver  the simplex cpu just returns
		unsigned int status = 5;
		return status;
	}
	__device__ void setConstraints(CudaMatrix *m, CudaVector *v, int sign) {
		setConstraints(m, v->data, sign);
	}
	__device__ void setConstraints(CudaMatrix *coeffMatrix_for_X, double* columVector_for_X, int UnitBall) { //New implementation
//		BoundValue_for_X = columVector_for_X;
//
//		bool single_bound_flag_for_X = false;
//		int count = 0;
//
//		for (int i = 0; i < coeffMatrix_for_X.rows; i++) {
//			single_bound_flag_for_X = false;
//			for (int j = 0; j < coeffMatrix_for_X.cols; j++) {
//				if (coeffMatrix_for_X.get(i, j) != 0) {
//					count++; //keeping track of the single variable
//				}
//			}
//			if (count == 1) {
//				single_bound_flag_for_X = true;
//				count = 0;
//			}
//		}
//
//		//single_bound_flag = false;
//		if (single_bound_flag_for_X) {
//			single_bound_flag_result = true;
//		} else {
//			single_bound_flag_result = false;
//			printf("\n CudaSolver->Model does not has Bounded Variable!!!\n"); //Todo::Take care of this condition
//		}

		setConstraints(coeffMatrix_for_X, columVector_for_X);
	}
	__device__ void setMin_Or_Max(int a) {

		//TODO not used in simplex cput ?
	}
	__device__ void setConstraints(CudaMatrix *A, double* B) {
//		unsigned int N_S = number_of_LPs;
		unsigned int N_S = 1;//like in simplex_cpu
		orig_CoefficientMatrix .assign(A);
		BoundValue = B;

		bool single_bound_flag = false;
		int count = 0;

		for (int i = 0; i < A->rows; i++) {
			single_bound_flag = false;
			for (int j = 0; j < A->cols; j++) {
				if (A->get(i, j) != 0) {
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

//			int No_O = A.size2();
			int No_O = A->cols;
//			int No_C = A.size1();
			int No_C = A->rows;
			int M = No_C + 1;
			int N = No_O + 3 + No_C;
			c = 1 + No_O;
			//MAT = (float *) calloc(N_S * M * N, sizeof(float));
			unsigned int memSize = N_S * M * N;
			cudaError_t err;
//			err = cudaMallocHost(&MAT, memSize); //Pinned memory Syntax:: cudaMallocHost(&h_ptr,bytes);
			MAT=new double[memSize]();
//			cudaMemset(MAT, 0, memSize); //initializing all elements to zero
			//std::cout << "Before OMP for MAT!!!" << std::endl;

			for (int s = 0; s < N_S; s++) {
				for (int i = 0; i < M - 1; i++) {
					for (int j = 0; j < N; j++) {
						if (j == 0) {
							MAT[(int) ((i + j * M) + (M * N * s))] = c + i;
						} else if (j > 1) {
							if (j < (No_O + 2)) {
								//Coefficient of A
								MAT[(int) ((i + j * M) + (M * N * s))] = (double) A->get(
										i, j - 2);
							} else if (j == N - 1) {
								MAT[(int) ((i + j * M) + (M * N * s))] =
										(double) B[i];
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
	}	//setting constraints of simplex

	__device__ double compute_LLP(double* C, int size) {
		//TODO frm gimplex
		return 1;

	}

};

class CudaPolytope {
public:
	unsigned int systemDimension;
	bool isEmpty;
	bool isUniverse;
	CudaMatrix coeffMatrix;
	int InEqualitySign = 1;
	int number_facets = 0;
	CudaVector columnVector;
	__device__ void setSystemDimension(unsigned int systemDimension) {
		this->systemDimension = systemDimension;
	}
	__host__ __device__ CudaPolytope() {
		// The default polytope inequality sign is <=
		InEqualitySign = 1;
		number_facets = 0;
		systemDimension = 0;
		// The default polytope inequality sign is <=
		InEqualitySign = 1;
		isUniverse = true; //It is a Universe polytope
		isEmpty = false;
	}
	__device__ __host__ void assign(CudaPolytope* p) {
		this->systemDimension = p->systemDimension;
		this->isEmpty = p->isEmpty;
		this->isUniverse = p->isUniverse;
		this->coeffMatrix.assign(&p->coeffMatrix);
		this->InEqualitySign = p->InEqualitySign;
		this->number_facets = p->number_facets;
		this->columnVector.assign(&p->columnVector);
	}
	__device__ CudaPolytope(CudaMatrix* coeffMatrix, CudaVector* columnVector, int InEqualitySign) {

//		this->setNumberFacets(coeffMatrix.size1());
		number_facets = coeffMatrix->rows;
//		this->setSystemDimension(coeffMatrix.size2());
		setSystemDimension(coeffMatrix->cols);
		this->coeffMatrix.assign(coeffMatrix);
//		this->columnVector.resize(this->number_facets);
		this->columnVector.assign(columnVector);
		this->InEqualitySign = InEqualitySign;
		//Since constraints are set So it is not empty and is Universe but 'bounded'
//		this->setIsEmpty(false); //Not an empty Polytope
//		this->setIsUniverse(false); //Not a universe Polytope and is now a 'Bounded' polytope
		isEmpty = false;
		isUniverse = false;

		//lp.setConstraints(coeffMatrix, columnVector,InEqualitySign);
	}

	__device__ static void post_assign_exact(CudaPolytope * newPolytope, CudaMatrix * R, CudaVector * w, CudaPolytope *ret) {
		CudaMatrix AA, A_dash;
		CudaMatrix R_inverse(R->rows, R->cols);
		CudaVector b_p, b_dash, term2;
		bool invertible;
		//polytope post_assign;
//		b_p = newPolytope->getColumnVector();
		b_p.assign(&newPolytope->columnVector);
		//	cout<<"bp size = "<<b_p.size()<<endl;
		//	cout<<"Testing 1\n";

		invertible = R->inverse(&R_inverse);	//Size of R_inverse has to be assigned otherwise error
		//	cout<<"R_inverse (rows,col) = ("<<R_inverse.size1()<<" , "<<R_inverse.size2()<<")"<<endl;

		if (invertible) {
			//cout<<"\nAmit check a\n";
			/*
			 cout << "\nAmit here\n";
			 for (int i = 0; i < R_inverse.size1(); i++) {
			 for (int j = 0; j < R_inverse.size2(); j++)
			 cout << R_inverse(i, j) << "\t";
			 cout << endl;
			 }
			 cout << "\nAmit there is Inverse Output\n";
			 */

			AA.assign(&newPolytope->coeffMatrix);
			//		cout<<"AA (rows,col) = ("<<AA.size1()<<" , "<<AA.size2()<<")"<<endl;
			/*		cout << "\nAmit here is AA\n";
			 for (int i = 0; i < AA.size1(); i++) {
			 for (int j = 0; j < AA.size2(); j++)
			 cout << AA(i, j) << "\t";
			 cout << endl;
			 }*/
			//AA.multiply(R_inverse, A_dash);

			CudaMatrix::times(&AA, &R_inverse, &A_dash);
			/*		cout << "\nAmit here is A_dash = AA * R_inverse \n";
			 for (int i = 0; i < A_dash.size1(); i++) {
			 for (int j = 0; j < A_dash.size2(); j++)
			 cout << A_dash(i, j) << "\t";
			 cout << endl;
			 }*/

//			A_dash.mult_vector(w, term2);
			A_dash.mult_vector(w, &term2);
			/*cout << "\nAmit here is A' x w \n";
			 for (int j = 0; j < term2.size(); j++)
			 cout << term2[j] << "\t";
			 */
			//		cout<<"Testing 2\n";
			//		cout<<"bp size = "<<b_p.size()<<endl;
			//		cout<<"term2 size = "<<term2.size()<<endl;
//			b_dash = vector_add(b_p, term2);

			CudaVector::vector_add(&b_p, &term2, &b_dash);
			/*
			 cout << "\nAmit here is bp + A' x w \n";
			 for (int j = 0; j < b_dash.size(); j++)
			 cout << b_dash[j] << "\t";
			 */

		} else {
//			std::cout << "\nThe Transition Dynamics Matrix is not Invertible!!!\n";
			//post_assign_approx(newPolytope, R, W,)//TODO WHAT ABOUT THIS ?
		}
		CudaPolytope r;
		CudaPolytope(&A_dash, &b_dash, 1);
		ret->assign(&r);
	}
	__device__ static void post_assign_approx_deterministic(CudaPolytope * newPolytope, CudaMatrix * R, CudaVector * w, CudaMatrix * Directions, CudaPolytope *res) {

		int max_or_min = 2;	//Maximizing
//		std::vector<double> b(Directions.size1()), each_direction(
//				Directions.size2()), direction_trans;
		CudaVector b(Directions->rows);
		CudaVector each_direction(Directions->cols);
		CudaVector direction_trans;
		CudaMatrix R_transpose;
		R->transpose(&R_transpose);
		//create glpk object to be used by the polytope

		CudaLpSolver lp;
		lp.setMin_Or_Max(max_or_min);
		lp.setConstraints(&newPolytope->coeffMatrix, &newPolytope->columnVector, newPolytope->InEqualitySign);
		for (unsigned int i = 0; i < Directions->rows; i++) {
			for (unsigned int j = 0; j < Directions->cols; j++)
//				each_direction[j] = Directions(i, j);
				each_direction.set(j, Directions->get(i, j));
			R_transpose.mult_vector(&each_direction, &direction_trans);

			b.set(i, newPolytope->computeSupportFunction(&direction_trans, &lp) + CudaVector::dot_product(&each_direction, w));
		}
		CudaPolytope r = CudaPolytope(Directions, &b, 1);
		res->assign(&r);
	}

	__device__ void GetPolytope_Intersection(CudaPolytope *P2, CudaPolytope *res) {

		CudaMatrix total_coeffMatrix, m1;
		m1.assign(&coeffMatrix); //assigning constant matrix to matrix m1 so that matrix_join function can be called
		m1.matrix_join(&P2->coeffMatrix, &total_coeffMatrix);
		CudaVector total_columnVector;
		CudaVector::vector_join(&this->columnVector, &P2->columnVector, &total_columnVector);

		CudaPolytope newp(&total_coeffMatrix, &total_columnVector, 1);
		res->assign(&newp);

//		return newp;
	}

	__device__ bool check_polytope_intersection(CudaPolytope * p2) {
		bool flag = false;
		/*
		 * Process: Add all constraints of P1(the calling polytope object) and P2 to form new constraints
		 * then run lp_solver's TestConstraints function to test if the constraints have No Feasible Solution,
		 *`
		 */
		int Min_Or_Max = 2;

		CudaLpSolver lp2;
		//lp2.setDefalultObject();		//executed by default now
		lp2.setMin_Or_Max(Min_Or_Max);

		CudaMatrix m1;
//		m1 = this->getCoeffMatrix(); //assigning constant matrix to matrix m1 so that matrix_join function can be called
		m1.assign(&coeffMatrix);
		CudaMatrix total_coeffMatrix;
		m1.matrix_join(&p2->coeffMatrix, &total_coeffMatrix);

		CudaVector total_columnVector;

		CudaVector::vector_join(&this->columnVector, &p2->columnVector, &total_columnVector);

		/*cout<<"Printing the total ColumnVector = "<<endl;
		 for(unsigned int i=0;i<total_columnVector.size();i++)
		 cout<<total_columnVector[i]<<"\t";
		 cout<<endl;*/
//		lp2.setConstraints(total_coeffMatrix, total_columnVector,
//					this->InEqualitySign);
		lp2.setConstraints(&total_coeffMatrix, &total_columnVector, this->InEqualitySign);

		unsigned int checkStatus = lp2.TestConstraints();
		//Later i may have to perform checking if its an empty polytope should not execute this.
		//cout << "\n\nResult of TestConstriants = " << checkStatus << endl;
		//4 for GLP_NOFEAS; 3 for GLP_INFEAS; 6 for solution is unbounded
		if (checkStatus == 4 || checkStatus == 3 || checkStatus == 6)
			flag = false;
		else
			flag = true;

		return flag;
	}
	__device__ void setMoreConstraints(CudaVector * coeff_constraint, double bound_value) {
//		this->setSystemDimension(coeff_constraint.size());	//or can be obtained from the map_size()
		setSystemDimension(coeff_constraint->size);
//		this->setIsUniverse(false); //Not a Universe Polytope and is now 'Bounded' polytope
		isUniverse = false;
		this->InEqualitySign = 1; // assuming that always less than ineq cons is added.
		// todo: make the impl to accept ineq sign as param or
		// change the func name to add_lt_inequalityCons().

		unsigned int row_size, col_size;
		row_size = coeffMatrix.rows;
		col_size = coeffMatrix.cols; //dimension of the polytope
		if (col_size == 0) // The poly is currently empty
			col_size = coeff_constraint->size;
		else
			assert(col_size == coeff_constraint->size);
		this->coeffMatrix.resize(row_size + 1, col_size, true); //adding one more constraint
		this->columnVector.resize(row_size + 1); //adding one more constraint's bound value
		for (unsigned int i = 0; i < col_size; i++) {
			this->coeffMatrix.set(row_size, i, coeff_constraint->get(i));
		}
		this->columnVector.set(row_size, bound_value);
	}

	__device__ void setMoreConstraints(CudaMatrix * coeff_constraints, CudaVector * bound_values) {
//		this->setSystemDimension(coeff_constraints.size2());
		setSystemDimension(coeff_constraints->cols);
//		this->setIsUniverse(false); //Not a Universe Polytope and is now 'Bounded' polytope
		isUniverse = false;

		unsigned int row_size, dim_size, rows_new;
		row_size = this->coeffMatrix.rows;
		dim_size = this->coeffMatrix.cols; //dimension of the polytope
		rows_new = coeff_constraints->rows;
		//dim_size =coeff_constraints.size2();
		unsigned int new_total_rows = row_size + rows_new;

		this->coeffMatrix.resize(new_total_rows, dim_size, true); //adding more constraints
		this->columnVector.resize(new_total_rows); //adding more constraint's bound values
		for (unsigned int i = 0; i < rows_new; i++) {
			for (unsigned int j = 0; j < dim_size; j++) {
				this->coeffMatrix.set(row_size + i, j, coeff_constraints->get(i, j));
			}
			this->columnVector.set(row_size + i, bound_values->get(i));
		}
	}
	__device__ double computeSupportFunction(CudaVector* direction, CudaLpSolver *lp) {

		assert(direction->size > 0);
		double sf;

		//	std::cout<<"Entered inside ComputeSupportFunction 1 !!\n";
		if (this->isEmpty) {
			//	throw std::runtime_error("\nCompute Support Function called for an Empty Polytope.\n");
			sf = 0; //returns zero for empty polytope
		} else if (this->isUniverse)
			printf("\n Cannot Compute Support Function of a Universe Polytope.\n");
		else {
			//		std::cout<<"Before Compute_LLP !!\n";
			sf = lp->compute_LLP(direction->data, direction->size); //since lp has already been created and set
		}								//with constraints at the time of creation

		return sf;
	}


	__device__ double max_norm(unsigned int dim_for_Max_Norm) {
		unsigned int dimension_size = dim_for_Max_Norm;
		double Max_A, sf, Max = 0.0;
		if (this->isEmpty)
			sf = 0; //returns zero for empty polytope
		else if (this->isUniverse)
			printf("\n ERROR Universe Unbounded Polytope!!!\n");
		else {
			//sf = lp.Compute_LLP(direction);	//since lp has already been created and set with constraints at the time of creation
			//std::vector < std::vector<double> > generator_directions; //this vector-vector is used only in this function not affecting the rest of the codes
			double **generator_directions = new double*[dimension_size * 2];
			for (int i = 0; i < dimension_size * 2; i++) {
				generator_directions[i] = new double[dimension_size];
			}

			//Generator for Positive Directions for example Right and Up
			for (unsigned int i = 0; i < dimension_size; i++) {
//							std::vector<double> directions(dimension_size, 0.0);
				for (unsigned int j = 0; j < dimension_size; j++) {
					generator_directions[i][j] = 0.0;
				}
				generator_directions[i][i] = 1; //Positive Generators

			}
			//Generator for Negative Directions for example Left and Down
			for (unsigned int i = 0; i < dimension_size; i++) {
				for (unsigned int j = 0; j < dimension_size; j++) {
					generator_directions[i + dimension_size][j] = 0.0;
				}
				generator_directions[i][i] = -1;
			}

			CudaLpSolver solver;
//						lp1.setMin_Or_Max(2); //Setting GLP_MAX
			CudaMatrix mat;
			mat.assign(&this->coeffMatrix);
			solver.setConstraints(&mat, this->columnVector.data);
			//Finding the maximum of all Direction : Returns the max element
			for (unsigned int i = 0; i < dimension_size * 2; i++) {
				double* each_generator;
				each_generator = generator_directions[i];
				//cout<<"Each Generator = (" << each_generator[0]<<" , "<<each_generator[1]<<") "<<endl;
				sf = solver.compute_LLP(each_generator, dimension_size);
				delete[] generator_directions[i];
				Max_A = (abs(sf));
				if (Max_A > Max)
					Max = Max_A;
			}
			delete[] generator_directions;
		}
		return Max;

	}
};

class CudaPolyhedra {
public:
	__device__ __host__ void assign(CudaPolyhedra *from) {
		Matrix_SupportFunction.assign(&from->Matrix_SupportFunction); //Note if it has invariants_dirs then Matrix_SupportFunction will also have bound_value
		template_Directions.assign(&from->template_Directions);	//only the template directions

		/* This  (Matrix_InvariantBound,invariant_Directions) will be replaced with my structure idea later */
		Matrix_InvariantBound.assign(&from->Matrix_InvariantBound);	//Note now Matrix_SupportFunction will NOT have the bound_value
		//	math::matrix<double> All_Directions;//Number of Rows or facets of the Convex Set/Polytoe including the invariants(template_dirs + invariant_dirs)
		invariant_Directions.assign(&from->invariant_Directions);	//only the invariant directions

		//unsigned int total_Directions;	//Number of rows or the number of faces/constraints of the Convex Set/Polytoe Omega's
		total_template_Directions = (from->total_template_Directions);	//total number of template directions
		total_invariant_Directions = (from->total_template_Directions);	//total number of invariant directions
		total_iterations = (from->total_template_Directions);	//Number of Columns or the number of Convex Set/Polytope Omega's

	}
	__device__ __host__ CudaPolyhedra(CudaMatrix * matrix_support_function, CudaMatrix * template_directions) {
		//	this->setTotalIterations(matrix_support_function.size2());
		this->setMatrixSupportFunction(matrix_support_function);
		//	this->setAllDirections(all_directions);
		this->setTemplateDirections(template_directions); //no more working with all_directions
		//	this->setTotalTemplateDirections(all_directions.size1());
		this->setTotalInvariantDirections(0);
	}

	__device__ void getPolytope(unsigned int Iterations_Number, CudaPolytope *res) {
		CudaPolytope p;

		//p->setCoeffMatrix(this->getTemplateDirections());//no more working with all_directions
		p.coeffMatrix.assign(&template_Directions);
		//	std::vector<double> boundValue(this->getTotalTemplateDirections());	//no more working with all_directions
		CudaVector vec(this->getTotalTemplateDirections());
		p.columnVector.assign( &vec);
		for (unsigned int i = 0; i < Matrix_SupportFunction.rows; i++) {
			p.columnVector.set(i, Matrix_SupportFunction.get(i, Iterations_Number));
		}

//		p->setInEqualitySign(1);	//Ax<=b		0 for = sign, 1 for <= sign and 2 for >= sign
		p.InEqualitySign = 1;
		res->assign(&p);
//		return p;
	}

	__device__ void polys_intersectionSequential_optimize(CudaPolytope * G, CudaVectorTemplate<CudaPair<unsigned int> >*res) {
		size_type row = 0;
		size_type col = 0;
		CudaMatrix mat_sf(row, col);
		bool is_intersected = false, intersection_started = false, immediate_false = false, intersection_ended = false;
		int foundIntersection = 0, intersection_start;
		//	std::list<template_polyhedra::ptr> intersected_region;

//		std::vector<bool> intersects(this->Matrix_SupportFunction.cols, false); //all false initially
		//bool intersects[this->Matrix_SupportFunction.cols] = { 0 };//
		CudaVectorTemplate<bool> intersects(this->Matrix_SupportFunction.cols, false);
		for (unsigned int i = 0; i < this->Matrix_SupportFunction.cols; i++) {
			//std::cout<<"\n Inner thread Template_polyhedra omp_get_num_threads() = "<< omp_get_num_threads()<<"\n";
			CudaPolytope p;
			this->getPolytope(i, &p);
			//CudaVector constraint_bound_values(this->getInvariantDirections().rows);
			CudaVector constraint_bound_values;
			this->getInvariantBoundValue(i, &constraint_bound_values);
//			p->setMoreConstraints(this->getInvariantDirections(), constraint_bound_values);
			CudaMatrix invDir;
			invDir.assign(this->getInvariantDirections());
			p.setMoreConstraints(&invDir, &constraint_bound_values);
			intersects.set(i, p.check_polytope_intersection(G)); //result of intersection
		} //end of parallel-loop :: we have the list of intersected polys

//		std::list<std::pair<unsigned int, unsigned int> > intersected_range;

//		CudaVectorTemplate<CudaPair <unsigned int> >intersected_range(this->Matrix_SupportFunction.cols+1);
		res->size = (this->Matrix_SupportFunction.cols + 1);
		res->deAlocData();
		res->alocData();

		//	cout << "Is This Big = " << this->Matrix_SupportFunction.size2() << "\n";
		int count = 0;
//		std::pair<unsigned int, unsigned int> inte_range;
		CudaPair<unsigned int> inte_range;
		for (unsigned int i = 0; i < this->Matrix_SupportFunction.cols; i++) { //sequential reading of an boolean_array that tells intersected polys
			is_intersected = intersects.get(i);

			if (is_intersected == true) { //if intersects create a new template_polyhedra subset
				intersection_started = true;
				if (foundIntersection == 0) { //Intersection start
					foundIntersection++; //intersection_start = i;	//Intersection started at this position of 'i'
					inte_range.data[0] = i; //	cout << "\nIntersection Found at = " << i << endl;
				}
			}
			if (intersection_started == true && is_intersected == false) {
				intersection_ended = true;
				inte_range.data[1] = i;
				foundIntersection = 0;
			}
			if (intersection_ended == true && intersection_started == true) {
//				intersected_range.push_back(inte_range);
				res->set(count, inte_range);
				count++;
				intersection_started = false;
				intersection_ended = false; //cout << "\nIntersection Ended at = " << i << "\n";
			}
		} //end of parallel for-loop

		if (intersection_started == true && intersection_ended == false) {
			inte_range.data[1] = this->Matrix_SupportFunction.cols - 1;
			//intersected_range.push_back(inte_range);  //	cout << "\nIntersection did not End = " << i2 << "\n";
			res->set(count, inte_range);
		}
		//intersected_range REMOVE UNUISED positions TODO  coount<size
//		return intersected_range;
	}

	__device__ void flowpipe_intersectionSequential(CudaPolytope * guard, CudaVectorTemplate<CudaPolytope*> * res) {

		CudaVectorTemplate<CudaPair<unsigned int> > range_list;
//		range_list =
		polys_intersectionSequential_optimize(guard, &range_list);
		//cout <<"range_list.size = "<<range_list.size();
//		std::list<polytope::ptr> polys;

//		CudaVectorTemplate<CudaPolytope>polys(range_list.size);
		res->deAlocData();
		res->size = range_list.size;
		res->alocData();

		unsigned int poly_dir_size = this->template_Directions.rows + this->invariant_Directions.rows;
		CudaVector colVector(poly_dir_size);
//		for (std::list<std::pair<unsigned int, unsigned int> >::iterator range_it = range_list.begin(); range_it != range_list.end(); range_it++) {
		for (int i = 0; i < range_list.size; i++) {

//			unsigned int start = (*range_it).first;
			unsigned int start = range_list.get(i).data[0];
//			unsigned int end = (*range_it).second;
			unsigned int end = range_list.get(i).data[1];
			//cout << "first = " << start << "  second = " << end << std::endl;
//			for (unsigned int eachTemplateDir = 0; eachTemplateDir < this->template_Directions.size1(); eachTemplateDir++) {
			for (unsigned int eachTemplateDir = 0; eachTemplateDir < this->template_Directions.rows; eachTemplateDir++) {
				double Max_sf = this->Matrix_SupportFunction.get(eachTemplateDir, start);
				for (int i = start + 1; i <= end; i++) {
					double sf = this->Matrix_SupportFunction.get(eachTemplateDir, i);
					if (sf > Max_sf)
						Max_sf = sf;
				}  //end of each intersected region
//				colVector[eachTemplateDir] = Max_sf;
				colVector.set(eachTemplateDir, Max_sf);
			}  //end of each template direction ALSO HAVE TO PERFORM INVARIANT DIRECTION
			unsigned int total_dir = this->template_Directions.rows;
//			for (unsigned int eachInvDir = 0; eachInvDir < this->invariant_Directions.size1(); eachInvDir++) {
			for (unsigned int eachInvDir = 0; eachInvDir < this->invariant_Directions.rows; eachInvDir++) {
				double Max_sf = this->Matrix_InvariantBound.get(eachInvDir, start);
				for (int i = start + 1; i <= end; i++) {
					double sf = this->Matrix_InvariantBound.get(eachInvDir, i);
					if (sf > Max_sf)
						Max_sf = sf;
				} //end of each intersected region
//				colVector[total_dir + eachInvDir] = Max_sf;
				colVector.set(total_dir + eachInvDir, Max_sf);
			}
			CudaMatrix allDirs;
			this->template_Directions.matrix_join(&this->invariant_Directions, &allDirs);
//			polytope::ptr p = polytope::ptr(new polytope(allDirs, colVector, 1));
			CudaPolytope* p =new  CudaPolytope(&allDirs, &colVector, 1);
//			polys.push_back(p);			this->template_Directions.matrix_join(&this->invariant_Directions,&);

			res->set(i, p);//TODO THIS CAUSES A MEMORY LEEK IF THERE IS NOC DESTRUCTOR FOR CUDAVECTOR TEMPLATE
		} //end of multiple intersected region
		  //cout<<"polys.size = "<<polys.size()<<"\n";
//		return polys;
	}

	__device__ __host__ CudaPolyhedra() {
		total_iterations = 0;
		total_template_Directions = 0;
		total_invariant_Directions = 0;
	}
	__device__ __host__ void setMatrixSupportFunction(CudaMatrix * matrixSupportFunction) {
		//	cout<<"Called\n";
		Matrix_SupportFunction.assign(matrixSupportFunction);
		//cout<<"Called 2\n";
		this->setTotalIterations(matrixSupportFunction->cols);
	}

	__device__ __host__ unsigned int getTotalIterations() const {
		return total_iterations;
	}

	__device__ __host__ void setTotalIterations(unsigned int totalIterations) {
		total_iterations = totalIterations;
	}

	__device__ __host__ void setMatrix_InvariantBound(CudaMatrix * matrix_invariantBound) {
		Matrix_InvariantBound.assign(matrix_invariantBound);
	}
	__device__ __host__ unsigned int getTotalTemplateDirections() const {
		return total_template_Directions;
	}
	__device__ __host__ void setTotalTemplateDirections(unsigned int total_template_directions) {
		total_template_Directions = total_template_directions;
	}

	__device__ __host__ unsigned int getTotalInvariantDirections() const {
		return total_invariant_Directions;
	}
	__device__ __host__ void setTotalInvariantDirections(unsigned int total_invariant_directions) {
		total_invariant_Directions = total_invariant_directions;
	}

	__device__ __host__ CudaMatrix* getInvariantDirections()  {
		return &invariant_Directions;
	}

	__device__ void getInvariantBoundValue(int Iterations_Number, CudaVector *res) {	//need testing
		res->deAlocData();
		res->size = (this->getInvariantDirections()->rows);
		res->alocData();
		for (unsigned int i = 0; i < this->invariant_Directions.rows; i++) {
//			bound_value[i] = this->Matrix_InvariantBound(i, Iterations_Number);
			res->set(i, this->Matrix_InvariantBound.get(i, Iterations_Number));
		}
//		return bound_value;
	}

	__device__ __host__ void setTemplateDirections(CudaMatrix * template_Directions) {
		this->template_Directions.assign(template_Directions);
		this->setTotalTemplateDirections(template_Directions->rows);
	}

	__device__ __host__ void setInvariantDirections(CudaMatrix * invariant_Directions) {
		this->invariant_Directions.assign(invariant_Directions);
		this->setTotalInvariantDirections(invariant_Directions->rows);
	}
	__device__ __host__ CudaPolyhedra(CudaMatrix * matrix_support_function, CudaMatrix * matrix_invariant_bounds, CudaMatrix * template_directions, CudaMatrix * invariant_directions) {
		//this->setTotalIterations(matrix_support_function.size2());
		this->setMatrixSupportFunction(matrix_support_function);

		//this->setTotalTemplateDirections(template_directions.size1());
		this->setMatrix_InvariantBound(matrix_invariant_bounds);

		this->setTemplateDirections(template_directions);
		this->setInvariantDirections(invariant_directions);
	}

private:
	CudaMatrix Matrix_SupportFunction; //Note if it has invariants_dirs then Matrix_SupportFunction will also have bound_value
	CudaMatrix template_Directions;	//only the template directions

	/* This  (Matrix_InvariantBound,invariant_Directions) will be replaced with my structure idea later */
	CudaMatrix Matrix_InvariantBound;	//Note now Matrix_SupportFunction will NOT have the bound_value
	//	math::matrix<double> All_Directions;//Number of Rows or facets of the Convex Set/Polytoe including the invariants(template_dirs + invariant_dirs)
	CudaMatrix invariant_Directions;	//only the invariant directions

	//unsigned int total_Directions;	//Number of rows or the number of faces/constraints of the Convex Set/Polytoe Omega's
	unsigned int total_template_Directions;	//total number of template directions
	unsigned int total_invariant_Directions;	//total number of invariant directions
	unsigned int total_iterations;	//Number of Columns or the number of Convex Set/Polytope Omega's
};

class CudaDynamics {
public:
	bool isEmptyMatrixA;	//True if empty otherwise False
	CudaMatrix matrixA;
	bool isEmptyMatrixB;	//True if empty otherwise False
	CudaMatrix matrixB;
	CudaPolytope U;
	CudaVector C;
	bool isEmptyC;
	bool calculatedbTrans = false;
	bool calculatedphiTrans = false;
	CudaMatrix B_trans;
	CudaMatrix phi_trans;

	__device__ CudaMatrix calculateBtrans() {
		// transpose to be done once and stored
		if (!calculatedbTrans) {
			calculatedbTrans = true;
			CudaMatrix transRes;
			matrixB.transpose(&transRes);
			B_trans.assign(&transRes);
		}
	}
};
class CudaAssign {
public:
	CudaMatrix Map;
	CudaVector b;

	__device__ __host__  void assign(CudaAssign *a) {
		Map.assign(&a->Map);
		b.assign(&a->b);
	}
};
class CudaLocation;
class CudaTransition {

public:
	int trans_id;
//	string label;TODO LABELS FOR TRANSITIONS
	CudaLocation * source_location;
	CudaLocation * destination_location;
	CudaPolytope Gaurd;
	CudaAssign assignT;__device__ __host__ CudaTransition() {

	}
	//__device__ __host__ CudaTransition(int trans_id, char label, int source_id, int destination_id, polytope::ptr Gaurd, Assign& Assign_T);

//	__device__ __host__ CudaTransition(int transition_id, char label_name, int source_id, int dest_loc_id, CudaPolytope gaurd) {
//		trans_id = transition_id;
////		label = label_name; TODO LABEL
//		source_location_id = source_id;
//		destination_location_id = dest_loc_id;
//		Gaurd = gaurd;
//		//Assign_T = assign_Trans;
//	}

//	__device__ __host__ void setDestinationLocation( int dest_loc_id){
//		destination_location_id = dest_loc_id;
//	}
	__device__ __host__ CudaAssign getAssignT() {
		return assignT;
	}

	__device__ __host__ void setAssignT(CudaAssign * assingt) {
		assignT.assign(assingt);
	}

	__device__ __host__ CudaPolytope* getGaurd() {
		return &Gaurd;
	}

	__device__ __host__ void setGaurd(CudaPolytope * gaurd) {
		Gaurd.assign(gaurd);
	}

//	__device__ __host__ int getSource_Location_Id() {
//		return source_location_id;
//	}

//	__device__ __host__ void setSource_Location_Id(int source_loc_id) {
//		source_location_id = source_loc_id;
//	}
	__device__ __host__ int getTransitionId() const {
		return trans_id;
	}

	__device__ __host__ void setTransitionId(int transId) {
		trans_id = transId;
	}
};

class CudaTransMinkPoly: public CudaPolytope {
public:
	/** Constructor to represent convex sets which is a linear transformation of another convex set only: C' = Trans. C */
	CudaPolytope X0;
	CudaPolytope U;
	CudaMatrix TRANS;
	CudaVector C;
	CudaMatrix B_TRANS;
	double time;
	double beta;
	bool Cempty;__device__ __host__ void assign(CudaTransMinkPoly* p) {
		this->systemDimension = p->systemDimension;
		this->isEmpty = p->isEmpty;
		this->isUniverse = p->isUniverse;
		this->coeffMatrix.assign(&p->coeffMatrix);
		this->InEqualitySign = p->InEqualitySign;
		this->number_facets = p->number_facets;
		this->columnVector.assign(&p->columnVector);
		X0.assign(&p->X0);
		U.assign(&p->U);
		TRANS.assign(&p->TRANS);
		C.assign(&p->C);
		B_TRANS.assign(&p->B_TRANS);
		time = p->time;
		beta = p->beta;
		Cempty = p->Cempty;
	}
	__host__ __device__ CudaTransMinkPoly(CudaPolytope * myX0, CudaPolytope * myU, CudaVector * c, CudaMatrix * myTRANS, CudaMatrix * myB_TRANS, double mytime, double mybeta) {
		X0.assign(myX0);
		U.assign(myU);
		TRANS.assign(myTRANS);
		B_TRANS.assign(myB_TRANS);
		time = mytime;
		beta = mybeta;
		Cempty = false;	//C is NOT empty here as it is supplied
		C.assign(c);
	}
	__device__ CudaTransMinkPoly(CudaPolytope * myX0, CudaPolytope * myU, CudaMatrix * myTRANS, CudaMatrix * myB_TRANS, double mytime, double mybeta) {
		X0.assign(myX0);
		U.assign(myU);
		TRANS.assign(myTRANS);
		B_TRANS.assign(myB_TRANS);
		time = mytime;
		beta = mybeta;
		Cempty = true;	//C is empty here as not supplied
	}
	__device__ CudaTransMinkPoly(CudaPolytope * myX0, CudaMatrix * myTRANS) {
		X0.assign(myX0);
		TRANS.assign(myTRANS);
		time = 0;
		beta = 0;
		Cempty = true;	//C is empty here as not supplied
	}
	__device__ CudaTransMinkPoly() {

		time = 0;
		beta = 0;
		Cempty = true;	//C is empty here as not supplied
	}

	__device__ double computeSupportFunction(CudaVector * direction, CudaLpSolver * lp) {
		//this function is also called from compute_beta, compute_alfa, etc
		CudaVector dprime;
		TRANS.mult_vector(direction, &dprime);
		//	cout << "\nCalling transMinkPoly ComputerSupportFunction\n";

		double res1 = 0;
		if (!X0.isEmpty) {
			res1 = X0.computeSupportFunction(&dprime, lp);
		}
		//	cout << "\t res1 = " << res1;
		double res2 = 0.0;
		if (!U.isEmpty) {
			CudaLpSolver lp_U;
			lp_U.setConstraints(&U.coeffMatrix, &U.columnVector, U.InEqualitySign);

			B_TRANS.mult_vector(direction, &dprime);
			res2 = U.computeSupportFunction(&dprime, &lp_U);
		}
		//	cout << "\t  res2 = " << res2;
		double res3 = 0.0;	//for C
		if (!Cempty) {
			//cout<<"B_Trans = "<<B_TRANS<<std::endl;

			B_TRANS.mult_vector(direction, &dprime);
			res3 = dot_product(&dprime, &C);
		}
		//cout<<"\t res3 = "<<res3<<std::endl;
		double res = res1 + time * res2 + time * res3;
		if (beta != 0) {
			double dir_norm = support_unitball_infnorm(direction);
			return res + beta * dir_norm;
		} else
			return res;
	}

	__device__ double dot_product(CudaVector * vector1, CudaVector * vector2) {
		double res = 0;
		assert(vector1->size == vector2->size);
		for (int i = 0; i < vector1->size; i++) {
			res = res + vector1->get(i) * vector2->get(i);
		}
		return res;
	}
	__device__ double support_unitball_infnorm(CudaVector * dir) {
		double sum = 0.0;
		for (unsigned int i = 0; i < dir->size; i++)
			sum += abs(dir->get(i));
		return sum;
	}
};

class CudaLocation {
public:
	static int const BAD=0;
	static int const GOOD=1;
	static int const FINAL=2;
	static int const  UNSAVE=3;


	CudaDynamics system_dynamics;
	bool invariantExists;
	CudaPolytope invariant;
	CudaTransition* Out_Going_Transitions;
	int Out_Going_Transitions_size = 0;
	int name=-1;
	int id;
};

class CudaState {
public:
	CudaPolytope initialSet;
	CudaTransition * transition;
	CudaLocation* location;
	CudaPolyhedra continuousSet;
	CudaDiscreteSet discreteSet;
	int location_id;
	int transition_id;__host__ __device__ CudaState() {

	}
	__host__ __device__ CudaState(int locationId, CudaPolytope * initialSet) {
		location_id = locationId;
		this->initialSet.assign(initialSet);
	}


};
class AGJHGPU: public reachability {
public:
	AGJHGPU();
	virtual ~AGJHGPU();
	std::list<symbolic_states::ptr> gpuReachability(std::list<abstractCE::ptr>& ce_candidates);
	void initReachabilityParameters();
	void loadLocation(CudaLocation *cLocation, location* loc);

	void copyDataToGpu();
	CudaState convertStateToCuda(initial_state::ptr state);
	CudaState *d_queue;
	CudaState *d_hashTable;

	CudaTransition *d_transitions;
	CudaLocation *d_location;
};
/*

 struct cSet{


 }

 struct dSet{


 }


 struct state {
 dSet discreteSet;
 cSet continiousSet;
 state parent;
 int transitionId;
 polytope::ptr initialContiniousSet
 }


 */

#endif /* AGJHGPU_H_ */
