#include <thrust/device_vector.h>
#include "core_system/math/GPU/matrix_vector.cuh"
#include "Utilities/GPU_utility/gpu_sf_utility.cuh"
#include "Utilities/Template_Polyhedra.h"
#include "Utilities/invariantBoundaryCheck.h"
#include "core_system/math/matrix.h"


//Mixing CPU with GPU for Matrix-Vector Multiplication


typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;

template_polyhedra::ptr reachabilitySequential_GPU_MatrixVector_Multiply(
		Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen) {

	int numVectors = ReachParameters.Directions.size1();
	int dimension = Initial->getSystemDimension();
	unsigned int shm_NewTotalIteration = ReachParameters.Iterations; //Shared Variable for resize iterations number on crossing with invariant
	int Min_Or_Max = 2;

	int phi_row, phi_col;
	phi_row = ReachParameters.phi_trans.size1(); //original row
	phi_col = ReachParameters.phi_trans.size2(); //original col
	thrust::device_vector<float> phi_trans(phi_row * phi_col);
	 convert_thrust_matrix(ReachParameters.phi_trans, phi_trans);

	int Btrans_row, Btrans_col;
	Btrans_row = ReachParameters.B_trans.size1(); //original row
	Btrans_col = ReachParameters.B_trans.size2(); //original col
	thrust::device_vector<float> Btrans(Btrans_row * Btrans_col);
	convert_thrust_matrix(ReachParameters.B_trans, Btrans);

	math::matrix<double> MatrixValue; //Shared Matrix for all child thread
	size_type row = numVectors, col = shm_NewTotalIteration;
	if (isInvariantExist == true) { //if invariant exist. Computing
	//todo :: fix this later
		//shm_NewTotalIteration = InvariantBoundaryCheck(SystemDynamics, Initial, ReachParameters, invariant, lp_solver_type_choosen);
	} //End of Invariant Directions

	if (shm_NewTotalIteration == 1) {
		template_polyhedra::ptr poly_emptyp;
		return poly_emptyp;
	}

	col = shm_NewTotalIteration; //if invariant exist col will be resized
	MatrixValue.resize(row, col);
	int solver_type = lp_solver_type_choosen;
	lp_solver s_per_thread_I(solver_type), s_per_thread_U(solver_type),
			s_per_thread_inv(solver_type);

	s_per_thread_I.setMin_Or_Max(2);
	if (!ReachParameters.X0->getIsEmpty()) //set glpk constraints If not an empty polytope
		s_per_thread_I.setConstraints(ReachParameters.X0->getCoeffMatrix(),
				ReachParameters.X0->getColumnVector(),
				ReachParameters.X0->getInEqualitySign());

	s_per_thread_U.setMin_Or_Max(2);
	if (SystemDynamics.U->getIsEmpty()) { //empty polytope
		//Polytope is empty so no glpk object constraints to be set
	} else {
		s_per_thread_U.setConstraints(SystemDynamics.U->getCoeffMatrix(),
				SystemDynamics.U->getColumnVector(),
				SystemDynamics.U->getInEqualitySign());
	}

	std::vector<double> r1Variable; //now single dimension
	r1Variable.resize(dimension);
	std::vector<double> rVariable;
	rVariable.resize(dimension);
	thrust::device_vector<float> d_rVariable(rVariable.size());
	thrust::device_vector<float> d_r1Variable(r1Variable.size());
	//thrust::device_vector<float> C(rowB * colA);

	cublasHandle_t handle;
	cublasCreate(&handle);

	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {

		double zIInitial = 0.0, zI = 0.0, zV = 0.0;
		double sVariable, s1Variable;
		for (int i = 0; i < dimension; i++) {
			rVariable[i] = ReachParameters.Directions(eachDirection, i);
			d_rVariable[i] = rVariable[i]; //copy from  host to device
		}
		unsigned int loopIteration = 0;
		sVariable = 0.0; //initialize s0

		zIInitial = GPU_Omega_Support(handle, phi_trans, phi_row, phi_col,
				Btrans, Btrans_row, Btrans_col, ReachParameters, d_rVariable,
				Initial, SystemDynamics, s_per_thread_I, s_per_thread_U,
				Min_Or_Max);
		//cout<<"Above Only Initial output\n";
		MatrixValue(eachDirection, loopIteration) = zIInitial;
		loopIteration++;
		for (; loopIteration < shm_NewTotalIteration;) { //Now stopping condition is only "shm_NewTotalIteration"
			double TempOmega;

			//use GPU for multiplication
			//d_rVariable.assign(rVariable.begin(), rVariable.end());		 // problem due to double-to-float

			/*for (int i = 0; i < rVariable.size(); i++)
			 d_rVariable[i] = rVariable[i];*/

			/*std::cout << "\nd_rVariable ";
			 for (int i = 0; i < d_rVariable.size(); i++)
			 std::cout << d_rVariable[i] << " ";*/

			GPU_Multiply_Matrix_Vector(handle, phi_trans,
					phi_col, phi_row, d_rVariable, 1, rVariable.size(), d_r1Variable);

			/*d_r1Variable = GPU_Multiply_Matrix_Vector(phi_trans, phi_col,
			 phi_row, d_rVariable, 1, rVariable.size());*/
			/*std::cout<<"\nd_r1Variable = ";
			 for (int i = 0; i < d_r1Variable.size(); i++)
			 std::cout << d_r1Variable[i] << " ";*/

			/*for (int i = 0; i < d_r1Variable.size(); i++)
			 r1Variable[i] = d_r1Variable[i];*/

			/*thrust::copy(d_r1Variable.begin(), d_r1Variable.end(),
			 r1Variable.begin()); //coping from device to r1Variable*/ //problem due to double-to-float
			//ReachParameters.phi_trans.mult_vector(rVariable, r1Variable);
			/*zV = W_Support(ReachParameters, SystemDynamics, rVariable,
			 s_per_thread_U, Min_Or_Max);*/

			zV = GPU_W_Support(handle, Btrans, Btrans_row, Btrans_col,
					ReachParameters, SystemDynamics, d_rVariable,
					s_per_thread_U, Min_Or_Max);
			s1Variable = sVariable + zV;

			/*zI = Omega_Support(ReachParameters, r1Variable, Initial,
			 SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);*/
			zI = GPU_Omega_Support(handle, phi_trans, phi_row, phi_col, Btrans,
					Btrans_row, Btrans_col, ReachParameters, d_r1Variable,
					Initial, SystemDynamics, s_per_thread_I, s_per_thread_U,
					Min_Or_Max);

			TempOmega = zI + s1Variable; //Y1
			MatrixValue(eachDirection, loopIteration) = TempOmega; //Y1

			//rVariable = CopyVector(r1Variable); //source to destination
			for (int i = 0; i < dimension; i++)
				d_rVariable[i] = d_r1Variable[i];	//device swapping between r and r1 variable

			sVariable = s1Variable;
			loopIteration++; //for the next Omega-iteration or Time-bound
		} //end of while for each vector

		//cout<<endl;
	}
	// Destroy the handle
	cublasDestroy(handle);

	if (isInvariantExist == true) { //if invariant exist. Computing
		math::matrix<double> inv_sfm;
		//ReachParameters.TotalDirections = ReachParameters.Directions;
		//int lastDirs = ReachParameters.TotalDirections.size1() - 1;	//total number of rows (0 to rows-1)
		int num_inv = invariant->getColumnVector().size(); //number of Invariant's constriants
		inv_sfm.resize(num_inv, shm_NewTotalIteration);

		//int newDirectionSize = lastDirs + 1 + num_inv;
		//ReachParameters.TotalDirections.resize(newDirectionSize, dimension,true);
		//MatrixValue.resize(newDirectionSize, shm_NewTotalIteration, true);//Matrix resized
		for (int eachInvDirection = 0; eachInvDirection < num_inv;
				eachInvDirection++) {
			//for (int i = 0; i < dimension; i++) {
			//	ReachParameters.TotalDirections(lastDirs + eachInvDirection + 1,
			//			i) = invariant->getCoeffMatrix()(eachInvDirection, i);
			//}
			for (unsigned int i = 0; i < shm_NewTotalIteration; i++) {
				//MatrixValue(lastDirs + eachInvDirection + 1, i) = invariant->getColumnVector()[eachInvDirection];
				inv_sfm(eachInvDirection, i) =
						invariant->getColumnVector()[eachInvDirection];
			}
		}
		//cout<<"\nAmit"<<MatrixValue.size2()<<"\n";
		return template_polyhedra::ptr( new template_polyhedra(MatrixValue, inv_sfm,
				ReachParameters.Directions, invariant->getCoeffMatrix()));
		//return template_polyhedra(MatrixValue, ReachParameters.TotalDirections,ReachParameters.Directions,invariant->getCoeffMatrix());
		//return template_polyhedra(MatrixValue, ReachParameters.TotalDirections);
	} else {
		return template_polyhedra::ptr( new  template_polyhedra(MatrixValue, ReachParameters.Directions));
	}
}

