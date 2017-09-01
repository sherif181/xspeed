#include "Utilities/GPU_utility/gpu_sf_utility.cuh"

double GPU_W_Support(cublasHandle_t handle,
		thrust::device_vector<float> &Btrans, int Btrans_row, int Btrans_col,
		const ReachabilityParameters& ReachParameters,
		Dynamics& system_dynamics, thrust::device_vector<float> &d_direction,
		lp_solver &lp, int Min_Or_Max) {

	int dir_size = d_direction.size(); //this size is same for direction, trans_dir, etc Accessing this will be faster than device.size()
	std::vector<double> direction(dir_size);
	for (int i = 0; i < dir_size; i++)
		direction[i] = d_direction[i]; //copy element device to host

	std::vector<double> trans_dir(dir_size);
	thrust::device_vector<float> d_trans_dir(dir_size);

	double result;
	/*math::matrix<double> B_trans;
	 B_trans = ReachParameters.B_trans;
	 B_trans.mult_vector(direction, trans_dir);*/
	GPU_Multiply_Matrix_Vector(handle, Btrans, Btrans_col,
			Btrans_row, d_direction, 1, dir_size, d_trans_dir); //GPU multiplication
	for (int i = 0; i < dir_size; i++)
		trans_dir[i] = d_trans_dir[i]; //copy element device to host

	double res1 = ReachParameters.time_step
			* system_dynamics.U->computeSupportFunction(trans_dir, lp);
	double beta = ReachParameters.result_beta;
	double res_beta = beta * (double) support_unitball_infnorm(direction); //direction used here
	result = res1 + res_beta;

	return result;
}

/* phi_trans and Btrans are transposed already but row/col are still original to be inverted when calling GPU_Multiply */
double GPU_Omega_Support(cublasHandle_t handle,
		thrust::device_vector<float> &phi_trans, int phi_row, int phi_col,
		thrust::device_vector<float> &Btrans, int Btrans_row, int Btrans_col,
		const ReachabilityParameters& ReachParameters,
		thrust::device_vector<float> &d_direction,
		supportFunctionProvider::ptr Initial_X0, Dynamics& system_dynamics,
		lp_solver &lp, lp_solver &lp_U, int Min_Or_Max) {
	int dir_size = d_direction.size(); //this size is same for direction, trans_dir,r, etc Accessing this will be faster than device.size()

	std::vector<double> direction(dir_size);

	for (int i = 0; i < dir_size; i++)
		direction[i] = d_direction[i]; //copy element device to host

	double res1;
	res1 = Initial_X0->computeSupportFunction(direction, lp); //direction used here

	std::vector<double> trans_dir(dir_size);
	thrust::device_vector<float> d_trans_dir(dir_size);

	//math::matrix<double> A, B, B_trans;
	math::matrix<double> phi_tau_Transpose;
	std::vector<double> r(dir_size);
	thrust::device_vector<float> d_r(dir_size);
	/*phi_tau_Transpose = ReachParameters.phi_trans;
	 phi_tau_Transpose.mult_vector(direction, r);*/

	GPU_Multiply_Matrix_Vector(handle, phi_trans, phi_col, phi_row,
			d_direction, 1, dir_size, d_r); //GPU multiplication
	for (int i = 0; i < dir_size; i++)
		r[i] = d_r[i]; //copy element device to host

	double term1, term2, term3, term3a, term3b, res2;
	term1 = Initial_X0->computeSupportFunction(r, lp); //r used here

	/*B_trans = ReachParameters.B_trans;
	 B_trans.mult_vector(direction, trans_dir);*/
	GPU_Multiply_Matrix_Vector(handle, Btrans, Btrans_col,
			Btrans_row, d_direction, 1, dir_size, d_trans_dir); //GPU multiplication
	for (int i = 0; i < dir_size; i++)
		trans_dir[i] = d_trans_dir[i]; //copy element device to host

	term2 = ReachParameters.time_step
			* system_dynamics.U->computeSupportFunction(trans_dir, lp_U); //trans_dir used here
	term3a = ReachParameters.result_alfa;
	term3b = support_unitball_infnorm(direction);
	term3 = term3a * term3b;
	res2 = term1 + term2 + term3;
	//return (res1 > res2 ? res1:res2);
	if (res1 > res2)
		return res1;
	else
		return res2;
}
