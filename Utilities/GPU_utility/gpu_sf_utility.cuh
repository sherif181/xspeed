/*
 * gpu_sf_utility.cuh
 *
 *  Created on: 14-Apr-2015
 *      Author: amit
 */

#ifndef GPU_SF_UTILITY_CUH_
#define GPU_SF_UTILITY_CUH_

#include <thrust/device_vector.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "application/DataStructureDirections.h"
#include "application/sf_utility.h"		//for support_unitball_infnorm()
#include "core_system/math/GPU/matrix_vector.cuh"

/*double GPU_Omega_Support(const ReachabilityParameters& ReachParameters,
		std::vector<double> direction, supportFunctionProvider::ptr Initial_X0,
		Dynamics& system_dynamics, lp_solver &lp, lp_solver &lp_U,
		int Min_Or_Max);*/
double GPU_Omega_Support(cublasHandle_t handle,
		thrust::device_vector<float> &phi_trans, int phi_row, int phi_col,
		thrust::device_vector<float> &Btrans, int Btrans_row, int Btrans_col,
		const ReachabilityParameters& ReachParameters,
		thrust::device_vector<float> &d_direction,
		supportFunctionProvider::ptr Initial_X0, Dynamics& system_dynamics,
		lp_solver &lp, lp_solver &lp_U, int Min_Or_Max);

double GPU_W_Support(cublasHandle_t handle,
		thrust::device_vector<float> &Btrans, int Btrans_row, int Btrans_col,
		const ReachabilityParameters& ReachParameters,
		Dynamics& system_dynamics, thrust::device_vector<float> &d_direction,
		lp_solver &lp, int Min_Or_Max);


#endif /* GPU_SF_UTILITY_CUH_ */
