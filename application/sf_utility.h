/*
 * SupportFunction_Utility.h
 *
 *  Created on: 04-Jul-2014
 *      Author: amit
 */

#ifndef SUPPORTFUNCTION_UTILITY_H_
#define SUPPORTFUNCTION_UTILITY_H_

#include "core_system/math/matrix.h"
#include "core_system/continuous/Polytope/Polytope.h"
#include "application/DataStructureDirections.h"
#include "core_system/math/basic_functions.h"
#include <vector>
#include "core_system/continuous/ConvexSet/transMinkPoly.h"
#include "core_system/math/lp_solver/lp_solver.h"

using namespace std;
using namespace math;

/* Change to remove matrix_exponentiation. Pass as Parameter structure */
double Omega_Support(const ReachabilityParameters& ReachParameters,
		std::vector<double> direction, supportFunctionProvider::ptr Initial_X0,
		Dynamics& system_dynamics, lp_solver &lp, lp_solver &lp_U,
		int Min_Or_Max);

/* This function is called by the algorithm which avoids support function computation for the same/similar directions  */
double Omega_Support_Similar_Direction(
		const ReachabilityParameters& ReachParameters,
		std::vector<double> direction, supportFunctionProvider::ptr Initial_X0,
		Dynamics& system_dynamics, lp_solver &lp, lp_solver &lp_U,
		int Min_Or_Max, bool same_dir, Optimize_Omega_Support& optimize_omega);

double W_Support(const ReachabilityParameters& ReachParameters,
		Dynamics& system_dynamics, std::vector<double> direction, lp_solver &lp,
		int Min_Or_Max);

//dot product of vector1 and vector2
template<typename scalar_type>
scalar_type dot_product(std::vector<scalar_type> vector1,
		std::vector<scalar_type> vector2) {
	scalar_type res = 0;
	assert(vector1.size()==vector2.size());
	for (int i = 0; i < vector1.size(); i++) {
		res = res + vector1[i] * vector2[i];
	}
	return res;
}

/**
 * Returns the support function of the unit ball with infinity norm in a passed direction
 */
template<typename scalar_type>
scalar_type support_unitball_infnorm(std::vector<scalar_type> dir) {
	scalar_type sum = 0.0;
	for (unsigned int i = 0; i < dir.size(); i++)
		sum += abs(dir[i]);
	return sum;
}

/**
 * get the dynamics matrix, time step tau and the input polytope V.
 * return the computed beta value to the reachability value.
 * The beta is a value multiplied with the unit ball.
 */

template<typename scalar_type>
scalar_type compute_beta(Dynamics& SysD, scalar_type& tau,
		int lp_solver_type_choosen) {
	scalar_type norm_A = 0.0, result;
	unsigned int dim_for_Max_norm;
	if (!SysD.isEmptyMatrixA){ //if Not Empty
		norm_A = SysD.MatrixA.norm_inf();
	}
	math::matrix<scalar_type> Btrans;
	double V_max_norm = 0.0;

	if (!SysD.isEmptyMatrixB) { //if NOT Empty
		SysD.MatrixB.transpose(Btrans);
		dim_for_Max_norm = SysD.MatrixB.size1();	//dimension for computing Max_Norm(V): V=(B)29x6 . (u)6x1 = (dim of V)29x1
		supportFunctionProvider::ptr Vptr = transMinkPoly::ptr(
				new transMinkPoly(SysD.U, Btrans));
		V_max_norm = Vptr->max_norm(lp_solver_type_choosen, dim_for_Max_norm);
	}

	if (SysD.isEmptyMatrixA){ //if A is Empty
		result = 0;	//norm_A will be zero and which is common term
	}else {
		result = (exp(tau * norm_A) - 1 - tau * norm_A) * (V_max_norm / norm_A);
	}
	//result = (exp(tau * norm_A) - 1 - tau * norm_A) * (V_max_norm / norm_A);
//	cout<<"\nBeta = "<<(double)result<<endl;
	return result;
}

/**
 * get the dynamics matrix, time step tau and the input polytope V and the Initial polytope I(X_o).
 * return the computed beta value to the reachability value.
 * The beta is a value multiplied with the unit ball.
 */

template<typename scalar_type>
scalar_type compute_alfa(scalar_type tau, Dynamics& system_dynamics,
		supportFunctionProvider::ptr I, int lp_solver_type_choosen) // polytope V
		{
	scalar_type norm_A = 0.0, result;
	unsigned int dim_for_Max_norm=0;
	double V_max_norm = 0.0, I_max_norm = 0.0;
	if (!system_dynamics.isEmptyMatrixA){ //if Not Empty
		norm_A = system_dynamics.MatrixA.norm_inf();
	}
	//cout<<"\nInside Testing Matrix A's infinity norm = "<<norm_A<<endl;
	dim_for_Max_norm = I->getSystemDimension();	//I is initial polytope
	I_max_norm = I->max_norm(lp_solver_type_choosen, dim_for_Max_norm); //R_X_o ie max_norm of the Initial polytope
	//cout << "\nInside Testing I.max_norm = " << I_max_norm << endl;

	math::matrix<scalar_type> Btrans;
	if (!system_dynamics.isEmptyMatrixB) { //if NOT Empty
		system_dynamics.MatrixB.transpose(Btrans);
	//	cout <<"test11111\n";
		supportFunctionProvider::ptr Vptr = transMinkPoly::ptr(
				new transMinkPoly(system_dynamics.U, Btrans));
	//	cout <<"test2222222\n";
		dim_for_Max_norm = system_dynamics.MatrixB.size1();	//dimension for computing Max_Norm(V): V=(B)29x6 . (u)6x1 = (dim of V)29x1
	//	cout <<"dim_for_Max_norm = "<<dim_for_Max_norm<<"\n";
		V_max_norm = Vptr->max_norm(lp_solver_type_choosen, dim_for_Max_norm);
	//	cout <<"test33333\n";
	}

	//double V_max_norm = system_dynamics.U->max_norm();	incorrect as V=B.U
//	cout<<"\nInside Testing V_max_norm = "<<V_max_norm <<endl;
	if (system_dynamics.isEmptyMatrixA){ //if A is Empty
		result = 0;	//norm_A will be zero and which is common term
	}else {
//		cout<<"exp(tau * norm_A) = " << exp(tau * norm_A)<<"\n";
		result = (exp(tau * norm_A) - 1 - tau * norm_A) * (I_max_norm + (V_max_norm / norm_A));
	}
//	cout<<"\nAlfa = "<<(double)result<<endl;
	return result;
}

template<typename scalar_type>
void get_ublas_matrix(std::vector<std::vector<scalar_type> > std_Matrix,
		math::matrix<scalar_type>& res) {
	int r = std_Matrix.size();
	int c = std_Matrix[0].size();
	typedef typename math::matrix<scalar_type>::size_type size_type;
	size_type row = r;
	size_type col = c;
	//math::matrix<scalar_type> m = math::matrix<scalar_type>(r, c, std_Matrix.data());
	math::matrix<scalar_type> m(row, col); //, std_Matrix.data());
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			m(i, j) = std_Matrix[i][j];

	res = math::matrix<scalar_type>(m.size1(), m.size2(), m.data());
}

#endif /* SUPPORTFUNCTION_UTILITY_H_ */
