/*
 * lp_gurobi_simplex.h
 *
 *  Created on: 03-Sep-2014
 *      Author: amit
 */

#ifndef LP_GUROBI_SIMPLEX_H_
#define LP_GUROBI_SIMPLEX_H_



#include <iostream>
#include <gurobi_c++.h>
#include "core_system/math/matrix.h"
#include <boost/shared_ptr.hpp>

class gurobi_lp_solver {
public:
	typedef boost::shared_ptr<gurobi_lp_solver> gurobi_ptr;

	//static GRBEnv env;
	void setConstraints(math::matrix<double> coeff_constraints, std::vector<double> bounds, int bound_signs);

	double Compute_LPP(std::vector<double> coeff_function);


//	void setObjective(std::vector<double> coeff_function);
//	double Solve_LP();


	void setIteration_Limit(int limits);		// :: Need to test the function
	void setMin_Or_Max(int Min_Or_Max);		//will be required at the time of Optimization unlike GLPK at the time of setting constraints
	int getMin_Or_Max();
	void setInitial_SimplexControlParameters();		//reset to default parameters :: Need to test the function
	unsigned int getStatus();					// :: Need to test the function
	unsigned int TestConstraints();	//set the Maximizing Objective to zero and compute simplex :: Need to test the function

	gurobi_lp_solver();
	~gurobi_lp_solver();
	static int gurobi_lp_count;
private:

	GRBEnv *env;
	GRBModel *model;		// = GRBModel(env);
	GRBVar  *variables;
	unsigned int number_of_variables;
	unsigned int number_of_constraints;
	int Minimize_Or_Maximize;
	std::vector<double> Maximizing_Variables;// values of the variable that maximizes the result

};




#endif /* LP_GUROBI_SIMPLEX_H_ */
