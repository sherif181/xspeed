/*
 * bulk_LP_Solver.h
 *
 *  Created on: 19-Apr-2015
 *      Author: amit
 */

#ifndef BULK_LP_SOLVER_H_
#define BULK_LP_SOLVER_H_

#include "core_system/math/lp_solver/lp_solver.h"
#include "core_system/math/matrix.h"
#include <vector>
#include <list>

class bulk_lp_solver {
private:
	lp_solver::lp_solver_ptr  lp_problem;
	std::vector<float> result;
	int solver_type;
	math::matrix<double> CoefficientMatrix;
	std::vector<double> BoundValue;
	unsigned int Max_Or_Min;
public:

	bulk_lp_solver(int lp_solver_type);

	//get the result of all simplex
	std::vector<float>& getResultAll();
	//get the status of all simplex
	std::vector<int> getStatusAll();
	void setConstratint(math::matrix<double> CoefficientMatrix,
			std::vector<double> boundValue);
	void setMaxMin(unsigned int Max_Or_Min);	//1 for Min and 2 for Max
	void ComputeLP(math::matrix<float>& List_of_ObjValue);

	//Only needed for comparison
	void ComputeLP_ListVector(std::list<std::vector<double> >& List_of_ObjValue);

};

#endif /* BULK_LP_SOLVER_H_ */
