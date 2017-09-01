#include "core_system/math/Bulk_LP_Solver/bulk_LP_Solver.h"
#include "application/CopyArray.h"


bulk_lp_solver::bulk_lp_solver(int lp_solver_type) {

	//lp_problem.lp_solver(lp_solver_type);
	lp_problem = lp_solver::lp_solver_ptr(new lp_solver(lp_solver_type));
	solver_type = lp_solver_type;
	this->BoundValue.resize(0);
	this->Max_Or_Min = 2;	//by default Maximum
}

std::vector<float>& bulk_lp_solver::getResultAll() {
	return result;
}
//get the status of all simplex
std::vector<int> bulk_lp_solver::getStatusAll() {
	;
}

void bulk_lp_solver::setMaxMin(unsigned int Max_Or_Min) {
	//1 for Min and 2 for Max
	this->Max_Or_Min = Max_Or_Min;
	lp_problem->setMin_Or_Max(Max_Or_Min);
}
void bulk_lp_solver::setConstratint(math::matrix<double> CoefficientMatrix,
		std::vector<double> boundValue) {
	lp_problem->setConstraints(CoefficientMatrix, boundValue, 1);
	this->CoefficientMatrix = CoefficientMatrix;
	this->BoundValue.resize(boundValue.size());
	this->BoundValue = CopyVector(boundValue);

}
void bulk_lp_solver::ComputeLP(math::matrix<float>& List_of_ObjValue) {	//16-July-2015: Modified to pass be reference
	unsigned int numLPs = List_of_ObjValue.size1();
	result.resize(numLPs);
	/* Complete new lp_solver object to avoid race condition
	 *
	 #pragma omp parallel for
	 for (unsigned int i = 0; i < numLPs; i++) {
	 lp_solver lp(solver_type);
	 //std::cout<<"solver_type = " <<solver_type<<"\n";
	 lp.setConstraints(this->CoefficientMatrix,this->BoundValue,1);
	 std::vector<double> dir(List_of_ObjValue.size2());
	 for (unsigned int x = 0; x < List_of_ObjValue.size2(); x++)
	 dir[x] = (double) List_of_ObjValue(i, x);
	 result[i] = (float) lp.Compute_LLP(dir);
	 }*/

	for (unsigned int i = 0; i < numLPs; i++) {
		std::vector<double> dir(List_of_ObjValue.size2());
		for (unsigned int x = 0; x < List_of_ObjValue.size2(); x++)
			dir[x] = (double) List_of_ObjValue(i, x);//Looks like this operation is expensive
		result[i] = (float) lp_problem->Compute_LLP(dir);
	}
	/*
	 * This will not work as lp_problem's(glpk's object does not exists in the thread memory)
	 #pragma omp parallel for
	 for (unsigned int i = 0; i < numLPs; i++) {
	 std::vector<double> dir(List_of_ObjValue.size2());
	 for (unsigned int x = 0; x < List_of_ObjValue.size2(); x++)
	 dir[x] = (double) List_of_ObjValue(i, x);
	 result[i] = (float) lp_problem->Compute_LLP(dir);
	 }
	 */
}

void bulk_lp_solver::ComputeLP_ListVector(
		std::list<std::vector<double> >& List_of_ObjValue) {//16-July-2015: Modified to pass be reference
	unsigned int numLPs = List_of_ObjValue.size();
	result.resize(numLPs);
//	std::cout<<"\nnumLPs = "<<numLPs<<std::endl;
//	std::cout<<"\nList_of_ObjValue.size = "<<List_of_ObjValue.size()<<std::endl;
	int i = 0;
	std::vector<double> dir;
	for (std::list<std::vector<double> >::iterator it =
			List_of_ObjValue.begin(); it != List_of_ObjValue.end(); it++) {
		//i = List_of_ObjValue.begin();
		dir = *it;
		result[i] = (float) lp_problem->Compute_LLP(dir);
		i++;
	}

	/*	for (unsigned int i = 0; i < numLPs; i++) {
	 std::vector<double> dir(List_of_ObjValue.size2());
	 for (unsigned int x = 0; x < List_of_ObjValue.size2(); x++)
	 dir[x] = (double) List_of_ObjValue(i, x);		//Looks like this operation is expensive
	 result[i] = (float) lp_problem->Compute_LLP(dir);
	 }*/
}

