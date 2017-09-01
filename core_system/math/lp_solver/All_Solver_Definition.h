/*
 * All_Solver_Definition.h
 *
 *  Created on: 11-Feb-2015
 *      Author: amit
 */

#ifndef ALL_SOLVER_DEFINITION_H_
#define ALL_SOLVER_DEFINITION_H_

#define GLPK_SOLVER			1	//GLPK SOLVER Algorithm
#define GUROBI_SOLVER		2	//GUROBI Algorithm
#define PAR_SIMPLEX	 		3	//Parallel SIMPLEX Implementation for Multi-Core CPUs
#define GIMPLEX	 			4	//Parallel SIMPLEX Implementation for GPUs

#define SIMPLEX_CPU_SOLVER		5	//Parallel SIMPLEX Implementation for GPUs


	/*
	 * Meaning					Common_Retun_Code 	GLPK_Code	Gurobi_Code
	 * Solution is Undefined			1				1
	 * Solution is feasible				2				2
	 * Solution is Infeasible			3				3			3
	 * no feasible Solution Exists		4				4
	 * Solution is Optimal				5				5			2
	 * Solution is unbounded			6				6			5
	 */




#endif /* ALL_SOLVER_DEFINITION_H_ */
