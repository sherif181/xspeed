/*
 * ImprovedReachSequentialWithDelete.h
 *
 *  Created on: 05-May-2014
 *      Author: gurung
 */

#ifndef IMPROVEDREACHSEQUENTIALWITHDELETE_H_
#define IMPROVEDREACHSEQUENTIALWITHDELETE_H_

#include "lp_solver.h"
//#include "matrixOperation.h"
#include "DataStructureDirections.h"
#include <fstream>
#include <omp.h>

#define UNSOLVED 0
#define SOLVED 1

using namespace std;

/*
 * coeffMatrix : All facets of the Polytope I
 * bMatrix : All Bound-Value of the Polytope I
 * boundBGreaterThenExpr : All Bound-Sign of the Polytope I  ... 1 for <= and 0 for >=
 *
 * VcoeffMatrix : All facets of the Polytope V
 * VbMatrix : All Bound-Value of the Polytope V
 * VboundBGreaterThenExpr : All Bound-Sign of the Polytope V  ... 1 for <= and 0 for >=
 *
 * AMatrix : is the flow matrix
 * AllDirections : is the new Data Structure with All the list of directions including the facets directions and
 * 				 : other Structure variables such as Row, Column and Flag(0:UNSOLVED & 1:SOLVED)
 * M_Matrix : is the matrix with the optimal reachability result
 *
 */
/*
std::vector<std::vector<double> > ImprovedReachabilitySequentialWithDelete(
		std::vector<std::vector<double> > coeffMatrix,
		std::vector<double> bMatrix, std::vector<int> boundBGreaterThenExpr,
		std::vector<std::vector<double> > VcoeffMatrix,
		std::vector<double> VbMatrix, std::vector<int> VboundBGreaterThenExpr,
		std::vector<std::vector<double> > AMatrix, int iterations,
		std::vector<D> AllDirections,
		std::vector<std::vector<double> > M_Matrix) {
	unsigned int dimension = coeffMatrix[0].size();		// =dimension2;
	//	unsigned int number_of_facets = coeffMatrix.size();	//number of facets of the Initial polytope I
	//unsigned int numVectors = AllDirections.size();		//total_iterated_Directions
	unsigned int number_of_directions = (AllDirections.size() / iterations);

	glpk_lp_solver lpI, lpV;
	lpI.setConstraints(coeffMatrix, bMatrix, boundBGreaterThenExpr);
	lpV.setConstraints(VcoeffMatrix, VbMatrix, VboundBGreaterThenExpr);

	std::vector<vector<double> > M1_for_I;////declaration of the two matrix M1 and M2
	std::vector<vector<double> > M2_for_V;
	M1_for_I.resize(number_of_directions);	//defining the size of the matrix M1
	for (unsigned int i = 0; i < number_of_directions; i++)
		M1_for_I[i].resize(iterations);
	M2_for_V.resize(number_of_directions);	//defining the size of the matrix M2
	for (unsigned int i = 0; i < number_of_directions; i++)
		M2_for_V[i].resize(iterations);

	std::vector<std::vector<double> > rVariable;//declaration of the intermediate vectors
	rVariable.resize(dimension);
	for (unsigned int i = 0; i < dimension; i++) {
		rVariable[i].resize(1);
	}
	unsigned int row, col;
	int INNER_LOOP;		//OUTER_LOOP,
	// While at least one direction UNSOLVED, do loop
	std::ofstream myfile;		//,myfile2;
	myfile.open("/home/amit/myfile.txt");

	//for(OUTER_LOOP=0;OUTER_LOOP<(int)AllDirections.size();OUTER_LOOP++)				//Step 1 all direction is in AllDirections[eachVector].v[0 to dimensions-1]
	int lp_solve_cnt_I = 0, lp_solve_cnt_V = 0; // FOr counting the number of lp solves
	do {
		int solved_I = 0;		//solved_I=1  if computed else 0
		int solved_V = 0;		//solved_V=1  if computed else 0

		for (unsigned int i = 0; i < dimension; i++) {
			rVariable[i][0] = AllDirections[0].v[i];//it is needed in all cases
		}
		row = AllDirections[0].R;
		col = AllDirections[0].C;

		if (AllDirections[0].flag_I == UNSOLVED) {
			AllDirections[0].flag_I = SOLVED;				// Step 3 (a)
			myfile << "Outer Loop = " << "(" << row << "," << col << ")\t";
			lpI.setInitial_SimplexControlParameters();//	this will/may re-initialize the computation process
			M1_for_I[row][col] = lpI.maximize(rVariable);// Step 3 (b)  Actually Computed
			solved_I = 1;
			// incr lp_solve counter here
			lp_solve_cnt_I++;
		}
		if (AllDirections[0].flag_V == UNSOLVED) {
			AllDirections[0].flag_V = SOLVED;				// Step 3 (a)
			myfile << "Outer Loop = " << "(" << row << "," << col << ")\n\n";
			lpV.setInitial_SimplexControlParameters();//	this will/may re-initialize the computation process
			M2_for_V[row][col] = lpV.maximize(rVariable);		// Step 3 (c)
			solved_V = 1;
			// incr lp_solve counter here
			lp_solve_cnt_V++;
		}
		if (AllDirections[0].flag_I && AllDirections[0].flag_V) {
			AllDirections.erase(AllDirections.begin() + 0);	//delete the first position
		}
		myfile << "Inner Loop \n";
		for (INNER_LOOP = 0;
				INNER_LOOP < (int) AllDirections.size()
						&& !AllDirections.empty(); INNER_LOOP++)// check for all unsolved directions if optimal is optimal
				{
			int row1, col1;		//,bothFlag=0;

			row1 = AllDirections[INNER_LOOP].R;
			col1 = AllDirections[INNER_LOOP].C;
			for (unsigned int i = 0; i < dimension; i++) {
				rVariable[i][0] = AllDirections[INNER_LOOP].v[i];
			}
			if (AllDirections[INNER_LOOP].flag_I == UNSOLVED && solved_I) {
				//	lpI.setInitial_SimplexControlParameters();
				lpI.setIteration_Limit(0);	//just testing
				double result1 = lpI.maximize(rVariable);
				if (lpI.getStatus() == 5) {		//5 is GLP_OPT
					AllDirections[INNER_LOOP].flag_I = SOLVED;
					//				bothFlag++;
					M1_for_I[row1][col1] = result1;		//M1_for_I[row][col];
					myfile << "Inner Loop-I= " << INNER_LOOP << "(" << row1
							<< "," << col1 << ")\t";
					//				myfile<<INNER_LOOP<<"\t";
				}
			}
			if (AllDirections[INNER_LOOP].flag_V == UNSOLVED && solved_V) {
				//	lpV.setInitial_SimplexControlParameters();
				lpV.setIteration_Limit(0);	//just testing
				double result2 = lpV.maximize(rVariable);
				if (lpV.getStatus() == 5) {		//5 is GLP_OPT
					AllDirections[INNER_LOOP].flag_V = SOLVED;
					//				bothFlag++;
					M2_for_V[row1][col1] = result2;		//M2_for_V[row][col];
					myfile << "Inner Loop-V= " << INNER_LOOP << "(" << row1
							<< "," << col1 << ")\n";
					//				myfile<<INNER_LOOP<<"\n";
				}
			}
			if (AllDirections[INNER_LOOP].flag_I == SOLVED
					&& AllDirections[INNER_LOOP].flag_V == SOLVED) {
				AllDirections.erase(AllDirections.begin() + INNER_LOOP);//delete the direction
			}
		}	//End of Inner Loop

	} while (!AllDirections.empty());
//	/ * for (unsigned int test=0;test < AllDirections.size();
//			test++
//			)
//			{
//				if (AllDirections[test].flag_I==UNSOLVED && AllDirections[test].flag_V==UNSOLVED)
//				cout<<"\nOh God how is it possible?\n";
//			}
//			* /
//			double sVariable = 0.0;		//initialize s0
			//populating the final matrix M of reachability
			for (unsigned int i = 0; i < number_of_directions; i++)	//rows of the matrix
					{
				sVariable = 0.0;
				for (int j = 0; j < iterations; j++) {	//columns of the matrix
					if (j == 0)	//	first column
							{
						M_Matrix[i][j] = M1_for_I[i][j];
						sVariable = M2_for_V[i][0];
					}

					else {
						// optimize the sum computation
						//sVariable = 0.0;
						//for (int temp=0;temp<=(j-1);temp++)
						M_Matrix[i][j] = M1_for_I[i][j] + sVariable;
						sVariable = sVariable + M2_for_V[i][j];
					}
				}
			}
			//myfile<<"\n\nNumber of times lp solver called = "<<lp_solve_cnt;
			myfile.close();
			return M_Matrix;
		}
*/
#endif /* IMPROVEDREACHSEQUENTIALWITHDELETE_H_ */
