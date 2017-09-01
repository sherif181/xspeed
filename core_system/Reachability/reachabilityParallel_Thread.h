/*
 * reachabilityParallel_Thread.h
 *
 *  Created on: 21-May-2014
 *      Author: gurung
 */

#ifndef REACHABILITYPARALLEL_THREAD_H_
#define REACHABILITYPARALLEL_THREAD_H_


/*
#include "lp_solver.h"
//#include "matrixOperation.h"
#include <fstream>

#include <pthread.h>



 using namespace std;


 std::vector<std::vector <double> > reachabilityParallel_PThread(int rows, int dimension1, std::vector<std::vector <double> >coeffMatrix,
 std::vector <double> bMatrix, std::vector<int> boundBGreaterThenExpr, int Vrows, int dimension2,
 std::vector<std::vector <double> >VcoeffMatrix, std::vector <double> VbMatrix,
 std::vector<int> VboundBGreaterThenExpr, std::vector<std::vector <double> > AMatrix,
 int numVectors, std::vector<std::vector <double> > Vector_R, int iterationNum,
 std::vector<std::vector <double> > zValue)
 {
 int dimension = dimension1;		// =dimension2;
 //double OutputMatrix[dimension],VOutputMatrix[dimension];

 / *
 if(pthread_create(&phil[i], 0, (void*(*)(void*))philosopher,(void *)i))
 return -1;

 * /




 #pragma omp parallel for
 for(int tid=0;tid<numVectors;tid++)
 //for(int eachVector=0;eachVector<numVectors;eachVector++)
 {
 double zIInitial=0.0,zI=0.0, zV=0.0;
 double sVariable, s1Variable;
 //int tid = omp_get_thread_num();		//tid is the number of Vectors/Directions

 std::vector<std::vector <double> > r1Variable;
 r1Variable.resize(dimension);
 for (int i=0;i<dimension;i++)
 r1Variable[i].resize(1);

 std::vector<std::vector <double> > rVariable;
 rVariable.resize(dimension);
 for (int i=0;i<dimension;i++){
 rVariable[i].resize(1);
 rVariable[i][0]=Vector_R[tid][i];
 }
 int loopIteration=0;
 sVariable = 0.0; 		//initialize s0
 //zIInitial= LPPResult("SampleEquation1", rows, dimension, coeffMatrix, bMatrix, boundBGreaterThenExpr, rVariable, OutputMatrix);
 glpk_lp_solver lpSingle;
 lpSingle.setConstraints(coeffMatrix,bMatrix,boundBGreaterThenExpr);
 zIInitial= lpSingle.maximize(rVariable);

 zValue[tid][loopIteration]=zIInitial;			//Y0 = pI(r0)
 loopIteration++;
 //just create 2 lp object and initialize with Polytope-I and Polytope-V
 //and inside the loopIteration pass only the new value or rVariable or r1Variable
 //GLP code Begins here
 glpk_lp_solver mylp1;
 glpk_lp_solver mylp2;

 mylp1.setInitial_SimplexControlParameters();
 mylp2.setInitial_SimplexControlParameters();

 mylp1.setConstraints(coeffMatrix, bMatrix, boundBGreaterThenExpr);
 mylp2.setConstraints(VcoeffMatrix, VbMatrix, VboundBGreaterThenExpr);
 //GLP code Ends here
 while(loopIteration < iterationNum){

 std::vector<std::vector <double> > AmatrixTranspose;

 AmatrixTranspose.resize(dimension);	//m=dimension by n=dimension is n=dimension by m=row
 for (int i=0;i<dimension;i++)
 AmatrixTranspose[i].resize(dimension);

 AmatrixTranspose = MatrixVectorTranspose(dimension, dimension, AMatrix);	//, AmatrixTranspose);

 r1Variable= MatrixVectorMultiply(AmatrixTranspose, rVariable);
 zV = mylp2.maximize(rVariable);

 //zV=LPPResult("SampleEquation2", Vrows, dimension, VcoeffMatrix, VbMatrix, VboundBGreaterThenExpr, rVariable, VOutputMatrix);
 s1Variable = sVariable + zV;

 zI = mylp1.maximize(r1Variable);

 //zI= LPPResult("SampleEquation1", rows, dimension, coeffMatrix, bMatrix, boundBGreaterThenExpr, r1Variable, OutputMatrix);

 double TempOmega;
 TempOmega = zI + s1Variable;	//Y1
 zValue[tid][loopIteration]=TempOmega ;	//Y1

 rVariable = MatrixDoubleCopy(r1Variable);	//source to destination

 sVariable = s1Variable;
 loopIteration++;	//for the next Omega-iteration or Time-bound

 }	//end of while for each vector
 mylp1.free_glpk_lp_solver();
 mylp2.free_glpk_lp_solver();
 //glpk_lp_solver::free_environment_glpk_lp_solver();

 //Completion of support function for different number of vectors
 //}	//end of for each Iterations or Time-Bound
 }	//end of pragma omp parallel for
 return zValue;
 }


 */

#endif /* REACHABILITYPARALLEL_THREAD_H_ */
