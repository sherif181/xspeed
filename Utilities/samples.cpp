/*
 * samples.cpp
 *
 *  Created on: 03-Apr-2014
 *      Author: gurung
 */
/*

 #include <iostream>
 #include <stdlib.h>
 #include <glpk.h>

 using namespace std;

 double LPPResult(char *problemName, int numConstrains, int nDimension, double **ConstrCoeffMatrix, double *BMatrix, int *boundSign, double *CoeffZ, double *OutputDimensionMatrix)
 {
 glp_prob *lp;
 int ia[1+1000], ja[1+1000];
 double ar[1+1000], z;
 lp = glp_create_prob();
 glp_set_prob_name(lp, problemName);	// eg "sample"
 glp_set_obj_dir(lp, GLP_MAX);

 glp_add_rows(lp, numConstrains);
 for (int i=1;i<= numConstrains;i++)
 {
 glp_set_row_name(lp, i, "p");
 if (boundSign[i]==1)
 glp_set_row_bnds(lp, i, GLP_UP, 0.0, BMatrix[i]);
 else
 glp_set_row_bnds(lp, i, GLP_LO, 0.0, BMatrix[i]);
 }

 glp_add_cols(lp, nDimension);
 for (int i=1;i<= nDimension;i++)
 {
 glp_set_col_name(lp, i, "xi");
 glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0);
 glp_set_obj_coef(lp, i, CoeffZ[i]);
 }

 int count=1;
 for (int i=0;i< numConstrains;i++)
 {
 for (int j=0;j< nDimension;j++)
 {
 ia[count]=i+1, ja[count]=j+1, ar[count]=ConstrCoeffMatrix[i][j];
 count++;
 }
 }
 count--;

 glp_load_matrix(lp, count, ia, ja, ar);
 glp_simplex(lp, NULL);

 z = glp_get_obj_val(lp);
 for (int i=1;i<=nDimension;i++)
 OutputDimensionMatrix[i] = glp_get_col_prim(lp, i);

 glp_delete_prob(lp);
 return z;

 }

 //	double LPPResult(char *problemName, int numConstrains, int nDimension, double **ConstrCoeffMatrix,
 //						double *BMatrix, int *boundSign, double *CoeffZ, double *OutputDimensionMatrix)

 int main(void)
 {
 double z=0.0;
 double  bMatrix[1001], CoeffZ[1001], OutputMatrix[1001];
 int boundBGreaterThenExpr[1001];
 double **coeffMatrix;
 int rows,dimension;
 cout<<"Enter the number of constraints  "<<endl;cin>>rows;
 cout<<"Enter the number of variables/dimensions of the constraints  "<<endl;cin>>dimension;

 coeffMatrix = new  double *[rows];
 for (int i=0;i<rows;i++)
 coeffMatrix[i] = new double [dimension];

 cout<<"Enter the Constraints Matrix "<<endl;
 for (int i=0;i< rows;i++)
 for (int j=0;j< dimension;j++)
 cin>>coeffMatrix[i][j];

 cout<<"Enter the \'B\'- Value "<<endl;
 for (int i=1;i<=rows;i++)
 cin>>bMatrix[i];
 //bMatrix[1]=24;bMatrix[2]=10;bMatrix[3]=0;bMatrix[4]=0;
 cout<<"Enter the Equality-Sign of the Constraints [1 for B-value >= Expression otherwise 0]"<<endl;
 for (int i=1;i<=rows;i++){
 cout<<"# " <<i<<". Constraints Sign = ";
 cin>>boundBGreaterThenExpr[i];
 }
 //boundBGreaterThenExpr[1]=1;boundBGreaterThenExpr[2]=1;		//if 1 then Expr <= b if -1  then Expr >= b
 //boundBGreaterThenExpr[3]=-1;boundBGreaterThenExpr[4]=-1;	// if 1 then GLP_UP   if -1 then GLP_LO

 cout<<"Enter the Coefficient of the Maximizing Expression "<<endl;
 for (int i=1;i<=dimension;i++)
 cin>>CoeffZ[i];
 //CoeffZ[1]=2;CoeffZ[2]=3;
 z = LPPResult("Sample Equation", rows, dimension, coeffMatrix, bMatrix, boundBGreaterThenExpr, CoeffZ, OutputMatrix);

 cout <<"\n\nResult of Maximization z = "<< z<<endl;
 cout <<"Output of the Variables are :"<<endl;
 for (int i=1;i<=dimension;i++)
 {
 cout<<"x"<<i<<" = " << OutputMatrix[i]<<"\t";
 }
 }
 */
