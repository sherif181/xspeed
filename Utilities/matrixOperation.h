/*
 * matrixOperation.h
 *
 *  Created on: 08-Apr-2014
 *      Author: gurung
 */

#ifndef MATRIXOPERATION_H_
#define MATRIXOPERATION_H_

#include <vector>
using namespace std;
/*
 inline std::vector<std::vector<double> > MatrixVectorTranspose(int m, int n,
 std::vector<std::vector<double> > matrixA)
 //,		std::vector<std::vector<double> >matrixATranspose)
 {
 std::vector<std::vector<double> > matrixATranspose;
 matrixATranspose.resize(n);
 for (int i = 0; i < n; i++)
 matrixATranspose[i].resize(m);

 for (int i = 0; i < m; i++)
 for (int j = 0; j < n; j++)
 matrixATranspose[j][i] = matrixA[i][j];

 return matrixATranspose;
 }*/

/*
 inline void MatrixIntegerMultiply(int m, int n, int p, int q, int **matrixA,
 int **matrixB, int **matrixC) {
 if (n == p) {
 //	cout<<"Matrix dimensions are not compatible";
 //	return -1;
 //}
 //else
 //{
 for (int i = 0; i < m; i++) {
 for (int j = 0; j < q; j++) {
 matrixC[i][j] = 0;
 for (int k = 0; k < n; k++) {
 matrixC[i][j] = matrixC[i][j]
 + matrixA[i][k] * matrixB[k][j];
 }
 }
 }
 }
 }

 inline void MatrixDoubleMultiply(int m, int n, int p, int q, double **matrixA,
 double **matrixB, double **matrixC) {
 if (n == p) {
 //	cout<<"Matrix dimensions are not compatible";
 //	return -1;
 //}
 //else
 //{
 for (int i = 0; i < m; i++) {
 for (int j = 0; j < q; j++) {
 matrixC[i][j] = 0;
 for (int k = 0; k < n; k++) {
 matrixC[i][j] = matrixC[i][j]
 + matrixA[i][k] * matrixB[k][j];
 }
 }
 }
 }
 }

 inline void MatrixIntegerTranspose(int m, int n, int **matrixA,
 int **matrixATranspose) {
 for (int i = 0; i < m; i++)
 for (int j = 0; j < n; j++)
 matrixATranspose[j][i] = matrixA[i][j];
 }

 inline void MatrixDoubleTranspose(int m, int n, double **matrixA,
 double **matrixATranspose) {
 for (int i = 0; i < m; i++)
 for (int j = 0; j < n; j++)
 matrixATranspose[j][i] = matrixA[i][j];
 }

 inline void MatrixDoubleCopy(int m, int n, double **matrixSrc,
 double **matrixDest) {
 for (int i = 0; i < m; i++)
 for (int j = 0; j < n; j++)
 matrixDest[i][j] = matrixSrc[i][j];
 }

 template<typename dataT>
 inline std::vector<std::vector<dataT> > MatrixDoubleCopy(
 std::vector<std::vector<dataT> > matrixSrc) {
 int m, n;
 m = matrixSrc.size();
 n = matrixSrc[0].size();
 std::vector<std::vector<dataT> > matrixDest;
 matrixDest.resize(m);
 for (int i = 0; i < m; i++)
 matrixDest[i].resize(n);

 for (int i = 0; i < m; i++)
 for (int j = 0; j < n; j++)
 matrixDest[i][j] = matrixSrc[i][j];
 return matrixDest;
 }
 //template <typename T>
 inline std::vector<std::vector<double> > MatrixVectorMultiply(
 std::vector<vector<double> > matrixA,
 std::vector<vector<double> > matrixB) {
 //I assume that resize for the matrixC have been done in the calling function

 int m, n, p, q;
 m = matrixA.size();
 n = matrixA.at(0).size();

 p = matrixB.size();
 q = matrixB.at(0).size();

 std::vector<std::vector<double> > matrixC;
 matrixC.resize(m);
 for (int i = 0; i < m; i++)
 matrixC[i].resize(q);

 if (n == p) {
 //	cout<<"Matrix dimensions are not compatible";
 //	return -1;	// not compatible
 //}
 //else
 //{
 for (int i = 0; i < m; i++) {
 for (int j = 0; j < q; j++) {
 matrixC[i][j] = 0;
 for (int k = 0; k < n; k++) {
 matrixC[i][j] = matrixC[i][j]
 + matrixA[i][k] * matrixB[k][j];
 }
 }
 }
 }
 return matrixC;
 //	return 1;	// Successful
 }
 //template <typename dataType>
 */

/*

 inline std::vector<double> MatrixAndVectorMultiply(
 std::vector<vector<double> > A, std::vector<double> B) {

 * This function can be used to multiply a Matrix with a Vector to return a Vector
 * Created on: 30-June-2014

 int m, n, p;

 m = A.size();
 n = A.at(0).size();

 p = B.size();

 std::vector<double> C;
 C.resize(m);

 if (n == p) {
 //	cout<<"Matrix and Vector dimensions are not compatible";
 //	return -1;	// not compatible
 //}
 //else
 //{
 for (int i = 0; i < m; i++) {
 C[i] = 0;
 for (int k = 0; k < n; k++) {
 C[i] = C[i] + A[i][k] * B[k];
 }
 }
 }
 return C;
 //	return 1;	// Successful
 }
 */

#endif /* MATRIXOPERATION_H_ */
