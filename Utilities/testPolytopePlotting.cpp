/*
 * testPolytopePlotting.cpp
 *
 *  Created on: 17-Nov-2014
 *      Author: amit
 */

#include "testPolytopePlotting.h"
#include <fstream>



void GeneratePolytopePlotter(polytope::ptr poly){
	/*
	 * Generate files for plotting in matlab
	 */

	std::ofstream MatLabfile, MatLabfile2;

	MatLabfile.open("/home/amit/matlabTest/ProjectOutput/Dirs_A_printing.txt");
	MatLabfile2.open("/home/amit/matlabTest/ProjectOutput/polytope_b_printing.txt");
math::matrix<double> A;
std::vector<double> b;
A = poly->getCoeffMatrix();
b = poly->getColumnVector();

		for (int i = 0; i < A.size1(); i++)
		{
			for (int j = 0; j < A.size2(); j++)
			{
				MatLabfile << A(i, j) << " ";
			}
			MatLabfile << std::endl;
			MatLabfile2 <<b[i]<<std::endl;
		}

	MatLabfile.close();
	MatLabfile2.close();
}

/*

void GenerateInitialPolytopePlotter(std::list<polytope::ptr> initial_polys_list){

	 * Generate files for plotting in matlab from a list of Polytopes


	std::ofstream MatLabfile, MatLabfile2;

	MatLabfile.open("/home/amit/matlabTest/ProjectOutput/Dirs_List_Initail_Poly.txt");
	MatLabfile2.open("/home/amit/matlabTest/ProjectOutput/b_List_Initail_Poly.txt");
math::matrix<double> A;
std::vector<double> b;



for ()
for (std::list<polytope::ptr>::iterator it=initial_polys_list.begin(); it!=initial_polys_list.end(); it++ )
{

	 b	 = (*it)->getColumnVector();
		MatLabfile2 <<b[i]<<" ";
	}
A= (*it)->getCoeffMatrix();
}


//Only for Direction
for (int i = 0; i < A.size1(); i++)
{
	for (int j = 0; j < A.size2(); j++)
		{
			MatLabfile << A(i, j) << " ";
		}
	MatLabfile << std::endl;
}
MatLabfile.close();
*/
