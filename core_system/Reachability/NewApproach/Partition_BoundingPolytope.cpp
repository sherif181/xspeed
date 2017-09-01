/*
 * Partition_BoundingPolytope.cpp
 *
 *  Created on: 03-Mar-2015
 *      Author: amit
 */

#include "core_system/Reachability/NewApproach/Partition_BoundingPolytope.h"
/*
 * Returns the list of polytopes: total number of partitions = (p)^nVar
 * unsigned int nVar : number of Variables to be partitioned
 * unsigned int p : number of partitions for each variables(eg 2 for 2 equal halfs or 3 equal halfs,etc
 */
std::list<polytope::ptr> Partition_BoundingPolytope(polytope::ptr S,
		unsigned int nVar, unsigned int p) {

	unsigned int dim = S->getSystemDimension();
	if (nVar > dim) {
		nVar = dim;	//all the variables will be partitioned
	}
	std::vector<std::vector<double> > newDirections;
	math::matrix<double> Real_Directions;
	std::list<polytope::ptr> polys, polys_temp;

	//Axis Directions
	newDirections = generate_axis_directions(dim);
	get_ublas_matrix(newDirections, Real_Directions); //it returns vector vector so need to do conversion here****** Temporary solution
//	std::cout << "Directions = " << Real_Directions;
//	math::matrix <double> coeff(2*dim,dim);

	lp_solver lp(1);	//1 indicate glpk solver used
	lp.setMin_Or_Max(2);
	lp.setConstraints(S->getCoeffMatrix(), S->getColumnVector(), 1);

	std::vector<double> eachDirs(dim);
	std::vector<double> boundvalue(2 * dim);
	for (unsigned int i = 0; i < Real_Directions.size1(); i++) {
		for (unsigned int j = 0; j < Real_Directions.size2(); j++) {
			eachDirs[j] = Real_Directions(i, j);
		}
		boundvalue[i] = S->computeSupportFunction(eachDirs, lp);
	}
	/*for (unsigned int i = 0; i < boundvalue.size(); i++) {
	 cout << boundvalue[i] << endl;
	 }*/

	polytope::ptr s1;//new over-approximated initial set created with bounding as box
	s1 = polytope::ptr(new polytope(Real_Directions, boundvalue, 1));//1 for EqualitySign
	// *************** New Box over-approximation is done **************************

	polys.push_back(s1);
	polys_temp.push_back(s1);
	//Now partition initial polytope s1 into only 2 partitions (later more partitions)
	for (unsigned int i = 0; i < nVar; i++) {
		polys_temp = polys;
		polys.clear();	//remove all items of the list
		for (std::list<polytope::ptr>::iterator it = polys_temp.begin();
				it != polys_temp.end(); it++) {
			std::vector<double> boundvalue1(2 * dim), boundvalue2(2 * dim);
			math::matrix<double> coeff1, coeff2;

			boundvalue1 = (*it)->getColumnVector();	//copying same value
			boundvalue2 = (*it)->getColumnVector();	//copying same value
			coeff1 = (*it)->getCoeffMatrix();	//copying same value = Real_Directions
			coeff2 = (*it)->getCoeffMatrix();	//copying same value

			double dist_apart, h;
			double b = boundvalue[2 * i];	//positive directions
			double a = boundvalue[2 * i + 1] * -1;	//negative directions

			if (a < 0 && b > 0) {	//condition with l+(a<0) and l-(b>0)
				dist_apart = ((-1 * a) + b) / 2.0;		//2 here is p
				h = a + dist_apart;
			} else	//for all other cases eg a>=0 and b>=0
			{
				h = (a + b) / 2.0;			//2 here is p
			}
//using 0 and 1 as indices represent l1(+ive) and l2(-ive):: Later can generalize
			boundvalue1[2 * i] = h;	//positive directions for partition1 <=h
			boundvalue2[2 * i + 1] = -h;//negative directions for partition2 >=h or -partition2 <= -h

			polytope::ptr set1, set2; //new partitioned initial sets
			set1 = polytope::ptr(new polytope(coeff1, boundvalue1, 1)); //1 for EqualitySign
			set2 = polytope::ptr(new polytope(coeff2, boundvalue2, 1)); //1 for EqualitySign
			polys.push_back(set1);
			polys.push_back(set2);

			/*cout << "\ncoeff1 = " << coeff1 << endl;
			 cout << "\nboundvalue1 = " << endl;
			 for (unsigned int i = 0; i < boundvalue1.size(); i++) {
			 cout << boundvalue1[i] << endl;
			 }
			 cout << "\ncoeff2 = " << coeff2 << endl;
			 cout << "\nboundvalue2 = " << endl;
			 for (unsigned int i = 0; i < boundvalue2.size(); i++) {
			 cout << boundvalue2[i] << endl;
			 }*/
		}
	}
	return polys;
}

