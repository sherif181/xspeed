/*
 * Template_Polytope.cpp
 *
 *  Created on: 10-Jul-2014
 *      Author: amit
 */

#include "Utilities/Template_Polyhedra.h"
#include "core_system/continuous/Polytope/Polytope.h"
#include <omp.h>
#include "Utilities/testPolytopePlotting.h"
#include "boost/timer/timer.hpp"

typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;

template_polyhedra::template_polyhedra() {
	total_iterations = 0;
	total_template_Directions = 0;
	total_invariant_Directions = 0;
}
template_polyhedra::template_polyhedra(
		math::matrix<double> matrix_support_function,
		math::matrix<double> template_directions) {
//	this->setTotalIterations(matrix_support_function.size2());
	this->setMatrixSupportFunction(matrix_support_function);
//	this->setAllDirections(all_directions);
	this->setTemplateDirections(template_directions); //no more working with all_directions
//	this->setTotalTemplateDirections(all_directions.size1());
	this->setTotalInvariantDirections(0);
}

/*template_polyhedra::template_polyhedra(
 math::matrix<double> matrix_support_function,
 math::matrix<double> all_directions,
 math::matrix<double> template_directions,
 math::matrix<double> invariant_directions) {
 this->setTotalIterations(matrix_support_function.size2());
 this->setMatrixSupportFunction(matrix_support_function);
 this->setAllDirections(all_directions);
 this->setInvariantDirections(invariant_directions);
 this->setTemplateDirections(template_directions);
 }*/

template_polyhedra::template_polyhedra(
		math::matrix<double> matrix_support_function,
		math::matrix<double> matrix_invariant_bounds,
		math::matrix<double> template_directions,
		math::matrix<double> invariant_directions) {
	//this->setTotalIterations(matrix_support_function.size2());
	this->setMatrixSupportFunction(matrix_support_function);

	//this->setTotalTemplateDirections(template_directions.size1());
	this->setMatrix_InvariantBound(matrix_invariant_bounds);

	this->setTemplateDirections(template_directions);
	this->setInvariantDirections(invariant_directions);
}

const math::matrix<double>& template_polyhedra::getTemplateDirections() const {
	return template_Directions;
}
void template_polyhedra::setTemplateDirections(
		math::matrix<double>& template_Directions) {
	this->template_Directions = template_Directions;
	this->setTotalTemplateDirections(template_Directions.size1());
}
const math::matrix<double>& template_polyhedra::getInvariantDirections() const {
	return invariant_Directions;
}
void template_polyhedra::setInvariantDirections(
		math::matrix<double>& invariant_Directions) {
	this->invariant_Directions = invariant_Directions;
	this->setTotalInvariantDirections(invariant_Directions.size1());
}

const math::matrix<double>& template_polyhedra::getMatrix_InvariantBound() const {
	return Matrix_InvariantBound;
}
void template_polyhedra::setMatrix_InvariantBound(
		math::matrix<double>& matrix_invariantBound) {
	Matrix_InvariantBound = matrix_invariantBound;
}

/*
 unsigned int template_polyhedra::getTotalDirections() const {
 return total_Directions;
 }
 void template_polyhedra::setTotalDirections(unsigned int totalDirections) {
 total_Directions = totalDirections;
 }
 */
unsigned int template_polyhedra::getTotalTemplateDirections() const {
	return total_template_Directions;
}
void template_polyhedra::setTotalTemplateDirections(
		unsigned int total_template_directions) {
	total_template_Directions = total_template_directions;
}

unsigned int template_polyhedra::getTotalInvariantDirections() const {
	return total_invariant_Directions;
}
void template_polyhedra::setTotalInvariantDirections(
		unsigned int total_invariant_directions) {
	total_invariant_Directions = total_invariant_directions;
}

std::vector<double> template_polyhedra::getInvariantBoundValue(
		int Iterations_Number) {	//need testing
	std::vector<double> bound_value(this->getInvariantDirections().size1());

	for (unsigned int i = 0; i < this->invariant_Directions.size1(); i++) {
		bound_value[i] = this->Matrix_InvariantBound(i, Iterations_Number);
	}
	return bound_value;
}

polytope::ptr template_polyhedra::getPolytope(unsigned int Iterations_Number) {
	polytope::ptr p;
	p = polytope::ptr(new polytope());
	p->setCoeffMatrix(this->getTemplateDirections());//no more working with all_directions
	std::vector<double> boundValue(this->getTotalTemplateDirections());	//no more working with all_directions
	for (unsigned int i = 0; i < Matrix_SupportFunction.size1(); i++) {
		boundValue[i] = Matrix_SupportFunction(i, Iterations_Number);
	}
	p->setColumnVector(boundValue);
	p->setInEqualitySign(1);	//Ax<=b		0 for = sign, 1 for <= sign and 2 for >= sign
	return p;
}

/*
 polytope template_polyhedra::getIntersectedRegion(polytope gaurd){

 }
 */

const std::list<template_polyhedra::ptr> template_polyhedra::polys_intersectionSequential(polytope::ptr G, int lp_solver_type_choosen) { //need testing due to modification
	size_type row = 0;
	size_type col = 0;
//	cout << "1 = Testing INSIDE \n";
	math::matrix<double> mat_sf(row, col);
	bool is_intersected = false, intersection_started = false, immediate_false =
			false, intersection_ended = false;
	int foundIntersection = 0, intersection_start;
	std::list<template_polyhedra::ptr> intersected_region;
	for (unsigned int i = 0; i < this->Matrix_SupportFunction.size2(); i++) {
		//		cout<<"3 INSIDE .....LOOP no = "<<i<<endl;
		polytope::ptr p;
		p = this->getPolytope(i);
		std::vector<double> constraint_bound_values(this->getInvariantDirections().size1());
		constraint_bound_values = this->getInvariantBoundValue(i);
		/* from the reach_region that is the calling template_polyhedra add the invariant constraint and bound_value to p*/
		p->setMoreConstraints(this->getInvariantDirections(), constraint_bound_values);
		is_intersected = p->check_polytope_intersection(G, lp_solver_type_choosen);	//result of intersection

		if (is_intersected == true) {//if intersects create a new template_polyhedra subset
//#pragma omp critical
				//		{
			intersection_started = true;
			/*if (foundIntersection == 0) {	//Intersection start
			 foundIntersection++;	//intersection_start = i;	//Intersection started at this position of 'i'
			 }*/
			col = col + 1;//		cout<<"5 ... INSIDE .... INTERSECTED .....\n";
			mat_sf.resize(p->getCoeffMatrix().size1(), col, true);//resizing to append the polytope p of the flowpipe that intersects
			//for (unsigned int j = 0; j < this->getTotalTemplateDirections();j++)
			for (unsigned int j = 0; j < p->getCoeffMatrix().size1(); j++) {
				mat_sf(j, col - 1) = p->getColumnVector()[j];//	mat_sf.addColumn(p.getColumnVector()); add column in matrix class
			}
			//		}
			//	cout << "\nIntersection Found at = " << i << endl;
		}

		if (is_intersected == false && intersection_started == true) {
//#pragma omp critical
			//	{
			intersection_ended = true;
			//	}
		}

		if (intersection_started == true && intersection_ended == true) {
//#pragma omp critical
			//		{
			intersected_region.push_back(template_polyhedra::ptr(new template_polyhedra(mat_sf, p->getCoeffMatrix())));
			mat_sf.clear();	//reset all values to zero
			mat_sf.resize(0, 0);
			col = 0;
			intersection_started = false;
			intersection_ended = false;
			//		} 		//cout << "\nIntersection Ended at = " << i << "\n";
		}
		//	i2 = i;
	}	//end of parallel for-loop
	if (intersection_started == true && intersection_ended == false) {
		math::matrix<double> directions, templateDirs, invDirs;
		templateDirs = this->getTemplateDirections();
		invDirs = this->getInvariantDirections();
		templateDirs.matrix_join(invDirs, directions);

		intersected_region.push_back(template_polyhedra::ptr(new template_polyhedra(mat_sf, directions)));//intersected_region.push_back(template_polyhedra(mat_sf, p.getCoeffMatrix()));
		//	cout << "\nIntersection did not End = " << i2 << "\n";
	}

	/*	if (foundIntersection == 0)
	 cout << "\nSorry Intersection Not Found!!!!!!\n";
	 else
	 cout << "\nIntersection Found at = " << intersection_start << endl;
	 */

	/*
	 * After checking all elements(Omegas) of the template_polyhedra of the reachable set
	 * the function returns the point_of_intersection_started with the guard G.
	 * And only returns the interescted region with the guard G as a new template_polyhedra:::with added invariant_directions
	 */
//	point_of_intersection = intersection_start;
	return intersected_region;
}

const std::list<template_polyhedra::ptr> template_polyhedra::polys_intersectionParallel(
		polytope::ptr G, int lp_solver_type_choosen) { //need testing due to modification

	size_type row = 0;
	size_type col = 0;
	math::matrix<double> mat_sf(row, col);
	bool is_intersected = false, intersection_started = false, immediate_false =false, intersection_ended = false;
	int foundIntersection = 0, intersection_start;
	std::list<template_polyhedra::ptr> intersected_region;
	//intersected_region = template_polyhedra::ptr(new template_polyhedra());

	std::vector<bool> intersects(this->Matrix_SupportFunction.size2(), false);//all false initially
//	omp_set_dynamic(0);	//handles dynamic adjustment of the number of threads within a team
	//omp_set_nested(3);	//enable nested parallelism
	//omp_set_num_threads(10);
	//omp_set_max_active_levels(3);

	boost::timer::cpu_timer t100;
	t100.start();

#pragma omp parallel for //num_threads(4)
	for (unsigned int i = 0; i < this->Matrix_SupportFunction.size2(); i++) {
		//std::cout<<"\n Inner thread Template_polyhedra omp_get_num_threads() = "<< omp_get_num_threads()<<"\n";
		polytope::ptr p;
		p = this->getPolytope(i);
		std::vector<double> constraint_bound_values(
				this->getInvariantDirections().size1());
		constraint_bound_values = this->getInvariantBoundValue(i);
		p->setMoreConstraints(this->getInvariantDirections(),
				constraint_bound_values);

		intersects[i] = p->check_polytope_intersection(G, lp_solver_type_choosen);	//result of intersection
	} //end of parallel-loop :: we have the list of intersected polys

	t100.stop();
	double clock100;
	clock100 = t100.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
	double return100 = clock100 / (double) 1000;
	std::cout << "\nOnly Parallel Intersection Time:Wall(Seconds) = " << return100 << std::endl;

	boost::timer::cpu_timer t101;
	t101.start();
	polytope::ptr p;
	cout<<"Is This Big = "<<this->Matrix_SupportFunction.size2()<<"\n";
	for (unsigned int i = 0; i < this->Matrix_SupportFunction.size2(); i++) {//sequential reading of an boolean_array that tells intersected polys
		is_intersected = intersects[i];
		if (is_intersected == true) {//if intersects create a new template_polyhedra subset
			intersection_started = true;
			/*if (foundIntersection == 0) {	//Intersection start
			 foundIntersection++;	//intersection_start = i;	//Intersection started at this position of 'i'
			 }*/
			p = this->getPolytope(i);
			std::vector<double> constraint_bound_values(this->getInvariantDirections().size1());
			constraint_bound_values = this->getInvariantBoundValue(i);
			p->setMoreConstraints(this->getInvariantDirections(), constraint_bound_values);

			col = col + 1;//		cout<<"5 ... INSIDE .... INTERSECTED .....\n";
			mat_sf.resize(p->getCoeffMatrix().size1(), col, true);//resizing to append the polytope p of the flowpipe that intersects
			//for (unsigned int j = 0; j < this->getTotalTemplateDirections();j++)
			for (unsigned int j = 0; j < p->getCoeffMatrix().size1(); j++) {
				mat_sf(j, col - 1) = p->getColumnVector()[j];//	mat_sf.addColumn(p.getColumnVector()); add column in matrix class
			}
			//	cout << "\nIntersection Found at = " << i << endl;
		}
		if (intersection_started == true && is_intersected == false) {
			intersection_ended = true;
		}
		if (intersection_ended == true && intersection_started == true) {
			intersected_region.push_back(template_polyhedra::ptr(new template_polyhedra(mat_sf, p->getCoeffMatrix())));
			cout<< "   polys = "<<mat_sf.size2() <<std::endl;
			mat_sf.clear();	//reset all values to zero
			mat_sf.resize(0, 0);
			col = 0;
			intersection_started = false;
			intersection_ended = false;
			//cout << "\nIntersection Ended at = " << i << "\n";
		}
		//	i2 = i;
	}	//end of parallel for-loop

	if (intersection_started == true && intersection_ended == false) {
		math::matrix<double> directions, templateDirs, invDirs;
		templateDirs = this->getTemplateDirections();
		invDirs = this->getInvariantDirections();
		templateDirs.matrix_join(invDirs, directions);
		intersected_region.push_back(template_polyhedra::ptr(new template_polyhedra(mat_sf, directions)));//intersected_region.push_back(template_polyhedra(mat_sf, p.getCoeffMatrix()));
		//cout<< "   polys = "<<mat_sf.size2() <<std::endl;
		//	cout << "\nIntersection did not End = " << i2 << "\n";
	}
	t101.stop();
	double clock101;
	clock101 = t101.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
	double return101 = clock101 / (double) 1000;
	std::cout << "\nSequential section Time:Wall(Seconds) = " << return101 <<std::endl;


	/*	if (foundIntersection == 0)
	 cout << "\nSorry Intersection Not Found!!!!!!\n";
	 else
	 cout << "\nIntersection Found at = " << intersection_start << endl;
	 */

	/*
	 * After checking all elements(Omegas) of the template_polyhedra of the reachable set
	 * the function returns the point_of_intersection_started with the guard G.
	 * And only returns the interescted region with the guard G as a new template_polyhedra:::with added invariant_directions
	 */
//	point_of_intersection = intersection_start;
	return intersected_region;
}

std::list<std::pair<unsigned int, unsigned int> > template_polyhedra::polys_intersectionSequential_optimize(polytope::ptr G, int lp_solver_type_choosen) {
	size_type row = 0;
	size_type col = 0;
	math::matrix<double> mat_sf(row, col);
	bool is_intersected = false, intersection_started = false, immediate_false = false, intersection_ended = false;
	int foundIntersection = 0, intersection_start;
	std::list<template_polyhedra::ptr> intersected_region;
	printf("suport function cols %d \n",this->Matrix_SupportFunction.size2());
	std::vector<bool> intersects(this->Matrix_SupportFunction.size2(), false); //all false initially
	for (unsigned int i = 0; i < this->Matrix_SupportFunction.size2(); i++) {
		//std::cout<<"\n Inner thread Template_polyhedra omp_get_num_threads() = "<< omp_get_num_threads()<<"\n";
		polytope::ptr p;
		p = this->getPolytope(i);
		std::vector<double> constraint_bound_values(this->getInvariantDirections().size1());
		constraint_bound_values = this->getInvariantBoundValue(i);
		p->setMoreConstraints(this->getInvariantDirections(), constraint_bound_values);

		intersects[i] = p->check_polytope_intersection(G, lp_solver_type_choosen); //result of intersection
		cout<<"intersects "<<intersects[i]<<"\n";
	} //end of parallel-loop :: we have the list of intersected polys

	std::list<std::pair<unsigned int, unsigned int> > intersected_range;
//	cout << "Is This Big = " << this->Matrix_SupportFunction.size2() << "\n";

	std::pair<unsigned int, unsigned int> inte_range;
	for (unsigned int i = 0; i < this->Matrix_SupportFunction.size2(); i++) { //sequential reading of an boolean_array that tells intersected polys
		is_intersected = intersects[i];
		if (is_intersected == true) { //if intersects create a new template_polyhedra subset
			intersection_started = true;
			if (foundIntersection == 0) { //Intersection start
				foundIntersection++; //intersection_start = i;	//Intersection started at this position of 'i'
				inte_range.first = i; //	cout << "\nIntersection Found at = " << i << endl;
			}
		}
		if (intersection_started == true && is_intersected == false) {
			intersection_ended = true;
			inte_range.second = i;
			foundIntersection = 0;
		}
		if (intersection_ended == true && intersection_started == true) {
			intersected_range.push_back(inte_range);
			intersection_started = false;
			intersection_ended = false; //cout << "\nIntersection Ended at = " << i << "\n";
		}
	} //end of parallel for-loop

	if (intersection_started == true && intersection_ended == false) {
		inte_range.second = this->Matrix_SupportFunction.size2() - 1;
		intersected_range.push_back(inte_range);  //	cout << "\nIntersection did not End = " << i2 << "\n";
	}
	return intersected_range;
}

std::list<polytope::ptr> template_polyhedra::flowpipe_intersectionSequential(polytope::ptr guard, int lp_solver_type_choosen){

	std::list<std::pair<unsigned int, unsigned int> > range_list;
	range_list = polys_intersectionSequential_optimize(guard,lp_solver_type_choosen);
	//cout <<"range_list.size = "<<range_list.size();
	std::list<polytope::ptr> polys;
	unsigned int poly_dir_size = this->template_Directions.size1() +this->invariant_Directions.size1();
	std::vector<double> colVector(poly_dir_size);
	printf("range list size %d",range_list.size());
	for (std::list<std::pair<unsigned int, unsigned int> >::iterator range_it = range_list.begin(); range_it != range_list.end();
			range_it++) {
		unsigned int start = (*range_it).first, end=(*range_it).second;
		//cout << "first = " << start << "  second = " << end << std::endl;
		for (unsigned int eachTemplateDir=0;eachTemplateDir<this->template_Directions.size1();eachTemplateDir++){

			double Max_sf=this->Matrix_SupportFunction(eachTemplateDir,start);
			for (int i = start+1; i <= end; i++) {
				double sf=this->Matrix_SupportFunction(eachTemplateDir,i);
				if (sf > Max_sf)
					Max_sf = sf;
			}//end of each intersected region
			colVector[eachTemplateDir] = Max_sf;
		}//end of each template direction ALSO HAVE TO PERFORM INVARIANT DIRECTION
		unsigned int total_dir = this->template_Directions.size1();
		for (unsigned int eachInvDir=0;eachInvDir<this->invariant_Directions.size1();eachInvDir++){
			double Max_sf=this->Matrix_InvariantBound(eachInvDir,start);
			for (int i = start + 1; i <= end; i++) {
				double sf = this->Matrix_InvariantBound(eachInvDir, i);
				if (sf > Max_sf)
					Max_sf = sf;
			} //end of each intersected region
			colVector[total_dir + eachInvDir] = Max_sf;
		}
		math::matrix<double> allDirs;
		this->template_Directions.matrix_join(this->invariant_Directions, allDirs);
		polytope::ptr p = polytope::ptr(new polytope(allDirs,colVector,1));
		polys.push_back(p);
	}//end of multiple intersected region
	//cout<<"polys.size = "<<polys.size()<<"\n";
	return polys;
}

const polytope::ptr template_polyhedra::getTemplate_approx(int lp_solver_type_choosen) {
	/*
	 * Here the calling template_polyhedra will have Matrix_SupportFunction and template_directions
	 * without any invariant_direction or Matrix_InvariantBound:::::: it ALSO has invariant_bounds included
	 * Get the number of directions as rows for each rows compute sf of all omega's and take the maximum of
	 these sf and add it as column vector then return the polytop p as p(direction, columnvector);
	 */
	unsigned int rows, max_or_min = 2;	//Maximizing
	double Max_sf, sf;
	int type = lp_solver_type_choosen;
	rows = this->getTotalTemplateDirections();//cout << "\nrows here = " << rows << endl;
	std::vector<double> columnVector(rows);
	std::vector<double> each_direction(this->template_Directions.size2());//can try putting it inside the loop
	for (unsigned int i = 0; i < rows; i++) {
		polytope::ptr p;
		lp_solver lp(type);
		p = this->getPolytope(0);	//First Omega's support function value
		lp.setMin_Or_Max(2);
		lp.setConstraints(p->getCoeffMatrix(), p->getColumnVector(), p->getInEqualitySign());

		for (unsigned int dim = 0; dim < this->template_Directions.size2();dim++)
			each_direction[dim] = this->getTemplateDirections()(i, dim);//get each direction

		Max_sf = p->computeSupportFunction(each_direction, lp);//First Omega's support function value
		cout<<"Max_sf  = "<<Max_sf <<std::endl;
		cout<<"This-value = "<<this->Matrix_SupportFunction(i,0)<<std::endl;

		for (unsigned int j = 1; j < this->total_iterations; j++) {
			lp_solver lp1(type);
			p = this->getPolytope(j);
			lp1.setMin_Or_Max(2);
			lp1.setConstraints(p->getCoeffMatrix(), p->getColumnVector(), p->getInEqualitySign());
			sf = p->computeSupportFunction(each_direction, lp1);

			if (sf > Max_sf)
				Max_sf = sf;
		}
		columnVector[i] = Max_sf;//Maximum of sf in one direction of all Omegas
	}
//std::cout<<"\nTotal_iterations = " <<this->total_iterations<<"\n";
	return polytope::ptr(new polytope(this->getTemplateDirections(), columnVector, 1));//1 as all signs are <=
}

/*const math::matrix<double>& template_polyhedra::getAllDirections() const {
 return All_Directions;
 }

 void template_polyhedra::setAllDirections(
 math::matrix<double>& allDirections) {
 this->All_Directions = allDirections;
 this->setTotalDirections(allDirections.size1());

 }*/

const math::matrix<double>& template_polyhedra::getMatrixSupportFunction() const {
	return Matrix_SupportFunction;
}

void template_polyhedra::setMatrixSupportFunction(
		math::matrix<double>& matrixSupportFunction) {
//	cout<<"Called\n";
	Matrix_SupportFunction = matrixSupportFunction;
	//cout<<"Called 2\n";
	this->setTotalIterations(matrixSupportFunction.size2());
}

unsigned int template_polyhedra::getTotalIterations() const {
	return total_iterations;
}

void template_polyhedra::setTotalIterations(unsigned int totalIterations) {
	total_iterations = totalIterations;
}

void template_polyhedra::resize_matrix_SupportFunction(int dir_nums,
		int iterations_before_intersection) {
	size_type row = dir_nums;
	size_type col = iterations_before_intersection;
	Matrix_SupportFunction.resize(row, col, true);
}

template_polyhedra::ptr template_polyhedra::union_TemplatePolytope(
		template_polyhedra::ptr& Tpoly) {

	if (this->total_iterations == 0)	//if the calling polyhedra is empty
		return template_polyhedra::ptr(
				new template_polyhedra(Tpoly->getMatrixSupportFunction(),
						Tpoly->getTemplateDirections()));

	//std::cout<<"\nEntered inside Union_templatePolytope\n";
	size_type rows = Tpoly->getMatrixSupportFunction().size1();	//rows will not change only column size will increase
	unsigned int k;
	size_type cols = Matrix_SupportFunction.size2()
			+ Tpoly->getMatrixSupportFunction().size2();
	//std::cout<<"\nRows = "<<rows<<"Cols = "<<cols<<"\n";

	math::matrix<double> new_SFMatrix;
	new_SFMatrix = Matrix_SupportFunction;//copy entire calling object's SFmatrix
	new_SFMatrix.resize(rows, cols, true);//size increased to append Tpoly's SFMatrix

	for (unsigned int i = 0; i < rows; i++) {
		k = Matrix_SupportFunction.size2();	//for each row k should begin from Maximum column size of calling Object's SFMatrix
		for (unsigned int j = 0; j < Tpoly->getMatrixSupportFunction().size2();
				j++) {
			new_SFMatrix(i, k) = Tpoly->getMatrixSupportFunction()(i, j);
			//	std::cout<<Tpoly.getMatrixSupportFunction()(i,j)<<"\t";
			k++;
		}
		//std::cout<<endl;
	}
	return template_polyhedra::ptr(
			new template_polyhedra(new_SFMatrix, this->getTemplateDirections()));
}
