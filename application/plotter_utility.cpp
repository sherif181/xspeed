/*
 * plotter_utility.cpp
 *
 *  Created on: 18-Nov-2015
 *      Author: amit
 */
#include "plotter_utility.h"

//std::list<template_polyhedra>
void plot_output(std::list<symbolic_states::ptr>& Symbolic_states) {
	/*Steps:
	 * 1) Generate a single polytope from the reachability_sfm(which is a list of template_polyhedra)
	 * 2) Call create constraint representation of the single polytope in the form of a string
	 *
	 * 3) Call the SpaceEx plotter_function(have to create an interface which takes string as input and
	 * 	  returns a matrix(n-dim) of vertices for plotting)
	 *
	 * 4) Plotting variables(x,y) can be selected from user and then only the selected columns from the matrix can be written on to a FILE.
	 *    This FILE can be plotted using the system(graph -T -B .... )utility
	 */

	std::list<symbolic_states::ptr>::iterator it;
	for (it = Symbolic_states.begin(); it != Symbolic_states.end(); it++) { //Each (*it) is an individual flowpipe

		math::matrix<double> sfm, template_direction, invariant_bound_values,
				invariant_direction;

		// Append template_direction and invariant_directions to form matrix "A"
		// Append sfm and invariant_bound_values to form Matrix of SupportFunction_Values including the Invariant_Bound (each column makes 1 polytope)

		math::matrix<double> SupportFunction_Values, All_Directions;

		sfm = (*it)->getContinuousSetptr()->getMatrixSupportFunction();	//MatrixSupportFunction
		template_direction = (*it)->getContinuousSetptr()->getTemplateDirections();	//Direction
		invariant_direction = (*it)->getContinuousSetptr()->getInvariantDirections(); //invariant_directions
		invariant_bound_values = (*it)->getContinuousSetptr()->getMatrix_InvariantBound(); //invariant_bound_matrix

		sfm.matrix_join(invariant_bound_values, SupportFunction_Values);
		template_direction.matrix_join(invariant_direction, All_Directions);

		std::list<polytope> polytope_list;
		std::vector<double> column_vector(SupportFunction_Values.size1()); //total rows

		for (int i = 0; i < SupportFunction_Values.size2(); i++) { //each column
			for (int j = 0; j < SupportFunction_Values.size1(); j++) {
				column_vector[j] = SupportFunction_Values(j, i); //retrives column by column
			}
			polytope_list.push_back(polytope(All_Directions, column_vector, 1)); //1 for <= sign
		} // A list of Polytope is created instead of calling the plotter with a single polytope

		std::string str_lin_constraints;
		str_lin_constraints = "";
		for (std::list<polytope>::iterator poly = polytope_list.begin();
				poly != polytope_list.end(); poly++) {
			math::matrix<double> A;
			std::vector<double> b;
			A = (*poly).getCoeffMatrix();
			b = (*poly).getColumnVector();

			for (int i = 0; i < A.size1(); i++) {
				if (i > 0)
					str_lin_constraints.append(" & ");
				bool found_var_before = false;
				for (int j = 0; j < A.size2(); j++) {

					if (A(i, j) == 1) {
						if (j == 0) {
							str_lin_constraints.append("x");
							int num = j + 1;
							str_lin_constraints.append(
									boost::lexical_cast<string>(num)); //x1, x2  etc
						} else {
							if (found_var_before) {
								str_lin_constraints.append("+x");
								int num = j + 1;
								str_lin_constraints.append(
										boost::lexical_cast<string>(num)); //x1, x2  etc
							} else {
								str_lin_constraints.append("x");	//x1, x2 etc have not been found yet. This is the 1st occurrence
								int num = j + 1;
								str_lin_constraints.append(
										boost::lexical_cast<string>(num)); //x1, x2  etc
							}
						}
						found_var_before = true;
					} else if (A(i, j) == -1) {
						str_lin_constraints.append("-x");
						int num = j + 1;
						str_lin_constraints.append(
								boost::lexical_cast<string>(num)); //-x1, -x2  etc
						found_var_before = true;
					} // == 0 do not add any constraint

				}//end of one constraint
				str_lin_constraints.append("<=");
				str_lin_constraints.append(boost::lexical_cast<string>(b[i]));
			}//end of the polytope
			std::cout<<"\nConstraint Representation  = "<<str_lin_constraints<<std::endl;

			//Now here i can call SpaceEx plotter


		}//end of the list of polytope

	} //End of a single FlowPipe

}

