/*
 * io_utility.h
 *
 *  Created on: 19-Jan-2016
 *      Author: amit
 */
#include "io_utility.h"

typedef std::vector<std::pair<double, double> > Intervals;
void Interval_Generator(std::list<symbolic_states::ptr>& Symbolic_states_list,
		std::list<std::pair<int, Intervals> > & location_interval_outputs,
		initial_state::ptr& init_state) {

	Intervals Interval_Outputs(
			init_state->getInitialSet()->getSystemDimension());

	double max_value;
	std::list<symbolic_states::ptr>::iterator SS;
	for (SS = Symbolic_states_list.begin(); SS != Symbolic_states_list.end();
			SS++) {
		//Each sysmbolic_state or each Location
		int locID;
		discrete_set ds;
		ds = (*SS)->getDiscreteSet();
		for (std::set<int>::iterator it = ds.getDiscreteElements().begin();
				it != ds.getDiscreteElements().end(); ++it) {
			locID = (*it); //Assuming only a single element exist in the discrete_set
		}
		std::pair<int, Intervals> loc_interval;
		loc_interval.first = locID;
		math::matrix<double> each_sfm;
		each_sfm = (*SS)->getContinuousSetptr()->getMatrixSupportFunction();

		int Totaldirs = each_sfm.size1();
		for (int i = 0; i < Totaldirs; i++) { //i==row_number
			for (unsigned int k = 0; k < each_sfm.size2(); k++) { //k==col_number
				double sfm_value = each_sfm(i, k);
				if (k == 0) {
					max_value = sfm_value;
				} else {
					if (sfm_value > max_value) {
						max_value = sfm_value;
					}
				}
			}	//getting the max_value for each row
			int index = i / (int) 2; //getting the variable_index
			if ((i % 2) == 0) { //even row is right_value of the interval(ie Max value/positive direction value)
				Interval_Outputs[index].second = max_value;
			} else { //left_value of the interval(ie Min value/negative direction value)
				Interval_Outputs[index].first = -1 * max_value;
			}
		} //end of sfm returns vector of all variables[min,max] intervals

		loc_interval.second = Interval_Outputs;
		location_interval_outputs.push_back(loc_interval);
	} // end-of-SS
}

//Creating a bounding_box polytope with constraints as (template_directions + invariant_directions)
// and bounds as maximum of (sfm + invariant_bounds)
polytope::ptr convertBounding_Box(template_polyhedra::ptr sfm) {
	polytope::ptr boundingPolytope;
	boundingPolytope = polytope::ptr(new polytope());

	math::matrix<double> directional_constraints, all_dirs;
	math::matrix<double> each_sfm;

	if (sfm->getInvariantDirections().size2()==0){	//no invariant exists
		all_dirs = sfm->getTemplateDirections();
	}else{	//invariants exists so join the directions
		//std::cout<<" called  ";
		directional_constraints = sfm->getTemplateDirections();
		directional_constraints.matrix_join(sfm->getInvariantDirections(), all_dirs);
	}
	each_sfm = sfm->getMatrixSupportFunction();
	//boundingPolytope->setCoeffMatrix(directional_constraints);
	boundingPolytope->setCoeffMatrix(all_dirs);

	double max_value;
	int Totaldirs = each_sfm.size1();
	std::vector<double> polytope_bounds(Totaldirs);

	for (int i = 0; i < Totaldirs; i++) { //i==row_number
		for (unsigned int k = 0; k < each_sfm.size2(); k++) { //k==col_number
			double sfm_value = each_sfm(i, k);
			if (k == 0) {
				max_value = sfm_value;
			} else {
				if (sfm_value > max_value) {
					max_value = sfm_value;
				}
			}
		}	//getting the max_value for each row
		polytope_bounds[i] = max_value;
	} //end of sfm returns vector of all variables[min,max] intervals

	std::vector<double> all_polytope_bounds;

	if (sfm->getInvariantDirections().size2()==0){	//no invariant exists
		all_polytope_bounds = polytope_bounds;
	}else{	////invariants exists so join the bounds
		std::vector<double> inv_bounds(sfm->getMatrix_InvariantBound().size1());
		for (unsigned int k = 0; k < sfm->getInvariantDirections().size1(); k++) { //k==rows or number of invariants
			inv_bounds[k] = sfm->getMatrix_InvariantBound()(k,0);
		}
		all_polytope_bounds = vector_join(polytope_bounds, inv_bounds);
	}

	boundingPolytope->setColumnVector(all_polytope_bounds);
	boundingPolytope->setInEqualitySign(1);//Indicating all <= sign

	return boundingPolytope;
}
