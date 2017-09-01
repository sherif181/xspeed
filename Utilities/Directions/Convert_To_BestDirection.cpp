/*
 * Convert_To_BestDirection.cpp
 *
 *  Created on: 07-Dec-2014
 *      Author: amit
 */

#include "Utilities/Directions/Convert_To_BestDirection.h"

std::vector<double> getBestDirection(std::vector<double> dir) {
	std::vector<double> r_dir;

	return r_dir;
}

std::vector<double> getBestDirection(polytope::ptr initial,
		std::vector<double> dir, int lp_solver_type_choosen) {
	std::vector<double> best_dir, anticlock_dir, clock_dir;
	double min_sf, sf_clockwise, sf_anticlockwise, sin_val, cos_val;
	//sin(radian) and 1 radian = 57.295 degrees; 1 degree = (1/57.295)radians

	int flag_clockwise = 1, flag_anticlockwise = 1;
	math::matrix<double> R_clock_wise(2, 2), R_anti_clock_wise(2, 2);
	int type = lp_solver_type_choosen;
	lp_solver lp(type), lp_dummy(type);
	lp.setMin_Or_Max(2);

	if (!initial->getIsEmpty()) //set glpk constraints If not an empty polytope
		lp.setConstraints(initial->getCoeffMatrix(), initial->getColumnVector(),
				initial->getInEqualitySign());
	clock_dir = dir;
	anticlock_dir = dir;
	best_dir = dir;
	min_sf = initial->computeSupportFunction(dir, lp);
	//Now rotate the direction in +/- 45 degrees(can be improved later) and compute sf and get the min of all sf
	for (int degree = 1; degree <= 45; degree += 3) {
		sin_val = sin(degree / 57.295); //sin(-theta) = - sin(+theta)
		cos_val = cos(degree / 57.295); //cos(+theta) = cos(-theta)
		if (flag_anticlockwise == 1) {
			R_anti_clock_wise(0, 0) = cos_val;
			R_anti_clock_wise(0, 1) = -1 * sin_val;
			R_anti_clock_wise(1, 0) = sin_val;
			R_anti_clock_wise(1, 1) = cos_val;
			R_anti_clock_wise.mult_vector(dir, anticlock_dir); //best_dir = R * best_dir;	//rotated/transformation of dir is computed
			sf_anticlockwise = initial->computeSupportFunction(anticlock_dir,
					lp);
		}
		if (flag_clockwise == 1) {
			R_clock_wise(0, 0) = cos_val;
			R_clock_wise(0, 1) = sin_val;
			R_clock_wise(1, 0) = -1 * sin_val;
			R_clock_wise(1, 1) = cos_val;
			R_clock_wise.mult_vector(dir, clock_dir); //best_dir = R * best_dir;	//rotated/transformation of dir is computed
			sf_clockwise = initial->computeSupportFunction(clock_dir, lp);
		}
		if (sf_anticlockwise < min_sf) {
			min_sf = sf_anticlockwise;
			best_dir = anticlock_dir;
			flag_clockwise = 0;
		}
		if (sf_clockwise < min_sf) {
			min_sf = sf_clockwise;
			best_dir = clock_dir;
			flag_anticlockwise = 0;
		}
	}
	return best_dir;
}

bool Issame_direction(std::vector<double> dir1, std::vector<double> dir2) {
	double epsilon = 0.2; //more appropriate function is needed
	std::vector<double> difference(dir1.size());
	for (unsigned int i = 0; i < dir1.size(); i++) {
		difference[i] = abs(dir1[i] - dir2[i]); //absolute(differece) <= epsilon
		if (difference[i] > epsilon)
			return false; //far apart vectors
	}
	return true;
}
bool IsDirectionExist(math::matrix<double> directions,
		std::vector<double> best_dir) {
	std::vector<double> dir(best_dir.size());
	for (unsigned int i = 0; i < directions.size1(); i++) {
		for (unsigned int j = 0; j < directions.size2(); j++) {
			dir[j] = directions(i, j);
		}
		if (Issame_direction(best_dir, dir)) {
			return true;
		}
	} //finished comparing with all the directions
	return false;
}

math::matrix<double> getBestTemplateDirections(polytope::ptr initial,
		math::matrix<double> directions, int lp_solver_type_choosen) {
	math::matrix<double> r_directions;
	int row = 0;
	std::vector<double> dir(directions.size2()), best_dir(directions.size2());
	for (unsigned int i = 0; i < directions.size1(); i++) {
		for (unsigned int j = 0; j < directions.size2(); j++) {
			dir[j] = directions(i, j);
		}
		best_dir = getBestDirection(initial, dir, lp_solver_type_choosen);
		if (i == 0) { //for the first direction just append
			r_directions.resize(row + 1, dir.size()); //only 1 directions in the r_directions matrix
			for (unsigned int col_index = 0; col_index < dir.size();
					col_index++)
				r_directions(row, col_index) = best_dir[col_index];
			row++; //for the next entry of r_directions
		} else {
			if (!IsDirectionExist(r_directions, best_dir)) { //append only the best_dir if it does not already exists in r_directions
				r_directions.resize(row + 1, dir.size(), true); // the r_directions matrix
				for (unsigned int col_index = 0; col_index < dir.size();
						col_index++)
					r_directions(row, col_index) = best_dir[col_index];
				row++; //for the next entry of r_directions
			}
		}
	}

	return r_directions;
}

math::matrix<double> getBestEqualTemplateDirections(polytope::ptr initial,
		math::matrix<double> directions, int lp_solver_type_choosen) {
	math::matrix<double> r_directions;
	int row = 0;
	std::vector<double> dir(directions.size2()), best_dir(directions.size2());
	for (unsigned int i = 0; i < directions.size1(); i++) {
		for (unsigned int j = 0; j < directions.size2(); j++) {
			dir[j] = directions(i, j);
		}
		best_dir = getBestDirection(initial, dir, lp_solver_type_choosen);
		if (i == 0) { //for the first direction just append
			r_directions.resize(row + 1, dir.size()); //only 1 directions in the r_directions matrix
			for (unsigned int col_index = 0; col_index < dir.size();
					col_index++)
				r_directions(row, col_index) = best_dir[col_index];
		} else {
			r_directions.resize(row + 1, dir.size(), true); // resize the r_directions matrix
			if (!IsDirectionExist(r_directions, best_dir)) { //append only the best_dir if it does not already exists in r_directions
				for (unsigned int col_index = 0; col_index < dir.size();
						col_index++)
					r_directions(row, col_index) = best_dir[col_index];
			} else { //otherwise append the same direction
				for (unsigned int col_index = 0; col_index < dir.size();
						col_index++)
					r_directions(row, col_index) = dir[col_index];
			}
		}
		row++; //for the next entry of r_directions
	}

	return r_directions;
}
