/*
 * nlpFunctions.h
 *
 *  Header file with the required functions to call NLP solver: NLOPT
 *  Created on: 25-Sep-2016
 *      Author: rajarshi
 */

#ifndef NLP_FUNC_H
#define NLP_FUNC_H

#include "core_system/math/matrix.h"

#define VALIDATION

struct polyConstraints {
	math::vector<double> a;
	double b;
	unsigned int sstate_index;
};

struct boundConstriant {
	double bound;
	unsigned int var_index;
	bool is_ge; // to mark if bound is a >= constraint
};

/**
 * Objective function for splicing with HA constraints only. (Sriram et. al.)
 */
double myobjfunc1(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);

/**
 * Objective function for splicing with Flowpipe constraints only. (Sergiy's Idea)
 */
double myobjfunc2(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);

/**
 * Objective function for splicing with mixed NLP-LP . (Goran's Idea)
 */
double myobjfunc3(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);


double myconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data);

double myBoundConstraint(const std::vector<double> &x, std::vector<double> &grad, void *data);


#endif /* nlpFunctions.h */
