/*
 * utility_functions.h
 *
 *  Created on: 19-Mar-2015
 *      Author: amit
 */

#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include <math.h>
#include <vector>
#include <assert.h>
/*
 * Returns the angle between the two vectors l1 and l2 iff (|l1|.|l2| != 0) otherwise returns -999
 * 999 indicates that angle between them can not be computed
 */
double compute_theta(std::vector<double> l1, std::vector<double> l2);


#endif /* UTILITY_FUNCTIONS_H_ */
