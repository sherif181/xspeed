/*
 * analyticODESol.h
 *
 *  Created on: 06-Sep-2016
 *      Author: rajarshi
 */

#ifndef ANALYTICODESOL_H_
#define ANALYTICODESOL_H_

#include <vector>
#include "application/DataStructureDirections.h"

/*
 * Returns the solution of an ODE of the form X' = A*X+b at some time
 */
std::vector<double> ODESol(std::vector<double> x0, const Dynamics& D, double time);


#endif /* ANALYTICODESOL_H_ */
