/*
 * DataStructureHybridAutomata.h
 *
 *  Created on: 09-Jul-2014
 *      Author: amit
 */

#ifndef DATASTRUCTUREHYBRIDAUTOMATA_H_
#define DATASTRUCTUREHYBRIDAUTOMATA_H_

#include "core_system/math/matrix.h"
/*
 * Assignment of X' = M x + b
 * X' = R X + w
 */
struct Assign {
	math::matrix<double> Map;
	std::vector<double> b;
};

#endif /* DATASTRUCTUREHYBRIDAUTOMATA_H_ */
