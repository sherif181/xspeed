/*
 * StandardVector.h
 *
 *  Created on: 12-Jul-2014
 *      Author: amit
 */

#ifndef STANDARDVECTOR_H_
#define STANDARDVECTOR_H_

#include <vector>
#include "assert.h"

using namespace std;

inline std::vector<double> vector_join(std::vector<double> v1,
		std::vector<double> v2) {

	std::vector<double> result;
	result = v1;
	unsigned int tot_size;
	tot_size = result.size() + v2.size();
	result.resize(tot_size);

	for (unsigned int i = v1.size(), j = 0; j < v2.size(); i++, j++)
		result[i] = v2[j];

	return result;
}

inline std::vector<double> vector_add(std::vector<double> v1,
		std::vector<double> v2) {

	std::vector<double> result;
	assert(v1.size() == v2.size());
	result.resize(v2.size());

	for (unsigned int i = 0; i < v2.size(); i++)
		result[i] = v1[i] + v2[i];

	return result;
}

#endif /* STANDARDVECTOR_H_ */
