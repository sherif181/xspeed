/*
 * utility_functions.cpp
 *
 *  Created on: 19-Mar-2015
 *      Author: amit
 */

#include <core_system/Reachability/NewApproach/utility_functions.h>

double compute_theta(std::vector<double> l1, std::vector<double> l2){

	assert(l1.size() == l2.size());

	double sum1=0.0, sum2=0.0, prod_sum=0.0;

	for (unsigned int i=0;i<l1.size();i++){
		sum1 = sum1 + l1[i] * l1[i];
		sum2 = sum2 + l2[i] * l2[i];
		prod_sum = prod_sum + l1[i] * l2[i];
	}

	double res_deno = sqrt(sum1) * sqrt(sum2);
	if (res_deno == 0)
		return 999;

	double res = prod_sum / res_deno;
	double theta = acos (res);
	return theta;
}


