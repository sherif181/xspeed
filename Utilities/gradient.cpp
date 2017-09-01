/*
 * gradient.cpp
 *
 *  Created on: 18-Aug-2016
 *      Author: rajarshi
 */

#include <Utilities/gradient.h>

/*
 * Computes derivative of point to polytope distance x to I w.r.t x.
 * The definition of point to polytope distance is as defined in the
 * polytope.cpp class implementation
 */
std::vector<double> dist_grad(std::vector<double> trace_end_pt, polytope::ptr I, std::vector<double> chain_mult)
{
	assert(chain_mult.size() == trace_end_pt.size());
	std::vector<double> grad(trace_end_pt.size(),0 );

	math::matrix<double> C = I->getCoeffMatrix();
	std::vector<double> b = I->getColumnVector();

	assert(trace_end_pt.size() == C.size2());
	assert(I->getInEqualitySign() == 1);

	double facet_distance = 0;
	double coef_sq_sum = 0;
	for(unsigned int i=0;i<C.size1();i++){
		for(unsigned int j=0;j<C.size2();j++){
			facet_distance += trace_end_pt[j]*C(i,j);
			coef_sq_sum += C(i,j)*C(i,j);
		}
		facet_distance -=b[i];
		if(facet_distance > 0){
			for(unsigned int k=0;k<grad.size();k++){
				grad[k] += C(i,k)/math::sqrt(coef_sq_sum);
			}
		}
		coef_sq_sum = 0;
		facet_distance = 0;
	}
	for(unsigned int k=0;k<trace_end_pt.size();k++)
		grad[k] = grad[k]*chain_mult[k];
	return grad;
}
