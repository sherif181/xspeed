/*
 * transMinkPoly.h
 *
 *  Created on: 16-Nov-2014
 *      Author: amit
 */

#include "core_system/continuous/ConvexSet/supportFunctionProvider.h"
#include "core_system/continuous/Polytope/Polytope.h"
#include <boost/shared_ptr.hpp>
#include <vector>

#ifndef TRANSMINKPOLY_H_
#define TRANSMINKPOLY_H_

/**
 * This class is to implement convex sets of the form:
 * C' = ATRANS. C \oplus t. BTRANS.U \oplus \oplus beta.B , where C, U are polytopes
 * and ATRANS, BTRANS are transformation matrices, Beta is a constant and B is a unit ball over the specified norm
 */

class transMinkPoly : public supportFunctionProvider
{
	polytope::ptr X0;
	polytope::ptr U;
	math::matrix<double> TRANS;
	math::matrix<double> B_TRANS;
	double beta;
	double time;

	bool Cempty;
	std::vector<double> C;

public:
	typedef boost::shared_ptr<transMinkPoly> ptr;

	/* most explicit constructor */
	transMinkPoly(polytope::ptr myX0, polytope::ptr myU, math::matrix<double> myTRANS, math::matrix<double> myB_TRANS, double mytime, double mybeta);

	/** Constructor to represent convex sets which is a linear transformation of another convex set only: C' = Trans. C */
	transMinkPoly(polytope::ptr myX0, math::matrix<double> myTRANS);

	/*
	 * constructor with x(t), u(t) and + C of the ODE x'(t) = Ax(t) + Bu(t) + C
	 */
	transMinkPoly(polytope::ptr myX0, polytope::ptr myU, std::vector<double> c,
			math::matrix<double> myTRANS, math::matrix<double> myB_TRANS, double mytime, double mybeta);

	unsigned int getSystemDimension() const;

	/*
	 * Returns TRUE when both X0 and U are empty polytopes
	 * returns FALSE if any one of the polytopes are empty
	 */
	bool getIsEmpty() const;

	double computeSupportFunction(std::vector<double> direction, lp_solver &lp);
	double max_norm(int lp_solver_type_choosen, unsigned int dim_for_Max_Norm);

};



#endif /* TRANSMINKPOLY_H_ */
