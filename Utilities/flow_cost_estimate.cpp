
/*
 * flow_cost_estimate.cpp
 *
 *  Created on: 18-Apr-2016
 *      Author: ray
 */

#include "flow_cost_estimate.h"

double flow_cost_estimate(polytope::ptr X0, polytope::ptr I, Dynamics d, double time_horizon, double fine_time_step)
{
	unsigned int dim = X0->getSystemDimension();
	assert(I->getSystemDimension() == dim);

	// Choose a coarse time-step for efficiency
	unsigned int granularity_index = 100;
	//double time_step = time_horizon/granularity_index;
	double time_step = fine_time_step;

	simulation sim(dim,time_step,d);
	std::vector<double> x0;

	//get a random point inside X0
	std::vector<double> obj(dim, 0);
	obj[0] = 1;
	lp_solver lp(GLPK_SOLVER);
	lp.setConstraints(X0->getCoeffMatrix(), X0->getColumnVector(),
			X0->getInEqualitySign());
	lp.Compute_LLP(obj);
	x0 = lp.get_sv();

	bound_sim simv; bool status;
	simv = sim.bounded_simulation(x0,time_horizon,I,status);
	double coarse_cross_time = simv.cross_over_time;

	//debug
	//std::cout << "Printing the coarse time returned :" << simv.cross_over_time << std::endl;
	//sim.set_time_step(fine_time_step);
	//--
	// re-simulate with a fine time-step for accuracy
	//simv = sim.bounded_simulation(simv.v, time_step, I);
	//return coarse_cross_time - time_step + simv.cross_over_time;
	return coarse_cross_time;
}



double flow_cost_estimate_invFace(polytope::ptr X0, polytope::ptr I, Dynamics d, double time_horizon, double fine_time_step)
{
	unsigned int dim = X0->getSystemDimension();
	assert(I->getSystemDimension() == dim);

	std::vector<double> allCrossingTime(I->getColumnVector().size());
/*
	std::set<std::vector<double> > pts;
	for (int i=0;i<->getCoeffMatrix().size1();i++) {
				for (int j = 0; j < I->getCoeffMatrix().size2(); j++) {
					obj[j] = I->getCoeffMatrix()(eachInv, j);
					pts.insert(obj);	//removes the redundant points
				}
			}*/

//for (int eachInv = 0;eachInv < I->getColumnVector().size();eachInv++){

	std::vector<double> x0;
	//get a random point inside X0
	std::vector<double> obj(dim, 0);

	//obj[0] = 1;

	lp_solver lp(GLPK_SOLVER);
	lp.setConstraints(X0->getCoeffMatrix(), X0->getColumnVector(),
			X0->getInEqualitySign());
	lp.Compute_LLP(obj);
	x0 = lp.get_sv();

	// Choose a coarse time-step for efficiency
	unsigned int granularity_index = 100;
	//double time_step = time_horizon/granularity_index;
	double time_step = fine_time_step;

	simulation sim(dim,time_step,d);


	bound_sim simv;bool status;
	simv = sim.bounded_simulation(x0,time_horizon,I,status);
	double coarse_cross_time = simv.cross_over_time;
//}
	//debug
	//std::cout << "Printing the coarse time returned :" << simv.cross_over_time << std::endl;
	//sim.set_time_step(fine_time_step);
	//--
	// re-simulate with a fine time-step for accuracy
	//simv = sim.bounded_simulation(simv.v, time_step, I);
	//return coarse_cross_time - time_step + simv.cross_over_time;
	return coarse_cross_time;
}
