
/*
 * flow_cost_estimate.h
 *
 *  Created on: 18-Apr-2016
 *      Author: ray
 */

#ifndef UTILITIES_FLOW_COST_ESTIMATE_H_
#define UTILITIES_FLOW_COST_ESTIMATE_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "application/DataStructureDirections.h"
#include "counterExample/simulation.h"

/**
 * Returns an estimated time after which the flow dynamics
 * takes the flow outside the invariant, starting from the
 * initial states X0
 */
double flow_cost_estimate(polytope::ptr X0, polytope::ptr I, Dynamics d, double time_horizon, double fine_time_step);

/*
 * Returns an estimated time after which the flow dynamics
 * takes the flow outside the invariant, starting from the
 * initial states X0.
 * It is more expensive than flow_cost_estimate(...) as more simulations are generated one for each
 * invariant face.
 */
double flow_cost_estimate_invFace(polytope::ptr X0, polytope::ptr I, Dynamics d, double time_horizon, double fine_time_step);

#endif /* UTILITIES_FLOW_COST_ESTIMATE_H_ */
