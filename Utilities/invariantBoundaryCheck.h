/*
 * invariantCheckWithReachability.h
 *
 *  Created on: 15-Oct-2014
 *      Author: amit
 */

#ifndef INVARIANTCHECKWITHREACHABILITY_H_
#define INVARIANTCHECKWITHREACHABILITY_H_

#include "application/DataStructureDirections.h"
#include "application/sf_utility.h"
#include "application/CopyArray.h"
#include <omp.h>
#include "core_system/continuous/ConvexSet/transMinkPoly.h"
#include "core_system/continuous/ConvexSet/supportFunctionProvider.h"
#include <stdlib.h>	//for using abs() absolute function

/*
 * Sequential Algorithm using support-function-algorithm to perform invariant crossing check
 * Returns the iterations(=Time-horizon/time-step), when the flowpipe(Omega) is completely outside
 * searches in-efficiently one after another invariant face
 * NOTE::
 * This is our Old IN-EFFICIENT implementation where positive and negative normals to invariant's faces are computed
 * to check the points of intersection START and END. But START point were redundant task and was the caused for XSpeed's slow-performance
 */
void InvariantBoundaryCheck1(Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,	int lp_solver_type_choosen, unsigned int &newTotIters);

/*
 * Sequential Algorithm using support-function-algorithm to perform invariant crossing check
 * searches EFFICIENTLY at one go to all invariant faces.
 * Also computes ONLY negative direction of the invariant's faces
 *
 * Returns the iterations(=Time-horizon/time-step), for which flowpipe's Omegas are completely inside the invariant
 * NOTE::
 * it returns based on "1-indexing" i.e., if Omega_0, Omega_1 and Omega_2 are completely inside it returns 3 as the
 * number of Omega's that need to be constructed in flowpipe is three.
 */
void InvariantBoundaryCheck(Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,	int lp_solver_type_choosen, unsigned int &newTotIters);

/*
 * Using Binary Search Algorithm to perform a Quick Check if an Omega crossing the invariant (uses actual set computation)
 * Returns the iterations(=Time-horizon/time-step), when the flowpipe is completely outside
 */
void quickInvariantBoundaryCheck(Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant, int lp_solver_type_choosen, unsigned int &newTotIters);

/*
 * Returns an initialSet at time equals to START_TIME uses the transMinkPoly class to get the transposed Set at time equals to START_TIME
 */
supportFunctionProvider::ptr getInitialSet(double START_TIME, ReachabilityParameters& ReachParameters, Dynamics& SystemDynamics);

/*
 * Uses the Binary Search algorithm to perform invariant boundary check
 * Also performs a quick invariant check with the Last Omega if included just returns the total
 * iterations that is supplied by the user otherwise, when the invariant is in-between then Binary Search Algo is performed
 */
unsigned int invariantCheck(std::vector<double>& pos_dir, std::vector<double>& neg_dir, double SearchKey, ReachabilityParameters& ReachParameters,
		Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial, int lp_solver_type_choosen);


//Invariant check using Slow-Start incremental time-step
//NOTE :: this approach of increasing the time-step does not work in support-function algorithm
// as change in time-step increases rate of bloating 'error'
void SlowStartInvariantBoundaryCheck(Dynamics& SystemDynamics,supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParametersOld, polytope::ptr invariant, int lp_solver_type_choosen, unsigned int &newTotIters);

/*
 * Uses the definition of able to compute a set at time instance 't' given the initial set.
 * Using COARSE time-step it Generate sets and test if it crosses the invariant, if so than
 * it starts sequential FINE search to pin point the exact time
 * Returns the iterations(=Time-horizon/time-step), for which flowpipe's Omegas are completely inside the invariant
 * NOTE::
 * it returns based on "1-indexing" i.e., if Omega_0, Omega_1 and Omega_2 are completely inside it returns 3 as the
 * number of Omega's that need to be constructed in flowpipe is three.
 */
void jumpInvariantBoundaryCheck(Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant, int lp_solver_type_choosen, unsigned int &newTotIters);

/*
 * Returns THE time when a convex set just crosses the invariant boundary one after another of a invariant polyhedra.
 * Note: This is an in-efficient method
 */
double invariantCrossingCheck(double START_TIME, double time_step, double time_horizon, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant, Dynamics& SystemDynamics, int lp_solver_type_choosen);

/*
 * Checking for each invariant direction individually
 * Returns the time-bound when it just crosses the invariant bound
 */
double invariantFaceCrossingCheck(std::vector<double>& neg_dir, double invBound, double START_TIME, double time_step, double time_horizon,
		supportFunctionProvider::ptr Initial, ReachabilityParameters& ReachParameters, Dynamics& SystemDynamics, int lp_solver_type_choosen);

/*
 * Returns the time-bound when a convex set represented by support-function (Initial), just crosses the invariant bound
 * efficient approach to detect flowpipe cost estimation.
 */
double invariantCrossingCheck1(double START_TIME, double time_step, double time_horizon, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant, Dynamics& SystemDynamics, int lp_solver_type_choosen, std::string fileName);

/*
 * Returns a bounded polytope from a convex set 'p' represented using support-function using templated directions in ReachParameters
 */
polytope::ptr getBoundedConvexSet(supportFunctionProvider::ptr& p, ReachabilityParameters& ReachParameters);

#endif /* INVARIANTCHECKWITHREACHABILITY_H_ */
