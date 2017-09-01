/*
 * sf_directions.h
 *
 *  Created on: 01-May-2014
 *      Author: gurung
 */

#ifndef SF_DIRECTIONS_H_
#define SF_DIRECTIONS_H_
#define UNSOLVED 0

#include <vector>
#include "application/DataStructureDirections.h"
#include "application/CopyArray.h"
#include <list>
#include "omp.h"
//#include "matrixOperation.h"

/**
 * Function to get the axis directions
 * N is the dimension of the system
 */
std::vector<std::vector<double> > generate_axis_directions(unsigned int N);

//typedef std::vector<std::vector<double> > direction_list;
std::vector<std::vector<double> > get_octagonal_directions(unsigned int dim);

/**
 * Create the list of all directions that the sf algorithm calls for lp solving
 * on initial set (I) and input set (V).
 *
 * facet_dirs: set of initial directions (including the Axis directions) in which to compute the support function of the reach set.
 * ITERS: Number of iterations in the reachability algorithm
 * A : The flow matrix
 *
 */
std::vector<D> get_directions(ReachabilityParameters reach_parameter);

/*
 * Creates all the list of directions for computing support function using GPU
 */
std::vector<AllDirection> get_DirectionList(
		ReachabilityParameters &ReachParameters, unsigned int newiters);

/*
 * Generate two lists ie list1 and list2 of matrix of directions for X0 and U for computing SF using GPU
 */
void getDirectionList_X0_and_U_OLD(std::vector<AllDirection> &directionList,
		unsigned int numDirections, unsigned int TotalIterations,
		math::matrix<float> &list1, math::matrix<float> &list2);

void getDirectionList_X0_and_U(int numCoresAvail, ReachabilityParameters &ReachParameters,
		unsigned int newiters, math::matrix<float> &list_X0,
		math::matrix<float> &list_U, bool U_empty, Dynamics& SystemDynamics);

void getDirectionList_X0_and_U_OnlyForGLPK(
		ReachabilityParameters &ReachParameters, unsigned int newiters,
		std::list<std::vector<double> > &list_X0,
		std::list<std::vector<double> > &list_U, bool U_empty, Dynamics& SystemDynamics);

#endif /* SF_DIRECTIONS_H_ */
