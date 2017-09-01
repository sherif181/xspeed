/*
 * AGJH.h
 *
 *  Created on: 28-Oct-2016
 *      Author: amit
 */

#ifndef AGJH_H_
#define AGJH_H_

#include "reachability.h"

class agjh : public reachability{
public:
	/**
	 * Adaption of GH algorithm for parallel Breadth first search
	 */
	std::list<symbolic_states::ptr> ParallelBFS_GH();
private:
	/**
	 * Computes postC from sym state s.
	 */
	template_polyhedra::ptr postC(initial_state::ptr s);

	/**
	 * computes postD from a template polyhedra computed by postC
	 */

	std::list<initial_state::ptr> postD(symbolic_states::ptr symb);
};


#endif /* AGJH_H_ */
