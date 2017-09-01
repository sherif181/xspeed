/*
 * abstract_symbolic_state.h
 *
 *  Created on: 19-Jan-2016
 *      Author: amit
 */

#ifndef ABSTRACT_SYMBOLIC_STATE_H_
#define ABSTRACT_SYMBOLIC_STATE_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "core_system/discrete/discrete_set/discrete_set.h"
#include "core_system/symbolic_states/initial_state.h"
#include "core_system/continuous/ConvexSet/supportFunctionProvider.h"

/**
 * Class to represent the abstract flowpipe. An abstract flowpipe
 * is a polyhedral approximation of the computed flowpipe from the
 * postC operator. The postC generates an over-approximate flowpipe
 * and this class stores a polyhedral(hyperbox) abstraction of the flowpipe.
 */

class abstract_symbolic_state {
public:
	typedef boost::shared_ptr<abstract_symbolic_state> ptr;

	abstract_symbolic_state();
	abstract_symbolic_state(discrete_set& discreteSet,
			polytope::ptr continuousSet, polytope::ptr polyI);

	/**
	 * Returns the polyhedral abstraction of the flowpipe
	 */
	polytope::ptr getContinuousSet() const;
	/**
	 * Returns the initial polytope of the abstract flowpipe.
	 */
	polytope::ptr getInitialSet(){
		return I;
	};
	void setContinuousSet(polytope::ptr continuousSet);
	/**
	 * Return the location id of this abstract symbolic state
	 */
	const discrete_set getDiscreteSet(){
		return DiscreteSet;
	}

	const discrete_set& getDiscreteSet() const;
	void setDiscreteSet(const discrete_set& discreteSet);
	void setInitialSet(polytope::ptr polyI){
		I=polyI;
	};

private:
	discrete_set DiscreteSet; //locationID
	polytope::ptr ContinuousSet;	// flowpipe abstraction as a polytope
	polytope::ptr I; // The class ptr with the exact initial set of the flow-pipe
};

#endif /* ABSTRACT_SYMBOLIC_STATE_H_ */
