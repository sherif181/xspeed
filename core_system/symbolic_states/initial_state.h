/*
 * initial_state.h
 *
 *  Created on: 03-Dec-2015
 *      Author: amit
 */

#ifndef INITIAL_STATE_H_
#define INITIAL_STATE_H_

#include "core_system/continuous/Polytope/Polytope.h"
#include "symbolic_states.h"

class initial_state {
public:
	typedef boost::shared_ptr<initial_state> ptr;
	initial_state();
	initial_state(int locationId, polytope::ptr initialSet);
	initial_state(int locationId, polytope::ptr initialSet,
			symbolic_states::ptr parentPtrSymbolicState, int transitionId);
	~initial_state();

	const polytope::ptr getInitialSet() const;
	void setInitialSet(const polytope::ptr initialSet);
	int getLocationId() const;
	void setLocationId(int locationId);

	int getTransitionId() const;
	void setTransitionId(int transitionId);

	const symbolic_states::ptr getParentPtrSymbolicState() const;
	void setParentPtrSymbolicState(
			const symbolic_states::ptr parentPtrSymbolicState);

private:
	int location_id;
	polytope::ptr initial_set;

	/*
	 * For ease of implementation transition_id  and parent_pointer_symbolic_states is kept here
	 *
	 * transition_id:: maintains the trans_id that generates this initial_state in any particular
	 * location, except for the originating location or the start location it can be set to 0;
	 *
	 * parentPtr_symbolic_state:: maintains pointer to parent symbolic_state
	 */
	int transition_id;
	symbolic_states::ptr parentPtr_symbolic_state;
};

#endif /* INITIAL_STATE_H_ */
