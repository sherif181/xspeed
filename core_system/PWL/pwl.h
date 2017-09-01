/*
 * pwl.h
 *
 *  Created on: 12-Nov-2014
 *      Author: amit
 */

#ifndef PWL_H_
#define PWL_H_
/*
 * Passed-Waiting-List of symbolic_states. A class that handles Queue data_structure of sysmbolic_states
 */

#include "core_system/symbolic_states/symbolic_states.h"
#include <list>

class pwl {
public:
	const std::list<symbolic_states>& getPassedList() const;
	void setPassedList(const std::list<symbolic_states>& passedList);
	const std::list<symbolic_states>& getWaitingList() const;
	void setWaitingList(const std::list<symbolic_states>& waitingList);

	//deletes a symbolic state from the front of the list and returns the deleted symbolic state
	symbolic_states WaitingList_delete_front();

	//inserts a symbolic state at the end of the passed_list
	void PassedList_insert(symbolic_states s);

	//inserts a symbolic state at the end of the waiting_list
	void WaitingList_insert(symbolic_states s);
	//return TRUE if the waiting_list is empty otherwise FALSE
	bool WaitingList_isEmpty();
	// Get the size of the waiting list, i.e. the number of elements in the list
	unsigned int get_waiting_list_size();

private:
	std::list<symbolic_states> waiting_list;
	std::list<symbolic_states> passed_list;
	unsigned int waiting_list_size;	//number of elements currently present in the waiting_list
};

#endif /* PWL_H_ */
