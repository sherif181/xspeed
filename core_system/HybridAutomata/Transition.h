/*
 * Transition.h
 *
 *  Created on: 09-Jul-2014
 *      Author: amit
 */

#ifndef TRANSITION_H_
#define TRANSITION_H_

#include "core_system/HybridAutomata/DataStructureHybridAutomata.h"
#include "core_system/continuous/Polytope/Polytope.h"
#include <boost/shared_ptr.hpp>

class transition {
	int trans_id;
	string label;
	int source_location_id;
	int destination_location_id;
	polytope::ptr Gaurd;
	Assign Assign_T;
public:
	typedef boost::shared_ptr<transition> ptr;
	transition();
	transition(int trans_id, string label, int source_id, int destination_id,
			polytope::ptr Gaurd, Assign& Assign_T);
	Assign& getAssignT();
	void setAssignT(Assign assignT);
	int getDestination_Location_Id();
	void setDestination_Location_Id(int dest_loc_id);
	polytope::ptr getGaurd();
	void setGaurd(polytope::ptr gaurd);
	const string& getLabel() const;
	void setLabel(const string& label);
	int getSource_Location_Id();
	void setSource_Location_Id(int source_loc_id);
	int getTransitionId() const;
	void setTransitionId(int transId);
};

#endif /* TRANSITION_H_ */
