/*
 * Location.cpp
 *
 *  Created on: 09-Jul-2014
 *      Author: amit
 */

#include "core_system/HybridAutomata/Location.h"
//#include "core_system/HybridAutomata/Transition.h"

using namespace std;

location::location() {
	Name = "";
}

location::location(int Loc_ID, string name, Dynamics system_dynamics, polytope::ptr invariant, bool inv_exists, std::list<transition::ptr> Out_Going_Trans) {
	loc_id = Loc_ID;
	Name = name;
	System_Dynamics = system_dynamics;
	Invariant = invariant;
	InvariantExists = inv_exists;
	Out_Going_Transitions = Out_Going_Trans;
}

Dynamics& location::getSystem_Dynamics() {
	return System_Dynamics;
}

void location::setSystem_Dynamics(const Dynamics& system_dynamics) {
	System_Dynamics = system_dynamics;
}

polytope::ptr location::getInvariant() {
	return Invariant;
}

void location::setInvariant(polytope::ptr invariant) {
	Invariant = invariant;
}

const string& location::getName() const {
	return Name;
}

void location::setName(const string& name) {
	this->Name = name;
}

int location::getLocId() const {
	return loc_id;
}

void location::setLocId(int locId) {
	loc_id = locId;
}

std::list<transition::ptr>& location::getOut_Going_Transitions(){
	return Out_Going_Transitions;
}
void location::add_Out_Going_Transition(transition::ptr t){
	this->Out_Going_Transitions.push_back(t);
		//Adj_Transitions.max_size()		//returns the size of the adjacent transitions/locations
}

transition::ptr location::getTransition(int trans_id){
	transition::ptr temp;
	std::list<transition::ptr>::iterator it;
	for (it=Out_Going_Transitions.begin(); it != Out_Going_Transitions.end();it++){
		//cout <<"(*it)->getTransitionId()= "<<(*it)->getTransitionId()<<endl;
		int transID = (*it)->getTransitionId();
		if (transID==trans_id){
			temp = (*it);
		}
	}
	return temp;
}

bool location::isInvariantExists() const {
	return InvariantExists;
}

void location::setInvariantExists(bool invariantExists) {
	InvariantExists = invariantExists;
}


