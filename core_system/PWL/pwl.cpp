#include "core_system/PWL/pwl.h"

const std::list<symbolic_states>& pwl::getPassedList() const {
	return passed_list;
}

void pwl::setPassedList(const std::list<symbolic_states>& passedList) {
	passed_list = passedList;
}

const std::list<symbolic_states>& pwl::getWaitingList() const {
	return waiting_list;
}

void pwl::setWaitingList(const std::list<symbolic_states>& waitingList) {
	waiting_list = waitingList;
}

symbolic_states pwl::WaitingList_delete_front() {
	symbolic_states s;

	s = waiting_list.front();
	waiting_list.pop_front();
	return s;
}
bool pwl::WaitingList_isEmpty() {
	if (waiting_list.empty())	//list is empty
		return true;
	else
		return false;
}

void pwl::PassedList_insert(symbolic_states s) {
	passed_list.push_back(s);
}

void pwl::WaitingList_insert(symbolic_states s) {
	waiting_list.push_back(s);
}

unsigned int pwl::get_waiting_list_size() {

	return waiting_list.size();
}
