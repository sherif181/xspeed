/*
 * HybridAutomata.cpp
 *
 *  Created on: 09-Jul-2014
 *      Author: amit
 */

#include "Hybrid_Automata.h"

using namespace std;

hybrid_automata::hybrid_automata() {
	dimension = 0;
}
hybrid_automata::hybrid_automata(std::map<int, location::ptr>& list_locs, location::ptr init_loc,
		int dim) {
	list_locations = list_locs;	//assigning a map to another map
	initial_loc = init_loc;
	dimension = dim;
}

location::ptr& hybrid_automata::getInitial_Location() {
	return initial_loc;
}

void hybrid_automata::addInitial_Location(location::ptr& initLoc) {
	initial_loc = initLoc;
}

location::ptr hybrid_automata::getLocation(int Loc_Id){
	//iterate through the transition and find for given input Name that match with the source_location
	assert(list_locations.count(Loc_Id)!=0);
	location::ptr l;
	l = list_locations[Loc_Id];
//	std::cout<<"\ntest 3\n";
	return l;
}

int hybrid_automata::getDimension() const {
	return dimension;
}

void hybrid_automata::setDimension(int dim) {
	this->dimension = dim;
}

void hybrid_automata::addMapped_Locations_List(std::map<int, location::ptr>& mapped_location_list){
	list_locations = mapped_location_list;
}
void hybrid_automata::addLocation(location::ptr& loc){
	int key = loc->getLocId();
	list_locations[key] = loc;	//storing the loc with the proper loc_id as the key
}
