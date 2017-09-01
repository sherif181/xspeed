/*
 * vartoindexmap.cpp
 *
 *  Created on: 07-Feb-2016
 *      Author: rajarshi
 */

#include <core_system/HybridAutomata/vartoindexmap.h>
#include <iostream>

var_to_index_map::map_ptr var_to_index_map::var_index_map_ptr =
		var_to_index_map::map_ptr(new std::map<std::string,unsigned int>());

var_to_index_map::var_to_index_map() {
	// TODO Auto-generated constructor stub
}

var_to_index_map::~var_to_index_map() {
	// TODO Auto-generated destructor stub
}

void var_to_index_map::print_var_index_map()
{
	unsigned int i = 0;
	std::cout << "The variable to index map is:\n";
	for(std::map<std::string, unsigned int>::iterator it = var_index_map_ptr->begin(); it!=var_index_map_ptr->end();it++){
		std::cout << "Variable = " << (*it).first;
		std::cout << " Value = " << (*it).second << std::endl;
		i++;
	}
}

/**
 * Returns the size of the map, i.e., the number of variables of the map
 */
unsigned int var_to_index_map::map_size()
{
	return var_index_map_ptr->size();
}
