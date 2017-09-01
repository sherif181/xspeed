/*
 * LockAvoidanceUtility.h
 *
 *  Created on: 24-Feb-2016
 *      Author: amit
 */

#ifndef LOCKAVOIDANCEUTILITY_H_
#define LOCKAVOIDANCEUTILITY_H_

#include <vector>
#include "core_system/PWL/pwlist.h"
#include "core_system/symbolic_states/initial_state.h"
#include "core_system/Reachability/reachDataStructure.h"

// QpwList[t][list] for read and QpwList[t][list] for write
//Returns true if the entire vector of pwlist is empty otherwise false
bool isEmpty_Qpw_list(std::vector<pwlist::ptr> Qpw_list);

//Returns the size of the Qpw_list
unsigned int getSize_Qpw_list(std::vector<pwlist::ptr> Qpw_list);

//Returns a vector of all the waiting intial_state to be processed for flowpipe computation
//std::vector<initial_state::ptr> getAllpw_list(std::vector<pwlist::ptr>& Qpw_list,unsigned int size, pwlist::ptr& allPassedList);
std::vector<initial_state::ptr> getAllpw_list(std::vector<std::vector<pwlist::ptr> > & Qpw_list, int t, unsigned int size,
		pwlist::ptr& allPassedList);

//Compute the total number of directions for polytope X0 and U in(including) all the symbolic states
void getCountTotal(std::vector<LoadBalanceData>& LoadBalanceDS,
		unsigned int& countTotal_X, unsigned int& countTotal_U);

//Given the index i, returns "the index of the Symbolic State" and
// "the index of the direction for listDir_X0 within the symbolic state" as a reference variable
void search_SymState_dirsX0Index(unsigned int i,
		std::vector<LoadBalanceData>& LoadBalanceDS, int& SymStateIndex,
		unsigned int& dirsIndex);


//Given the index i, returns "the index of the SFM" and
// "the index of the SFM's column within each SFM" as a reference variable
void search_sfmIndex_colIndex(unsigned int i,
		std::vector<LoadBalanceData_PostD>& LD_post_D, int& sfmIndex,
		unsigned int& colIndex);

//Given the index i, returns "the index of the Symbolic State" and
// "the index of the direction for listDir_U within the symbolic state" as a reference variable
void search_SymState_dirsUIndex(unsigned int i,
		std::vector<LoadBalanceData>& LoadBalanceDS, int& SymStateIndex,
		unsigned int& dirsIndex);

#endif /* LOCKAVOIDANCEUTILITY_H_ */
