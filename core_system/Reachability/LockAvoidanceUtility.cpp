/*
 * LockAvoidanceUtility.cpp
 *
 *  Created on: 24-Feb-2016
 *      Author: amit
 */

#include "LockAvoidanceUtility.h"
#include "assert.h"

using namespace std;

// QpwList[t][list] for read and QpwList[t][list] for write
//Returns true if the entire vector of pwlist is empty otherwise false
bool isEmpty_Qpw_list(std::vector<pwlist::ptr> Qpw_list) {
	bool emptyFlag = true;
	for (int i = 0; i < Qpw_list.size(); i++) {
		emptyFlag = Qpw_list[i]->isEmpty_WaitingList();
		if (!emptyFlag)
			break;
	}
	return emptyFlag;
}

unsigned int getSize_Qpw_list(std::vector<pwlist::ptr> Qpw_list) {
	unsigned int list_size=0, total_size = 0;
	for (int i = 0; i < Qpw_list.size(); i++) {
		list_size = Qpw_list[i]->getWaitingListSize();
	//	cout<<"list_size = " << list_size<<std::endl;
		total_size = total_size + list_size;
	}
	return total_size;
}



std::vector<initial_state::ptr> getAllpw_list(std::vector<std::vector<pwlist::ptr> >& Qpw_list, int t, unsigned int size,
		pwlist::ptr& allPassedList) {
//std::vector<initial_state::ptr> getAllpw_list(std::vector<pwlist::ptr>& Qpw_list, unsigned int size, pwlist::ptr& allPassedList) {
	std::vector<initial_state::ptr> allList(size); //size is required for static allocate the allList vector size
	unsigned int vectorSize = Qpw_list[t].size();
	unsigned int index = 0;

	for (int i = 0; i < vectorSize; i++) { //verify by printing  Qpw_list.size()
		//each Qpw_list[i] is a pwlist ... so have to iterate
		unsigned int listSize = Qpw_list[t][i]->getWaitingListSize();
		for (int j = 0; j < listSize; j++) {
			allList[index] = Qpw_list[t][i]->WaitingList_delete_front();//POP operation
			allPassedList->PassedList_insert(allList[index]);	//PUSH operation
			index++;
		}
	}
	//Qpw_list[t].clear();
	assert(size == index);
	return allList;
}

void getCountTotal(std::vector<LoadBalanceData>& LoadBalanceDS,
		unsigned int& countTotal_X, unsigned int& countTotal_U) {
	unsigned int countTot = 0, countTot2 = 0;
	for (int i = 0; i < LoadBalanceDS.size(); i++) { //for each symbolic-states
		countTot = countTot + LoadBalanceDS[i].List_dir_X0.size1();
		countTot2 = countTot2 + LoadBalanceDS[i].List_dir_U.size1();
	}
	countTotal_X = countTot;
	countTotal_U = countTot2;
}

//Running time complexity is order of the "number of symbolic states per breadth-level
void search_SymState_dirsX0Index(unsigned int i,
		std::vector<LoadBalanceData>& LoadBalanceDS, int& SymStateIndex,
		unsigned int& dirsIndex) {
	unsigned int tot = 0;
	for (int id = 0; id < LoadBalanceDS.size(); id++) {
		unsigned int siz = LoadBalanceDS[id].List_dir_X0.size1(); //returns the number
	//	cout<<"siz = "<<siz<<"\t";
		tot = tot + siz;
		if (i < tot) { //Note <= is incorrect as ==tot will fall in the next symbolic_state, due to 0-indexing
			SymStateIndex = id;
			unsigned int term =(i - (tot - siz)); //(i - tot - siz - 1);
			//cout<<"i  = "<<i <<"  tot = " << tot<<"  siz = "<<siz <<"  term = "<<term<<"\n";
			dirsIndex = term;
			/*if (term == 0)
				dirsIndex = 0;
			else
				dirsIndex = term - 1;	//-1 due to 0-indexing*/
			//cout<<"i  = "<<i <<"  tot = " << tot<<"  siz = "<<siz <<"  term = "<<term<<"   dirsIndex = "<<dirsIndex<<"\n";
			break; //when found stop searching the rest of the symbolic state
		}
	}
}


//Running time complexity is order of the "number of symbolic states per breadth-level
void search_sfmIndex_colIndex(unsigned int i,
		std::vector<LoadBalanceData_PostD>& LD_post_D, int& sfmIndex,
		unsigned int& colIndex) {
	unsigned int tot = 0;
	for (int id = 0; id < LD_post_D.size(); id++) {
		unsigned int siz = LD_post_D[id].sfm->getTotalIterations();
		tot = tot + siz;
		if (i < tot) { //Note <= is incorrect as ==tot will fall in the next symbolic_state, due to 0-indexing
			sfmIndex = id;
			unsigned int term =(i - (tot - siz));
			colIndex = term;
			break; //when found stop searching the rest of the symbolic state
		}
	}
}



void search_SymState_dirsUIndex(unsigned int i,
		std::vector<LoadBalanceData>& LoadBalanceDS, int& SymStateIndex,
		unsigned int& dirsIndex) {
	unsigned int tot = 0;
	for (int id = 0; id < LoadBalanceDS.size(); id++) {
		unsigned int siz = LoadBalanceDS[id].List_dir_U.size1(); //returns the number
		tot = tot + siz;
		if (i < tot) { //Note <= is incorrect as ==tot will fall in the next symbolic_state, due to 0-indexing
			SymStateIndex = id;
			unsigned int term = (i - (tot - siz)); //(i - tot - siz - 1);
			if (term == 0)
				dirsIndex = 0;
			else
				dirsIndex = term - 1;	//-1 due to 0-indexing
			///dirsIndex = (i - tot - siz - 1); //-1 due to 0-indexing
			break;//when found stop searching the rest of the symbolic state
		}
	}
}
