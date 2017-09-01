/*
 * simplex1.hpp
 *
 *  Created/Modified to #define-header on: 05-Jan-2016
 *      Author: amit
 */

#ifndef SIMPLEX1_HPP_
#define SIMPLEX1_HPP_

//#include "simplex.h"
#include <climits>
#include <iomanip>
#include <set>

/*
 bool maybe_equal(double x, double y){
 double tol = 1e6; // precision op to 6 decimal places
 double xt= 1= x*tol;
 double yt = y*tol;
 } */

template<typename T>
simplex<T>::simplex() {
	dim = 0;
}
template<typename T>
simplex<T>::~simplex() {
	//dtor
}
template<typename T>
void simplex<T>::pivot(int e, int lv) {
//	std::cout << "pivot called with entering var =" << e << "leaving var = "<< lv << std::endl;
	// Compute the coefficeints of the equation for new basic variables x_e
	math::matrix<T> Anew(As.size1(), As.size2());
	std::vector<T> bnew(bs.size(), 0);
	//std::cout << "As(lv,e) = " << As(lv,e) << std::endl;

	assert(As(lv, e) != 0);
	bnew[e] = bs[lv] / As(lv, e);
//	std::cout << " bs[lv] = " << bs[lv] << " As[lv,e] = " << As(lv, e)<< "bnew[e] = " << bnew[e] << std::endl;

	typename std::set<int>::iterator it, it1;
	for (it = N.begin(); it != N.end(); it++) {
		if (*it == e)
			continue;
		int j = *it;
		Anew(e, j) = As(lv, j) / As(lv, e);
		//	std::cout<<" Anew(e, j) = " << Anew(e, j) <<std::endl;
	}
	Anew(e, lv) = 1 / As(lv, e);

	// Compute the coefficients of the remaining constraints
	for (it = B.begin(); it != B.end(); it++) {
		int i = *it;
		if (i == lv) {
			for (int entire_row = 0; entire_row <= B.size() + N.size();
					entire_row++)	// Amit
				Anew(lv, entire_row) = 0;//Amit entire row for the leaving variable coefficient set to 0
			continue;
		}
		bnew[i] = bs[i] - As(i, e) * bnew[e];
		for (it1 = N.begin(); it1 != N.end(); it1++) {
			int j = *it1;
			if (j == e) {
				Anew(i, e) = 0;	//AMIT for the same entering variable, coefficient is set to 0
				continue;
			}
			Anew(i, j) = As(i, j) - As(i, e) * Anew(e, j);
		//	std::cout << " Anew(i, j) = " << Anew(i, j) << std::endl;
		}
		Anew(i, lv) = -1 * As(i, e) * Anew(e, lv);
		Anew(i, i) = 0;	//AMIT for the same basic variable, coefficient is set to 0
	//	std::cout << " Anew(i, lv) AMIT = " << Anew(i, lv) << std::endl;
	}

	/*std::cout << "state of the constraints matrix Amit Check here\n";
	 for (int i = 1; i < Anew.size1(); i++) {
	 for (int j = 1; j < Anew.size2(); j++)
	 std::cout << Anew(i, j) << " ";
	 std::cout << std::endl;
	 }*/

	// Compute the objective function
	obj_val = obj_val + cs[e] * bnew[e];
	//std::vector<T> cnew(dim + 1, 0);
	std::vector<T> cnew(cs.size(), 0);	//Code tried by Amit
	for (it1 = N.begin(); it1 != N.end(); it1++) {
		int j = *it1;
		if (j == e)
			continue;
		/*std::cout << "cs at " << j << " = " << cs[j] << std::endl;
		 std::cout << "cs at e =" << cs[e] << std::endl;
		 std::cout << "As at e,j = " << As(e,j) << std::endl; */
		cnew[j] = cs[j] - cs[e] * Anew(e, j);
	}
	cnew[lv] = -1 * cs[e] * Anew(e, lv);
	// Compute new sets of basic and nonbasic variables
	cs = cnew;
	As = Anew;
	bs = bnew;
	N.erase(e);
	N.insert(lv);
	B.erase(lv);
	B.insert(e);
}
/**
 * Converts the LP problem from Standard form to Slack form, initialises the data structures
 * to be used by the simplex algorithm. The slack form is such that the initial basic solution
 * is feasible. The function stops if the lp is infeasible or unbounded with a proper message.
 */
template<typename T>
void simplex<T>::initialize() {
	obj_val = 0;
	int M = A.size1(); // Number of constrains in the LP.
	int D = A.size2(); // Number of variables of the LP.
	dim = M + D;
	// set variable 1 to N as the non-basic variables
	for (unsigned int i = 1; i <= D; i++)
		N.insert(i);
	//set variables N+1 to M as the basic variables (the slack variables)
	for (unsigned int i = 1; i <= M; i++)
		B.insert(D + i);
	// A larger coefficient matrix, constants vector, objective vector defined
	As = math::matrix<T>(dim + 1, dim + 1);
	bs = std::vector<T>(dim + 1, 0);
	cs = std::vector<T>(dim + 1, 0);
	// populate the members of the As matrix, bs vector
	typename std::set<int>::iterator it;
	assert(B.size() == b.size());
	// As,bs,cs are initialised as 1 index arrays from A,b,c
	int i = 0;
	for (it = B.begin(); it != B.end(); it++, i++) {
		int m = *it;
		for (unsigned int j = 0; j < D; j++) {
			As(m, j + 1) = A(i, j);
		}
		bs[m] = b[i];
	}
	for (it = N.begin(); it != N.end(); it++) {
		cs[*it] = c[*it - 1];
	}
}
template<typename T>
void simplex<T>::display_state() {
	std::cout << "state of the simplex :\n";
	std::set<int>::iterator it1;
	std::cout << "Non basic variables:\n";
	for (it1 = N.begin(); it1 != N.end(); it1++) {
		std::cout << *it1 << std::endl;
	}
	std::cout << "Basic variables:\n";
	for (it1 = B.begin(); it1 != B.end(); it1++) {
		std::cout << *it1 << std::endl;
	}
	std::cout << "\nSize of the Total Variable : " << cs.size() << "\n";
	std::cout << "objective function coefficients:\n";
	for (int i = 1; i < cs.size(); i++)
		std::cout << cs[i] << std::endl;

	std::cout << "constants:\n";
	for (int i = 1; i < bs.size(); i++)
		std::cout << bs[i] << std::endl;

	std::cout << "state of the constraints matrix\n";
	for (int i = 1; i < As.size1(); i++) {
		for (int j = 1; j < As.size2(); j++)
			std::cout << As(i, j) << " ";
		std::cout << std::endl;
	}
	std::cout << "Objective function value at current iteration = " << obj_val
			<< std::endl;
	std::cout
			<< "============================================================= \n";
}

template<typename T>
std::vector<T> simplex<T>::solve() {
	int n = N.size() + B.size();
	std::vector<T> x(n, 0);
	/*	//Amit:: Debug
	 for (unsigned int i = 1; i <= n; i++) {
	 std::cout << "Debug Amit :: delta at Creation " << i << " = " << delta[i] << std::endl;
	 }*/
	int e, l, iters = 0;
	//typename std::set<int>::iterator it, it1,it2;

	// start solving with simplex algorithm
	while (iters <= max_iters) {
		std::vector<T> delta(n + 1, INT_MAX);//AMIT : this declaration should be inside the loop to see the INT_MAX changes reflect
		e = -1;
		for (std::set<int>::iterator it = N.begin(); it != N.end(); it++) {
			if (cs[*it] <= 0)
				continue;
			e = *it;
			//std::cout << "Adding entering var =" << e << ", as cs[e] = " <<  cs[e] << std::endl;
			break;// after choosing the entering variable, find the leaving variable
		}
		if (e == -1) {
			/*
			 for(IT = B.begin(); IT!=B.end();IT++){
			 x[*IT] = bs[*IT];
			 }*/
			//		std::cout << "Landed inside exit condition\n";
			break;// break from the while loop;
		} else {
		//	std::cout << "The entering variable is " << e << std::endl;
			for (std::set<int>::iterator it2 = B.begin(); it2 != B.end();
					it2++) {
				int i = *it2;
				if (As(i, e) > 0) {
					delta[i] = bs[i] / As(i, e);
					/*std::cout << "delta at var " << i << " = " << delta[i]
					 << std::endl;
					 std::cout << "bs at " << i << " = " << bs[i] << std::endl;
					 std::cout << "As at i,e" << "= " << As(i, e) << std::endl;*/
				}
				//std::cout << "Amit: Coefficient at " << i << " = " << As(i, e)<<"\n";
			}

			double min = INT_MAX;
			for (unsigned int i = 1; i <= n; i++) {
				//Amit:: Debug
				//	std::cout << "Debug Amit :: delta at var " << i << " = " << delta[i] << std::endl;

				if (delta[i] < min) {
					min = delta[i];
					l = i;
				}
			}
		//	std::cout << "min delta = " << min << std::endl;
		//	std::cout << "leaving variable = " << l << std::endl;
			if (min == INT_MAX) {
				std::cout << "Unbounded LP";
				exit(0);
				//throw 10;
				//return;
			} else {
			//	std::cout << "simplex state BEFORE calling pivot from solve:\n";
				//	display_state();

				pivot(e, l);
			//	std::cout << "simplex state AFTER calling pivot from solve:\n";
				//	display_state();
			}
			/*
			 std::cout << "simplex state after pivoting\n"  << std::endl;
			 display_state();
			 */
			iters++; // one simplex iteration complete
		}
	} // end of while
	return x;
}

template<typename T>
void simplex<T>::process_lp() {
	unsigned int k, m, n;
	double min = INT_MAX;
	for (unsigned int i = 0; i < b.size(); i++) {
		if (b[i] < min) {
			min = b[i];
			k = i;
		}
	}
	if (b[k] >= 0) {
		initialize();
		//	std::cout << "LP processed\n";
		return;
	}
	k = k + 1; // since we start indexing arrays at 1 in our simplex implementation
	m = A.size1(); // no. of constraints of the LP
	n = A.size2(); // no. of variables in the LP
//	std::cout << "k= " << k << " m =" << m << ", n = " << n << std::endl;

	std::vector<T> x(m + n, 0);
	bool flag = false;
	math::matrix<T> Anew(m, n + 1); // no. of vars is with the added extra variable of the auxiliary LP.
	c.resize(n + 1);
	for (unsigned int i = 0; i < m; i++) {
		for (unsigned int j = 0; j < n; j++) {
			Anew(i, j) = A(i, j);
			if (!flag) // if c not yet changed to the obj of L_aux, then change it
				c[j] = 0;
		}
		Anew(i, n) = -1;
		flag = true;
	}
	c[n] = -1;
	A = Anew;
	/* debug purpose
	 std::cout << "Anew matrix\n" << std::endl;
	 for (unsigned int i = 0; i < A.size1(); i++) {
	 for (unsigned int j = 0; j < A.size2(); j++) {
	 std::cout << A(i, j) << " ";
	 }
	 std::cout << std::endl;
	 }
	 std::cout << "elements of cnew\n" << std::endl;
	 for (unsigned int i = 0; i < c.size(); i++) {
	 std::cout << c[i] << std::endl;
	 }
	 std::cout << "elements of bnew\n" << std::endl;
	 for (unsigned int i = 0; i < b.size(); i++) {
	 std::cout << b[i] << std::endl;
	 }
	 // */

	initialize(); // get the slack form
//	std::cout << "Printing After Initialize\n" << std::endl;

	int lv = N.size() + k; // set the leaving variable
	int e = N.size(); // the entering variable is the var added in the Auxiliary LP.
/*
	std::cout << "leaving  = " << lv << ", entering = " << e << std::endl;
	std::cout << "Before pivot Call\n";
	display_state();
*/
	pivot(e, lv);	//Data Mismatch for matrix A detected here
	//display_state();
//	std::cout << "After first pivot successful\n";
//	display_state();

	solve(); // run simplex on the previously obtained slack form. it is guaranteed that the initial basic solution is feasible.

//	std::cout << "process_lp: reached after call to solve auxiliary LP"<< std::endl;
//	display_state();

	if (get_obj_val() == 0) { // To be changed to check the value of x[e], e being the index of the auxiliary variable
		//if(x[lv] == 0) // the aux variable of the lp_aux has value 0 in the optimal vector
		int en;
		//std::cout<<"Amit:: x0 is lv =  "<< lv <<" or e = "<< e <<std::endl;
		for (std::set<int>::iterator it = B.begin(); it != B.end(); it++) {
			if (*it == e) { // aux variable as basic variable
				for (std::set<int>::iterator it1 = N.begin(); it1 != N.end();
						it1++) {
					if (As(e, *it1) != 0) {
						en = *it1;
						break;
					}
				}
				pivot(en, e);
			}
		}

	//	std::cout << "STATE BEFORE RESTORING TO ORIGINAL OBJECTIVE FUNCTION\n";
		//	display_state();

		//Anew = math::matrix<T>(As.size1() - 1, As.size2() - 1);
		Anew.resize(As.size1() - 1, As.size2() - 1); //AMIT resizing the array

		std::vector<T> bnew(bs.size() - 1, 0);
		std::vector<T> cnew(cs.size() - 1, 0);
		// remove the row and column for the auxiliary variable from As
		int row = 0;
		for (unsigned int i = 1; i < As.size1(); i++) {
			if (i == e) {
				continue;
			}
			row++;
			int k = 1;
			for (unsigned int j = 1; j < As.size2(); j++) {
				if (j != e) {
					Anew(row, k) = As(i, j);
					k++;
				}
			}
		}

		/*
		 //AMIT debug
		 std::cout<<"debugging \n";
		 for (int x=1;x<Anew.size1();x++){
		 for (int y=1;y<Anew.size2();y++)
		 std::cout<<"Anew(x,y) = "<<Anew(x,y)<<"\t";
		 std::cout<<std::endl;
		 }
		 */

		int k = 1;
		for (unsigned int j = 1; j < As.size2(); j++) {
			if (j != e) {
				bnew[k] = bs[j];
				k++;
			}
		}
		As = Anew;
		bs = bnew;
		T v = T(0);
		T alpha;

		for (unsigned int i = 0; i < c_orig.size(); i++) {
			cnew[i + 1] = c_orig[i];
		}

	//	std::cout << "STATE AFTER RESTORING THE ORIGINAL OBJECTIVE FUNCTION\n";
		//	display_state();

		//debug
//		std::cout << "values at cnew\n";
//		for(unsigned int i=1;i<cnew.size();i++)
//		{
//				std::cout << cnew[i] << std::endl;
//		}
		//---
		//std::cout << "Value of new objective functions = " << std::endl;
		for (unsigned int i = 1; i < cnew.size(); i++) {
			if (cnew[i] != 0) {
				for (typename std::set<int>::iterator iter = B.begin();
						iter != B.end(); iter++) {
					if (*iter == i) {
						alpha = cnew[i];
						cnew[i] = 0;
						for (unsigned int j = 1; j < As.size2(); j++) {
							cnew[j] += -As(i, j) * alpha;
							//std::cout << "Amit::" << cnew[j] << std::endl;
						}
						//std::cout << "end of new obj function\n" << std::endl;
						v += bs[i] * alpha;
					}
					iter = B.end();
					break;
				}
			}
		}

		//AMIT debug
		/*std::cout << "values at cnew\n";
		 for (unsigned int i = 1; i < cnew.size(); i++) {
		 std::cout << cnew[i] << std::endl;
		 }*/
		cs = cnew;
		/*  Debug::
		 //AMIT: Now resize cs as size - 1 for removing auxiliary variable
		 std::cout<<"Amit Gurung :: cs Size BEFORE = "<<cs.size()<<"\n";
		 cs.resize(cs.size() - 1);
		 std::cout<<"Amit Gurung :: cs Size AFTER = "<<cs.size()<<"\n";
		 */

		/*
		 //AMIT debug
		 std::cout << "values at cs after copying\n";
		 for(unsigned int i=1;i<cs.size();i++)
		 {
		 std::cout << cs[i] << std::endl;
		 }
		 */
		//std::cout << "Value of new objective function = " << v << std::endl;
		obj_val = v;

		/*
		 N.clear();
		 B.clear();
		 int D = c_orig.size();
		 for (unsigned int i = 1; i <= D; i++)
		 N.insert(i);

		 for (unsigned int i = 1; i <= A_orig.size1(); i++)
		 B.insert(D + i);
		 */

		//Code tried by :: AMIT
		std::set<int> tempN, tempB; // .
		tempN = N;
		tempB = B;
		N.clear();
		B.clear();

		for (typename std::set<int>::iterator iter = tempN.begin();
				iter != tempN.end(); iter++) {
			if (*iter == e) {	//e is the x0 variable or the auxiliary variable
				continue;	//skiping or removing the auxiliary variable
			}
			//Inserting/Renaming the remaining variable but if the variables are greater than auxiliary variable name/index subtracted by 1
			if (*iter > e)
				N.insert(*iter - 1);		//requires variable renaming
			else
				N.insert(*iter);		//Does Not requires variable renaming
		}

		for (typename std::set<int>::iterator iter = tempB.begin();
				iter != tempB.end(); iter++) {
			//Inserting/Renaming the variable but if the variables are greater than auxiliary variable name/index subtracted by 1
			if (*iter > e)
				B.insert(*iter - 1);		//requires variable renaming
			else
				B.insert(*iter);		//Does Not requires variable renaming
		}

	//	std::cout << "LP processed\n";
		//	display_state();
	//	std::cout<< "\n===============Initialize Simplex Call Over=====================\n";
	} else {
		std::cout << "Simplex::process_lp: LP has infeasible solution\n";
		exit(0);
	}
}




#endif /* SIMPLEX1_HPP_ */
