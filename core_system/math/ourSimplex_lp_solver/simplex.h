#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <vector>
#include <set>
#include "core_system/math/matrix.h"
#include <climits>

/**
 A templated class simplex which implements the basic simplex algorithm as defined on the
 Introduction to algorithms book by Cormen et. al.
 The template type is kept for possible accurate computation with Rational type

 Usages sequence :
 	 	 simplex<double> s;
  	  	 s = simplex<double>(A, b, obj_dir);
		 s.process_lp();
		 r = s.solve();
		 std::cout << "Output from Simplex = " << s.get_obj_val();

 */

template<typename T>
class simplex {
public:
	simplex();
	/**
	 * Ap*x <= b provides the constraint space in standard form and c is the objective function direction.
	 * It is assumed in the implementation that the Ap,bp,cp are 0 indexed arrays.
	 */
	simplex(const math::matrix<T> Ap, const std::vector<T> bp,
			const std::vector<T> cp) {
		A = math::matrix<T>(Ap);
		b = std::vector<T>(bp);
		c = std::vector<T>(cp);

		assert(A.size2() == c.size());
		assert(b.size() == A.size1());

		N = std::set<int>(); // empty set
		B = std::set<int>(); // empty set
		dim = N.size() + B.size();
		max_iters = INT_MAX;
		// save the original LP copy
		A_orig = math::matrix<T>(Ap);
		b_orig = std::vector<T>(bp);
		c_orig = std::vector<T>(cp);

	}
	;
	/**
	 * converts the LP into a slack form such that the initial basic solution is feasible. This function also
	 * notifies if the given LP is infeasible or unbounded and terminates.
	 */
	void process_lp();
	/**
	 * Generates a slack form assuming that the LP has the initial basic solution as feasible.
	 */
	void initialize();
	/** modifies the A,B,N,c,b once the entering and the leaving variable is passed to the function */
	void pivot(int e, int lv);
	/** the simplex algorithm. Returns the optimal vector. The optimal value of the objective
	 *  function is stored in var obj_val.
	 *  solve() runs assuming that the slack form of the LP has initial basic solution to be feasible.
	 *  For the above assumption, a call to solve must be prepended with a call to  process_lp*/
	std::vector<T> solve();
	/**
	 * Return the optimal objective function value
	 */
	T get_obj_val() {
		return obj_val;
	}
	/**
	 * Displays the state of the simplex data structures
	 */
	void display_state();
	virtual ~simplex();
protected:
private:
	int dim; // the number of variables of the linear programming problem, including the slack variables.
	int max_iters; // stops the simplex iteration after max_iters times.
	std::vector<T> c; // c is the vector of objective function coefficients.
	std::vector<T> b; // b is the vector of constants in the linear equalities in the standard form.
	T obj_val; // Keeps track of the current objective function value
	std::set<int> N, B; // N, B  is the set of indices of the non-basic and basic variables respectively.
	math::matrix<T> A, As, A_orig; // matrix A stores the LP constraint coefficients in standard form.
	// matrix As stores the intermediate matrix used by the simplex algorithm.
	// matrix A_orig storing the original constraints of the LP.

	std::vector<T> bs, b_orig; // bs stores the intermediate vector of constant terms used
							   //  b_orig stores the original LP constraint constant terms

							   // by the simplex algorithm
	std::vector<T> c_orig, cs; //  cs stores the intermediate vector of objective function
	// c_orig stores the original objective function vector

};

#include "simplex.hpp"
#endif // SIMPLEX_H
