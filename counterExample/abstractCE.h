/*
 * abstractCE.h
 *
 *  Created on: 12-Jan-2016
 *      Author: rajarshi
 */

#ifndef ABSTRACTCE_H_
#define ABSTRACTCE_H_


#include "core_system/continuous/ConvexSet/supportFunctionProvider.h"

#include "core_system/HybridAutomata/Transition.h"
#include "core_system/HybridAutomata/Hybrid_Automata.h"
#include <list>
#include <boost/shared_ptr.hpp>
#include "nlopt.hpp"
#include "counterExample/abstract_symbolic_state.h"
#include "counterExample/concreteCE.h"
#include "nlpFunctions.h"


#include <fstream>
#include <string>



/**
 * This class is a data structure to store the abstract counter-example generated
 * by XSpeed when given an unsafe symbolic state. An abstract counter example
 * is an ordered list of symbolic states starting from an initial and ending in an
 * unsafe symbolic symbolic state.
 *
 * @author: Rajarshi
 */
//using NEWMAT::ColumnVector;
extern unsigned int N;
extern unsigned int dim;
extern hybrid_automata::ptr HA;
extern std::vector<int> locIdList;
extern std::list<transition::ptr> transList;
extern polytope::ptr bad_poly;
extern std::list<refinement_point> ref_pts; // a list of invariant violating points to refine the search and obtained a validated trajectory
extern std::vector<std::vector<double> > X0; // list of start point of the trajectory segments. Used only in the NLP-LP mixed program
extern std::list<abstract_symbolic_state::ptr> ce_sym_states; // list of CE abstract sym states. Used only in the NLP-LP mixed problem

class abstractCE
{
public:
	typedef boost::shared_ptr<abstractCE> ptr;

	/**empty constructor */
	abstractCE() {
	}
	;
	/* another constructor */
	abstractCE(std::list<abstract_symbolic_state::ptr> s_states,
			std::list<transition::ptr> ts, hybrid_automata::ptr h, polytope::ptr fpoly);
	/* destructor */
	~abstractCE() {
	}
	;
	const std::list<abstract_symbolic_state::ptr> get_CE_sym_states() const {
		return sym_states;
	}
	const std::list<transition::ptr> get_CE_transitions() const {
		return trans;
	}

	const abstract_symbolic_state::ptr get_first_symbolic_state() const;

	/**
	 * The semantics assumes that the last abstract_symbolic_state in the list contains the
	 * unsafe polytope
	 */
	const abstract_symbolic_state::ptr get_last_symbolic_state() const;

	/**
	 * Returns the forbidden polytope
	 */
	const polytope::ptr get_forbidden_poly(){
		return forbid_poly;
	}
	/**
	 * Returns the i-th symbolic state from the CE
	 */
	abstract_symbolic_state::ptr get_symbolic_state(unsigned int i) const;

	const unsigned int get_length() const {
		return length;
	}

	void set_length(unsigned int len) {
		length = len;
	}

	void set_sym_states(std::list<abstract_symbolic_state::ptr> sym);

	void set_transitions(std::list<transition::ptr> transitions) {
		trans = transitions;
	}
	/**
	 * Sets the reference to the hybrid automaton to which this CE refers.
	 */
	void set_automaton(hybrid_automata::ptr h){
		H = h;
	}
	/**
	 * Sets the forbidden polytope of this abstract counter example
	 */
	void set_forbid_poly(polytope::ptr fpoly){
		forbid_poly = fpoly;
	}
	hybrid_automata::ptr get_automaton(){
		return H;
	}
	/**
	 * returns a validated counterexample trace, a trace that satisfies the invariant
	 */
	concreteCE::ptr get_validated_CE(double tolerance);

	/**
	 * Plot the counter example projected along dimensions passed
	 * as parameters
	 */
	void plot(unsigned int i, unsigned int j);

private:
	/**
	 * The first symbolic state is the initial symbolic state and the last one
	 * is the unsafe symbolic state
	 */
	std::list<abstract_symbolic_state::ptr> sym_states;

	/**
	 * The list of transitions taken from the initial abstract_symbolic_state to the
	 * final abstract_symbolic_state.
	 */
	std::list<transition::ptr> trans;
	/**
	 * length shows how many discrete transitions are present in the abstract counter
	 * example.
	 */
	unsigned int length;

	/**
	 * The reference to the automaton to which this is a counter example
	 */
	hybrid_automata::ptr H;
	/**
	 * The reference to the forbidden polytope given by the user
	 */
	polytope::ptr forbid_poly;

	/**
	 * Returns an instance of the concrete counter-example from the abstract using NLP and flowpipe constraints
	 */
	concreteCE::ptr gen_concreteCE(double tolerance, const std::list<refinement_point>& refinements);

	/**
	 * Returns an instance of the concrete counter-example, if it exists, using mixed NLP-LP
	 */
	concreteCE::ptr gen_concreteCE_NLP_LP(double tolerance, const std::list<refinement_point>& refinements);
	/**
	 * Returns an instance of the concrete counter-example, if it exists, using NLP and limited HA constraints
	 */
	concreteCE::ptr gen_concreteCE_NLP_HA(double tolerance, const std::list<refinement_point>& refinements);
};


std::vector<double> simulate_trajectory(const std::vector<double>& x0,
		Dynamics& D, const double& time, double& distance, polytope::ptr I, std::vector<double>&);

#endif /* ABSTRACTCE_H_ */
