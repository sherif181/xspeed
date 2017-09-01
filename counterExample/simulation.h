/*
 * simulation.h
 *
 *  Created on: 12-Feb-2016
 *      Author: rajarshi
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_
#include <vector>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include "application/DataStructureDirections.h"
#include "Utilities/gradient.h"

/*
 * A structure to return the last point in Invariant and the
 * time of crossing the invariant
 */
struct bound_sim
{
	std::vector<double> v;
	double cross_over_time;
};

/**
 * This class provides methods to simulate an ODE equation for a given
 * initial value.
 */
class simulation : public var_to_index_map {

	/** The number of discrete samples to be computed
	 * in solving the ODE, in order to get the simulation. */

	unsigned int N;

	/**
	 * The dimension of the ODE system.
	 */
	unsigned int dimension;
	Dynamics D;
	double reltol;
	double abstol;
	string filename;
	unsigned int x1; // the first output dimension for plotting.
	unsigned int x2; // the second output dimension for plotting.

public:
	typedef boost::shared_ptr<simulation> ptr;
	/** To store the end point of the simulation trajectory and its distance
	 *  from a given polytope.
	 */
//	typedef std::pair<std::vector<double>, double> simD;

	simulation();
	simulation(unsigned int dim, unsigned int steps, Dynamics Dyn, double rel_tol=1e-8, double abs_tol=1e-8){
		dimension = dim;
		N = steps;
		reltol = rel_tol;
		abstol = abs_tol;
		D = Dyn;
		filename=std::string();
		// default ploting dimension

		x1 = 0; // The default plotting of this dimension
		x2 = 1; // The default plotting of this dimension

	}
	virtual ~simulation();
	/**
	 * Sets the name of the output file to the parameter string.
	 * The parameter string should be the absolute path.
	 * The simulation shall be printed to this file if filename is not
	 * empty
	 */
	void set_outfile(std::string s){
		filename = s;
	}
	/**
	 * Returns the dimension of the simulation object.
	 */
	unsigned int get_system_dimension()
	{
		return dimension;
	}
	/**
	 * sets the projection dimensions to output the simulation points
	 * in a file
	 */
	void set_out_dimension(unsigned int i, unsigned int j){
		x1 = i;
		x2 = j;
	}
	/*
	 * Set the number of simulation samplings or steps
	 */
	void set_steps(unsigned int n)
	{
		N = n;
	}
	/**
	 * Generates a simulation trace for time duration, starting at start_time.
	 * The initial state is given by the first parameter
	 */
	std::vector<double> simulate(std::vector<double>, double time);

	/**
	 * Generates a simulation trace for time duration, starting at start_time.
	 * The time instant, within the simulation time, when the polytope I is
	 * violated by the trace is returned and with the first lsimulation point
	 * that violated I, as a struct object. status is set to false if invariant
	 * is violated.
	 */
	bound_sim bounded_simulation(std::vector<double>, double time, polytope::ptr I, bool &status);

	/**
	 * Simulate and also compute the distance of the trajectory with a polytope,
	 * gradient of the distance of the trace to the given polytope invariant I w.r.t
	 * the trace end point and the gradient of the distance of the trace
	 * to the invariant w.r.t time.
	 */
	std::vector<double> metric_simulate(std::vector<double> x, double time, double& distance, polytope::ptr Inv,
			std::vector<double>& grad);

};

#endif /* SIMULATION_H_ */
