/*
 * userOptions.h
 *
 *  Created on: 09-May-2016
 *      Author: rajarshi
 */

#ifndef USEROPTIONS_H_
#define USEROPTIONS_H_

#include <cstring>
#include <iostream>

class userOptions {
	std::string model_filename; // filename of the automata model .xml file
	std::string config_filename; // filename of the configuration file .cfg
	std::string forbidden_state; // the string of forbidden state description
	unsigned int output_var_X; // first  dimension for plotting
	unsigned int output_var_Y; // second dimension for plotting
	unsigned int output_var_Z; // third dimension for plotting

	unsigned int model; // name of the pre-defined model to run for reachability.
	unsigned int direction_template; // template used for approximating support functions
	unsigned int time_horizon; // time horizon for reachability
	double time_step; // the time step of the support function algorithm
	unsigned int level; // the breadth level in bfs to stop reachability
//	unsigned int flow_algorithm; // Choice of the reachability algorithm
//	unsigned int automata_exploration_algorithm; // choice of algorithm for exploration of the graph

	unsigned int algo;	//Common arg for all types of Algorithm
	unsigned int total_slice_size;	//total number of partition-size or number of slices
	unsigned int stream_size;	//total number of streams selected for GPU streaming

public:
	userOptions();
	virtual ~userOptions();

	std::string get_modelFile();
	void set_modelFile(std::string modefile);
	std::string get_configFile();
	void set_configFile(std::string configfile);
	unsigned int get_first_plot_dimension();
	void set_first_plot_dimension(unsigned int outdim);
	unsigned int get_second_plot_dimension();
	void set_second_plot_dimension(unsigned int outdim);
	unsigned int get_third_plot_dimension();
	void set_third_plot_dimension(unsigned int outdim);

	double get_timeStep();
	void set_timeStep(double t);
	double get_timeHorizon();
	void set_timeHorizon(double timeHorizon);
	unsigned int get_model();
	void set_model(unsigned int m);
	unsigned int get_directionTemplate();
	void set_directionTemplate(unsigned int d);
	unsigned int get_bfs_level();
	void set_bfs_level(unsigned int l);
/*	unsigned int get_flow_algorithm();
	void set_flow_algorithm(unsigned int alg);
	unsigned int get_automata_exploration_algorithm();
	void set_automata_exploration_algorithm(unsigned int exp_alg);*/
	std::string get_forbidden_state();
	void set_forbidden_state(std::string);

	unsigned int get_algorithm();	//returns the selected Algorithm
	void set_algorithm(unsigned int alg);	//assigns the Algorithm selected by the user
	unsigned int getStreamSize() const;
	void setStreamSize(unsigned int streamSize);
	unsigned int getTotalSliceSize() const;
	void setTotalSliceSize(unsigned int totalSliceSize);
};

#endif /* USEROPTIONS_H_ */
