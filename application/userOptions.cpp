/*
 * userOptions.cpp
 *
 *  Created on: 09-May-2016
 *      Author: rajarshi
 */

#include <application/userOptions.h>

userOptions::userOptions() {
	// TODO Auto-generated constructor stub
	model = 0; //default model, bouncing ball
	direction_template = 0; //default directions, box template
	output_var_X = 0; // default first dimension of plot
	output_var_Y = 1; // default second dimension of plot
	output_var_Z = 0; //default third plot dimension

//	automata_exploration_algorithm = 12; // sequential BFS
//	flow_algorithm = 1;	// SEQ

	stream_size =1;	//default set to 1 streams
}
std::string userOptions::get_modelFile()
{
	return model_filename;
}
void userOptions::set_modelFile(std::string modelfile)
{
	model_filename = modelfile;
}
std::string userOptions::get_configFile()
{
	return config_filename;
}
void userOptions::set_configFile(std::string configfile)
{
	config_filename = configfile;
}
unsigned int userOptions::get_first_plot_dimension()
{
	return output_var_X;
}
void userOptions::set_first_plot_dimension(unsigned int outdim)
{
	output_var_X = outdim;
}
unsigned int userOptions::get_second_plot_dimension()
{
	return output_var_Y;
}
void userOptions::set_second_plot_dimension(unsigned int outdim)
{
	output_var_Y = outdim;
}

unsigned int userOptions::get_third_plot_dimension()
{
	return output_var_Z;
}
void userOptions::set_third_plot_dimension(unsigned int outdim)
{
	output_var_Z = outdim;
}

double userOptions::get_timeStep()
{
	return time_step;
}
void userOptions::set_timeStep(double t)
{
	time_step = t;
}
double userOptions::get_timeHorizon()
{
	return time_horizon;
}
void userOptions::set_timeHorizon(double timeHorizon)
{
	time_horizon = timeHorizon;
}
unsigned int userOptions::get_model()
{
	return model;
}
void userOptions::set_model(unsigned int m)
{
	model = m;
}
unsigned int userOptions::get_directionTemplate()
{
	return direction_template;
}
void userOptions::set_directionTemplate(unsigned int d)
{
	 direction_template = d;
}
unsigned int userOptions::get_bfs_level()
{
	return level;
}
void userOptions::set_bfs_level(unsigned int l)
{
	level = l;
}
/*unsigned int userOptions::get_flow_algorithm()
{
	return flow_algorithm;
}
void userOptions::set_flow_algorithm(unsigned int alg)
{
	flow_algorithm = alg;
}
unsigned int userOptions::get_automata_exploration_algorithm()
{
	return automata_exploration_algorithm;
}
void userOptions::set_automata_exploration_algorithm(unsigned int exp_alg)
{
	automata_exploration_algorithm = exp_alg;
}*/

std::string userOptions::get_forbidden_state()
{
	return forbidden_state;
}
//void userOptions::set_forbidden_state(std::__cxx11::string forbid_s){ //creates some other error if use -std=c++11 options
void userOptions::set_forbidden_state(std::string forbid_s){

	forbidden_state = forbid_s;
}
userOptions::~userOptions() {
	// TODO Auto-generated destructor stub
}


unsigned int userOptions::get_algorithm() {
	return algo;	//returns the selected Algorithm
}

void userOptions::set_algorithm(unsigned int alg) {
	algo =alg;	//assigns the Algorithm selected by the user
}

unsigned int userOptions::getStreamSize() const {
	return stream_size;
}

void userOptions::setStreamSize(unsigned int streamSize) {
	stream_size = streamSize;
}

unsigned int userOptions::getTotalSliceSize() const {
	return total_slice_size;
}

void userOptions::setTotalSliceSize(unsigned int totalSliceSize) {
	total_slice_size = totalSliceSize;
}


