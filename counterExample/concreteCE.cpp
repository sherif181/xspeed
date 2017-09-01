/*
 * concreteCE.cpp
 *
 *  Created on: 20-Jan-2016
 *      Author: rajarshi
 */

#include <counterExample/concreteCE.h>
#include <counterExample/abstractCE.h>
#include <cassert>

concreteCE::concreteCE() {
	// TODO Auto-generated constructor stub

}

concreteCE::~concreteCE() {
	// TODO Auto-generated destructor stub
}

const concreteCE::traj_segment concreteCE::get_segment(unsigned int i) const {
	assert(i < T.size());
	std::list<traj_segment>::const_iterator it = T.begin();
	for (unsigned int k = 0; k != i; k++) {
		it++;
	}
	return *it;

}
const concreteCE::trans_id concreteCE::getTransitionId(unsigned int i) const {
	assert(i < jumps.size());
	std::list<trans_id>::const_iterator it = jumps.begin();
	for (unsigned int k = 0; k != i; k++) {
		it++;
	}
	return *it;
}
void concreteCE::plot_ce(std::string filename, unsigned int x1, unsigned int x2)
{
	traj_segment seg;
	unsigned int locId;
	sample simulation_sample;
	simulation::ptr sim;
	std::cout << "Inside concrete trace plotter, no. of trajectory segs:" << T.size() << std::endl;
	double steps = 100;
	//debug
	//printing the trace start point in a different file for plotting with distinguishible symbol
	std::ofstream myfile;
	// truncate the contents of the bad_trace.o file
	myfile.open(filename.c_str(),ofstream::trunc);
	myfile.close();

	myfile.open("./start_pts");

	for(trajectory::iterator it = T.begin(); it!=T.end();it++){
		seg = *it;
		locId = seg.first;
		std::cout << "plotCE: loc id = " << locId << std::endl;
		simulation_sample = seg.second;
		sim = simulation::ptr(new simulation(simulation_sample.first.size(),steps,ha->getLocation(locId)->getSystem_Dynamics()));
		sim->set_outfile(filename);
		sim->set_out_dimension(x1, x2);
	//	for(unsigned int i=0;i<sim->get_system_dimension();i++)
		myfile << simulation_sample.first[x1] << " " << simulation_sample.first[x2];

	//  std::cout << "simulation: dwell time:" << simulation_sample.second;
		sim->simulate(simulation_sample.first, simulation_sample.second);
	}
	myfile.close();

}

bool concreteCE::valid(struct refinement_point& ref_pt)
{
	traj_segment seg;
	unsigned int locId;
	sample simulation_sample;
	simulation::ptr sim;

	double steps = 100; // defines the precision of validation
	unsigned int seq_no = 0;

	// debug
	bool flag = true;
	//--
	for(trajectory::iterator it = T.begin(); it!=T.end();it++){
		seg = *it;
		locId = seg.first;
		simulation_sample = seg.second;
		//debug

//		std::vector<double> start = simulation_sample.first;
//		std::cout << "The start point of the first segment:\n";
//		if(flag)
//		{
//			for(unsigned int i=0;i<start.size();i++)
//			{
//				std::cout << start[i] << " " << std::endl ;
//			}
//			flag = false;
//		}
		//---
		polytope::ptr Inv;
		Inv = ha->getLocation(locId)->getInvariant();

		sim = simulation::ptr(new simulation(simulation_sample.first.size(),steps,ha->getLocation(locId)->getSystem_Dynamics()));
		bound_sim b; bool status=true;
		b =  sim->bounded_simulation(simulation_sample.first,simulation_sample.second,Inv, status);

		if(!status){
			ref_pt.seq_no = seq_no;
			ref_pt.violating_pt = b.v;
			ref_pt.time = b.cross_over_time;
			return false;
		}
		seq_no++;
	}
	return true;
}
