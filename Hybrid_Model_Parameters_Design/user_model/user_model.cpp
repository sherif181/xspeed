// Created by Hyst v1.3
// Hybrid Automaton in XSpeed
// Converted from file: /home/amit/cuda-workspace/XSpeed/Release/bball_timed.xml
// Command Line arguments: -tool xspeed "" -verbose -output /home/amit/cuda-workspace/XSpeed/Release/nv09test.cpp -input /home/amit/cuda-workspace/XSpeed/Release/bball_timed.xml /home/amit/cuda-workspace/XSpeed/Release/bball.cfg


#include "Hybrid_Model_Parameters_Design/user_model/user_model.h"
void user_model(hybrid_automata& Hybrid_Automata, std::list<initial_state::ptr>& init_state_list, ReachabilityParameters& reach_parameters, userOptions& op) {


typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;

// ************* Section required for Forbidden Inputs *************
unsigned int Directions_Type = 1;
unsigned int iter_max = 6;
double time_horizon = 10.0; 
double sampling_time = 1.0E-4;

op.set_timeStep(sampling_time);
op.set_timeHorizon(time_horizon);
op.set_bfs_level(iter_max);
op.set_directionTemplate(Directions_Type);


unsigned int dim;
size_type row, col;

polytope::ptr initial_polytope_I0, forbid_polytope;
polytope::ptr invariant0;

polytope::ptr guard_polytope0;

Dynamics system_dynamics;

math::matrix<double> ConstraintsMatrixI , ConstraintsMatrixV, invariantConstraintsMatrix , guardConstraintsMatrix , Amatrix , Bmatrix,forbiddenMatrixI;

std::vector<double> boundValueI,boundValueV , C , invariantBoundValue , guardBoundValue, boundValueF;

int boundSignI=1, invariantBoundSign=1, guardBoundSign=1, boundSignV=1;


// The mode name is  always_running

row = 3;
col = 3;
Amatrix.resize(row, col);
Amatrix.clear();
Amatrix(0 , 1) = 1.0;
system_dynamics.isEmptyMatrixA = false;
system_dynamics.MatrixA = Amatrix;

system_dynamics.isEmptyMatrixB = true;

C.resize(row );
C.assign(row,0);
C[1] = -1.0;
C[2] = 1.0;
system_dynamics.isEmptyC = false;
system_dynamics.C = C;


row = 3;
col = 3;
invariantConstraintsMatrix.resize(row, col);
invariantConstraintsMatrix.clear();
invariantConstraintsMatrix(0,0)= -1.0;

invariantBoundValue.resize(row);
invariantBoundValue.assign(row,0);
invariant0 = polytope::ptr(new polytope(invariantConstraintsMatrix, invariantBoundValue,invariantBoundSign));

system_dynamics.U = polytope::ptr(new polytope(true));


std::list<transition::ptr> Out_Going_Trans_fromalways_running;

// The transition label ishop

// Original guard: x <= 0.1 & v < 0

row = 2;
col = 3;

guardConstraintsMatrix.resize(row, col);
guardConstraintsMatrix.clear();
guardConstraintsMatrix(0,0) = 1.0;
guardConstraintsMatrix(1,1) = 1.0;

guardBoundValue.resize(row);
guardBoundValue.assign(row,0);
guardBoundValue[0] = 0.1;
guard_polytope0 = polytope::ptr(new polytope(guardConstraintsMatrix, guardBoundValue, guardBoundSign));


math::matrix<double> R;
row = 3;
col = 3;
R.resize(row, col);
R.clear();
R(0,0) =  1.0;
R(1,1) =  -0.75;
R(2,2) =  1.0;
std::vector<double> w(row);
w.assign(row,0);


Assign assignment;
assignment.Map = R;
assignment.b = w;

transition::ptr t1 = transition::ptr(new transition(1,"hop",1,1,guard_polytope0,assignment));

Out_Going_Trans_fromalways_running.push_back(t1);
location::ptr l1 = location::ptr(new location(1, "always_running", system_dynamics, invariant0, true, Out_Going_Trans_fromalways_running));

Hybrid_Automata.addInitial_Location(l1);
Hybrid_Automata.addLocation(l1);


row = 6;
col = 3;
ConstraintsMatrixI.resize(row, col);
ConstraintsMatrixI.clear();
ConstraintsMatrixI(0 , 0) = 1;
ConstraintsMatrixI(1 , 0) = -1;
ConstraintsMatrixI(2 , 1) = 1;
ConstraintsMatrixI(3 , 1) = -1;
ConstraintsMatrixI(4 , 2) = 1;
ConstraintsMatrixI(5 , 2) = -1;
boundValueI.resize(row );
boundValueI.assign(row,0);
boundValueI[0]=10.2;
boundValueI[1]=-10;


initial_polytope_I0 = polytope::ptr(new polytope(ConstraintsMatrixI, boundValueI, boundSignI));

dim = initial_polytope_I0->getSystemDimension();
int transition_id = 0;
unsigned int initial_location_id =1;

symbolic_states::ptr S0;

initial_state::ptr I0 = initial_state::ptr(new initial_state(initial_location_id, initial_polytope_I0, S0, transition_id));

init_state_list.push_back(I0);
Hybrid_Automata.setDimension(dim);



Hybrid_Automata.insert_to_map("x",0);
Hybrid_Automata.insert_to_map("v",1);
Hybrid_Automata.insert_to_map("t",2);



}
