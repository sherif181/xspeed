/*
 * simulation.cpp
 *
 *  Created on: 12-Feb-2016
 *      Author: Rajarshi
 */

#include "counterExample/simulation.h"
#include "core_system/continuous/Polytope/Polytope.h"
#include <fstream>

simulation::simulation() {
	// TODO Auto-generated constructor stub

}

simulation::~simulation() {
	// TODO Auto-generated destructor stub
}

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
	Dynamics* D = (Dynamics *)(f_data);
	math::matrix<double> A;

	if(D->isEmptyMatrixA)
	{
		assert(!D->isEmptyC);
		unsigned int dim = D->C.size();
		A = math::matrix<double>(D->C.size(),D->C.size());
		for(unsigned int i=0;i<dim;i++)
			for(unsigned int j=0;j<dim;j++)
				A(i,j)=0;
	}
	else
		A = math::matrix<double>(D->MatrixA);

	std::vector<double> C(A.size1());

	assert(A.size1() == A.size2());
	C = D->C;

	double sum;
	// Input set not included in the dynamics.

	for(unsigned int i=0;i<A.size1();i++){
		sum=0;
		for(unsigned int j=0;j<A.size2();j++){
			sum += A(i,j)*NV_Ith_S(y,j);
		}
		NV_Ith_S(ydot,i) = sum + C[i];
	}
	return 0;
}

std::vector<double> lin_interpolate(std::vector<double> x_l2, std::vector<double> x_l1, double factor)
{
	assert(x_l1.size() == x_l2.size());
	std::vector<double> res(x_l1.size());
	for(unsigned int i=0;i<x_l1.size();i++)
		res[i] = factor*x_l2[i] + x_l1[i]*(1-factor);
	return res;
}
std::vector<double> simulation::simulate(std::vector<double> x, double time)
{
	int flag;
	realtype T0 = 0;
	realtype Tfinal = time;
	realtype t=0;
	Dynamics *data = &D;


	N_Vector y = NULL;
	N_Vector u = NULL;

	assert(x.size() == dimension);
	u = N_VNew_Serial(dimension);

	// return the same initial point if the simulation time is smaller than the simulation time_step
	double time_step = Tfinal/N;

	for(unsigned int i=0;i<dimension;i++)
		NV_Ith_S(u,i) = x[i];

	void *cvode_mem;
	cvode_mem = NULL;
	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

	if( check_flag((void *)cvode_mem, "CVodeCreate", 0))
	{
		throw std::runtime_error("CVODE failed\n");
	}

	/** Input user data */
	CVodeSetUserData(cvode_mem, (void *)data);
	/* Call CVodeInit to initialize the integrator memory and specify the
	* user's right hand side function in u'=f(t,u), the inital time T0, and
	* the initial dependent variable vector u. */

	flag = CVodeInit(cvode_mem, f, T0, u);

	if(check_flag(&flag, "CVodeInit", 1)){
		throw std::runtime_error("CVODE failed\n");
	}

	flag = CVDense(cvode_mem, dimension);
	if (check_flag(&flag, "CVDense", 1))
	{
		throw std::runtime_error("CVODE failed\n");
	}

	 /* Call CVodeSStolerances to specify the scalar relative tolerance
	  * and scalar absolute tolerance */
	flag = CVodeSStolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSStolerances", 1))
	{
		throw std::runtime_error("CVODE failed\n");
	}

	/* In loop over output points: call CVode, print results, test for errors */

	//printing simulation trace in a file for debug purpose, in the plot_dim dimension

	double time_offset = x[get_index("t")];

	bool print_flag = false;
	std::ofstream myfile;
	if(!filename.empty()){
		myfile.open(this->filename.c_str(),std::fstream::app);
		print_flag = true;
	}

	std::vector<double> last(dimension);

	if(print_flag){
		// We plot the initial point also
		myfile << x[this->x1] << "  " << x[this->x2];
		myfile << "\n";

		for(unsigned int k=1;k<=N;k++) {
			double tout = k*time_step;
			// remember this point in the last vector
			for(unsigned int i=0;i<dimension;i++)
				last[i] = NV_Ith_S(u,i);
			flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
			if(check_flag(&flag, "CVode", 1)) break;
			//myfile << NV_Ith_S(u,this->x) << "  " << NV_Ith_S(u,this->y);
			//myfile << time_offset + t << "  " << NV_Ith_S(u,this->x);
			myfile << NV_Ith_S(u,this->x1) << "  " << NV_Ith_S(u,this->x2);
			myfile << "\n";
		}
		myfile << "\n";
		myfile.close();
	}
	else { // no printing the simulation points to file
		for(unsigned int k=1;k<=N;k++) {
			double tout = (k*Tfinal)/N;
			// remember this point in the last vector
			for(unsigned int i=0;i<dimension;i++)
				last[i] = NV_Ith_S(u,i);
			flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
			if(check_flag(&flag, "CVode", 1)) break;
		}
	}

	std::vector<double> res(dimension);
	for(unsigned int i=0;i<dimension;i++)
	{
		res[i] = NV_Ith_S(u,i);
	}
	N_VDestroy_Serial(u); /* Free u vector */
	CVodeFree(&cvode_mem); /* Free integrator memory */

	return res;
}

/**
 * Bounded simulation. Simulation bounded within a polytope.
 * Sets status to false if invariant violated. Otherwise,
 * status is set to true.
 *
 */

bound_sim simulation::bounded_simulation(std::vector<double> x, double time, polytope::ptr I, bool& status)
{
	int flag;
	realtype T0 = 0;
	realtype Tfinal = time;
	realtype t=0;
	Dynamics *data = &D;
	bound_sim b;
	N_Vector u = NULL;

	assert(x.size() == dimension);
	u = N_VNew_Serial(dimension);

	double dist = math::abs(I->point_distance(x));
	if(dist > 1e-5 ){
		throw std::runtime_error("bounded simulation: initial point outside invariant. NLP problem constrains not set correctly\n");
	}
	for(unsigned int i=0;i<dimension;i++)
		NV_Ith_S(u,i) = x[i];

	void *cvode_mem;
	cvode_mem = NULL;
	// Call CVodeCreate to create the solver memory and specify the
	// Backward Differentiation Formula and the use of a Newton iteration

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

	if( check_flag((void *)cvode_mem, "CVodeCreate", 0))
	{
		throw std::runtime_error("CVODE failed\n");
	}

	// Input user data
	CVodeSetUserData(cvode_mem, (void *)data);

	// Call CVodeInit to initialize the integrator memory and specify the
	//ser's right hand side function in u'=f(t,u), the inital time T0, and
	// the initial dependent variable vector u.

	flag = CVodeInit(cvode_mem, f, T0, u);

	if(check_flag(&flag, "CVodeInit", 1)){
		throw std::runtime_error("CVODE failed\n");
	}

	flag = CVDense(cvode_mem, dimension);
	if (check_flag(&flag, "CVDense", 1))
	{
		throw std::runtime_error("CVODE failed\n");
	}

	flag = CVodeSStolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSStolerances", 1))
	{
		throw std::runtime_error("CVODE failed\n");
	}

	//printing simulation trace in a file for debug purpose, in the plot_dim dimension


	bool print_flag = false;

	std::ofstream myfile;
	if(!filename.empty()){
		myfile.open(this->filename.c_str(),std::fstream::app);
		print_flag = true;
	}

	std::vector<double> v(dimension),prev_v(dimension);
	prev_v = x;

	double time_step = Tfinal/N;
//	std::cout << "Total simulation time:" << time << std::endl;
	status = true;
	if(print_flag){
		// We plot the initial point also
		myfile << x[this->x1] << "  " << x[this->x2];
		myfile << "\n";

		for(unsigned int k=1;k<=N;k++) {
			double tout = k*time_step;
			flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
			if(check_flag(&flag, "CVode", 1)) break;
			myfile << NV_Ith_S(u,this->x1) << "  " << NV_Ith_S(u,this->x2);
			myfile << "\n";

		}
		myfile << "\n";
		myfile.close();
	}
	else{ // no printing the simulation points to file
		for(unsigned int k=1;k<=N;k++) {
			double tout = (k*Tfinal)/N;
			flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
			if(check_flag(&flag, "CVode", 1)) break;
			// check polytope satisfaction
			for(unsigned int i=0;i<dimension;i++)
				v[i] = NV_Ith_S(u,i);
			double dist = I->point_distance(v);

			if(math::abs(dist) > 1e-5){
//				std::cout << "time:" << t << " ";
//				std::cout << "distance = " << dist << std::endl;
				status = false;
				break;
			}
			else
				prev_v = v;
		}
	}
	b.v = v;
	b.cross_over_time = t;
	N_VDestroy_Serial(u);
	CVodeFree(&cvode_mem);
	return b;
}

std::vector<double> simulation::metric_simulate(std::vector<double> x, double time,
		double& distance, polytope::ptr I, std::vector<double>& grad)
{
	int flag;
	realtype T0 = 0;
	realtype Tfinal = time;
	realtype t=0;
	Dynamics *data = &D;

	N_Vector u = NULL;

	assert(x.size() == dimension);
	u = N_VNew_Serial(dimension);

	assert(grad.size() == dimension);
	// initialize

	for(unsigned int i=0;i<dimension;i++){
		NV_Ith_S(u,i) = x[i];
		grad[i] = 0;
	}
	void *cvode_mem;
	cvode_mem = NULL;
	// Call CVodeCreate to create the solver memory and specify the
	// Backward Differentiation Formula and the use of a Newton iteration

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

	if( check_flag((void *)cvode_mem, "CVodeCreate", 0))
	{
		throw std::runtime_error("CVODE failed\n");
	}

	// Input user data
	CVodeSetUserData(cvode_mem, (void *)data);

	// Call CVodeInit to initialize the integrator memory and specify the
	//ser's right hand side function in u'=f(t,u), the inital time T0, and
	// the initial dependent variable vector u.

	flag = CVodeInit(cvode_mem, f, T0, u);

	if(check_flag(&flag, "CVodeInit", 1)){
		throw std::runtime_error("CVODE failed\n");
	}

	flag = CVDense(cvode_mem, dimension);
	if (check_flag(&flag, "CVDense", 1))
	{
		throw std::runtime_error("CVODE failed\n");
	}

	flag = CVodeSStolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSStolerances", 1))
	{
		throw std::runtime_error("CVODE failed\n");
	}

	//printing simulation trace in a file for debug purpose, in the plot_dim dimension

	std::vector<double> v(dimension);

	math::matrix<double> expAt;
	std::vector<double> g(dimension);
	distance=0;
	for(unsigned int k=1;k<=N;k++) {
		double tout = (k*Tfinal)/N;
		flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
		if(check_flag(&flag, "CVode", 1)) break;
		// check polytope satisfaction
		for(unsigned int i=0;i<dimension;i++)
			v[i] = NV_Ith_S(u,i);
		distance += I->point_distance(v);
	}
	N_VDestroy_Serial(u);
	CVodeFree(&cvode_mem);
	return v;
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */

  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */

  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }
  }

  /* Check if function returned NULL pointer - no memory allocated */

  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  }

  return(0);
}
