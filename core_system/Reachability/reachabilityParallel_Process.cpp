/*
 * reachabilityParallel_Process.cpp
 *
 *  Created on: 24-Nov-2014
 *      Author: amit
 */

#include "core_system/Reachability/reachabilityParallel_Process.h"

void reachFunction(unsigned int eachDirection, Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		lp_solver &s_per_thread_I, lp_solver &s_per_thread_U, double* sf_vals) {
//cout<<"\n\nI am Called \n";
	int dimension = Initial->getSystemDimension();

	//double invariant_SupportFunction = invariant->getColumnVector()[0]; //only one constraint say, lamda

	double zIInitial = 0.0, zI = 0.0, zV = 0.0;
	double sVariable, s1Variable;

	std::vector<double> r1Variable; //now single dimension
	r1Variable.resize(dimension);

	// ****** SharedMemory Requirement Code	Begins :: For Client or Child Process to access or modify SharedMemory ********* /
	int shmid;
	key_t key;
	int *shm_NewTotalIteration;
	// *** We need to get the segment named "5678", created by the server. ******  /
	key = 5678;
	// ******** Locate the segment. ********* /
	if ((shmid = shmget(key, SHMSZ, 0666)) < 0) {
		perror("Error Locating the Segment ::  shmget");
		exit(1);
	}

	// ********Now we attach the segment to our data space. ********** /
	if ((shm_NewTotalIteration = (int *) shmat(shmid, NULL, 0)) == (int *) -1) {
		perror("Error Attaching the Segment to our Data Variable::  shmat");
		exit(1);
	}
	//std::cout<<*shm_NewTotalIteration<<"\n";
	// ****** Now read what the server put in the memory, i.e., initially the Actual TotalIteration  or if modified by child process with vector==1 *** /
	//reading is done in the while loop
	// ****** SharedMemory Requirement Code Ends	 ********* /

	std::vector<double> rVariable;
	rVariable.resize(dimension);
	for (int i = 0; i < dimension; i++) {
		rVariable[i] = ReachParameters.Directions(eachDirection, i);
	}

	sf_vals[0] = eachDirection; //First Element is the Direction/Vector Number
	int loopIteration = 1; //So loop can begin from 1
	sVariable = 0.0; //initialize s0
	int Min_Or_Max = 2;

	zIInitial = Omega_Support(ReachParameters, rVariable, Initial,
			SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);

	sf_vals[loopIteration] = zIInitial; //Y0 = pI(r0)
	loopIteration++;
	//int conditionMet = 0;	//for Optimization of First Vector from computing after crossing invariant's boundary

	//int found = 0;	//has found the crossing point with the Invariant

	for (; loopIteration <= *shm_NewTotalIteration;) { //shm_NewTotalIteration is reading the shared variable
		double TempOmega;
		ReachParameters.phi_trans.mult_vector(rVariable, r1Variable);
		/** Precompute beta and send it as parameter */
		zV = W_Support(ReachParameters, SystemDynamics, rVariable,
				s_per_thread_U, Min_Or_Max);
		s1Variable = sVariable + zV;
		zI = Omega_Support(ReachParameters, r1Variable, Initial, SystemDynamics,
				s_per_thread_I, s_per_thread_U, Min_Or_Max);
		TempOmega = zI + s1Variable; //Y1
		sf_vals[loopIteration] = TempOmega;
		rVariable = CopyVector(r1Variable); //source to destination
		sVariable = s1Variable;
		loopIteration++; //for the next Omega-iteration or Time-bound

	} //end of while for each vector
	//	cout<<"\nChild Function has Ended\n"<<eachDirection<<endl<<loopIteration<<endl<<found;
}

//Implementation of reachability parallel using PROCESS FORKING approach

const template_polyhedra::ptr reachabilityParallel_Process(
		Dynamics& SystemDynamics, supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen) {

	int numVectors = ReachParameters.Directions.size1();
	int iterationNum = ReachParameters.Iterations; //Initial number of iterations this variable should be made sharedVariable

	// ****** SharedMemory Requirement Code	Begins ********* /
	int shmid;
	key_t key;
	int *shm_number;
	// *** We'll name our shared memory segment "5678". ******  /
	key = 5678;
	// ******** Create the segment. ********* /
	if ((shmid = shmget(key, SHMSZ, IPC_CREAT | 0666)) < 0) {
		perror("Error Creating the Segment ::  shmget");
		exit(1);
	}

	// ********Now we attach the segment to our data space. ********** /
	if ((shm_number = (int *) shmat(shmid, NULL, 0)) == (int *) -1) {
		perror("Error Attaching the Segment to our Data Variable::  shmat");
		exit(1);
	}
	// ************ Now put some things into the memory for the other process to read.************ /

	if (isInvariantExist == true) { //if invariant exist. Computing
		/*
		 * Computing support function for polytope for the pair of invairant's direction
		 * to determine the iteration's number at which the polytope is completely outside the invariant's region
		 * and setting the newIteration number as the value obtained here.
		 */
		unsigned int newIter;
		InvariantBoundaryCheck(SystemDynamics, Initial,
				ReachParameters, invariant, lp_solver_type_choosen, newIter);
		*shm_number = newIter;
		//	cout <<"\nInvariant Exists!!!\n";
	} //End of Invariant Directions
	else {
		*shm_number = iterationNum; //initial Actual Total iterations
	}

//	cout << "\nTotal Initial Iteration before = " << *shm_number;

	// ****** SharedMemory Requirement Code Ends	 ********* /

	// IPC with pipes to get computed results from the child processes
	//int fd[2];	//single pipe for all child process
	int fd[numVectors][2]; //as many as number of vectors/Process so as the number of Reading pipes

//	//
//	 if (pipe(fd) < 0) {
//	 std::cout << "reachabilityParallel_Process: Error in creating pipes\n";
//	 exit(0);
//	 }
	//cout << "\nRunning Prallel Processes\n";
	pid_t pID, returnStatus, arr[numVectors];
//	std::cout << "Parent process ID = " << getpid() << std::endl;
	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {

		if (pipe(fd[eachDirection]) < 0) {
			std::cout
					<< "reachabilityParallel_Process: Error in creating pipes\n";
			std::cout << "Error : " << errno << "\n";
			exit(0);
		}

		if ((pID = fork()) == -1) {
			cout << "Error in FORKing child processes" << endl;
			exit(1);
		}
		if (pID == 0) // child
				{
			// define glpk object here. We have one glpk object per process local.
			// polytope.set_glpk_obj(glpk_obj* obj);

			int type = lp_solver_type_choosen;

			lp_solver s_per_thread_I(type), s_per_thread_U(type);
			s_per_thread_I.setMin_Or_Max(2);
			if (!Initial->getIsEmpty()) //set glpk constraints If not an empty polytope
				s_per_thread_I.setConstraints(
						ReachParameters.X0->getCoeffMatrix(),
						ReachParameters.X0->getColumnVector(),
						ReachParameters.X0->getInEqualitySign());

			s_per_thread_U.setMin_Or_Max(2);
			if (SystemDynamics.U->getIsEmpty()) { //empty polytope
				//Polytope is empty so no glpk object constraints to be set
			} else {
				s_per_thread_U.setConstraints(
						SystemDynamics.U->getCoeffMatrix(),
						SystemDynamics.U->getColumnVector(),
						SystemDynamics.U->getInEqualitySign());
			}

			// Code only executed by child process
			close(fd[eachDirection][0]); // reading end of the pipe closed by the child process
//
//			  Array use in pipe communication for passing the values of each row of SupportFunction Matrix(SFM) generated by
//			  the called processes along with the last index as the row number of the SFM.
			double *R = new double[iterationNum + 1];

			// Calling the reachFucntion to compute reachability state for each Direction/Vector
			reachFunction(eachDirection, SystemDynamics, Initial,
					ReachParameters, invariant, s_per_thread_I, s_per_thread_U,
					R);

			//Accessing the SharedMemory to know the size of the last Iterations
			//Avoid inserting direction here in last element as array size varies
			//	R[*shm_number] = eachDirection; // last entry to contain the direction/vector.
			//R[iterationNum] = eachDirection; // last entry to contain the direction/vector.
			if (write(fd[eachDirection][1], R,
					sizeof(double) * (iterationNum + 1)) < 0) { //Accessing the SharedMemory to know the size of the last Iterations
				cout
						<< "reachabilityParallel_Process: Error in writing into the pipe by the child\n";
				printf(
						"Error Writing into the Pipe:Error Code : %d  Description: %s\n",
						errno, strerror(errno));
			}
			exit(1); //important to stop the child from forking again
		} else if (pID < 0) // failed to fork
				{
			std::cout << "reachabilityParallel_Process: Failed to fork" << endl;
			//printf ("Error Writing into the Pipe Description: %s\n", strerror( errno ) );
			exit(1);
		} else { // Code only executed by Parent process
			//	std::cout << "Child process created with Process ID = " << pID << std::endl;
			arr[eachDirection] = pID; //List of process id's
		}
	}
//	close(fd[1]); // write end of the pipe closed at the parent side
	double *BUF = new double[iterationNum + 1];
	unsigned int n;

	math::matrix<double> MatrixValue;
	size_type row = numVectors, col = iterationNum;
	MatrixValue.resize(row, col);

//#pragma omp parallel for
	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {
		close(fd[eachDirection][1]); // write end of the pipe closed at the parent side for eachDirection
		//if (pID == arr[x]) {
		n = read(fd[eachDirection][0], BUF,
				sizeof(double) * (iterationNum + 1)); //Reading all buffer that has been written in the pipe
		if (n < 0) {
			std::cout << "reachabilityParallel_Process: Write to pipe failed"
					<< std::endl;
			exit(0);
		} else {
			//int d = BUF[0];		//First value i.e., eachDirection number
			for (int i = 0; i < iterationNum; i++)
				MatrixValue(eachDirection, i) = BUF[i + 1];
			//MatrixValue(d, i) = BUF[i + 1];
		}
		//}
	} //end of pragma omp parallel for

	////// *
	// while ((n = read(fd[0], BUF, sizeof(double) * (iterationNum + 1))) > 0) {	//Reading all buffer that has been written in the pipe
	//	while ((n = read(fd[0], BUF, sizeof(double) * (*shm_number + 1))) > 0) {//Reading all buffer that has been written in the pipe

	//int d = BUF[0];		//First value i.e., eachDirection number
	//for (int i = 0; i < iterationNum; i++)
	//MatrixValue(d, i) = BUF[i + 1];
	//}
	//if (n < 0) {
	//std::cout << "reachabilityParallel_Process: Write to pipe failed"
	//<< std::endl;
	//exit(0);
	//}

	// * /////

	for (int i = 0; i < numVectors; i++) {
		waitpid(arr[i], &returnStatus, 0); // Parent process waits here for all child to terminate.
	}

	for (int i = 0; i < numVectors; i++)
		close(fd[i][0]);

	//cout<<"Outside Parallel\n";
	/*
	 * Appending invariant directions and invariant constraints/bounds(alfa)
	 * Goal : To truncate the reachable region within the Invariant region
	 */
	if (isInvariantExist == true) { //if invariant exist. Computing
		math::matrix<double> inv_sfm;
		int num_inv = invariant->getColumnVector().size(); //number of Invariant's constriants
		inv_sfm.resize(num_inv, *shm_number);
		for (int eachInvDirection = 0; eachInvDirection < num_inv;
				eachInvDirection++) {
			for (unsigned int i = 0; i < *shm_number; i++) {
				inv_sfm(eachInvDirection, i) =
						invariant->getColumnVector()[eachInvDirection];
			}
		}

		/*	ReachParameters.TotalDirections = ReachParameters.Directions;
		 int lastDirs = ReachParameters.TotalDirections.size1() - 1;	//total number of rows (0 to rows-1)
		 int num_inv = invariant->getColumnVector().size();//number of Invariant's constriants
		 int newDirectionSize = lastDirs + 1 + num_inv;
		 int dimension = ReachParameters.X0->getCoeffMatrix().size2();
		 ReachParameters.TotalDirections.resize(newDirectionSize, dimension,
		 true);
		 MatrixValue.resize(newDirectionSize, *shm_number, true);//Matrix resized
		 for (int eachInvDirection = 0; eachInvDirection < num_inv;
		 eachInvDirection++) {
		 for (int i = 0; i < dimension; i++) {
		 ReachParameters.TotalDirections(lastDirs + eachInvDirection + 1,
		 i) = invariant->getCoeffMatrix()(eachInvDirection, i);
		 }
		 for (int i = 0; i < *shm_number; i++) {
		 MatrixValue(lastDirs + eachInvDirection + 1, i) =
		 invariant->getColumnVector()[eachInvDirection];
		 }
		 }
		 */
		shmdt(shm_number);
		shmctl(shmid, IPC_RMID, NULL);
		// ***************************
		/*
		 * Now free the environment for all glpk objects created inside the pthread when all threads completed
		 */
		//return template_polyhedra(MatrixValue, ReachParameters.TotalDirections,ReachParameters.Directions, invariant->getCoeffMatrix());
		return template_polyhedra::ptr(
				new template_polyhedra(MatrixValue, inv_sfm,
						ReachParameters.Directions,
						invariant->getCoeffMatrix()));
		//return template_polyhedra(MatrixValue, ReachParameters.TotalDirections);
	} else {
		/*
		 * Now free the environment for all glpk objects created inside the pthread when all threads completed
		 */
		MatrixValue.resize(numVectors, *shm_number, true);
		shmdt(shm_number);
		shmctl(shmid, IPC_RMID, NULL);
		// ***************************
		//but writing or resizing only upto the NewTotalIteration
		return template_polyhedra::ptr(
				new template_polyhedra(MatrixValue, ReachParameters.Directions));
	}
}

