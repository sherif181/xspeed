/*
 * reachabilityParallel_Iterations.cpp
 *
 *  Created on: 24-Nov-2014
 *      Author: amit
 */

#include "core_system/Reachability/reachabilityParallel_Iterations.h"
/*

const template_polyhedra reachabilityParallelIters(Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial, ReachabilityParameters& ReachParameters,
		polytope::ptr invariant, bool isInvariantExist) {

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;

	std::vector<D>AllDirections = ReachParameters.AdvanceTransposeDirection;

	//omp_set_num_threads(numVectors);
	//#pragma omp parallel
	//int tid;
	int numVectors = ReachParameters.Directions.size1();
	int dimension = Initial->getSystemDimension();
	unsigned int shm_NewTotalIteration = ReachParameters.Iterations;//Shared Variable for resize iterations number on crossing with invariant
	int Min_Or_Max = 2;

	math::matrix<double> MatrixValue;	//Shared Matrix for all child thread
	size_type row = numVectors, col = shm_NewTotalIteration;
	MatrixValue.resize(row, col);

	if (isInvariantExist == true) {	//if invariant exist. Computing
		/ *
		 * Computing support function for polytope for the pair of invairant's direction
		 * to determine the iteration's number at which the polytope is completely outside the invariant's region
		 * and setting the newIteration number as the value obtained here.
		 * /
		shm_NewTotalIteration = InvariantBoundaryCheck(SystemDynamics, Initial,
				ReachParameters, invariant);
		//	cout <<"\nInvariant Exists!!!\n";
	}	//End of Invariant Directions
//	cout <<"\nInvariant DOES NOT Exists!!!\n";
	int flag = 0;
	cout<<"\nNumber of Iterations to be executed = "<<shm_NewTotalIteration<<std::endl;

	int NCores = 8;		//number of cores
	int PartSize = shm_NewTotalIteration / NCores;		//partition size to be handled in each cores
	cout<<"\nPartSize  = "<<PartSize <<std::endl;
	//Becarefull about the fractional partsize

//#pragma omp parallel for		//can try later
	for (int eachDirection = 0; eachDirection < numVectors; eachDirection++) {
		//int Start_zIInitial=0;
//		omp_set_num_threads(NCores);
//#pragma omp parallel for
		for (int iterLoop= 0;iterLoop < NCores; iterLoop++)
			{


				//	std::cout<<"\nThread ID / Each Direction = "<<eachDirection<<endl;
				//int tid = omp_get_thread_num();		//tid is the number of Vectors/Directions
				glpk_lp_solver s_per_thread_I, s_per_thread_U, s_per_thread_inv;

				s_per_thread_I.setMin_Or_Max(2);
				if (!Initial->getIsEmpty())//set glpk constraints If not an empty polytope
					s_per_thread_I.setConstraints(ReachParameters.X0->getCoeffMatrix(),
							ReachParameters.X0->getColumnVector(), ReachParameters.X0->getInEqualitySign());

				s_per_thread_U.setMin_Or_Max(2);
				if (SystemDynamics.U->getIsEmpty()) {	//empty polytope
					//Polytope is empty so no glpk object constraints to be set
				} else {
					s_per_thread_U.setConstraints(SystemDynamics.U->getCoeffMatrix(),
							SystemDynamics.U->getColumnVector(),
							SystemDynamics.U->getInEqualitySign());
				}
				double zIInitial = 0.0, zI = 0.0, zV = 0.0;
				double sVariable, s1Variable;
				//	std::cout<<"\nCheck 1"<<endl;

				std::vector<double> r1Variable;	//now single dimension
				r1Variable.resize(dimension);
				std::vector<double> rVariable;
				rVariable.resize(dimension);
				unsigned int loopIteration = 0;
				sVariable = 0.0; 		//initialize s0

				int li = iterLoop * PartSize;
		//	cout<<"\nOmega_Support(EachDirection) = "<<eachDirection<<endl;

				//rVariable = getTransposedVector(eachDirection,li);
				int position = eachDirection * shm_NewTotalIteration + li;
				//int position = eachDirection * shm_NewTotalIteration + 0;
				if (AllDirections[position].R == eachDirection && AllDirections[position].C == li)
					std::cout<<"\nPosition of the Data Structure of All Directions Correctly mapped!!!\n";
				else
					std::cout<<"\nPosition of the Data Structure of All Directions NOT mapped Correctly ??????XXXXX\n";

				//rVariable = AllDirections[position].v;	//this might work
				for (int i = 0; i < dimension; i++) {
					rVariable[i] = AllDirections[position].v[i];
				}

				if (iterLoop == 0){	//for first partition initial state = Initial
					zIInitial = Omega_Support(ReachParameters, rVariable, Initial,
							SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
					//Start_zIInitial = zIInitial;
					MatrixValue(eachDirection, li) = zIInitial;
				}
				else	// Computing initial state for the remaining  partition
				{

				//	#pragma omp critical
					{
						double original_tau = ReachParameters.time_step;
						double new_delta_tau = ReachParameters.time_step * li;	//shifted dt
						ReachParameters.time_step = new_delta_tau;
						double InitialOmega;
						double temp_zV=0.0, temp_zI=0.0;

						int same_position = eachDirection * shm_NewTotalIteration + li;
						for (int i = 0; i < dimension; i++) {
							r1Variable[i] = AllDirections[same_position].v[i];		//transposed direction
						}
						/ * //ReachParameters.phi_trans.mult_vector(rVariable, r1Variable);
						int start_position = eachDirection * shm_NewTotalIteration + 0;
						for (int i = 0; i < dimension; i++) {
							//should be the Initial direction
							rVariable[i] = AllDirections[start_position].v[i];
						}* /

						temp_zV = W_Support(ReachParameters, SystemDynamics, rVariable,
								s_per_thread_U, Min_Or_Max);
						s1Variable = sVariable + temp_zV;
						temp_zI = Omega_Support(ReachParameters, r1Variable, Initial,
								SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
						InitialOmega = temp_zI + s1Variable; 		//Y1
						MatrixValue(eachDirection, li) = InitialOmega; 		//Y1

						std::cout<<"\n Transposed Initial Omega = "<<InitialOmega<<"\t";
						//now restore back the original time step
						ReachParameters.time_step = original_tau;

						//copying the intermediate value for the next omega computation
						rVariable = CopyVector(r1Variable);		//source to destination
						sVariable = s1Variable;

					}
					std::cout<<"\n Good Very Good!!!!\t";
					std::cout<<"\n sVariable  = "<<sVariable <<"\t";
				}
					std::cout<<"\n SFM = "<<zIInitial<<"\t";
				//std::cout<<"\nOmega Support Matrix Value= "<<MatrixValue(eachDirection, loopIteration)<<endl;

				for (int i=1; i<= PartSize; i++){

					int iter = iterLoop * PartSize + i; //each index of the particular direction
					if (iter == shm_NewTotalIteration){
						break;
					}
					int index = eachDirection * shm_NewTotalIteration + iter;
					//ReachParameters.phi_trans.mult_vector(rVariable, r1Variable);
					for (int x=0;x<dimension;x++)
						r1Variable[x]  = AllDirections[index].v[x];

					if (AllDirections[index].R == eachDirection && AllDirections[index].C == iter)
						std::cout<<"\nPosition Correctly mapped!!!\n";
					else
						std::cout<<"\nPosition NOT Correctly mapped ??????XXXXX\n";


					double TempOmega;
					//		std::cout<<"\tHello = "<<loopIteration;

					/ ** Precompute beta and send it as parameter * /
					zV = W_Support(ReachParameters, SystemDynamics, rVariable,
							s_per_thread_U, Min_Or_Max);
					//cout<<"zV= "<<zV<<"\t";
					s1Variable = sVariable + zV;
					zI = Omega_Support(ReachParameters, r1Variable, Initial,
							SystemDynamics, s_per_thread_I, s_per_thread_U, Min_Or_Max);
					//	cout<<"zi= "<<zI<<"\t";
					TempOmega = zI + s1Variable; 		//Y1
					MatrixValue(eachDirection, iter) = TempOmega; 		//Y1
					//	std::cout<<TempOmega<<"\t";
					rVariable = CopyVector(r1Variable);		//source to destination
					sVariable = s1Variable;
		//amit			loopIteration++;	//for the next Omega-iteration or Time-bound
					loopIteration++;	//for the next Omega-iteration or Time-bound
				}	//end of while for each vector
					//cout<<endl;
			}	//end of partition loop
	}	//end of pragma omp parallel for

	cout<<"Outside Parallel\n";
	/ *
	 * Appending invariant directions and invariant constraints/bounds(alfa)
	 * Goal : To truncate the reachable region within the Invariant region
	 * /
	if (isInvariantExist == true) {	//if invariant exist. Computing
		ReachParameters.TotalDirections = ReachParameters.Directions;
		//	cout << "Rows = " << ReachParameters.TotalDirections.size1() << endl;
		//	cout << "Cols = " << ReachParameters.TotalDirections.size2() << endl;
		int lastDirs = ReachParameters.TotalDirections.size1() - 1;	//total number of rows (0 to rows-1)
		int num_inv = invariant->getColumnVector().size();	//number of Invariant's constriants
		//	cout << "num_inv = "<<num_inv <<endl;
		//	cout << "dimension = "<<dimension <<endl;
		//	cout << "lastDirs = "<<lastDirs<<endl;
		int newDirectionSize = lastDirs + 1 + num_inv;
		ReachParameters.TotalDirections.resize(newDirectionSize, dimension,
				true);
		MatrixValue.resize(newDirectionSize, shm_NewTotalIteration, true);//Matrix resized
		//	cout<<"OK 2\n";
		//cout<<"ReachParameters.TotalDirections.size1() = "<<ReachParameters.TotalDirections.size1()<<endl;
		for (int eachInvDirection = 0; eachInvDirection < num_inv;
				eachInvDirection++) {
			for (int i = 0; i < dimension; i++) {
				/ *ReachParameters.TotalDirections(lastDirs + eachInvDirection + 1,
						i) = ReachParameters.InvariantDirections(
						2 * eachInvDirection, i); * /
				ReachParameters.TotalDirections(lastDirs + eachInvDirection + 1, i) = invariant->getCoeffMatrix()(eachInvDirection, i);
			}
			for (int i = 0; i < shm_NewTotalIteration; i++) {
				MatrixValue(lastDirs + eachInvDirection + 1, i) =
						invariant->getColumnVector()[eachInvDirection];
			}
		}
		//	cout<<"OK 3\n";
		return template_polyhedra(MatrixValue, ReachParameters.TotalDirections);
	} else {
		MatrixValue.resize(numVectors, shm_NewTotalIteration, true);//but writing or resizing only upto the NewTotalIteration
		return template_polyhedra(MatrixValue, ReachParameters.Directions);
	}

}

*/


