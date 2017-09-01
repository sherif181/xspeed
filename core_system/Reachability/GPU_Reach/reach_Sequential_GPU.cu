/*
 * reach_Sequential_GPU.cpp
 *
 *  Created on: 18-April-2015
 *      Author: amit
 */
#include "core_system/Reachability/GPU_Reach/reach_Sequential_GPU.cuh"
#include "core_system/math/Gimplex/simplex.cuh"
#include "core_system/math/Bulk_LP_Solver/bulk_LP_Solver.h"
#include "boost/timer/timer.hpp"
#include <list>

//Correct implementation for All sizes of Block_LPs
void bulk_Solver(math::matrix<double> constraint_matrix,
		std::vector<double> boundValue, math::matrix<float> list_obj_funs,
		unsigned int number_of_streams, unsigned int no_lps_possible,
		std::vector<float> &res) {

	unsigned int tot_lp = list_obj_funs.size1();
	std::cout << "Total LPs " << tot_lp << std::endl;

	//unsigned int lp_block_size = 28672;	//4004; //183500; //input how many LPs you want to solve at a time ??????
	//unsigned int lp_block_size = tot_lp; //4004; //183500; //input how many LPs you want to solve at a time ??????
	unsigned int lp_block_size = no_lps_possible;

	unsigned int number_of_blocks;
	bool equalBlockSize = true;
	if (tot_lp % lp_block_size == 0) {
		number_of_blocks = tot_lp / lp_block_size;
		equalBlockSize = true;
	} else {
		number_of_blocks = (tot_lp / lp_block_size) + 1; //This Last Block with LESS LPs must be taken care
		equalBlockSize = false;
	}
	std::cout << "Total Blocks " << number_of_blocks << std::endl;

	std::list<block_lp> bulk_lps; //list of sub-division of LPs
	struct block_lp myLPList;
	myLPList.block_obj_coeff.resize(lp_block_size, list_obj_funs.size2());

	if (equalBlockSize == true) { //equal block size :: All blocks are of equal size
		//Equal-Block-Size so making it fit for OMP-Parallel
		for (int i = 0; i < number_of_blocks; i++) {
			bulk_lps.push_back(myLPList);
		} //iterator is ready now
		  //Did not get to test this parallelizing
#pragma omp parallel for
		for (unsigned int lp_number = 0; lp_number < lp_block_size;
				lp_number++) {
			int index = 0;
			for (std::list<block_lp>::iterator it = bulk_lps.begin();
					it != bulk_lps.end(); it++) {
				for (unsigned int i = 0; i < list_obj_funs.size2(); i++) {
					(*it).block_obj_coeff(lp_number, i) = list_obj_funs(
							lp_number + index * lp_block_size, i);
				}
				index++; //now index is just the number of partitions 2 or 3 or 4 may be
			}
		} //end of all LPs
		  //std::cout << "\nEqual Block Size\n";
	} //end of equal-block-size

	if (equalBlockSize == false) { //unequal block size :: Last block has less LPs so solving separately
		//std::cout << "\nUnEqual Block Size!!!!\n";
		int count = 0;
		//Equal-Block-Size so making it fit for OMP-Parallel for 1st part
		for (int i = 0; i < (number_of_blocks - 1); i++) {
			bulk_lps.push_back(myLPList);
		} //iterator is ready now	for all equal-block-size expect the last block

#pragma omp parallel for
		for (unsigned int lp_number = 0; lp_number < lp_block_size;
				lp_number++) {
			int index = 0;
			for (std::list<block_lp>::iterator it = bulk_lps.begin();
					it != bulk_lps.end(); it++) { //will have one less block(i.e., Skips the last block)
				for (unsigned int i = 0; i < list_obj_funs.size2(); i++) {
					(*it).block_obj_coeff(lp_number, i) = list_obj_funs(
							lp_number + index * lp_block_size, i);
				}
				index++;
			}
		} //end of all LPs

		//Taking care about the Last Block with LESS number of LPs which will be SKIPPED in the above if-statement
		//Reading the remaining LAST BLOCK of LPs which was SKIPPED in the above if-statement
		struct block_lp myLPList2;
		int last_block_size = tot_lp - (tot_lp / lp_block_size) * lp_block_size;
		myLPList2.block_obj_coeff.resize(last_block_size,
				list_obj_funs.size2());
		unsigned int index = 0;
		unsigned int lp_number;
		lp_number = (number_of_blocks - 1) * lp_block_size; //starting of Last Block
#pragma omp parallel for
		for (int lp_left = lp_number; lp_left < tot_lp; lp_left++) {
			index = lp_left - lp_number; //index starts from 0 to last LP
			for (unsigned int i = 0; i < list_obj_funs.size2(); i++) {
				myLPList2.block_obj_coeff(index, i) = list_obj_funs(lp_left, i);
			}
		}
		bulk_lps.push_back(myLPList2); //pushing the Last block
	} //end of unequal_block_size

	std::list<block_lp_result> bulk_result;
	struct block_lp_result eachBlock;
	for (std::list<block_lp>::iterator it = bulk_lps.begin();
			it != bulk_lps.end(); it++) {
		unsigned int each_bulk_size = (*it).block_obj_coeff.size1();
		//std::cout << "each_bulk_size = " << each_bulk_size << std::endl;
		eachBlock.results.resize(each_bulk_size);

		Simplex lp_problem(each_bulk_size); //GPU computation
		lp_problem.setConstratint(constraint_matrix, boundValue);
		lp_problem.ComputeLP((*it).block_obj_coeff, number_of_streams); //actual GPU computation

		eachBlock.results = lp_problem.getResultAll();
		bulk_result.push_back(eachBlock);
		// std::cout<<"Result Computed\n";
	}

	res.resize(tot_lp);
	unsigned int index_res = 0;
	for (std::list<block_lp_result>::iterator it = bulk_result.begin();
			it != bulk_result.end(); it++) {
		unsigned int block_result_size = (*it).results.size();
		for (unsigned int i = 0; i < block_result_size; i++) {
			res[index_res] = (*it).results[i];
			index_res++;
		}
	}
//std::cout << "Result size = " << res.size() << std::endl;

}

void bulk_Solver_With_UnitBall(int UnitBall,
		math::matrix<double> constraint_matrix, std::vector<double> boundValue,
		math::matrix<float> list_obj_funs, unsigned int number_of_streams,
		unsigned int no_lps_possible, std::vector<float> &result_X,
		std::vector<float> &result_UnitBall) {

	unsigned int tot_lp = list_obj_funs.size1();
	std::cout << "Total LPs " << tot_lp << std::endl;

	//unsigned int lp_block_size = 28672;	//4004; //183500; //input how many LPs you want to solve at a time ??????
	//unsigned int lp_block_size = tot_lp; //4004; //183500; //input how many LPs you want to solve at a time ??????
	unsigned int lp_block_size = no_lps_possible;

	unsigned int number_of_blocks;
	bool equalBlockSize = true;
	if (tot_lp % lp_block_size == 0) {
		number_of_blocks = tot_lp / lp_block_size;
		equalBlockSize = true;
	} else {
		number_of_blocks = (tot_lp / lp_block_size) + 1; //This Last Block with LESS LPs must be taken care
		equalBlockSize = false;
	}
	std::cout << "Total Blocks " << number_of_blocks << std::endl;

	std::list<block_lp> bulk_lps; //list of sub-division of LPs
	struct block_lp myLPList;
	myLPList.block_obj_coeff.resize(lp_block_size, list_obj_funs.size2());

	if (equalBlockSize == true) { //equal block size :: All blocks are of equal size
		//Equal-Block-Size so making it fit for OMP-Parallel
		for (int i = 0; i < number_of_blocks; i++) {
			bulk_lps.push_back(myLPList);
		} //iterator is ready now
		  //Did not get to test this parallelizing
#pragma omp parallel for
		for (unsigned int lp_number = 0; lp_number < lp_block_size;
				lp_number++) {
			int index = 0;
			for (std::list<block_lp>::iterator it = bulk_lps.begin();
					it != bulk_lps.end(); it++) {
				for (unsigned int i = 0; i < list_obj_funs.size2(); i++) {
					(*it).block_obj_coeff(lp_number, i) = list_obj_funs(
							lp_number + index * lp_block_size, i);
				}
				index++; //now index is just the number of partitions 2 or 3 or 4 may be
			}
		} //end of all LPs
		  //std::cout << "\nEqual Block Size\n";
	} //end of equal-block-size

	if (equalBlockSize == false) { //unequal block size :: Last block has less LPs so solving separately
		//std::cout << "\nUnEqual Block Size!!!!\n";
		int count = 0;
		//Equal-Block-Size so making it fit for OMP-Parallel for 1st part
		for (int i = 0; i < (number_of_blocks - 1); i++) {
			bulk_lps.push_back(myLPList);
		} //iterator is ready now	for all equal-block-size expect the last block

#pragma omp parallel for
		for (unsigned int lp_number = 0; lp_number < lp_block_size;
				lp_number++) {
			int index = 0;
			for (std::list<block_lp>::iterator it = bulk_lps.begin();
					it != bulk_lps.end(); it++) { //will have one less block(i.e., Skips the last block)
				for (unsigned int i = 0; i < list_obj_funs.size2(); i++) {
					(*it).block_obj_coeff(lp_number, i) = list_obj_funs(
							lp_number + index * lp_block_size, i);
				}
				index++;
			}
		} //end of all LPs

		//Taking care about the Last Block with LESS number of LPs which will be SKIPPED in the above if-statement
		//Reading the remaining LAST BLOCK of LPs which was SKIPPED in the above if-statement
		struct block_lp myLPList2;
		int last_block_size = tot_lp - (tot_lp / lp_block_size) * lp_block_size;
		myLPList2.block_obj_coeff.resize(last_block_size,
				list_obj_funs.size2());
		unsigned int index = 0;
		unsigned int lp_number;
		lp_number = (number_of_blocks - 1) * lp_block_size; //starting of Last Block
#pragma omp parallel for
		for (int lp_left = lp_number; lp_left < tot_lp; lp_left++) {
			index = lp_left - lp_number; //index starts from 0 to last LP
			for (unsigned int i = 0; i < list_obj_funs.size2(); i++) {
				myLPList2.block_obj_coeff(index, i) = list_obj_funs(lp_left, i);
			}
		}
		bulk_lps.push_back(myLPList2); //pushing the Last block
	} //end of unequal_block_size

	std::list<block_lp_result> bulk_result;
	struct block_lp_result eachBlock;
	for (std::list<block_lp>::iterator it = bulk_lps.begin();
			it != bulk_lps.end(); it++) {
		unsigned int each_bulk_size = (*it).block_obj_coeff.size1();
		//std::cout << "each_bulk_size = " << each_bulk_size << std::endl;
		eachBlock.results.resize(each_bulk_size);

		Simplex Solver(UnitBall, each_bulk_size); //GPU computation
		Solver.setConstratint(constraint_matrix, boundValue, UnitBall);
		Solver.ComputeLP((*it).block_obj_coeff, UnitBall, number_of_streams); //actual GPU computation

		Solver.getResult_X(eachBlock.results);
		Solver.getResult_UnitBall(eachBlock.results_UnitBall);
		bulk_result.push_back(eachBlock);
		// std::cout<<"Result Computed\n";
	}

	//std::vector<float> res(tot_lp);
	result_X.resize(tot_lp);
	result_UnitBall.resize(tot_lp);

	unsigned int index_res = 0;
	for (std::list<block_lp_result>::iterator it = bulk_result.begin();
			it != bulk_result.end(); it++) {
		unsigned int block_result_size = (*it).results.size();
		for (unsigned int i = 0; i < block_result_size; i++) {
			//res[index_res] = (*it).results[i];
			result_X[index_res] = (*it).results[i];
			result_UnitBall[index_res] = (*it).results_UnitBall[i];
			index_res++;
		}
	}
}

/*
 * After optimising the duplicate Support Function computation
 */

void reachabilitySequential_GPU(unsigned int boundedTotIteration, Dynamics& SystemDynamics,
		supportFunctionProvider::ptr Initial,
		ReachabilityParameters& ReachParameters, polytope::ptr invariant,
		bool isInvariantExist, int lp_solver_type_choosen,
		unsigned int number_of_streams, int Solver_GLPK_Gurobi_GPU,
		template_polyhedra::ptr & reachableRegion) {
	//template_polyhedra::ptr reachRegion;

	typedef typename boost::numeric::ublas::matrix<double>::size_type size_type;
	unsigned int NewTotalIteration = ReachParameters.Iterations;
	bool U_empty = false;
	int num_inv = invariant->getColumnVector().size(); //number of Invariant's constriants
	std::vector<double> inv_bounds;
	inv_bounds = invariant->getColumnVector();
	math::matrix<double> inv_directions;
	inv_directions = invariant->getCoeffMatrix();
//	cout<<"Inv_directions are :: "<<inv_directions<<std::endl;
//	cout<<"inv_bounds are :: ";
//	for (int i=0;i<num_inv;i++){
//		cout<<inv_bounds[i]<<"\t";
//	}
//	std::cout << "num_inv = "<<num_inv<<"\n";
	if (isInvariantExist == true) { //if invariant exist. Computing
		std::cout << "Yes Invariant Exist!!!";
		NewTotalIteration = boundedTotIteration;
		std::cout << "NewTotalIteration = " << NewTotalIteration << std::endl;
	} //End of Invariant Directions
	if (NewTotalIteration <= 1) {
		template_polyhedra::ptr poly_empty;
		//return poly_empty; //NO need to proceed Algorithm further
		reachableRegion = poly_empty;
	}
	//reachableRegion = template_polyhedra::ptr(new template_polyhedra());
	//cout<<"reachableRegion->getTotalIterations() 1 = " <<reachableRegion->getTotalIterations()<<std::endl;
	if (SystemDynamics.U->getIsEmpty()) { //polytope U can be empty set
		U_empty = true;
		//std::cout<<"U is Empty!!!!";
	}

	int Solver = Solver_GLPK_Gurobi_GPU; //1 for CPU solver(GLPK); //2 for CPU solver(Gurobi); //3 for GPU solver(Gimplex)
//  ************* Generation of Directions Begins ***************
//std::vector<AllDirection> Direction_List;
	unsigned int numVectors = ReachParameters.Directions.size1();

	unsigned int totalDirList1 = numVectors * (NewTotalIteration + 1); //1 extra for loop1
	math::matrix<float> List_for_X0(totalDirList1,
			ReachParameters.Directions.size2());
	unsigned int totalDirList2 = numVectors * NewTotalIteration; //'n' dirs for each 'n' loops
	math::matrix<float> List_for_U(totalDirList2,
			ReachParameters.Directions.size2());

	/*std::cout << "\nNumber of Directions/LPs for X0 = " << totalDirList1;
	 if (!U_empty) {
	 std::cout << "\nNumber of Directions/LPs for U = " << totalDirList2;
	 }*/

	std::list<std::vector<double> > List_X0;
	std::list<std::vector<double> > List_U;
	//cout<<"reachableRegion->getTotalIterations() 2 = " <<reachableRegion->getTotalIterations()<<std::endl;
	boost::timer::cpu_timer DirectionsGenerate_time;
	DirectionsGenerate_time.start();
	if (Solver == 3) {
		//for OMP --cuda not supporting OMP-- so added library "lgomp"; build-stage -Xcompiler -fopenmp
		int numCoresAvail = omp_get_num_procs(); //get the number of cores
		getDirectionList_X0_and_U(numCoresAvail, ReachParameters, NewTotalIteration,
				List_for_X0, List_for_U, U_empty, SystemDynamics); //Optimized into a single function the 2 Tasks
	} else {
		//Only for profiling GLPK solver Time for comparison with boundary value implementation
		getDirectionList_X0_and_U_OnlyForGLPK(ReachParameters,
				NewTotalIteration, List_X0, List_U, U_empty, SystemDynamics); //Optimized into a single function the 2 Tasks
	}
	DirectionsGenerate_time.stop();
	double wall_clock1;
	wall_clock1 = DirectionsGenerate_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
	double return_Time1 = wall_clock1 / (double) 1000;
	std::cout
			<< "\nDirections Generation(parallel): Boost Time taken:Wall  (in Seconds) = "
			<< return_Time1 << std::endl;
//  ************* Generation of Directions Ends ***************

	int dimension = Initial->getSystemDimension();
	int Min_Or_Max = 2;
	size_type row = numVectors, col = NewTotalIteration;
	math::matrix<double> MatrixValue(row, col);

	std::vector<float> supp_func_X0, supp_func_U, supp_func_UnitBall,
			result_dotProduct;
	//cout<<"reachableRegion->getTotalIterations() 3  = " <<reachableRegion->getTotalIterations()<<std::endl;
	int device;
	cudaDeviceProp props;
	cudaGetDevice(&device);
	cudaGetDeviceProperties(&props, device);
	double MemorySize = props.totalGlobalMem / 1024; //converting into KiloBytes
//		std::cout << "\nMemorySize = " << MemorySize;
//		std::cout << "Each LP Size (KiloBytes) = " << eachLP_Size << std::endl;
	unsigned int no_lps_possible;
	//cout<<"reachableRegion->getTotalIterations() 3b  = " <<reachableRegion->getTotalIterations()<<std::endl;
	if (Solver == 3) { // ************ GRAPHIC PROCESSING UNIT Computations  ********
		//result obtained from GPU interface // OMP can also be tried and compared

		bool IsBoundedLPSolver = true; //Remember to make this false when testing General Gimplex
		double eachLP_Size;
		if (IsBoundedLPSolver) {
			/*
			 * ToDo::recompute the LP size based on NEW Implementataion of Single Kernel
			 * For bounded LP Solver:: each lp size is just objective function size * number of LPs to be solved
			 * and a constant factor:: size of the bound vectors
			 */
			eachLP_Size = ReachParameters.Directions.size2() * sizeof(float);
			eachLP_Size = eachLP_Size / (double) 1024; //converted into KiloBytes
			unsigned int boundValueSize =
					ReachParameters.X0->getColumnVector().size()
							* sizeof(float);

			no_lps_possible = (MemorySize - boundValueSize) / eachLP_Size; //Taking less by integer_casting NO PROBLEM
			std::cout << "Number of LP per bulk Possible is " << no_lps_possible
					<< std::endl;
			//	cout<<"reachableRegion->getTotalIterations() 3c  = " <<reachableRegion->getTotalIterations()<<std::endl;
		} else {
			/*
			 * lp_block_size can be computed from the Constraint_matrix Model
			 MaxSize per LP = row * (col + row + row_artificial + 4) * sizeof(float)
			 */
			eachLP_Size = (ReachParameters.X0->getCoeffMatrix().size1() + 1)
					* (2 * ReachParameters.X0->getCoeffMatrix().size2()
							+ ReachParameters.X0->getCoeffMatrix().size1() + 4)
					* sizeof(float);
			eachLP_Size = eachLP_Size / (double) 1024; //converted into KiloBytes
			no_lps_possible = MemorySize / eachLP_Size; //Taking less by integer_casting NO PROBLEM
			std::cout << "Number of LP per bulk Possible = " << no_lps_possible
					<< std::endl;
		}
		std::cout << "totalDirList1 LP = " << totalDirList1 << std::endl;
		bool single_bulk = true; //false -- if multiple bulk is required to be processed/solved
		if (totalDirList1 > no_lps_possible) { //totalDirList1 is the Maximum
			single_bulk = false;
		}
		boost::timer::cpu_timer onlyGimplex_time;
		onlyGimplex_time.start();
		if (single_bulk) { //If single bulk is solved this code-section takes less time than else-part
			if (IsBoundedLPSolver) {
				//:: Solve in a single call for both X0 and supFunUnitBall to same multiple memory transfer and separately for U
				if (!U_empty) { //polytope U can be empty set
					std::cout << "totalDirList2 LP = " << totalDirList2
							<< std::endl;
					Simplex simplex_for_U(totalDirList2);
					simplex_for_U.setConstratint(
							SystemDynamics.U->getCoeffMatrix(),
							SystemDynamics.U->getColumnVector());
					simplex_for_U.ComputeLP(List_for_U, number_of_streams);
					supp_func_U = simplex_for_U.getResultAll();
				} //working
				  //compute only X0 and supp_unitBall in one kernel
				int UnitBall = 1; //just to have different signature for the overloaded functions of class Simplex
				Simplex solver(UnitBall, totalDirList1);
				solver.setConstratint(ReachParameters.X0->getCoeffMatrix(),
						ReachParameters.X0->getColumnVector(), UnitBall);
				std::cout << "New imple with dotProduction:: Started\n";
				//	cout<<"reachableRegion->getTotalIterations() 3d  = " <<reachableRegion->getTotalIterations()<<std::endl;
				//solver.ComputeLP(List_for_X0, UnitBall, number_of_streams); OLD method
				solver.ComputeLP(List_for_X0, UnitBall, number_of_streams,
						SystemDynamics.C);
				//todo:: some memory issue exist here

				//	cout<<"reachableRegion->getTotalIterations() 3e  = " <<reachableRegion->getTotalIterations()<<std::endl;

				solver.getResult_X(supp_func_X0); //return the result as argument
				solver.getResult_UnitBall(supp_func_UnitBall);
				//std::cout<<"supp_func_UnitBall.size = "<<supp_func_UnitBall.size();
				solver.getResult_dotProduct(result_dotProduct);
				//std::cout<<"result_dotProduct.size = "<<result_dotProduct.size();
//				std::cout<<"result_dotProduct Values :: ";
//				for (int i=0;i<result_dotProduct.size();i++){
//						cout<<result_dotProduct[i]<<"\t";
//					}
				std::cout << "New imple with dotProduction:: Ended\n";
			} else { //OLD implementation with single call of kernel for X0 and U
					 //TODO::But now missing supp_fun_UnitBall_infinity_norm
				Simplex simplex_for_X0(totalDirList1), simplex_for_U(
						totalDirList2);
				simplex_for_X0.setConstratint(
						ReachParameters.X0->getCoeffMatrix(),
						ReachParameters.X0->getColumnVector());
				if (!U_empty) { //polytope U can be empty set
					simplex_for_U.setConstratint(
							SystemDynamics.U->getCoeffMatrix(),
							SystemDynamics.U->getColumnVector());
				} //else { U_empty = true;}
				simplex_for_X0.ComputeLP(List_for_X0, number_of_streams);
				supp_func_X0 = simplex_for_X0.getResultAll();
				if (!U_empty) {
					simplex_for_U.ComputeLP(List_for_U, number_of_streams);
					supp_func_U = simplex_for_U.getResultAll();
				}
			}
			//	cout<<"reachableRegion->getTotalIterations() 4 = " <<reachableRegion->getTotalIterations()<<std::endl;
			std::cout << "Single Bulk";
		} else { //IF SOLVING BY DIVISION IS REQUIRED
			if (IsBoundedLPSolver) {
				//:: Solve in a single call for both X0 and supFunUnitBall to same multiple memory transfer and separately for U

				int UnitBall = 1; //just to have different signature for the overloaded functions of class Simplex
				bulk_Solver_With_UnitBall(UnitBall,
						ReachParameters.X0->getCoeffMatrix(),
						ReachParameters.X0->getColumnVector(), List_for_X0,
						number_of_streams, no_lps_possible, supp_func_X0,
						supp_func_UnitBall); //ONLY UnitBall result will extra

				if (!U_empty) { //polytope U can be empty set	//NO CHANGE REQUIRED HERE IN BULK SOLVER
					bulk_Solver(SystemDynamics.U->getCoeffMatrix(),
							SystemDynamics.U->getColumnVector(), List_for_U,
							number_of_streams, no_lps_possible, supp_func_U);
				}
			} else {
				//OLD implementation with single call of kernel for X0 and U
				//TODO::But now missing supp_fun_UnitBall_infinity_norm
			}
		}
		onlyGimplex_time.stop();
		double wall_clock, user_clock, system_clock;
		wall_clock = onlyGimplex_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
		//user_clock = onlyGimplex_time.elapsed().user / 1000000;
		//system_clock = onlyGimplex_time.elapsed().system / 1000000;
		double return_Time = wall_clock / (double) 1000;
		std::cout
				<< "\nGPU-simplex Solver: Boost Time taken:Wall  (in Seconds) = "
				<< return_Time << std::endl;
		//std::cout << "\nGPU-simplex Boost Time taken:User  (in Seconds) = " << user_clock / (double) 1000 << std::endl;
		//std::cout << "\nGPU-simplex Boost Time taken:System  (in Seconds) = " << system_clock / (double) 1000 << std::endl;

		// *********************** GPU computation Over  *********************
	}
//	cout<<"reachableRegion->getTotalIterations() 5 = " <<reachableRegion->getTotalIterations()<<std::endl;
	if (Solver >= 1 && Solver < 3) { // ************ CPU Solver ****************

	//Todo::Similarly i have to implement Solver for UnitBall_infinity_norm if i have to use Final loop for reachAlgorithm
		boost::timer::cpu_timer onlyGimplex_time;
		onlyGimplex_time.start();
		bulk_lp_solver simplex_for_X0(Solver), simplex_for_U(Solver); //Solver = 1 for GLPK; = 2 for Gurobi
		bool U_empty = false;
		simplex_for_X0.setMaxMin(2); //2 for Maximum
		simplex_for_X0.setConstratint(ReachParameters.X0->getCoeffMatrix(),
				ReachParameters.X0->getColumnVector());
		if (!SystemDynamics.U->getIsEmpty()) { //polytope U can be empty set
			simplex_for_U.setMaxMin(2); //2 for Maximum
			simplex_for_U.setConstratint(SystemDynamics.U->getCoeffMatrix(),
					SystemDynamics.U->getColumnVector());
		} else {
			U_empty = true;
		}
		//simplex_for_X0.ComputeLP(List_for_X0);
		simplex_for_X0.ComputeLP_ListVector(List_X0); //only for GLPK comparison

		supp_func_X0 = simplex_for_X0.getResultAll();

		if (!U_empty) {
			//simplex_for_U.ComputeLP(List_for_U);
			simplex_for_U.ComputeLP_ListVector(List_U); //only for GLPK comparison
			supp_func_U = simplex_for_U.getResultAll();
		}
		onlyGimplex_time.stop();
		double wall_clock, user_clock, system_clock;
		wall_clock = onlyGimplex_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
		double return_Time = wall_clock / (double) 1000;
		std::cout
				<< "\nCPU(GLPK/Gurobi) Solver: Boost Time taken:Wall  (in Seconds) = "
				<< return_Time << std::endl;

	} // ************ CPU Solver Over ****************

//unsigned int index = 0, index_X0 = 0, index_U = 0;//indices upto totalDirList1 and totalDirList2
	//cout<<"reachableRegion->getTotalIterations() 6 = " <<reachableRegion->getTotalIterations()<<std::endl;

	std::cout << "\n Before Final Reach Algorithm ";
	std::cout << std::endl;
	//Breaking here	for TESTING/Reading LP_Solver
//	return template_polyhedra(MatrixValue, ReachParameters.Directions);

	boost::timer::cpu_timer reachLoop_time;
	reachLoop_time.start();

#pragma omp parallel for
	for (unsigned int eachDirection = 0; eachDirection < numVectors;
			eachDirection++) {
		unsigned int index_X0, index_U; //making the index suitable for parallelizing
		//unsigned int index; 	//index = eachDirection * NewTotalIteration;
		//here i have a list of result of Supp_fun_Of_UnitBall_infinity_norm
		if (eachDirection == 0) { //only starting loop begins with 0
			index_X0 = eachDirection * NewTotalIteration;
		} else { //
			index_X0 = eachDirection * NewTotalIteration + eachDirection; //only X0(list_X0) has 2 directions for first-iteration
		}
		if (!U_empty) {
			index_U = eachDirection * NewTotalIteration;
		}
		double res1;
		double term1, term2, term3, term3a, term3b, res2, term3c = 0.0;
		double zIInitial = 0.0, zI = 0.0, zV = 0.0;
		double sVariable = 0.0, s1Variable; //initialize s0
		std::vector<double> rVariable(dimension), r1Variable(dimension);
		unsigned int loopIteration = 0;
		//	std::cout<<"Testing 1\n";
		//  ************** Omega Function   ********************
		res1 = supp_func_X0[index_X0]; //X0->SF(direction)			//	0
		//	std::cout<<"Testing 2\n";
		//term3b = support_unitball_infnorm(Direction_List[index].direction);
		term3b = (double) supp_func_UnitBall[index_X0]; //  needed  0
		if (!SystemDynamics.isEmptyC) {
			term3c = ReachParameters.time_step * result_dotProduct[index_X0];
		}
		//	std::cout<<"Testing 3\n";
		index_X0++; //	made 1
		term1 = supp_func_X0[index_X0]; //X0->SF(phi_trans_dir)		//  1
		//	std::cout<<"Testing 4\n";
		index_X0++; //	made 2
		if (!U_empty) {
			term2 = ReachParameters.time_step * supp_func_U[index_U]; //U->SF(Btrans_dir)
			index_U++;
		} else
			term2 = 0;

		term3a = ReachParameters.result_alfa; //compute_alfa(ReachParameters.time_step,system_dynamics,Initial_X0);
		term3 = term3a * term3b;
		res2 = term1 + term2 + term3 + term3c; //term3c Added
		//zIInitial = (res1 > res2 ? res1:res2);
		if (res1 > res2)
			zIInitial = res1;
		else
			zIInitial = res2;
		//  ************** Omega Function   ********************
		MatrixValue(eachDirection, loopIteration) = zIInitial; //		index++;
		loopIteration++;
		for (; loopIteration < NewTotalIteration;) { //Now stopping condition is only "shm_NewTotalIteration"
			//  ************** W_Support Function   ********************
			double result, res_sup;
			if (!U_empty) {
				res1 = ReachParameters.time_step * supp_func_U[index_U - 1]; //replace previous value
				//index_U++;
			} else {
				res1 = 0;
			}
			double beta = ReachParameters.result_beta;
			//res_sup = (double) support_unitball_infnorm(Direction_List[index].direction);
			//res_sup = (double) supp_func_UnitBall[d_index];  d_index++;	//Should replace from previous computation
			//res_sup = term3b; //replaced from previous steps
			/*if (loopIteration == 1) // needed  0 again here
				res_sup = supp_func_UnitBall[index_X0 - 2]; //Should replace from previous computation*/

			//double res_beta = beta * res_sup;
			double res_beta = beta * term3b;

			result = res1 + res_beta + term3c; //Added term3c
			zV = result;
			//  ************** W_Support Function   ********************
			s1Variable = sVariable + zV;
			//  ************** Omega Function   ********************
			//double res1;
			res1 = supp_func_X0[index_X0 - 1]; ////replace previous value....  X0->SF(direction)		//	(2 -1)=1
			double term1, term2, term3, term3a, res2;
			term1 = supp_func_X0[index_X0]; //X0->SF(phi_trans_dir)		//	2
			index_X0++; // 	made 3
			if (!U_empty) {
				term2 = ReachParameters.time_step * supp_func_U[index_U]; //U->SF(Btrans_dir)
				index_U++;
			} else {
				term2 = 0;
			}
			term3a = ReachParameters.result_alfa; //compute_alfa(ReachParameters.time_step,system_dynamics,Initial_X0);
			//term3b = support_unitball_infnorm(Direction_List[index - 1].Phi_trans_dir);
			//term3b = support_unitball_infnorm(Direction_List[index].direction1);
			if (loopIteration == 1) {
				term3b = (double) supp_func_UnitBall[index_X0 - 2]; //Compute here		//needed 1
				if (!SystemDynamics.isEmptyC) {
					term3c = ReachParameters.time_step
							* result_dotProduct[index_X0 - 2];
				}
			} else {
				term3b = (double) supp_func_UnitBall[index_X0 - 1];
				if (!SystemDynamics.isEmptyC) {
					term3c = ReachParameters.time_step
							* result_dotProduct[index_X0 - 1];
				}
			}

			term3 = term3a * term3b;
			res2 = term1 + term2 + term3 + term3c;
			//zIInitial = (res1 > res2 ? res1:res2);
			if (res1 > res2)
				zI = res1;
			else
				zI = res2;
			//  ************** Omega Function   ********************
			double TempOmega;
			TempOmega = zI + s1Variable; //Y1
			MatrixValue(eachDirection, loopIteration) = TempOmega; //Y1
			sVariable = s1Variable; //	index++;
			loopIteration++; //for the next Omega-iteration or Time-bound
		} //end of all Iterations of each vector/direction
	} //end of for each vector/directions
	  //cout<<"reachableRegion->getTotalIterations() = 6" <<reachableRegion->getTotalIterations()<<std::endl;
	reachLoop_time.stop();
	double wall_clock;
	wall_clock = reachLoop_time.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
	double reach_Time = wall_clock / (double) 1000;
	std::cout << "\nFinal Reach Loop Time:Wall  (in Seconds) = " << reach_Time
			<< std::endl;

	/** Appending invariant directions and invariant constraints/bounds(alfa)
	 ** Goal : To truncate the reachable region within the Invariant region	 */
	//int num_inv = invariant->getColumnVector().size(); //number of Invariant's constriants
	//std::cout<<"working 2ab"<<std::endl;
	if (isInvariantExist == true) { //if invariant exist. Computing
		//std::cout<<"Test inside 1a\n";
		math::matrix<double> inv_sfm;
//		std::cout<<"num_inv = "<<num_inv <<"\n";
//		std::cout<<"working"<<std::endl;
		inv_sfm.resize(num_inv, NewTotalIteration);
		for (int eachInvDirection = 0; eachInvDirection < num_inv;
				eachInvDirection++) {
			//std::cout<<"working"<<std::endl;
			for (unsigned int i = 0; i < NewTotalIteration; i++) {
				//inv_sfm(eachInvDirection, i) = invariant->getColumnVector()[eachInvDirection];
				inv_sfm(eachInvDirection, i) = inv_bounds[eachInvDirection];
			}
		}
		//	std::cout<<"MatrixValue is ::"<<MatrixValue<<std::endl;
//		std::cout<<"inv_sfm is ::"<<inv_sfm<<std::endl;
//		std::cout<<"inv_directions is ::"<<inv_directions<<std::endl;
//		std::cout<<"ReachParameters.Directions is ::"<<ReachParameters.Directions<<std::endl;
		//return template_polyhedra::ptr(new template_polyhedra(MatrixValue, inv_sfm, ReachParameters.Directions, invariant->getCoeffMatrix()));
		//	std::cout<<"working a2"<<std::endl;

		reachableRegion = template_polyhedra::ptr(new template_polyhedra());
		//	std::cout<<"reachRegion size = "<<reachableRegion->getTotalIterations()<<std::endl;
		//	std::cout<<"working 2b"<<std::endl;
		reachableRegion->setTemplateDirections(ReachParameters.Directions);
		//std::cout<<"working 2b"<<std::endl;
		reachableRegion->setMatrix_InvariantBound(inv_sfm);
		//std::cout<<"working 3"<<std::endl;
		reachableRegion->setInvariantDirections(inv_directions);
		//std::cout<<"working 4"<<std::endl;
		reachableRegion->setMatrixSupportFunction(MatrixValue);
		//return template_polyhedra::ptr( new template_polyhedra(MatrixValue, inv_sfm, ReachParameters.Directions, inv_directions));
	} else {
		reachableRegion = template_polyhedra::ptr(new template_polyhedra());
		reachableRegion->setMatrixSupportFunction(MatrixValue);
		reachableRegion->setTemplateDirections(ReachParameters.Directions);
		//	return template_polyhedra::ptr(	new template_polyhedra(MatrixValue, ReachParameters.Directions));
	}
	//std::cout<<"working 5"<<std::endl;
	//return reachRegion;
}

