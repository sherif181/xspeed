/*
 * All_PP_Definition.h
 *
 *  Created on: 11-Dec-2014
 *      Author: amit
 */

#ifndef ALL_PP_DEFINITION_H_
#define ALL_PP_DEFINITION_H_

// *********** Various Algorithms ******************
#define SEQ_SF			1	//Sequential Reachability Algorithm
#define PAR_SF		 	2	//Parallel Reachability Algorithm using OpenMP Thread Creation, using parallelizing over Directions
#define TIME_SLICE	 	3	//Parallel Reachability Algorithm using OpenMP Thread Creation, using parallelizing over Iterations
#define AGJH_BFS		4	//Parallel Breadth First Search Algorithm using Gerard J. Holzmann Algorithm
#define TPBFS			5	//Parallel Breadth First Search Algorithm using Load Balancing Algorithm
#define GPU_SF			6	//PostC using GPU with Sequential BFS



//#define PAR_PROCESS	3	//Parallel Reachability Algorithm using Process Creation, using parallelizing over Directions
//#define PAR_ITER_DIR	5	//Parallel Reachability Algorithm using parallelizing over Iterations and Directions (using OpenMP Thread Creation)
//#define PAR_BY_PARTS	6	//Parallel Reachability Algorithm parallelizing over Initial Set Partitioning (using OpenMP Thread Creation)
//#define PAR_BY_PARTS_ITERS	7	//Parallel Reachability Algorithm parallelizing over multiple Initial Set ALSO Partitioning iterations.
//#define SAME_DIRS		8	//If dir1 = phi_transpose x dir then no need to compute support function for dir1
//#define ALL_DIRS		9	//compute all the list of transpose directions
//#define GPU_MULTI_SEQ	10	//Matrix-Vector GPU multiplication
//#define BFS			12	//Sequential Breadth First Search Algorithm
// *********** Various Algorithms ******************

// *********** Various Hybrid System Model ******************
#define BBALL				1		//Bouncing Ball
#define TBBALL				2		//Timed Bouncing Ball
#define HELICOPTER			3		//HELICOPTER Model
#define FIVEDIMSYS			4		//Model of A Five Dimensional System
#define NAVIGATION_1		5		//NAVIGATION Benchmark model - NAV01 (3 x 3)
#define NAVIGATION_2		6		//NAVIGATION Benchmark model - NAV02 (3 x 3)
#define NAVIGATION_3		7		//NAVIGATION Benchmark model - NAV03 (3 x 3)
#define NAVIGATION_4		8		//NAVIGATION Benchmark model - NAV04 (5 x 5)
#define NAVIGATION_5		9		//NAVIGATION Benchmark model - NAV05 (9 x 9)
#define CIRCLE_ONE_LOC		10		//CIRCLE model with one location
#define CIRCLE_TWO_LOC		11		//CIRCLE model with two locations
#define CIRCLE_FOUR_LOC		12		//CIRCLE model with four locations
#define OSCILLATOR			13		//OSCILLATOR model without any Filter

//**************** Hybrid Automata Definition ***********************

// *********** Different Template Directions ******************
#define BOX			1
#define OCT			2
//#define UNIFORM		//values > 3


#define CORE 	8	//Total number of existing cores in the working PC; Although we detect/read from program in Load Balance Algorithm
#define MAX_ALGO 	7	//Total number of existing Algorithm. This value need to be increased for new Algorithm when implemented


#endif /* ALL_PP_DEFINITION_H_ */
