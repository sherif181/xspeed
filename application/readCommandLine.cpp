/*
 * readCommandLine.cpp
 *
 *  Created on: 28-Oct-2016
 *      Author: hazel
 */

#include <application/readCommandLine.h>
namespace po = boost::program_options;

int readCommandLine(int argc, char *argv[], userOptions& user_options,
		hybrid_automata& Hybrid_Automata,
		std::list<initial_state::ptr>& init_state,
		ReachabilityParameters& reach_parameters) {

	std::pair<int, polytope::ptr> forbidden_set; //(locID1,Polytope1)}
	//	std::string bad_state; // string to capture the bad state description given by the user
	//	std::list<abstractCE::ptr> ce_candidates; //object of class counter_example
	//userOptions user_options;

	//int number_of_times = 1; //Make this 1 for Memory Profiling
	unsigned int number_of_streams = 1;
	//	int lp_solver_type_choosen = 1; //	1 for GLPK and 2 for Gurobi
	//	int Solver_GLPK_Gurobi_GPU = 3; //if Algorithm == 11 then (1 for GLPK; 2 for Gurobi; 3 for GPU)
	unsigned int Total_Partition; //for Parallel Iterations Algorithm :: number of partitions/threads

	bool isConfigFileAssigned = false, isModelParsed = false;

	po::options_description desc("XSpeed options");
	po::variables_map vm;

	if (argc == 1) { //No argument:: When Running directly from the Eclipse Editor
		std::cout << "Missing arguments!\n";
		std::cout << "Try XSpeed --help to see the command-line options\n";
		//exit(0);
		return 0;
	}

	desc.add_options()
	("help", "produce help message")
	("model", po::value<int>()->default_value(1), "set model for reachability analysis\n"
					"1.  Bouncing Ball Model: Variables{x,v} (Set to default)\n"
					"2.  Timed Bouncing Ball Model: Variables{x,v,t}\n"
					"3.  28-Dimensional Helicopter Controller Model: Variables{x1..x28}\n"
					"4.  Five dimensional Benchmark Model: Variables{x1..x5} \n"
					"5.  Navigation Benchmark Model-NAV01 (3 X 3): Variables{x1,x2,v1,v2}\n"
					"6.  Navigation Benchmark Model-NAV02 (3 X 3): Variables{x1,x2,v1,v2}\n"
					"7.  Navigation Benchmark Model-NAV03 (3 X 3): Variables{x1,x2,v1,v2}\n"
					"8.  Navigation Benchmark Model-NAV04 (5 X 5): Variables{x1,x2,v1,v2}\n"
					"9.  Navigation Benchmark Model-NAV05 (9 X 9): Variables{x1,x2,v1,v2}\n"
					"10. Circle with only ONE location model: Variables{x,y} \n"
					"11. Circle with TWO locations model: Variables{x,y} \n"
					"12. Circle with FOUR locations model: Variables{x,y} \n"
					"13. Oscillator model without any filters: Variables{x,y}\n"
					"14. Testing Model: Variables{depends on the model in test}\n")
	("directions", po::value<int>()->default_value(1), "Set the directions for template polyhedra:\n"
					"1. Box Directions (Set to default)\n"
					"2. Octagonal Directions \n"
					"n. 'n' uniform Directions \n")
	("time-horizon", po::value<double>(), "Set the Time horizon for the flowpipe computation per Location(Local time).")
	("time-step", po::value<double>(), "Set the sampling time for the flowpipe computation.")
	("transition-size", po::value<int>(), "Set the maximum number of Jumps(0 for no jump(breadth=1), 1 for first jump(breadth=2).")

	("algo,a",po::value<int>()->default_value(1), "Set the algorithm\n"
			"1/seq-SF -- Sequential Algorithm (both PostC and PostD are sequential) (Set to default)\n"
			"2/par-SF -- Lazy evaluation algorithm (Parallel PostC but Sequential PostD/BFS)\n"
			"3/time-slice-SF -- Time slicing algorithm (Parallel PostC but Sequential PostD/BFS)\n"
			"4/AGJH -- Adaptation of Gerard J. Holzmann\n"
			"5/TPBFS -- Load Balancing Algorithm\n"
			"6/gpu-postc -- Bounded input sets (PostC in GPU but Sequential PostD/BFS)\n")
	("number-of-streams", po::value<int>()->default_value(1), "Set the maximum number of GPU-streams (Set to 1 by default).")
	("time-slice", po::value<int>(), "Set the maximum number of Time Sliced(or partitions)for algo=time-slice-SF")
	("internal", "called internally when running hyst-xspeed model")
	("forbidden,F", po::value<std::string>(), "forbidden location_ID and forbidden set/region within that location") //better to be handled by hyst
	("include-path,I", po::value<std::string>(), "include file path")
	("model-file,m", po::value<std::string>(), "include model file")
	("config-file,c", po::value<std::string>(), "include configuration file")
	("output-file,o", po::value<std::string>(), "output file name for redirecting the outputs")
	("output-variable,v", po::value<std::string>(), "projecting variables for e.g., 'x,v' for Bouncing Ball") //better to be handled by hyst
;

	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	std::vector<std::string> output_vars(3); //stores the output/plotting variables

	//std::ofstream outFile;
	//math::matrix<double> vertices_list;
	const char *stFileNameWithPath;
	std::string input;
	for (int i = 1; i < argc; i++) {
		if (((std::string(argv[i]).find("-m") != string::npos)
				|| (std::string(argv[i]).find("-model-file") != string::npos)))
			i++;
		else if (((std::string(argv[i]).find("-c") != string::npos)
				|| (std::string(argv[i]).find("-config-file") != string::npos)))
			i++;
		else if (std::string(argv[i]).find("-o") != string::npos)
			i++;
		//if((!((std::string(argv[i]).find("-m") !=string::npos) ||(std::string(argv[i]).find("-model-file") !=string::npos))
		//		|| !((std::string(argv[i]).find("-c") !=string::npos) ||(std::string(argv[i]).find("-config-file") !=string::npos))) )
		else {
			input.append(argv[i]);
			input.append(" ");
		}
	}
	//std::cout << "The input is:" << input << std::endl;
	//std::cout << "argc is:" << argc << std::endl;
	if (argc > 1) { // Boost Options to be filled-up
		if (vm.count("help")) {
			cout << desc << "\n";
			exit(0);
			//return 0; //Just displaying the Help options and terminating the Project
		}
		std::string include_path = "", model_filename = "",
				config_filename = ""; //default set to empty
		if (vm.count("include-path")) {
			include_path = vm["include-path"].as<std::string>();
			//std::cout << "Include Path is: " << include_path << "\n";
		}

		if (vm.count("config-file")) {

			user_options.set_configFile(vm["config-file"].as<std::string>());
			//std::cout << "Configuration file is: " << user_options.get_configFile() << "\n";
		}
		if (vm.count("model-file")) {

			user_options.set_modelFile(vm["model-file"].as<std::string>());
			//std::cout << "Model file is: " << user_options.get_modelFile() << "\n";
		}
		// ********************** Setting for Output file **********************************
		std::string fileName, fullPath, fileWithPath;
		//std::string fileName, fullPath;
		//std::string fullPath,fileWithPath;
		//cout << endl << "Working here 2\n";
		if (vm.count("include-path")) {
			fullPath = vm["include-path"].as<std::string>();
			//std::cout << "Include Path is: " << fullPath << "\n";
		} else {
			fullPath = "./"; //default file path
		}
		fileWithPath.append(fullPath);
		if (vm.count("output-file")) {
			fileName = vm["output-file"].as<std::string>();
			//std::cout << "fileName is: " << fileName << "\n";
		} else {
			fileName = "out.txt";
		}
		fileWithPath.append(fileName);
		//std::cout << "fileWithPath is: " << fileWithPath << "\n";
		stFileNameWithPath = fileWithPath.c_str();
		//std::cout << "fileWithPath is: " << fileWithPath << "\n";
		//std::cout << "stFileNameWithPath = " << stFileNameWithPath << "\n";
		// ********************** Setting for Output file Done **********************************

		if (vm.count("model-file") && vm.count("config-file")) {
			std::cout << "Translating user model with Hyst\n";

			std::string cmdStr, replacingFile, SingleSpace = " ", projLocation,
					java_exeFile;

			//todo:: proper path to be handled from the relative/current installed location of the software
			replacingFile = "./user_model.cpp";

			java_exeFile = "java -jar";

			/*Old Implementation
			 *
			 cmdStr.append(java_exeFile);
			 cmdStr.append(SingleSpace);
			 cmdStr.append("./bin/Hyst-XSpeed.jar -xspeed -o");	//OLD hyst-Repo
			 cmdStr.append(SingleSpace);
			 cmdStr.append(replacingFile);
			 cmdStr.append(SingleSpace);
			 cmdStr.append(include_path);
			 cmdStr.append(user_options.get_modelFile());
			 cmdStr.append(SingleSpace);
			 cmdStr.append(include_path);
			 cmdStr.append(user_options.get_configFile());
			 const char *st, *st2, *st3, *st4, *st5;
			 st = cmdStr.c_str();
			 //std::cout <<"st = "<<st<<std::endl;
			 system(st); //calling hyst interface to generate the XSpeed model file
			 system("g++ -c -I./include/ user_model.cpp -o user_model.o");
			 system(
			 "g++ -L./lib/ user_model.o -lXSpeed -lboost_timer -lboost_system -lboost_chrono -lboost_program_options -lgomp -lglpk -lsundials_cvode -lsundials_nvecserial -lnlopt -o ./XSpeed.o");
			 //std::cout<<"file with path =" <<stFileNameWithPath<<std::endl;

			 */

			cmdStr.append(java_exeFile);
			cmdStr.append(SingleSpace);
			cmdStr.append("./bin/Hyst-XSpeed.jar -t xspeed \"\" -o");
			cmdStr.append(SingleSpace);
			cmdStr.append(replacingFile);
			cmdStr.append(" -i ");
			cmdStr.append(user_options.get_modelFile());//cmdStr.append(modelFile);
			cmdStr.append(" ");
			cmdStr.append(user_options.get_configFile());//cmdStr.append(configFile);
			const char *st, *st2, *st3, *st4, *st5;
			st = cmdStr.c_str();
			system(st); //calling hyst interface to generate the XSpeed model file
			system("g++ -c -I./include/ user_model.cpp -o user_model.o");
			system(
					"g++ -L./lib/ user_model.o -lXSpeed -lboost_timer -lboost_system -lboost_chrono -lboost_program_options -lgomp -lglpk -lsundials_cvode -lsundials_nvecserial -lnlopt -o ./XSpeed.o");

			string cmdStr1;
			cmdStr1.append("./XSpeed.o --internal -o");
			cmdStr1.append(SingleSpace);
			cmdStr1.append(stFileNameWithPath);
			cmdStr1.append(SingleSpace);
			cmdStr1.append(input);
			//std::cout<<"Command = "<<cmdStr1.c_str()<<"\n";
			//system("./XSpeed.o --internal");
			system(cmdStr1.c_str());
			exit(0);
		}

		if (vm.count("internal")) {
			// calls the hyst-xspeed generated model
			//int iters = (int) user_options.get_timeHorizon()/ user_options.get_timeStep();

			//Todo:: have to take inputs for options such as user_options.get_flow_algorithm() == Algorithm_Type,
			//		time-slice, user_options.get_automata_exploration_algorithm()==jumps
			//     all related inputs for gpu such as number-of-streams

			//user_model(Hybrid_Automata, init_state, reach_parameters, forbidden_set, user_options);

			user_model(Hybrid_Automata, init_state, reach_parameters,
					user_options);

			/*unsigned int x1 = Hybrid_Automata.get_index(output_vars[0]);
			 unsigned int x2 = Hybrid_Automata.get_index(output_vars[1]);
			 user_options.set_first_plot_dimension(x1);
			 user_options.set_second_plot_dimension(x2);*/ //Todo have to take care of this as well
			isConfigFileAssigned = true;
			isModelParsed = true;
		}

		//if (vm.count("directions") && isConfigFileAssigned == false) { //Compulsory Options but set to 1 by default
		if (vm.count("directions")) { //Compulsory Options but set to 1 by default
			user_options.set_directionTemplate(vm["directions"].as<int>());
			if (user_options.get_directionTemplate() <= 0) {
				std::cout << "Invalid Directions option specified\n";
				return 0;
			}
		}
		if (vm.count("output-variable")) {
			std::string VarStr;
			VarStr = vm["output-variable"].as<std::string>();

			typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
			boost::char_separator<char> sep(", ");
			tokenizer tokens(VarStr, sep);
			int index = 0;
			for (tokenizer::iterator tok_iter = tokens.begin();
					tok_iter != tokens.end(); ++tok_iter) {
				output_vars[index] = (std::string) (*tok_iter);
				//std::cout << "Output Variable = " << output_vars[index] << "\n";
				index++;
			}
		}

		//if (vm.count("forbidden") && isConfigFileAssigned == false) { //Compulsory Options but set to 1 by default
		if (vm.count("forbidden")) { //Compulsory Options but set to 1 by default
			user_options.set_forbidden_state(vm["forbidden"].as<std::string>());
		}

		//if (vm.count("model") && isConfigFileAssigned == false) { //Compulsory Options but set to 1 by default
		if (vm.count("model")) { //Compulsory Options but set to 1 by default
			user_options.set_model(vm["model"].as<int>());
			//if (user_options.get_model() < 1 || user_options.get_model() > 13) {
			//For testing model==14
			if (user_options.get_model() < 1 || user_options.get_model() > 14) {
				std::cout << "Invalid Model option specified\n";
				return 0;
			}
		}

		//if (vm.count("time-horizon") && isConfigFileAssigned == false) { //Compulsory Options
		if (vm.count("time-horizon")) { //Compulsory Options
			user_options.set_timeHorizon(vm["time-horizon"].as<double>());
			if (user_options.get_timeHorizon() <= 0) { //for 0 or negative time-bound
				std::cout
						<< "Invalid time-horizon option specified, A positive non zero bound expected\n";
				return 0;
			}
		} else if (isConfigFileAssigned == false) {
			std::cout << "Missing time-horizon option\n";
			return 0;
		}
		//if (vm.count("time-step") && isConfigFileAssigned == false) { //Compulsory Options
		if (vm.count("time-step")) { //Compulsory Options
			user_options.set_timeStep(vm["time-step"].as<double>());
			if (user_options.get_timeStep() > 0) {
				//std::cout << "\niterations_size = " << iterations_size;
			} else { //for 0 or negative sampling-time
				std::cout << "Invalid time-step option specified\n";
				return 0;
			}
		} else if (isConfigFileAssigned == false) {
			std::cout << "Missing time-step option\n";
			return 0;
		}
		//if (vm.count("transition-size") && isConfigFileAssigned == false) { //Compulsory Options
		if (vm.count("transition-size")) { //Compulsory Options
			user_options.set_bfs_level(vm["transition-size"].as<int>());
			if (user_options.get_bfs_level() < 0) {
				std::cout
						<< "Invalid bfs level specified, a positive number expected\n";
				return 0;
			}
		} else if (isConfigFileAssigned == false) {
			std::cout << "Missing transition-size option\n";
			return 0;
		}

		//Algorithm Preference
		if (vm.count("algo")) {
			user_options.set_algorithm(vm["algo"].as<int>());
			//std::cout << "algo = " << user_options.get_algorithm() << std::endl;
			if (user_options.get_algorithm() < 0
					|| user_options.get_algorithm() > MAX_ALGO) {
				std::cout << "Invalid algorithm option specified\n";
				return 0;
			}
		} else {
			std::cout << "Missing algo option\n";
			return 0;
		}
		if (user_options.get_algorithm() == 3) { //this argument will be set only if algorithm==time-slice or PAR_ITER
			if (vm.count("time-slice")) { //Compulsory Options if algorithm-type==Time-Slice(4)
				int partition_size = vm["time-slice"].as<int>();
				if (partition_size > 0) {
					//Total_Partition = partition_size;
					user_options.setTotalSliceSize(partition_size);
				} else { //for 0 or negative time-slice
					std::cout << "Invalid time-slice option specified\n";
					return 0;
				}
			} else {
				std::cout << "Missing time-slice option\n";
				return 0;
			}
		}
		if (user_options.get_algorithm() == 6) { //if gpu enabled then
			if (vm.count("number-of-streams")) { //Compulsory Options but set 1 by default
				int no_streams = vm["number-of-streams"].as<int>(); //default value ==1
				if (no_streams >= 1) {
					//number_of_streams = no_streams; //Number of GPU-Streams selected
					user_options.setStreamSize(no_streams);
				} else {
					std::cout << "Invalid number_of_streams option specified\n";
					return 0;
				}
			}
		}
	} //ALL COMMAND-LINE OPTIONS are set completely

	// Initialize the model with the parameters given by the user
	if (!isModelParsed) { //all command line options has been supplied
		load_model(init_state, Hybrid_Automata, user_options, reach_parameters,
				forbidden_set);

		unsigned int x1 = Hybrid_Automata.get_index(output_vars[0]);
		unsigned int x2 = Hybrid_Automata.get_index(output_vars[1]);
		if (!(output_vars[2].empty())) {
			unsigned int x3 = Hybrid_Automata.get_index(output_vars[2]);
			user_options.set_third_plot_dimension(x3);
		}

		user_options.set_first_plot_dimension(x1);
		user_options.set_second_plot_dimension(x2);

		if (!user_options.get_forbidden_state().empty())
			forbidden_set.second->print2file("./bad_poly", x1, x2);

	}
	return 1;
}

