
#include "args.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string.h>

Args::Args(const int &argc, char const *argv[]) {
	// Checks for a '-config' flag
	std::string config_file = "";
	for(int i = 1; i < argc; ++i) {
		if(strcmp(argv[i], "-config") == 0 && i + 1 < argc) {
			config_file = argv[i + 1];
			break;
		}
	}

	std::vector<std::string> cmd_args;
	if(config_file != "") {
		std::ifstream ifstr(config_file, std::ifstream::in);
		std::string token = "";

		while(ifstr >> token) {
			cmd_args.push_back(token);
		}
	}
	for(int i = 1; i < argc; ++i) {
		if(strcmp(argv[i], "-config") == 0 && i + 1 < argc) {
			++i;
			continue;
		}
		cmd_args.push_back(argv[i]);
	}

	// Initialize arguments and their flags
	args = {Arg_Val(&rand_seed, 0, "-rand_seed", "RAND_SEED", "Value to seed RNG with."),
			Arg_Val(&brief_out, "-brief_out", "If given, only outputs minimal info."),
			Arg_Val(&skip_s2pc, "-skip_s2pc", "If given, skips the secure algorithm and only runs the plaintext algo."),
			Arg_Val(&central_entities, "two_sms", "-ces", "CENTRAL_ENTITIES", "The central entities that will run the secure algorithm. Must be 'two_sms', or 'sm_ks'"),
			Arg_Val(&use_gt_rp_at_ss_from_pu, "-use_gt_rp", "If given, then the plaintext algorithm will use the real RP at each SS from PU instead of splitting the total using IDW."),
			Arg_Val(&num_pu, 1, "-npu", "NUM_PU", "Number of Primary Users to generate."),
			Arg_Val(&num_ss, 1, "-nss", "NUM_SS", "Number of Spectrum Sensors to generate."),
			Arg_Val(&num_su, 1, "-nsu", "NUM_SU", "Number of Secondary Users to generate."),
			Arg_Val(&location_range, 100.0, "-lr", "LOC_RANGE", "Range of locations entities will be placed"),
			Arg_Val(&num_pr_per_pu, 1, "-npr", "NUM_PR_PER_PU", "Number of PRs to generate per PU."),
			Arg_Val(&pr_range, 0.0, "-prr", "PR_RANGE", "Range to place PRs around PUs."),
			Arg_Val(&unit_type, "abs", "-ut", "UNIT_TYPE", "Unit type to use for calculations. Must be either \"abs\" or \"db\""),
			Arg_Val(&out_filename, "", "-out", "FILENAME", "File to write generated data to."),
			Arg_Val(&in_filename, "", "-inf", "FILENAME", "File to read generated data from."),
			Arg_Val(&propagation_model, "log_distance", "-pm", "PROPAGATION_MODEL", "Ground truth Propagation Model to use"),
			Arg_Val(&ld_path_loss0, 1.0, "-ld_pl0", "REFERENCE_PATH_LOSS", "For Log-Distance PM, the reference path loss (in dBm)"),
			Arg_Val(&ld_dist0, 5.0, "-ld_d0", "REFERENCE_DISTANCE", "For Log-Distance PM, the reference distance."),
			Arg_Val(&ld_gamma, 1.0, "-ld_g", "LD_GAMMA", "For Log-Distance PM, the gamma parameter."),
			Arg_Val(&splat_cmd, "splat", "-splat_cmd", "SPLAT_CMD", "Command to run SPLAT!. Must be either 'splat' or 'splat-hd'"),
			Arg_Val(&ref_lat, 40.75, "-ref_lat", "LATTITUDE", "Reference lattitude used to generate locations of entitites."),
			Arg_Val(&ref_long, 73.25, "-ref_long", "LONGITUDE", "Reference longitude used to generate locations of entitites."),
			Arg_Val(&splat_dir, "splat/", "-splat_dir", "DIR", "Directory of Splat data files"),
			Arg_Val(&sdf_dir, "sdf/", "-sdf_dir", "DIR", "Directory (relative to splat_dir) of sdf data."),
			Arg_Val(&return_dir, "../", "-return_dir", "DIR", "Directory (relative to splat_dir) of current location."),
			Arg_Val(&use_itwom_pl, "-itwom", "If given uses the itwom path loss from SPLAT's output"),
			Arg_Val(&num_ss_selection, 0, "-nss_s", "NUM_SS_SELECTION", "For each SU request, uses only the specified number of nearby SS. A value of 0 uses all SS."),
			Arg_Val(&num_pu_selection, 0, "-npu_s", "NUM_PU_SELECTION", "For each SU request, uses only the specified number of nearby PU. A value of 0 uses all SS."),
			Arg_Val(&do_plaintext_split, "-do_pt_split", "If given, will do preprocessing split of SS receive power in plaintext."),
			Arg_Val(&no_pr_thresh_update, "-no_pr_up", "If given, then pr thresholds aren't updated."),
			Arg_Val(&selection_algo, "sort", "-sel_algo", "ALGO", "Selection Algo used for SS. Must be \"none\", \"sort\", or \"random\""),
			Arg_Val(&secure_write_algo, "proposed", "-sec_write_algo", "ALGO", "ALgorithm for secure write. Must be either \"proposed\" or \"spc\""),
			Arg_Val(&grid_num_x, 0, "-grid_x", "NUM_X", "The number of times to divide the x dimension"),
			Arg_Val(&grid_num_y, 0, "-grid_y", "NUM_Y", "The number of times to divide the y dimension"),
			Arg_Val(&ss_receive_power_alpha, 1.0, "-rpa", "RECEIVE_POWER_ALPHA", "Parameter used when estimated the received power from PUs at SSs"),
			Arg_Val(&path_loss_alpha, 1.0, "-pla", "PATH_LOSS_ALPHA", "Parameter used when estimating the path loss between PRs and SUs."),
			Arg_Val(&num_float_bits, 16, "-float_bits", "NUM_FLOAT_BITS", "Number of bits of precision to use during the S2-PC calculations."),
			Arg_Val(&s2_pc_bit_count, 64, "-bit_count", "BIT_COUNT", "Number of total bits to use in S2-PC calculations."),
			Arg_Val(&pl_est_gamma, 1.0, "-plg", "PATH_LOSS_GAMMA", "Gmma used in path loss estimation.")
	};

	bool err = false;
	for(unsigned int i = 0; i < cmd_args.size(); ++i) {
		if(cmd_args[i] == "-h" || cmd_args[i] == "-help") {
			printHelp();
		}
		bool match_found = false;
		for(unsigned int j = 0; j < args.size(); ++j) {
			if(args[j].isMatch(cmd_args, i)) {
				args[j].setVal(cmd_args, i);
				match_found = true;
				break;
			}
		}
		// IF NONE FOUND PRINT AN ERROR MESSAGE
		if(!match_found) {
			std::cerr << "INVALID COMMAND LINE ARG: " << cmd_args[i] << std::endl;
			err = true;
		}
	}
	if(err) {
		for(unsigned int i = 0; i < cmd_args.size(); ++i) {
			std::cerr << cmd_args[i];
			if(i < cmd_args.size() - 1) {
				std::cerr << " ";
			}
		}
		std::cerr << std::endl;
		exit(0);
	}
}

void Args::printHelp() const {
	for(unsigned int i = 0; i < args.size(); ++i) {
		std::cout << " " << args[i].flag_str << " " << args[i].input_name << std::endl;
		std::cout << "    " << args[i].help_text << std::endl << std::endl;
	}

	// Prints the -config flag
	std::cout << " -config CONFIG_FILE" << std::endl;
	std::cout << "    Config file containing flags. Command line flags override values in config file." << std::endl;
	exit(0);
}

bool Args::Arg_Val::isMatch(const std::vector<std::string> &cmd_args, unsigned int &i) {
	if(cmd_args[i] == flag_str) {
		if(val_type == "f" || val_type == "i" || val_type == "s") {
			if(i >= cmd_args.size() - 1) {
				return false;
			}
			i++;
		}
		return true;
	}
	return false;
}

void Args::Arg_Val::setVal(const std::vector<std::string> &cmd_args, unsigned int &i) {
	if(val_type == "f") {
		*((float*) val_ptr) = atof(cmd_args[i].c_str());
	}
	if(val_type == "i") {
		*((int*) val_ptr) = atoi(cmd_args[i].c_str());
	}
	if(val_type == "s") {
		*((std::string*) val_ptr) = cmd_args[i];
	}
	if(val_type == "b") {
		*((bool*) val_ptr) = true;
	}
}
