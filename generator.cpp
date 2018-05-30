
#include "generator.h"

#include "itm.hpp"

#include "location.h"
#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_sensor.h"
#include "utils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <unistd.h>

// Earth's radius in meters
#define EARTH_RADIUS 6378000.0

void Generator::generateEntities(
		int num_pu, int num_ss, int num_su,
		std::vector<PU>* pus, std::vector<SS>* sss, std::vector<SU>* sus) const {
	// PU
	for(int j = 0; j < num_pu; ++j) {
		float pu_x = utils::randomFloat(0.0, location_range);
		float pu_y = utils::randomFloat(0.0, location_range);
	
		// Both in mW.
		float transmit_power = 2.0;
		float threshold = utils::fromdBm(utils::randomFloat(-130.0, -80.0));
		float pu_height = 100.0;

		if(utils::unit_type == utils::UnitType::ABS) {
			// Do Nothing
		} else if(utils::unit_type == utils::UnitType::DB) {
			transmit_power = utils::todBm(transmit_power);
			threshold = utils::todBm(threshold);
		} else {
			std::cerr << "Unsupported unit_type" << std::endl;
			exit(1);
		}

		PU pu(Location(pu_x, pu_y, pu_height), transmit_power, threshold);
		(*pus).push_back(pu);
	}

	// SS
	for(int i = 0; i < num_ss; ++i) {
		float x = utils::randomFloat(0.0, location_range);
		float y = utils::randomFloat(0.0, location_range);
		float ss_height = 10.0;
		Location loc(x, y, ss_height);

		float rp = 0.0;
		bool first = true;
		for(int j = 0; j < num_pu; ++j) {
			float path_loss = pm->getPathLoss((*pus)[j].loc, loc);
			if(utils::unit_type == utils::UnitType::ABS) {
				rp += path_loss * (*pus)[j].transmit_power;
			} else if(utils::unit_type == utils::UnitType::DB) {
				if(first) {
					first = false;
					rp = path_loss + (*pus)[j].transmit_power;
				} else {
					rp = utils::todBm(utils::fromdBm(rp) + utils::fromdBm(path_loss + (*pus)[j].transmit_power));
				}
			}
		}

		SS ss(loc, rp);
		(*sss).push_back(ss);
	}

	// SU
	for(int i = 0; i < num_su; ++i) {
		float x = utils::randomFloat(location_range / 5.0, 4.0 * location_range / 5.0);
		float y = utils::randomFloat(0.0, location_range);
		float su_height = 10.0;
		SU su(Location(x, y, su_height));
		(*sus).push_back(su);
	}
}

float Generator::computeGroundTruth(const SU& su, const std::vector<PU>& pus, float* ground_truth_path_loss) const {
	float tp = 0.0;
	bool first = true;
	for(unsigned int i = 0; i < pus.size(); ++i) {
		for(unsigned int j = 0; j < pus[i].prs.size(); ++j) {
			float this_tp = 0.0;
			float this_pl = pm->getPathLoss(su.loc, pus[i].prs[j].loc);;
			if(utils::unit_type == utils::UnitType::ABS) {
				this_tp = pus[i].prs[j].threshold / this_pl;
			} else if(utils::unit_type == utils::UnitType::DB) {
				this_tp = pus[i].prs[j].threshold - this_pl;
			} else {
				std::cerr << "Unsupported unit_type" << std::endl;
				exit(1);
			}
			if(first) {
				tp = this_tp;
				first = false;

				if(ground_truth_path_loss != nullptr) {
					*ground_truth_path_loss = this_pl;
				}
			} else if(this_tp < tp) {
				tp = this_tp;

				if(ground_truth_path_loss != nullptr) {
					*ground_truth_path_loss = this_pl;
				}
			}
		}
	}
	return tp;
}

float LogDistancePM::getPathLoss(const Location& loc1, const Location& loc2) const {
	float pl_dBm = pl0 + 10 * gamma * log10(loc1.dist(loc2) / d0);
	
	if(utils::unit_type == utils::UnitType::ABS) {
		return utils::fromdBm(-pl_dBm);
	} else if(utils::unit_type == utils::UnitType::DB) {
		return -pl_dBm;
	} else {
		std::cerr << "Unsupported unit_type" << std::endl;
		exit(1);
		return 0.0;
	}
}

float LongleyRicePM::getPathLoss(const Location& loc1, const Location& loc2) const {
	// Changes directory to splat folder
	chdir(splat_dir.c_str());

	// SDF dir

	// Convert location to longitude and lattitude
	float loc1_lat  = ref_lat  + (loc1.y / EARTH_RADIUS) * (180.0 / M_PI);
	float loc1_long = ref_long + (loc1.x / EARTH_RADIUS) * (180.0 / M_PI) / cos(ref_lat * M_PI / 180.0);

	float loc2_lat  = ref_lat  + (loc2.y / EARTH_RADIUS) * (180.0 / M_PI);
	float loc2_long = ref_long + (loc2.x / EARTH_RADIUS) * (180.0 / M_PI) / cos(ref_lat * M_PI / 180.0);

	// Create the needed files
	std::string loc1_fname = "loc1.qth";
	std::string loc2_fname = "loc2.qth";

	std::ofstream loc1_qth(loc1_fname);
	loc1_qth << "Loc1" << std::endl << loc1_lat << std::endl << loc1_long << std::endl << loc1.z << " m" << std::endl;
	loc1_qth.close();

	std::ofstream loc2_qth(loc2_fname);
	loc2_qth << "Loc2" << std::endl << loc2_lat << std::endl << loc2_long << std::endl << loc2.z << " m" << std::endl;
	loc2_qth.close();

	// Run the splat command
	std::stringstream sstr;
	float timeout_seconds = 300;
	sstr << "timeout " << timeout_seconds << "s splat-hd -t " << loc1_fname << " -r " << loc2_fname << " -d " << sdf_dir << " -L -metric -R 10 > /dev/null 2>&1";

	int sys_code = 1;
	int num_tries = 2;
	for(int i = 0; sys_code != 0 && i < num_tries; ++i) {
		sys_code = system(sstr.str().c_str());
	}

	if(sys_code != 0) {
		std::cerr << "Unable to run Splat-hd" << std::endl;
		exit(1);
	}

	// Read in the output file and get the path loss
	std::ifstream res("Loc1-to-Loc2.txt");

	std::string line = "";
	std::string free_space_prefix = "Free space path loss: ";
	std::string itwom_prefix = "ITWOM Version 3.0 path loss: ";
	float free_space_loss = 0.0;
	float itwom_loss = 0.0;
	while(std::getline(res, line)) {
		if(line.compare(0, free_space_prefix.size(), free_space_prefix) == 0) {
			free_space_loss = atof(line.substr(free_space_prefix.size(), line.size() - free_space_prefix.size() - 3).c_str());
		}
		if(line.compare(0, itwom_prefix.size(), itwom_prefix) == 0) {
			itwom_loss = atof(line.substr(itwom_prefix.size(), line.size() - itwom_prefix.size() - 3).c_str());
		}
	}
	res.close();
	chdir(return_dir.c_str());

	// Returns the greater of the two values.
	float rv = free_space_loss;
	if(itwom_loss > rv) {
		rv = itwom_loss;
	}

	if(utils::unit_type == utils::UnitType::ABS) {
		return utils::fromdBm(-rv);
	} else if(utils::unit_type == utils::UnitType::DB) {
		return -rv;
	} else {
		std::cerr << "Unsupported unit_type" << std::endl;
		exit(1);
		return 0.0;
	}
}
