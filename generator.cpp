
#include "generator.h"

#include "timer.h"
#include "location.h"
#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_sensor.h"
#include "utils.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <utility>

// Earth's radius in meters
#define EARTH_RADIUS 6378000.0

const float ss_height = 10.0;
const float su_height = 10.0;
const float pu_height = 50.0;

void Generator::generateEntities(
		int num_pu, int num_ss, int num_su, int num_pr_per_pu, float pr_range, const std::string & out_filename,
		std::vector<PU>* pus, std::vector<SS>* sss, std::vector<SU>* sus) const {
	const float pr_height = (num_pr_per_pu == 1 ? pu_height : ss_height);
	// SS
	for(int i = 0; i < num_ss; ++i) {
		if(out_filename != "") {
			std::cout << "Generating SS " << i + 1 << " of " << num_ss << std::endl;
		}

		float x = utils::randomFloat(0.0, location_range);
		float y = utils::randomFloat(0.0, location_range);

		Location loc(x, y, ss_height);

		SS ss(loc, 0.0);
		(*sss).push_back(ss);
	}

	// PU
	for(int j = 0; j < num_pu; ++j) {
		if(out_filename != "") {
			std::cout << "Generating PU " << j + 1 << " of " << num_pu << std::endl;
		}

		float pu_x = utils::randomFloat(0.0, location_range);
		float pu_y = utils::randomFloat(0.0, location_range);
	
		// Both in mW.
		float transmit_power = 2.0;

		if(utils::unit_type == utils::UnitType::ABS) {
			// Do Nothing
		} else if(utils::unit_type == utils::UnitType::DB) {
			transmit_power = utils::todBm(transmit_power);
		} else {
			std::cerr << "Unsupported unit_type" << std::endl;
			exit(1);
		}

		// PRs
		PU pu(Location(pu_x, pu_y, pu_height), transmit_power);
		for(int x = 0; x < num_pr_per_pu; ++x) {
			Location pr_loc;
			if(num_pr_per_pu > 1) {
				float pr_x = utils::randomFloat(-pr_range, pr_range) + pu.loc.x;
				float pr_y = utils::randomFloat(-pr_range, pr_range) + pu.loc.y;
				pr_loc = Location(pr_x, pr_y, pr_height);
				while(pr_loc.dist(pu.loc) > pr_range + 0.001) {
					pr_x = utils::randomFloat(-pr_range, pr_range) + pu.loc.x;
					pr_y = utils::randomFloat(-pr_range, pr_range) + pu.loc.y;
					pr_loc = Location(pr_x, pr_y, pr_height);
				}
			} else if(num_pr_per_pu == 1) {
				pr_loc = pu.loc;
			} else {
				std::cerr << "Invalid valid for num_pr_per_pu: " << num_pr_per_pu << std::endl;
				exit(1);
			}

			float threshold = utils::fromdBm(utils::randomFloat(-130.0, -80.0));

			if(utils::unit_type == utils::UnitType::ABS) {
				// Do Nothing
			} else if(utils::unit_type == utils::UnitType::DB) {
				threshold = utils::todBm(threshold);
			} else {
				std::cerr << "Unsupported unit_type" << std::endl;
				exit(1);
			}

			PR pr(pr_loc, threshold);

			pu.prs.push_back(pr);
		}

		// Run single and get data
		if(out_filename == "") {
			pm->preprocessPathLoss(&pu, ss_height, j);
			pm->loadANOFile(pu);
		}
		
		(*pus).push_back(pu);

		// Loop over SS and get the path loss between the PU and SS
		for(int i = 0; i < num_ss; ++i) {
			float path_loss = 1.0;
			if(out_filename == "") {
				path_loss = pm->getPathLoss((*pus)[j].loc, (*sss)[i].loc);
			}
			if(utils::unit_type == utils::UnitType::ABS) {
				(*sss)[i].received_power += path_loss * (*pus)[j].transmit_power;
			} else if(utils::unit_type == utils::UnitType::DB) {
				if(j == 0) {
					(*sss)[i].received_power = path_loss + (*pus)[j].transmit_power;
				} else {
					(*sss)[i].received_power = utils::todBm(utils::fromdBm((*sss)[i].received_power) + utils::fromdBm(path_loss + (*pus)[j].transmit_power));
				}
			}
		}
	}

	// SU
	float su_buffer = 0.0;
	for(int i = 0; i < num_su; ++i) {
		if(out_filename != "") {
			std::cout << "Generating SU " << i + 1 << " of " << num_su << std::endl;
		}

		float x = utils::randomFloat(su_buffer, location_range - su_buffer);
		float y = utils::randomFloat(su_buffer, location_range - su_buffer);
		

		Location su_loc(x, y, su_height);

		SU su(i, Location(x, y, su_height));
		(*sus).push_back(su);
	}

	if(out_filename != "") {
		outputEntities(out_filename, *pus, *sss, *sus);
		exit(0);
	}
}

void Generator::outputEntities(const std::string& out_filename, std::vector<PU>& pus, const std::vector<SS>& sss, const std::vector<SU>& sus) const {
	std::ofstream out(out_filename);

	// Write out Locations
	// PU + PRs
	for(unsigned int j = 0; j < pus.size(); ++j) {
		std::cout << "Outputting PU " << j + 1 << " of " << pus.size() << std::endl;
		out << "PU " << j << " " << pus[j].loc.x << " " << pus[j].loc.y << " " << pus[j].loc.z << " " << pus[j].transmit_power << std::endl;
		for(unsigned int x = 0; x < pus[j].prs.size(); ++x) {
			out << "PR " << j << " " << x << " " << pus[j].prs[x].loc.x << " " << pus[j].prs[x].loc.y << " " << pus[j].prs[x].loc.z << " " << pus[j].prs[x].threshold << std::endl;
		}
	}

	// SS
	for(unsigned int i = 0; i < sss.size(); ++i) {
		std::cout << "Outputting SS " << i + 1 << " of " << sss.size() << std::endl;
		out << "SS " << i << " " << sss[i].loc.x << " " << sss[i].loc.y << " " << sss[i].loc.z << std::endl;
	}

	// SU
	for(unsigned int i = 0; i < sus.size(); ++i) {
		std::cout << "Outputting SU " << i + 1 << " of " << sus.size() << std::endl;
		out << "SU " << i << " " << sus[i].loc.x << " " << sus[i].loc.y << " " << sus[i].loc.z << std::endl;
	}

	// Path losses
	// PU -> SS
	for(unsigned int j = 0; j < pus.size(); ++j) {
		std::cout << "Outputting PU_PL " << j + 1 << " of " << pus.size() << std::endl;
		pm->preprocessPathLoss(&(pus[j]), ss_height, j);
		pm->loadANOFile(pus[j]);
		for(unsigned int i = 0; i < sss.size(); ++i) {
			out << "PU_PL " << j << " " << i << " " << pm->getPathLoss(pus[j].loc, sss[i].loc) << std::endl;
		}
		// Delete file
		utils::deleteFile(pus[j].splat_ano_filename);
	}

	// PR -> SU
	int num_pr_pl = 0;
	for(unsigned int j = 0; j < pus.size(); ++j) {
		num_pr_pl += pus[j].prs.size();
	}

	int this_pr_pl = 1;
	for(unsigned int j = 0; j < pus.size(); ++j) {
		for(unsigned int x = 0; x < pus[j].prs.size(); ++x) {
			std::cout << "Outputting PR_PL " << this_pr_pl << " of " << num_pr_pl << std::endl;
			pm->preprocessPathLoss(&(pus[j].prs[x]), su_height, j, x);
			pm->loadANOFile(pus[j].prs[x]);
			for(unsigned int i = 0; i < sus.size(); ++i) {
				out << "PR_PL " << j << " " << x << " " << i << " " << pm->getPathLoss(pus[j].prs[x].loc, sus[i].loc) << std::endl;
			}
			// Delete PR file
			utils::deleteFile(pus[j].prs[x].splat_ano_filename);
			this_pr_pl++;
		}
	}
	out.close();
}

std::vector<float> Generator::computeGroundTruth(const std::vector<SU>& sus, const std::vector<PU>& pus, PathLossTable* path_loss_table) const {
	std::vector<std::vector<float> > path_losses;

	for(unsigned int j = 0; j < pus.size(); ++j) {
		pm->loadANOFile(pus[j]);
		path_losses.push_back(std::vector<float>());
		for(unsigned int i = 0; i < sus.size(); ++i) {
			float v = pm->getPathLoss(pus[j].loc, sus[i].loc);

			path_losses[j].push_back(v);
		}
	}

	std::vector<float> tps;
	for(unsigned int i = 0; i < sus.size(); ++i) {
		float tp = 0.0;
		bool first = true;
		for(unsigned int j = 0; j < pus.size(); ++j) {
			for(unsigned int x = 0; x < pus[j].prs.size(); ++x) {
				float this_tp = 0.0;

				float this_pl = path_losses[j][i];
				path_loss_table->addGroundTruthPathLoss(i, j, x, this_pl);
				if(pus[j].loc.dist(pus[j].prs[x].loc) > 0.01) {
					std::cerr << "Only supports RADAR system for now" << std::endl;
					exit(1);
				}
				
				if(utils::unit_type == utils::UnitType::ABS) {
					this_tp = pus[j].prs[x].threshold / this_pl;
				} else if(utils::unit_type == utils::UnitType::DB) {
					this_tp = pus[j].prs[x].threshold - this_pl;
				} else {
					std::cerr << "Unsupported unit_type" << std::endl;
					exit(1);
				}
				if(first) {
					tp = this_tp;
					first = false;
				} else if(this_tp < tp) {
					tp = this_tp;
				}
			}
		}
		tps.push_back(tp);
	}
	return tps;
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

	// Convert location to longitude and lattitude
	float loc1_lat  = ref_lat  + (loc1.y / EARTH_RADIUS) * (180.0 / M_PI);
	float loc1_long = ref_long + (loc1.x / EARTH_RADIUS) * (180.0 / M_PI) / cos(loc1_lat * M_PI / 180.0);

	float loc2_lat  = ref_lat  + (loc2.y / EARTH_RADIUS) * (180.0 / M_PI);
	float loc2_long = ref_long + (loc2.x / EARTH_RADIUS) * (180.0 / M_PI) / cos(loc2_lat * M_PI / 180.0);

	// Create the needed files
	std::string loc1_fname = "loc1.qth";
	std::string loc2_fname = "loc2.qth";

	std::ofstream loc1_qth(loc1_fname);
	loc1_qth << "Loc1" << std::endl << loc1_lat << std::endl << loc1_long << std::endl << loc1.z << "m" << std::endl;
	loc1_qth.close();

	std::ofstream loc2_qth(loc2_fname);
	loc2_qth << "Loc2" << std::endl << loc2_lat << std::endl << loc2_long << std::endl << loc2.z << "m" << std::endl;
	loc2_qth.close();

	// Run the splat command
	std::stringstream sstr;
	float timeout_seconds = 15;
	sstr << "timeout " << timeout_seconds << "s " << splat_cmd << " -t " << loc1_fname << " -r " << loc2_fname << " -d " << sdf_dir << " -L -metric -R 10 > /dev/null 2>&1";

	int sys_code = 1;
	int num_tries = 2;
	for(int i = 0; sys_code != 0 && i < num_tries; ++i) {
		sys_code = system(sstr.str().c_str());
	}

	if(sys_code != 0) {
		std::cerr << "Unable to run Splat, sys-code: " << sys_code << std::endl;
		std::cerr << "Command used: " << sstr.str() << std::endl;
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
	float path_loss = free_space_loss;
	if(use_itwom_pl && itwom_loss > path_loss) {
		path_loss = itwom_loss;
	}

	float low_pl_thresh = 55.0;
	if(path_loss <= low_pl_thresh) {
		path_loss = low_pl_thresh;
	}

	if(utils::unit_type == utils::UnitType::ABS) {
		return utils::fromdBm(-path_loss);
	} else if(utils::unit_type == utils::UnitType::DB) {
		return -path_loss;
	} else {
		std::cerr << "Unsupported unit_type" << std::endl;
		exit(1);
		return 0.0;
	}
}

void TransmitterOnlySplatPM::preprocessPathLoss(PU* pu, float ss_height, int pu_id) {
	// Changes directory to splat folder
	chdir(splat_dir.c_str());

	{
		std::stringstream sstr;
		sstr << splat_dir << "ano_files/pu_" << pu_id << ".ano";
		pu->splat_ano_filename = sstr.str();
	}

	// Convert location to longitude and lattitude
	float loc_lat  = ref_lat  + (pu->loc.y / EARTH_RADIUS) * (180.0 / M_PI);
	float loc_long = ref_long + (pu->loc.x / EARTH_RADIUS) * (180.0 / M_PI) / cos(loc_lat * M_PI / 180.0);

	std::string pu_fname = "PU.qth";
	
	std::ofstream pu_qth(pu_fname);
	pu_qth << "PU" << std::endl << loc_lat << std::endl << loc_long << std::endl << pu->loc.z << "m" << std::endl;
	pu_qth.close();

	std::stringstream sstr;
	float timeout_seconds = 30;
	sstr << "timeout " << timeout_seconds << "s " << splat_cmd << " -t " << pu_fname << " -d " << sdf_dir << " -L " << ss_height <<
			" -metric -R 20 -ano " << pu->splat_ano_filename << " > /dev/null 2>&1";

	int sys_code = 1;
	int num_tries = 2;
	for(int i = 0; sys_code != 0 && i < num_tries; ++i) {
		sys_code = system(sstr.str().c_str());
	}

	if(sys_code != 0) {
		std::cerr << "Unable to run Splat, sys-code: " << sys_code << std::endl;
		std::cerr << "Command used: " << sstr.str() << std::endl;
		exit(1);
	}
	chdir(return_dir.c_str());
	
	// loadANOFile(*pu);
}

void TransmitterOnlySplatPM::loadANOFile(const PU& pu) {
	if(pu.splat_ano_filename == "") {
		std::cerr << "Invalid filename" << std::endl;
		exit(1);
	}

	cur_loc = pu.loc;

	std::ifstream pu_ano(pu.splat_ano_filename);

	int line_num = 1;
	std::string line = "";
	cur_ano_data.clear();
	while(std::getline(pu_ano, line)) {
		if(line_num > 2) {
			// Remove last two characters if last character is '*'
			if(line.at(line.length() - 1) == '*') {
				line = line.substr(0, line.length() - 2);
			}

			// Split line by ', '
			std::string vs[5];
			std::stringstream sstr(line);
			sstr >> vs[0] >> vs[1] >> vs[2] >> vs[3] >> vs[4];

			for(int i = 0; i < 5; ++i) {
				if(vs[i].at(vs[i].length() - 1) == ',') {
					vs[i] = vs[i].substr(0, vs[i].length() - 1);
				}
			}

			// Get first two values and last value (lat, lon, path loss)
			float lat = atof(vs[0].c_str());
			float lon = atof(vs[1].c_str());
			float pl = -atof(vs[4].c_str());

			float y = (lat - ref_lat)  * EARTH_RADIUS * M_PI / 180.0;
			float x = (lon - ref_long) * EARTH_RADIUS * M_PI / 180.0 * cos(lat * M_PI / 180.0);
			Location data_loc(x, y, 0);
		
			cur_ano_data.push_back(std::make_pair(data_loc, pl));
		}
		line_num++;
	}
	pu_ano.close();
}

void TransmitterOnlySplatPM::preprocessPathLoss(PR* pr, float su_height, int pu_id, int pr_id) {
	// Changes directory to splat folder
	chdir(splat_dir.c_str());

	{
		std::stringstream sstr;
		sstr << splat_dir << "ano_files/pr_" << pu_id << "_" << pr_id << ".ano";
		pr->splat_ano_filename = sstr.str();
	}

	// Convert location to longitude and lattitude
	float loc_lat  = ref_lat  + (pr->loc.y / EARTH_RADIUS) * (180.0 / M_PI);
	float loc_long = ref_long + (pr->loc.x / EARTH_RADIUS) * (180.0 / M_PI) / cos(loc_lat * M_PI / 180.0);

	std::string pr_fname = "PR.qth";
	
	std::ofstream pr_qth(pr_fname);
	pr_qth << "PR" << std::endl << loc_lat << std::endl << loc_long << std::endl << pr->loc.z << "m" << std::endl;
	pr_qth.close();

	std::stringstream sstr;
	float timeout_seconds = 30;
	sstr << "timeout " << timeout_seconds << "s " << splat_cmd << " -t " << pr_fname << " -d " << sdf_dir << " -L " << su_height <<
			" -metric -R 20 -ano " << pr->splat_ano_filename << " > /dev/null 2>&1";

	int sys_code = 1;
	int num_tries = 2;
	for(int i = 0; sys_code != 0 && i < num_tries; ++i) {
		sys_code = system(sstr.str().c_str());
	}

	if(sys_code != 0) {
		std::cerr << "Unable to run Splat, sys-code: " << sys_code << std::endl;
		std::cerr << "Command used: " << sstr.str() << std::endl;
		exit(1);
	}
	chdir(return_dir.c_str());
	
	// loadANOFile(*pr);
}

void TransmitterOnlySplatPM::loadANOFile(const PR& pr) {
	if(pr.splat_ano_filename == "") {
		std::cerr << "Invalid filename" << std::endl;
		exit(1);
	}

	cur_loc = pr.loc;

	std::ifstream pr_ano(pr.splat_ano_filename);

	int line_num = 1;
	std::string line = "";
	cur_ano_data.clear();
	while(std::getline(pr_ano, line)) {
		if(line_num > 2) {
			// Remove last two characters if last character is '*'
			if(line.at(line.length() - 1) == '*') {
				line = line.substr(0, line.length() - 2);
			}

			// Split line by ', '
			std::string vs[5];
			std::stringstream sstr(line);
			sstr >> vs[0] >> vs[1] >> vs[2] >> vs[3] >> vs[4];

			for(int i = 0; i < 5; ++i) {
				if(vs[i].at(vs[i].length() - 1) == ',') {
					vs[i] = vs[i].substr(0, vs[i].length() - 1);
				}
			}

			// Get first two values and last value (lat, lon, path loss)
			float lat = atof(vs[0].c_str());
			float lon = atof(vs[1].c_str());
			float pl = -atof(vs[4].c_str());

			float y = (lat - ref_lat)  * EARTH_RADIUS * M_PI / 180.0;
			float x = (lon - ref_long) * EARTH_RADIUS * M_PI / 180.0 * cos(lat * M_PI / 180.0);
			Location data_loc(x, y, 0);
		
			cur_ano_data.push_back(std::make_pair(data_loc, pl));
		}
		line_num++;
	}
	pr_ano.close();
}

float TransmitterOnlySplatPM::getPathLoss(const Location& loc1, const Location& loc2) const {
	typedef struct{
		float ref_dist;
		float path_loss;
		float pu_dist;
	} PLDataType;

	struct comp {
		bool operator()(const PLDataType& v1, const PLDataType& v2) const {
			return v1.ref_dist < v2.ref_dist;
		}
	};

	if(loc1.dist(cur_loc) > 0.01) {
		std::cerr << "Mismatching locations" << std::endl;
		exit(1);
	}

	std::vector<PLDataType> selected_refs;
	for(unsigned int i = 0; i < cur_ano_data.size(); ++i) {
		PLDataType v;
		v.ref_dist = cur_ano_data[i].first.dist(loc2);
		v.path_loss = cur_ano_data[i].second;
		v.pu_dist = cur_ano_data[i].first.dist(loc1);

		if((int) selected_refs.size() < this->k || v.ref_dist < selected_refs.front().ref_dist) {
			if((int) selected_refs.size() >= this->k) {
				std::pop_heap(selected_refs.begin(), selected_refs.end(), comp());
				selected_refs.pop_back();
			}

			selected_refs.push_back(v);
			std::push_heap(selected_refs.begin(), selected_refs.end(), comp());
		}
	}

	float sum_weights = 0.0;
	std::vector<float> weights;
	for(unsigned int i = 0; i < selected_refs.size(); ++i) {
		weights.push_back(1.0 / pow(selected_refs[i].ref_dist, alpha));
		sum_weights += weights[i];
	}

	float path_loss = 0.0;
	for(unsigned int i = 0; i < selected_refs.size(); ++i) {
		float this_estimate = 0.0;
		if(gamma > 0.0) {
			// Log Distance
			this_estimate = selected_refs[i].path_loss - 10.0 * gamma * log10(loc1.dist(loc2) / selected_refs[i].pu_dist);
		} else {
			this_estimate = selected_refs[i].path_loss;
		}

		path_loss += weights[i] * this_estimate / sum_weights;
	}

	if(utils::unit_type == utils::UnitType::ABS) {
		return utils::fromdBm(path_loss);
	} else if(utils::unit_type == utils::UnitType::DB) {
		return path_loss;
	} else {
		std::cerr << "Unsupported unit_type" << std::endl;
		exit(1);
		return 0.0;
	}
}
