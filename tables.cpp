
#include "tables.h"

#include "primary_user.h"
#include "spectrum_sensor.h"

#include <iostream>
#include <map>
#include <vector>

void buildTables(const std::map<int, std::vector<const SSint*> >& ss_groups, const std::map<int, std::vector<const PUint*> >& pu_groups, const std::vector<PUint>& pus, int grid_size,
					GridTable* grid_table, PUTable* pu_table) {
	
	pu_table->num_pr_per_pu = -1;

	int ss_per_entry = -1;
	int pu_per_entry = -1;

	bool tmp_ss_set = false, tmp_pu_set = false;
	SSint tmp_ss;
	int tmp_pu_ind;

	for(auto ss_itr = ss_groups.begin(); ss_itr != ss_groups.end(); ++ss_itr) {
		auto pu_itr = pu_groups.find(ss_itr->first);

		if(pu_itr == pu_groups.end()) {
			std::cerr << "Invalid groups supplied to buildTables" << std::endl;
			exit(1);
		}

		if(ss_per_entry == -1) {
			ss_per_entry = int(ss_itr->second.size());
		}

		for(unsigned int i = 0; i < ss_itr->second.size(); ++i) {
			SSint ss_copy(*(ss_itr->second[i]));

			if(!tmp_ss_set) {
				tmp_ss = ss_copy;
				tmp_ss_set = true;
			}

			std::vector<int> new_rp;
			for(unsigned int j = 0; j < pu_itr->second.size(); ++j) {
				new_rp.push_back(ss_copy.received_power_from_pu[pu_itr->second[j]->index]);
			}

			ss_copy.received_power_from_pu = new_rp;

			grid_table->sss[ss_itr->first].push_back(ss_copy);
		}


		if(pu_per_entry == -1) {
			pu_per_entry = int(pu_itr->second.size());
		}

		for(unsigned int j = 0; j < pu_itr->second.size(); ++j) {
			if(pu_table->num_pr_per_pu == -1) {
				pu_table->num_pr_per_pu = pu_itr->second[j]->prs.size();
			} else if(pu_table->num_pr_per_pu != int(pu_itr->second[j]->prs.size())) {
				std::cerr << "Inconsistent number of PUs" << std::endl;
				exit(1);
			}

			if(!tmp_pu_set) {
				tmp_pu_ind = pu_itr->second[j]->index;
				tmp_pu_set = true;
			}

			grid_table->pu_refs[pu_itr->first].push_back(pu_itr->second[j]->index);

			if(pu_table->pus.count(pu_itr->second[j]->index) == 0) {
				pu_table->pus[pu_itr->second[j]->index] = *(pu_itr->second[j]);
			}
		}
	}


	for(unsigned int i = 0; i < pus.size(); ++i) {
		auto pu_itr = pu_table->pus.find(i);
		if(pu_itr == pu_table->pus.end()) {
			pu_table->pus[i] = pus[i];
		}
	}

	if(grid_size > 0) {
		for(int i = 0; i < grid_size; ++i) {
			auto ss_itr = grid_table->sss.find(i);
			if(ss_itr == grid_table->sss.end()) {
				grid_table->sss[i] = std::vector<SSint>(ss_per_entry, tmp_ss);
			}

			auto pu_itr = grid_table->pu_refs.find(i);
			if(pu_itr == grid_table->pu_refs.end()) {
				grid_table->pu_refs[i] = std::vector<int>(pu_per_entry, tmp_pu_ind);
			}
		}
	}

	// Check
/*	for(auto ss_itr = grid_table->sss.begin(); ss_itr != grid_table->sss.end(); ++ss_itr) {
		auto pu_itr = grid_table->pu_refs.find(ss_itr->first);
		if(pu_itr == grid_table->pu_refs.end()) {
			std::cerr << "AAAAAA" << std::endl;
			exit(1);
		}

		auto check_ss_itr = ss_groups.find(ss_itr->first);
		if(check_ss_itr == ss_groups.end()) {
			std::cerr << "BBBBBB" << std::endl;
			exit(1);
		}

		if(check_ss_itr->second.size() != ss_itr->second.size()) {
			std::cerr << "Mismatch SS size" << std::endl;
			exit(1);
		}

		for(unsigned int i = 0; i < ss_itr->second.size(); ++i) {
			if(ss_itr->second[i].loc.x != check_ss_itr->second[i]->loc.x ||
					ss_itr->second[i].loc.y != check_ss_itr->second[i]->loc.y ||
					ss_itr->second[i].received_power != check_ss_itr->second[i]->received_power) {
				std::cerr << "Basic info mismatch" << std::endl;
				exit(1);
			}

			for(unsigned int j = 0; j < ss_itr->second[i].received_power_from_pu.size(); ++j) {
				if(ss_itr->second[i].received_power_from_pu[j] != check_ss_itr->second[i]->received_power_from_pu[pu_itr->second[j]]) {
					std::cerr << "RP mismatch" << std::endl;
					exit(1);
				}
			}
		}

		auto check_pu_itr = pu_groups.find(pu_itr->first);
		if(check_pu_itr == pu_groups.end()) {
			std::cerr << "CCCCCCCCCCCCCCCCCC" << std::endl;
			exit(1);
		}

		if(check_pu_itr->second.size() != pu_itr->second.size()) {
			std::cerr << "Mismatch PU size" << std::endl;
			exit(1);
		}

		for(unsigned int i = 0; i < pu_itr->second.size(); ++i) {
			auto pu_int_itr = pu_table->pus.find(pu_itr->second[i]);
			if(pu_int_itr == pu_table->pus.end()) {
				std::cerr << "DDDDDDDDDDDDDDD" << std::endl;
				exit(1);
			}

			if(pu_int_itr->second.loc.x != check_pu_itr->second[i]->loc.x ||
					pu_int_itr->second.loc.y != check_pu_itr->second[i]->loc.y ||
					pu_int_itr->second.transmit_power != check_pu_itr->second[i]->transmit_power ||
					pu_int_itr->second.prs[0].threshold != check_pu_itr->second[i]->prs[0].threshold) {
				std::cerr << "EEEEEEEEEEEEEEEE" << std::endl;
				exit(1);
			}
		}
	}*/
}
