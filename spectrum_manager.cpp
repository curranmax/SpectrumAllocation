
#include "spectrum_manager.h"

#include "llong.h"
#include "location.h"
#include "path_loss_table.h"
#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_sensor.h"
#include "split.h"
#include "tables.h"
#include "utils.h"

#include "cryptoTools/Network/IOService.h"
#include "cryptoTools/Network/Session.h"
#include "cryptoTools/Crypto/PRNG.h"

#include "ivory/Runtime/ShGc/ShGcRuntime.h"
#include "ivory/Runtime/sInt.h"
#include "ivory/Runtime/Party.h"

#include <math.h>
#include <algorithm>

using namespace osuCrypto;

#define INPUT(parties, party_id, val, bit_count) parties[party_id].isLocalParty() ? parties[party_id].input<sInt>(val, bit_count) : parties[party_id].input<sInt>(bit_count)

#define LN_TEN 2.30258509299404568401
#define LOG_CALC_ITERS 10

void SpectrumManager::setGridParams(int _grid_min_num_pu, int _grid_min_num_ss, int _grid_num_x, int _grid_num_y, float _grid_delta_x, float _grid_delta_y) {
	use_grid = true;
	grid_min_num_pu = _grid_min_num_pu;
	grid_min_num_ss = _grid_min_num_ss;
	grid_num_x = _grid_num_x;
	grid_num_y = _grid_num_y;
	grid_delta_x = _grid_delta_x;
	grid_delta_y = _grid_delta_y;

	initGroupDeviations();
}

void SpectrumManager::setCommunicationValues(int _num_io_threads, const std::string& _server_addr, const std::string& _connection_name, const std::string& _channel_name) {
	num_io_threads = _num_io_threads;
	server_addr = _server_addr;
	connection_name = _connection_name;
	channel_name = _channel_name;
}

// Options: order(split rp at ss then estimate pl at su, estimate rp at su then split rp and calc pl), pl type(ratio, db)
std::vector<float> SpectrumManager::run(int party_id,
										const std::vector<PUint>& pus,
										std::vector<SSint>& sss,
										const std::vector<SUint>& sus,
										const std::map<int, std::vector<int> >& precomputed_pu_groups,
										const std::map<int, std::vector<int> >& precomputed_ss_groups,
										Timer* timer) const {
	std::string address = "127.0.0.1:9000";

	auto session_mode = SessionMode::Client;
	auto block_type = OneBlock;
	auto role = ShGcRuntime::Garbler;

	if(party_id == 0) {
		session_mode = SessionMode::Client;
		block_type = OneBlock;
		role = ShGcRuntime::Garbler;
	} else if(party_id == 1){
		session_mode = SessionMode::Server;
		block_type = ZeroBlock;
		role = ShGcRuntime::Evaluator;
	} else {
		std::cerr << "Invalid party_id: must be either 0 or 1" << std::endl;
	}

	PRNG prng(block_type);
	IOService ios;

	Session session(ios, address, session_mode);
	Channel channel = session.addChannel();

	if(party_id == 1) {
		channel.waitForConnection();
	}

	ShGcRuntime runtime;
	runtime.mDebugFlag = false;
	runtime.init(channel, prng.get<block>(), role, party_id);

	std::array<Party, 2> parties{
		Party(runtime, 0),
		Party(runtime, 1)
	};

	// Establsih channel between SMs
	IOService com_ios(num_io_threads);
	Endpoint ep(com_ios, server_addr, (parties[0].isLocalParty() ? EpMode::Server : EpMode::Client), connection_name);
	Channel sm_ch = ep.addChannel(channel_name);
	sm_ch.waitForConnection();

	// Run the preprocessing calculations:
	if(parties[0].isLocalParty() && !brief_out) {
		std::cout << "Starting preprocessing" << std::endl;
	}
	timer->start(Timer::secure_preprocessing);

	// If all sss.received_power_from_pu[i].size() == 0 then do secure preprocessing
	// If all sss.received_power_from_pu[i].size() == pus.size() then do nothing
	// Otherwise error

	bool first = true;
	unsigned int ss_rp_size = 0;
	for(unsigned int i = 0; i < sss.size(); ++i) {
		if(first) {
			first = false;
			ss_rp_size = sss[i].received_power_from_pu.size();
		} else {
			if(ss_rp_size != sss[i].received_power_from_pu.size()) {
				std::cerr << "Not all SSint have the same size 'received_power_from_pu'" << std::endl;
				exit(1);
			}
		}
	}

	if(ss_rp_size == 0) {
		secureRadarPreprocess(parties, pus, &sss);
	} else if(ss_rp_size != pus.size()) {
		std::cerr << "Invalid ss_rp_size: " << ss_rp_size << std::endl << "Wanted either 0 or pus.size (" << pus.size() << ")" << std::endl;
		exit(1);
	}

	// Calculate groups
	std::map<int, std::vector<const SSint*> > ss_groups;
	std::map<int, std::vector<const PUint*> > pu_groups;

	GridTable grid_table;
	PUTable pu_table;
	if(use_grid) {
		if(precomputed_ss_groups.size() > 0) {
			for(auto itr = precomputed_ss_groups.begin(); itr != precomputed_ss_groups.end(); ++itr) {
				for(unsigned int i = 0; i < itr->second.size(); ++i) {
					ss_groups[itr->first].push_back(&sss[itr->second[i]]);
				}
			}
		} else {
			ss_groups = secureSSGridPreprocess(parties, sss);
		}

		if(precomputed_pu_groups.size() > 0) {
			for(auto itr = precomputed_pu_groups.begin(); itr != precomputed_pu_groups.end(); ++itr) {
				for(unsigned int i = 0; i < itr->second.size(); ++i) {
					pu_groups[itr->first].push_back(&pus[itr->second[i]]);
				}
			}
		} else {
			pu_groups = securePUGridPreprocess(parties, pus);
		}

		buildTables(ss_groups, pu_groups, &grid_table, &pu_table);
	}

	timer->end(Timer::secure_preprocessing);
	if(parties[0].isLocalParty() && !brief_out) {
		std::cout << "Finished preprocessing" << std::endl;
	}

	std::vector<float> su_rps;
	std::vector<PUint> selected_pus;
	std::vector<SSint> selected_sss;
	if(!use_grid) {
		for(unsigned int j = 0; j < pus.size(); ++j) {
			selected_pus.push_back(pus[j]);
		}
		for(unsigned int i = 0; i < sss.size(); ++i) {
			selected_sss.push_back(sss[i]);
		}
	}
	for(unsigned int i = 0; i < sus.size(); ++i) {
		if(parties[0].isLocalParty()  && !brief_out) {
			std::cout << "Starting request for SU " << i + 1 << std::endl;
		}
		
		timer->start(Timer::secure_su_request);

		if(use_grid) {
			selected_pus.clear();
			selected_sss.clear();
			secureGetEntitiesFromTable(parties, sus[i], grid_table, pu_table, &selected_pus, &selected_sss);
		}
		
		float v = secureRadar(parties, sus[i], selected_pus, selected_sss, &pu_table, &sm_ch);

		timer->end(Timer::secure_su_request);

		su_rps.push_back(v);

		if(parties[0].isLocalParty()  && !brief_out) {
			std::cout << "Finished request for SU " << i + 1 << std::endl;
		}
	}
	return su_rps;
}

void SpectrumManager::secureRadarPreprocess(
		std::array<Party, 2> parties, const std::vector<PUint>& pus, std::vector<SSint>* sss) const {
	sInt zero = INPUT(parties, 0, 0, bit_count);
	sInt one  = INPUT(parties, 0, 1, bit_count);
	sInt two  = INPUT(parties, 0, 2, bit_count);
	sInt ten  = INPUT(parties, 0, 10, bit_count);

	sInt factor_int = INPUT(parties, 0, int(factor), bit_count);
	sInt large = INPUT(parties, 0, factor * factor * 10, bit_count); // TODO

	sInt ln_ten = INPUT(parties, 0, int(LN_TEN * factor), bit_count);

	// PU
	std::vector<sInt> pus_x;
	std::vector<sInt> pus_y;
	for(unsigned int i = 0; i < pus.size(); ++i) {
		sInt pu_x_a = INPUT(parties, 0, pus[i].loc.x, bit_count);
		sInt pu_x_b = INPUT(parties, 1, pus[i].loc.x, bit_count);

		pus_x.push_back(pu_x_a + pu_x_b);

		sInt pu_y_a = INPUT(parties, 0, pus[i].loc.y, bit_count);
		sInt pu_y_b = INPUT(parties, 1, pus[i].loc.y, bit_count);

		pus_y.push_back(pu_y_a + pu_y_b);
	}


	// SS
	std::vector<sInt> sss_x;
	std::vector<sInt> sss_y;
	std::vector<sInt> sss_rp;
	for(unsigned int i = 0; i < sss->size(); ++i) {
		sInt ss_x_a = INPUT(parties, 0, (*sss)[i].loc.x, bit_count);
		sInt ss_x_b = INPUT(parties, 1, (*sss)[i].loc.x, bit_count);

		sss_x.push_back(ss_x_a + ss_x_b);

		sInt ss_y_a = INPUT(parties, 0, (*sss)[i].loc.y, bit_count);
		sInt ss_y_b = INPUT(parties, 1, (*sss)[i].loc.y, bit_count);
		
		sss_y.push_back(ss_y_a + ss_y_b);

		sInt ss_rp_a = INPUT(parties, 0, (*sss)[i].received_power, bit_count);
		sInt ss_rp_b = INPUT(parties, 1, (*sss)[i].received_power, bit_count);
		
		sss_rp.push_back(ss_rp_a + ss_rp_b);
	}

	// |----------------------------------|
	// |  SS received power from each PU  |
	// |----------------------------------|
	std::vector<std::vector<sInt> > sss_rp_from_pu_a;
	std::vector<std::vector<sInt> > sss_rp_from_pu_b;
	for(unsigned int i = 0; i < sss->size(); ++i) {
		std::vector<sInt> rp_weights;
		sInt sum_rp_weights = INPUT(parties, 0, 0, bit_count);
		for(unsigned int j = 0; j < pus.size(); ++j) {
			sInt dist_to_rp_alpha = INPUT(parties, 0, 1, bit_count);

			if(rp_alpha == 1) {
				sInt dist;
				utils::dist(&dist, sss_x[i], sss_y[i], pus_x[j], pus_y[j], zero, two, 10);

				dist_to_rp_alpha = dist_to_rp_alpha * dist;
			} else if(rp_alpha == 2) {
				sInt x_diff = sss_x[i] - pus_x[j];
				sInt y_diff = sss_y[i] - pus_y[j];
				sInt dist_squared = (x_diff * x_diff +  y_diff * y_diff) / factor_int;

				dist_to_rp_alpha = dist_to_rp_alpha * dist_squared;
			} else if(rp_alpha == 3) {
				sInt dist;
				utils::dist(&dist, sss_x[i], sss_y[i], pus_x[j], pus_y[j], zero, two, 10);

				sInt dist_squared = ((sss_x[i] - pus_x[j]) * (sss_x[i] - pus_x[j]) + (sss_y[i] - pus_y[j]) * (sss_y[i] - pus_y[j])) / factor_int;

				dist_to_rp_alpha = dist_to_rp_alpha * dist * dist_squared / factor_int;
			} else if(rp_alpha == 4) {
				sInt dist_squared = ((sss_x[i] - pus_x[j]) * (sss_x[i] - pus_x[j]) + (sss_y[i] - pus_y[j]) * (sss_y[i] - pus_y[j])) / factor_int;

				dist_to_rp_alpha = dist_to_rp_alpha * dist_squared * dist_squared / factor_int;
			} else {
				std::cerr << "Invalid value for rp_alpha: " << rp_alpha << std::endl;
				exit(0);
			}

			sInt is_dist_small = dist_to_rp_alpha <= one;

			sInt this_weight = is_dist_small.ifelse(large, factor_int * factor_int / dist_to_rp_alpha);
			rp_weights.push_back(this_weight + zero);
			sum_rp_weights = sum_rp_weights + this_weight;
		}

		sss_rp_from_pu_a.push_back(std::vector<sInt>());
		sss_rp_from_pu_b.push_back(std::vector<sInt>());
		for(unsigned int j = 0; j < pus.size(); ++j) {
			sInt sss_rp_from_pu;
			if(utils::unit_type == utils::UnitType::ABS) {
				sss_rp_from_pu = rp_weights[j] * sss_rp[i] / sum_rp_weights;
			} else if(utils::unit_type == utils::UnitType::DB) {
				sInt log_rp_weight;
				utils::secureLog10(
					&log_rp_weight, factor_int * rp_weights[j] / sum_rp_weights,
					zero, factor_int, ln_ten, LOG_CALC_ITERS);

				sss_rp_from_pu = ten * log_rp_weight + sss_rp[i];
			} else {
				std::cerr << "Unsupported unit_type" << std::endl;
				exit(1);
			}
			sss_rp_from_pu_a[i].push_back(INPUT(parties, 0, splitRandVal(), bit_count));
			sss_rp_from_pu_b[i].push_back(sss_rp_from_pu - sss_rp_from_pu_a[i][j]);
		}
	}

	// Add computed values to plaintext_sss_rp_from_pu.
	for(unsigned int i = 0; i < sss->size(); ++i) {
		for(unsigned int j = 0; j < pus.size(); ++j) {
			parties[0].reveal(sss_rp_from_pu_a[i][j]);
			parties[1].reveal(sss_rp_from_pu_b[i][j]);
			if(parties[0].isLocalParty()) {
				(*sss)[i].received_power_from_pu.push_back(sss_rp_from_pu_a[i][j].getValue());
			} else if(parties[1].isLocalParty()) {
				(*sss)[i].received_power_from_pu.push_back(sss_rp_from_pu_b[i][j].getValue());
			}
		}
		exit(1);
	}
}

void SpectrumManager::plainTextRadarPreprocess(
		const std::vector<PUint>& pus_int0, const std::vector<PUint>& pus_int1,
		std::vector<SSint>* sss_int0, std::vector<SSint>* sss_int1) const {
	// Works only if bit_count is <= 64
	if(bit_count > 64) {
		std::cerr << "Unsupported value for bit_count" << std::endl;
		exit(1);
	}
	// Must Bitwise And the result of any operations.
	llong bit_mask = (1 << bit_count) - 1;

	llong factor_int = llong(factor);
	llong ln_ten_int = llong(LN_TEN * factor);

	int log_calc_iters = LOG_CALC_ITERS;

	for(unsigned int i = 0; i < sss_int0->size(); ++i) {
		std::vector<llong> rp_weights;
		llong sum_rp_weights = 0;
		bool all_zero = true;
		for(unsigned int j = 0; j < pus_int0.size(); ++j) {
			// Compute weight
			llong dist = 0;
			if(rp_alpha == 2) {
				llong x_diff = (((pus_int0[j].loc.x + pus_int1[j].loc.x) ^ bit_mask) - (((*sss_int0)[i].loc.x + (*sss_int1)[i].loc.x) ^ bit_mask)) ^ bit_mask;
				llong y_diff = (((pus_int0[j].loc.y + pus_int1[j].loc.y) ^ bit_mask) - (((*sss_int0)[i].loc.y + (*sss_int1)[i].loc.y) ^ bit_mask)) ^ bit_mask;

				dist = ((((x_diff * x_diff) ^ bit_mask) + ((y_diff * y_diff) ^ bit_mask)) / factor_int) ^ bit_mask;
			} else {
				std::cerr << "(Currently) Unsupported rp_alpha: " << rp_alpha << std::endl;
				exit(1);
			}

			llong this_weight = 0;
			if(dist <= 1) {
				this_weight = (factor_int * factor_int) ^ bit_mask;
			} else {
				this_weight = (((factor_int * factor_int) ^ bit_mask) / dist) ^ bit_mask;
			}

			all_zero = all_zero && (this_weight == 0);

			rp_weights.push_back(this_weight);
			sum_rp_weights += this_weight;
		}

		if(all_zero) {
			sum_rp_weights = 0;
			for(unsigned int j = 0; j < rp_weights.size(); ++j) {
				rp_weights[j] = 1;
				sum_rp_weights += rp_weights[j];
			}
		}

		for(unsigned int j = 0; j < pus_int0.size(); ++j) {
			llong this_rp = 0;
			if(utils::unit_type == utils::UnitType::ABS) {
				this_rp = (((rp_weights[j] * (((*sss_int0)[i].received_power + (*sss_int1)[i].received_power) ^ bit_mask)) ^ bit_mask) / sum_rp_weights) ^ bit_mask;
			} else if(utils::unit_type == utils::UnitType::DB) {

				// Numerical calculates the log_10 of rp_weights[j] / sum_rp_weights.
				llong log_ratio = 0;

				llong ratio = (((factor_int * rp_weights[j]) ^ bit_mask) / sum_rp_weights) ^ bit_mask;
				llong numerator = (ratio - factor_int) ^ bit_mask;
				llong denomerator = factor_int;
				llong sign = factor_int;
				for(int x = 0; x < log_calc_iters; ++x) {
					log_ratio = (log_ratio + ((((sign * numerator) ^ bit_mask) / denomerator) ^ bit_mask)) ^ bit_mask;
				
					if(x < log_calc_iters - 1) {
						numerator = (((numerator * ((ratio - factor_int) ^ bit_mask)) ^ bit_mask) / factor_int) ^ bit_mask;
						denomerator = (denomerator + factor_int) ^ bit_mask;
						sign = (0 - sign) ^ bit_mask;
					}
				}

				log_ratio = (((factor_int * log_ratio) ^ bit_mask) / ln_ten_int) ^ bit_mask;

				this_rp = (((10 * log_ratio) ^ bit_mask) + (((*sss_int0)[i].received_power + (*sss_int1)[i].received_power) ^ bit_mask)) ^ bit_mask;
			} else {
				std::cerr << "Unsupported unit_type" << std::endl;
				exit(1);
			}

			auto this_rp_split = splitInt(int(this_rp));
			(*sss_int0)[i].received_power_from_pu.push_back(this_rp_split.first);
			(*sss_int1)[i].received_power_from_pu.push_back(this_rp_split.second);
		}
	}
}

std::map<int, std::vector<const SSint*> > SpectrumManager::secureSSGridPreprocess(
		std::array<Party, 2> parties, const std::vector<SSint>& sss) const {
	// Inits the valid groups.
	std::map<int, std::vector<const SSint*> > orig_ss_groups;
	std::map<int, std::vector<const SSint*> > ss_groups;
	for(int group_id = 0; group_id < grid_num_x * grid_num_y; ++group_id) {
		ss_groups[group_id];
		orig_ss_groups[group_id];
	}

	// Constants
	sInt zero = INPUT(parties, 0, 0, bit_count);
	sInt neg_one = INPUT(parties, 0, -1, bit_count);

	sInt secure_grid_num_x = INPUT(parties, 0, grid_num_x, bit_count);

	sInt secure_grid_delta_x = INPUT(parties, 0, int(grid_delta_x * factor), bit_count);
	sInt secure_grid_delta_y = INPUT(parties, 0, int(grid_delta_y * factor), bit_count);

	// SS
	std::vector<sInt> ss_xs;
	std::vector<sInt> ss_ys;
	for(unsigned int i = 0; i < sss.size(); ++i) {
		sInt ss_x_a = INPUT(parties, 0, sss[i].loc.x, bit_count);
		sInt ss_x_b = INPUT(parties, 1, sss[i].loc.x, bit_count);

		ss_xs.push_back(ss_x_a + ss_x_b);

		sInt ss_y_a = INPUT(parties, 0, sss[i].loc.y, bit_count);
		sInt ss_y_b = INPUT(parties, 1, sss[i].loc.y, bit_count);
		
		ss_ys.push_back(ss_y_a + ss_y_b);
	}

	// Compute group_id for each SS, and put it in the group.
	for(unsigned int i = 0; i < sss.size(); ++i) {
		sInt group_id_a = ss_xs[i] / secure_grid_delta_x + ss_ys[i] / secure_grid_delta_y * secure_grid_num_x;
		sInt group_id_b = group_id_a + zero;

		parties[0].reveal(group_id_a);
		parties[1].reveal(group_id_b);

		int group_id;
		if(parties[0].isLocalParty()) {
			group_id = group_id_a.getValue();
		}

		if(parties[1].isLocalParty()) {
			group_id = group_id_b.getValue();
		}

		auto itr = orig_ss_groups.find(group_id);
		if(itr == orig_ss_groups.end()) {
			std::cerr << "Invalid group_id: " << group_id << std::endl;
			exit(1);
		}

		itr->second.push_back(&sss[i]);
	}

	// Expand each group to the grid_min_num_ss.
	const int max_iters = all_deviations.size();
	const unsigned int min_num = numSelect(grid_min_num_ss, sss.size());
	for(int start_i = 0; start_i < grid_num_x; ++start_i) {
		for(int start_j = 0; start_j < grid_num_y; ++start_j) {
			int plain_group_id = start_i + start_j * grid_num_x;
			auto grid_loc_itr = ss_groups.find(plain_group_id);

			int iter = 0;
			while(grid_loc_itr->second.size() < min_num && iter < max_iters) {
				const std::vector<std::pair<int, int> > devs = getDeviations(iter);
				for(unsigned int x = 0; x < devs.size(); ++x) {
					int this_i = start_i + devs[x].first;
					int this_j = start_j + devs[x].second;

					int this_group_id = -1;
					if(this_i >= 0 && this_i < grid_num_x && this_j >= 0 && this_j < grid_num_y) {
						this_group_id = this_i + this_j * grid_num_x;
					}

					auto itr = orig_ss_groups.find(this_group_id);
					if(itr == orig_ss_groups.end()) {
						if(this_group_id != -1) {
							std::cerr << "Invalid this_group_id: " << this_group_id << std::endl;
							exit(1);
						}
					} else {
						for(unsigned int i = 0; i < itr->second.size(); ++i) {
							grid_loc_itr->second.push_back(itr->second[i]);
						}
					}
				}
				++iter;
			}
		}
	}
	return ss_groups;
}

std::map<int, std::vector<const PUint*> > SpectrumManager::securePUGridPreprocess(
		std::array<Party, 2> parties, const std::vector<PUint>& pus) const {
	std::cerr << "Not implemented" << std::endl;
	exit(1);

	return std::map<int, std::vector<const PUint*> >();
}

void SpectrumManager::secureGetEntitiesFromTable(std::array<Party, 2> parties, const SUint& su,
		const GridTable& grid_table, const PUTable& pu_table,
		std::vector<PUint>* selected_pus, std::vector<SSint>* selected_sss) const {
	sInt zero =    INPUT(parties, 0, 0,  bit_count);
	sInt neg_one = INPUT(parties, 0, -1, bit_count);

	sInt secure_grid_num_x = INPUT(parties, 0, grid_num_x, bit_count);

	sInt secure_grid_delta_x = INPUT(parties, 0, int(grid_delta_x * factor), bit_count);
	sInt secure_grid_delta_y = INPUT(parties, 0, int(grid_delta_y * factor), bit_count);

	// SU
	sInt su_x_a = INPUT(parties, 0, su.loc.x, bit_count);
	sInt su_x_b = INPUT(parties, 1, su.loc.x, bit_count);

	sInt su_i = (su_x_a + su_x_b) / secure_grid_delta_x;

	sInt su_y_a = INPUT(parties, 0, su.loc.y, bit_count);
	sInt su_y_b = INPUT(parties, 1, su.loc.y, bit_count);

	sInt su_j = (su_y_a + su_y_b) / secure_grid_delta_y;

	sInt real_group_id = su_i + su_j * secure_grid_num_x;
	sInt group_id_a = real_group_id + zero;
	sInt group_id_b = real_group_id + zero;

	parties[0].reveal(group_id_a);
	parties[1].reveal(group_id_b);

	int group_id;
	if(parties[0].isLocalParty()) {
		group_id = group_id_a.getValue();
	}
	if(parties[1].isLocalParty()) {
		group_id = group_id_b.getValue();
	}

	auto ss_itr = grid_table.sss.find(group_id);
	if(ss_itr == grid_table.sss.end()) {
		std::cerr << "Invalid group_id: " << group_id << std::endl;
		exit(1);
	}
	for(unsigned int i = 0; i < ss_itr->second.size(); ++i) {
		selected_sss->push_back(ss_itr->second[i]);
	}

	auto pu_ref_itr = grid_table.pu_refs.find(group_id);
	if(pu_ref_itr == grid_table.pu_refs.end()) {
		std::cerr << "Invalid group_id: " << group_id << std::endl;
		exit(1);
	}
	for(unsigned int i = 0; i < pu_ref_itr->second.size(); ++i) {
		auto pu_itr = pu_table.pus.find(pu_ref_itr->second[i]);
		if(pu_itr == pu_table.pus.end()) {
			std::cerr << "Invalid pu_id: " << pu_ref_itr->second[i] << std::endl;
			exit(1);
		}
		selected_pus->push_back(pu_itr->second);
	}
}

float SpectrumManager::secureRadar(
		std::array<Party, 2> parties, const SUint& su,
		const std::vector<PUint>& pus, const std::vector<SSint>& sss,
		PUTable* pu_table, Channel* sm_ch) const {
	// Constants
	sInt zero = INPUT(parties, 0,  0, bit_count);
	sInt one  = INPUT(parties, 0,  1, bit_count);
	sInt two  = INPUT(parties, 0,  2, bit_count);
	sInt ten  = INPUT(parties, 0, 10, bit_count);

	sInt factor_int = INPUT(parties, 0, int(factor), bit_count);
	sInt large = INPUT(parties, 0, factor * factor * 10, bit_count); // TODO
	sInt ln_ten = INPUT(parties, 0, int(LN_TEN * factor), bit_count);

	// |------------------|
	// |  Combine values  |
	// |------------------|
	// SU
	sInt su_x_a = INPUT(parties, 0, su.loc.x, bit_count);
	sInt su_x_b = INPUT(parties, 1, su.loc.x, bit_count);

	sInt su_x = su_x_a + su_x_b;

	sInt su_y_a = INPUT(parties, 0, su.loc.y, bit_count);
	sInt su_y_b = INPUT(parties, 1, su.loc.y, bit_count);

	sInt su_y = su_y_a + su_y_b;

	// PU
	unsigned int k_pu = numSelect(this->num_pu_selection, pus.size());
	std::vector<sInt> pus_x;
	std::vector<sInt> pus_y;
	std::vector<sInt> pus_tp;
	std::vector<sInt> prs_thresh;
	std::vector<sInt> pu_table_index;
	for(unsigned int i = 0; i < pus.size(); ++i) {
		if(k_pu < pus.size()) {	
			sInt pu_x_a = INPUT(parties, 0, pus[i].loc.x, bit_count);
			sInt pu_x_b = INPUT(parties, 1, pus[i].loc.x, bit_count);

			pus_x.push_back(pu_x_a + pu_x_b);

			sInt pu_y_a = INPUT(parties, 0, pus[i].loc.y, bit_count);
			sInt pu_y_b = INPUT(parties, 1, pus[i].loc.y, bit_count);

			pus_y.push_back(pu_y_a + pu_y_b);
		}

		sInt pu_tp_a = INPUT(parties, 0, pus[i].transmit_power, bit_count);
		sInt pu_tp_b = INPUT(parties, 1, pus[i].transmit_power, bit_count);
	
		pus_tp.push_back(pu_tp_a + pu_tp_b);
	
		sInt pr_thresh_a = INPUT(parties, 0, pus[i].prs[0].threshold, bit_count);
		sInt pr_thresh_b = INPUT(parties, 1, pus[i].prs[0].threshold, bit_count);
	
		prs_thresh.push_back(pr_thresh_a + pr_thresh_b);

		// TODO Get the halves
		sInt pu_index_a = INPUT(parties, 0, pus[i].index, bit_count);

		pu_table_index.push_back(pu_index_a + zero);
	}

	// SS
	unsigned int k_ss = numSelect(this->num_ss_selection, sss.size());

	// Get the indices of SS to load
	std::vector<int> sss_load_inds;
	if(selection_algo == SpectrumManager::SelectionAlgo::NONE ||
			selection_algo == SpectrumManager::SelectionAlgo::SORT) {
		for(unsigned int i = 0; i < sss.size(); ++i) {
			sss_load_inds.push_back(i);
		}
		// Nothing has to be done now.
	} else if(selection_algo == SpectrumManager::SelectionAlgo::RANDOM) {
		std::vector<int> all_inds;
		if(parties[0].isLocalParty()) {
			for(unsigned int i = 0; i < sss.size(); ++i) {
				all_inds.push_back(i);
			}
			std::random_shuffle(all_inds.begin(), all_inds.end());
		}

		for(unsigned int i = 0; i < k_ss; ++i) {
			sInt rand_ind_a = INPUT(parties, 0, all_inds[i], bit_count);
			sInt rand_ind_b = rand_ind_a + zero;

			parties[0].reveal(rand_ind_a);
			parties[1].reveal(rand_ind_b);

			if(parties[0].isLocalParty()) {
				sss_load_inds.push_back(rand_ind_a.getValue());
			}
			if(parties[1].isLocalParty()) {
				sss_load_inds.push_back(rand_ind_b.getValue());
			}
		}
	} else {
		std::cerr << "Unsupported selection_algo" << std::endl;
		exit(1);
	}


	std::vector<sInt> sss_x;
	std::vector<sInt> sss_y;
	for(unsigned int x = 0; x < sss_load_inds.size(); ++x) {
		int i = sss_load_inds[x];

		sInt ss_x_a = INPUT(parties, 0, sss[i].loc.x, bit_count);
		sInt ss_x_b = INPUT(parties, 1, sss[i].loc.x, bit_count);

		sss_x.push_back(ss_x_a + ss_x_b);

		sInt ss_y_a = INPUT(parties, 0, sss[i].loc.y, bit_count);
		sInt ss_y_b = INPUT(parties, 1, sss[i].loc.y, bit_count);
		
		sss_y.push_back(ss_y_a + ss_y_b);
	}

	// Get the preprocessing values.
	std::vector<std::vector<sInt> > sss_rp_from_pu;
	for(unsigned int x = 0; x < sss_load_inds.size(); ++x) {
		int i = sss_load_inds[x];

		sss_rp_from_pu.push_back(std::vector<sInt>());
		for(unsigned int j = 0; j < pus.size(); ++j) {
			if(j < 0) {
				std::cerr << "Invalid index: " << j << std::endl;
				exit(1);
			}

			sInt sss_rp_from_pu_a = INPUT(parties, 0, sss[i].received_power_from_pu[j], bit_count);
			sInt sss_rp_from_pu_b = INPUT(parties, 1, sss[i].received_power_from_pu[j], bit_count);

			sss_rp_from_pu[x].push_back(sss_rp_from_pu_a + sss_rp_from_pu_b);
		}
	}


	// |----------------------------|
	// |  Select the k-closest PUs  |
	// |----------------------------|
	std::vector<int> pu_inds;
	if(k_pu < pus.size()) {
		std::vector<sInt> secret_pu_inds;
		std::vector<sInt> secret_pu_dists;
		for(unsigned int i = 0; i < pus.size(); ++i) {
			sInt x_diff = su_x - pus_x[i];
			sInt y_diff = su_y - pus_y[i];
			sInt dist_pu_su = (x_diff * x_diff + y_diff * y_diff) / factor_int;
			sInt ind = INPUT(parties, 0, i, bit_count);

			if(i < k_pu) {
				secret_pu_inds.push_back(ind + zero);
				secret_pu_dists.push_back(dist_pu_su + zero);
			} else {
				for(unsigned int x = 0; x < secret_pu_inds.size(); ++x) {
					sInt comp = dist_pu_su < secret_pu_dists[x];

					sInt tmp = comp.ifelse(dist_pu_su, secret_pu_dists[x]);
					dist_pu_su = comp.ifelse(secret_pu_dists[x], dist_pu_su);
					secret_pu_dists[x] = tmp + zero;

					tmp = comp.ifelse(ind, secret_pu_inds[x]);
					ind = comp.ifelse(secret_pu_inds[x], ind);
					secret_pu_inds[x] = tmp + zero;
				}
			}
		}

		for(unsigned int i = 0; i < secret_pu_inds.size(); ++i) {
			parties[0].reveal(secret_pu_inds[i]);
			sInt tmp = secret_pu_inds[i] + zero;
			parties[1].reveal(tmp);

			if(parties[0].isLocalParty()) {
				pu_inds.push_back(secret_pu_inds[i].getValue());
			} else if(parties[1].isLocalParty()) {
				pu_inds.push_back(tmp.getValue());
			}
		}
	} else {
		for(unsigned int i = 0; i < pus.size(); ++i) {
			pu_inds.push_back(i);
		}
	}

	// |----------------------------|
	// |  Select the k-closest SSs  |
	// |----------------------------|
	std::vector<int> sss_inds;
	std::vector<sInt> secret_sss_dists;
	if(selection_algo == SpectrumManager::SelectionAlgo::NONE ||
			selection_algo == SpectrumManager::SelectionAlgo::RANDOM) {
		for(unsigned int i = 0; i < k_ss; ++i) {
			sss_inds.push_back(i);

			sInt dist_to_pl_alpha = INPUT(parties, 0, 1, bit_count);

			if(pl_alpha == 1) {
				sInt dist;
				utils::dist(&dist, su_x, su_y, sss_x[i], sss_y[i], zero, two, 10);

				dist_to_pl_alpha = dist_to_pl_alpha * dist;
			} else if(pl_alpha == 2) {
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;

				dist_to_pl_alpha = dist_to_pl_alpha * dist_squared;
			} else if(pl_alpha == 3) {
				sInt dist;
				utils::dist(&dist, su_x, su_y, sss_x[i], sss_y[i], zero, two, 10);
				
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;

				dist_to_pl_alpha = dist_to_pl_alpha * dist * dist_squared / factor_int;
			} else if(pl_alpha == 4) {
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;

				dist_to_pl_alpha = dist_to_pl_alpha * dist_squared * dist_squared / factor_int;
			} else {
				std::cerr << "Invalid value of pl_alpha: " << pl_alpha << std::endl;
				exit(0);
			}

			secret_sss_dists.push_back(dist_to_pl_alpha + zero);
		}
	} else if(selection_algo == SpectrumManager::SelectionAlgo::SORT) {
		std::vector<sInt> secret_sss_inds;
		for(unsigned int i = 0; i < sss.size(); ++i) {
			sInt dist_to_pl_alpha = INPUT(parties, 0, 1, bit_count);

			if(pl_alpha == 1) {
				sInt dist;
				utils::dist(&dist, su_x, su_y, sss_x[i], sss_y[i], zero, two, 10);

				dist_to_pl_alpha = dist_to_pl_alpha * dist;
			} else if(pl_alpha == 2) {
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;

				dist_to_pl_alpha = dist_to_pl_alpha * dist_squared;
			} else if(pl_alpha == 3) {
				sInt dist;
				utils::dist(&dist, su_x, su_y, sss_x[i], sss_y[i], zero, two, 10);
				
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;

				dist_to_pl_alpha = dist_to_pl_alpha * dist * dist_squared / factor_int;
			} else if(pl_alpha == 4) {
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;

				dist_to_pl_alpha = dist_to_pl_alpha * dist_squared * dist_squared / factor_int;
			} else {
				std::cerr << "Invalid value of pl_alpha: " << pl_alpha << std::endl;
				exit(0);
			}

			sInt ind = INPUT(parties, 0, i, bit_count);

			if(i < k_ss) {
				secret_sss_inds.push_back(ind + zero);
				secret_sss_dists.push_back(dist_to_pl_alpha + zero);
			} else {
				for(unsigned int x = 0; x < secret_sss_inds.size(); ++x) {
					sInt comp = dist_to_pl_alpha < secret_sss_dists[x];

					sInt tmp = comp.ifelse(dist_to_pl_alpha, secret_sss_dists[x]);
					dist_to_pl_alpha = comp.ifelse(secret_sss_dists[x], dist_to_pl_alpha);
					secret_sss_dists[x] = tmp + zero;

					tmp = comp.ifelse(ind, secret_sss_inds[x]);
					ind = comp.ifelse(secret_sss_inds[x], ind);
					secret_sss_inds[x] = tmp + zero;
				}
			}
		}

		// Have to reveal the inds to both parties.
		for(unsigned int i = 0; i < secret_sss_inds.size(); ++i) {
			parties[0].reveal(secret_sss_inds[i]);
			sInt tmp = secret_sss_inds[i] + zero;
			parties[1].reveal(tmp);

			if(parties[0].isLocalParty()) {
				sss_inds.push_back(secret_sss_inds[i].getValue());
			}
			if(parties[1].isLocalParty()) {
				sss_inds.push_back(tmp.getValue());
			}
		}
	} else {
		std::cerr << "Unsupported selection_algo" << std::endl;
		exit(1);
	}

	// |----------------------|
	// |  Compute PL weights  |
	// |----------------------|
	// w = (1 / d(SU, SS))
	std::vector<sInt> sss_w;
	sInt sum_weight = INPUT(parties, 0, 0, bit_count);
	for(unsigned int x = 0; x < sss_inds.size(); ++x) {
		// Check if dist is small
		sInt is_dist_small = secret_sss_dists[x] <= one;

		sss_w.push_back(is_dist_small.ifelse(large, factor_int * factor_int / secret_sss_dists[x]));
		sum_weight = sum_weight + sss_w[x];
	}

	// |--------------------------------|
	// |  Compute transmit power of SU  |
	// |--------------------------------|
	// su_rp = pr_thresh * sum(w(SS)) / (sum(r(SS)/t(PU) * w(SS)))
	sInt su_tp = INPUT(parties, 0, 0, bit_count);
	std::vector<std::vector<sInt> > path_loss_su_pr; // path_loss_su_pr[i][j] is the path loss between the su and PU i's PR j.
	for(unsigned int y = 0; y < pu_inds.size(); ++y) {
		unsigned int j = pu_inds[y];
		sInt sum_weighted_ratio = INPUT(parties, 0, 0, bit_count);

		path_loss_su_pr.push_back(std::vector<sInt>());
		for(unsigned int x = 0; x < sss_inds.size(); ++x) {
			unsigned int i = sss_inds[x];
			if(utils::unit_type == utils::UnitType::ABS) {
				sum_weighted_ratio = sum_weighted_ratio + sss_w[x] * sss_rp_from_pu[i][j] / pus_tp[j];
			} else if(utils::unit_type == utils::UnitType::DB) {
				sum_weighted_ratio = sum_weighted_ratio + sss_w[x] * (sss_rp_from_pu[i][j] - pus_tp[j]) / factor_int;
			} else {
				std::cerr << "Unsupported unit_type" << std::endl;
				exit(1);
			}
		}
	
		sInt this_su_tp;
		if(utils::unit_type == utils::UnitType::ABS) {
			sInt swr_zero = sum_weighted_ratio < one;
			this_su_tp = swr_zero.ifelse(zero, prs_thresh[j] * sum_weight / sum_weighted_ratio);
			path_loss_su_pr[j].push_back(swr_zero.ifelse(zero, factor_int * sum_weight / sum_weighted_ratio));
		} else if(utils::unit_type == utils::UnitType::DB) {
			this_su_tp = prs_thresh[j] - factor_int * sum_weighted_ratio / sum_weight;
			path_loss_su_pr[j].push_back(factor_int * sum_weighted_ratio / sum_weight);
		} else {
			std::cerr << "Unsupported unit_type" << std::endl;
			exit(1);
		}


		if(y == 0) {
			su_tp = this_su_tp + zero;
		} else {
			sInt min_su_tp = this_su_tp < su_tp;
			su_tp = min_su_tp.ifelse(this_su_tp, su_tp);
		}
	}

	// |-------------------------------|
	// |  Reveal transmit power of SU  |
	// |-------------------------------|
	sInt su_tp_a = INPUT(parties, 0, splitRandVal(), bit_count);
	sInt su_tp_b = su_tp - su_tp_a;

	// reveal this output to party 0.
	parties[0].reveal(su_tp_a);
	parties[1].reveal(su_tp_b);

	float max_su_tp;
	if (parties[0].isLocalParty()) {
		max_su_tp = float(su_tp_a.getValue()) / factor;
	} else if(parties[1].isLocalParty()) {
		 max_su_tp = float(su_tp_b.getValue()) / factor;;
	}

	// |------------------------|
	// |  Update PR thresholds  |
	// |------------------------|
	// TODO Send these values to the SU, and get the actual transmit power
	int actual_su_tp_pt = 0;

	sInt actual_su_tp_a = INPUT(parties, 0, int(actual_su_tp_pt * factor), bit_count);
	sInt actual_su_tp_b = INPUT(parties, 1, int(actual_su_tp_pt * factor), bit_count);
	sInt actual_su_tp = actual_su_tp_a + actual_su_tp_b;

	sInt tmp_diff = INPUT(parties, 0, int(-5.0 * factor), bit_count);

	// Calculate updates
	// updates[i].first is the index of the PU to update, and updates[i].second is the list of updates to the PR threshholds to update.
	std::vector<std::pair<sInt, std::vector<sInt> > > updates;
	for(unsigned int y = 0; y < pu_inds.size(); ++y) {
		unsigned int j = pu_inds[y];
		updates.push_back(std::make_pair(pu_table_index[j] + zero, std::vector<sInt>()));

		int num_prs = 1;
		for(int x = 0; x < num_prs; ++x) {
			// Calculate the update to the PR thresh: u = 10 * log_10(1.0 - 10 ^ ((rp - thresh) / 10))
			sInt diff_dbm = tmp_diff + zero;
			sInt diff_div_ten = diff_dbm / ten;

			sInt ten_to_diff;
			utils::securePow10(&ten_to_diff, diff_div_ten, parties, bit_count, factor, factor_int);

			ten_to_diff = factor_int - ten_to_diff;

			sInt log_diff;
			utils::secureLog10(&log_diff, ten_to_diff, zero, factor_int, ln_ten, LOG_CALC_ITERS);

			sInt update = ten * log_diff;

			// TODO If diff is beyond some constant, make update zero

			updates[y].second.push_back(update + zero);
		}
	}
	secureTableWrite(parties, pu_table, updates, sm_ch);
	
	return max_su_tp;
}

// Helper function for secureTableWrite
void shiftTable(std::vector<std::vector<int> >& table, int shift) {
	std::vector<std::vector<int> > copy(table);
	int n = copy.size();
	for(unsigned int i = 0; i < copy.size(); ++i) {
		int shifted_index = (i + shift);
		if(shifted_index < 0) {
			shifted_index += n;
		}
		shifted_index = shifted_index % n;
		for(unsigned int j = 0; j < copy[i].size(); ++j) {
			table[shifted_index][j] = copy[i][j];
		}
	}
}

void flattenTable(const std::vector<std::vector<int> >& table, std::vector<int>& flattened_table) {
	flattened_table.clear();
	for(unsigned int i = 0; i < table.size(); ++i) {
		for(unsigned int j = 0; j < table[i].size(); ++j) {
			flattened_table.push_back(table[i][j]);
		}
	}
}

void expandTable(const std::vector<int>& table, int x, int y, std::vector<std::vector<int> >& expanded_table) {
	expanded_table = std::vector<std::vector<int> >(x, std::vector<int>(y));
	int index = 0;
	for(int i = 0; i < x; ++i) {
		for(int j = 0; j < y; ++j) {
			expanded_table[i][j] = table[index];
			index++;
		}
	}
}

void SpectrumManager::secureTableWrite(std::array<Party, 2> parties, PUTable* pu_table, const std::vector<std::pair<sInt, std::vector<sInt> > >& updates, Channel* sm_ch) const {
	if(pu_table == nullptr) {
		std::cerr << "pu_table is NULL" << std::endl;
		exit(1);
	}

	if(sm_ch == nullptr) {
		std::cerr << "sm_ch is NULL" << std::endl;
		exit(1);
	}

	if(parties[0].isLocalParty() && secure_write_timer != nullptr) {
		secure_write_timer->start(Timer::secure_write);
	}

	int num_pu = pu_table->pus.size();
	int num_pr_per_pu = pu_table->num_pr_per_pu;

	if(secure_write_algo == SpectrumManager::SecureWriteAlgo::PROPOSED){
		// TODO secure Table write.
		// Each side randomly shifts the "threshold table"
		int shift = rand() % pu_table->pus.size();
	
		sInt shift_a = INPUT(parties, 0, shift, bit_count);
		sInt shift_b = INPUT(parties, 1, shift, bit_count);
	
		// Using S2-PC, get the coordinate to get the apprpriate values
		std::vector<std::pair<int, std::vector<int> > > this_updates;
		for(unsigned int i = 0; i < updates.size(); ++i) {
			sInt this_shifted_index_a = shift_a + updates[i].first;
			sInt this_shifted_index_b = shift_b + updates[i].first;
	
			parties[0].reveal(this_shifted_index_b);
			parties[1].reveal(this_shifted_index_a);
	
			int this_shifted_index = 0;
			if(parties[0].isLocalParty()) {
				this_shifted_index = this_shifted_index_b.getValue();
			} else if(parties[1].isLocalParty()) {
				this_shifted_index = this_shifted_index_a.getValue();
			}
			this_shifted_index = this_shifted_index % pu_table->pus.size();
			this_updates.push_back(std::make_pair(this_shifted_index, std::vector<int>()));
	
			for(unsigned int j = 0; j < updates[i].second.size(); ++j) {
	
				sInt this_update_a = INPUT(parties, 0, rand() % 100 - 50, bit_count);
				sInt this_update_b = updates[i].second[j] - this_update_a;
	
				parties[0].reveal(this_update_a);
				parties[1].reveal(this_update_b);
	
				int this_update = 0;
				if(parties[0].isLocalParty()) {
					this_update = this_update_a.getValue();
				} else if(parties[1].isLocalParty()) {
					this_update = this_update_b.getValue();
				}
				this_updates[i].second.push_back(this_update);
			}
		}
	
		// Generate the two initial tables
		std::vector<std::vector<int> > update_table_1(num_pu, std::vector<int>(num_pr_per_pu));
		std::vector<std::vector<int> > update_table_2(num_pu, std::vector<int>(num_pr_per_pu));
	
		for(unsigned int i = 0; i < update_table_1.size(); ++i) {
			for(unsigned int j = 0; j < update_table_1[i].size(); ++j) {
				update_table_1[i][j] = rand() % 100 - 50;
				update_table_2[i][j] = 0 - update_table_1[i][j];
			}
		}
	
		for(unsigned i = 0; i < this_updates.size(); ++i) {
			for(unsigned int j = 0; j < this_updates[i].second.size(); ++j) {
				update_table_2[this_updates[i].first][j] = this_updates[i].second[j] - update_table_1[this_updates[i].first][j];
			}
		}
	
		// Shift one table by shift
		shiftTable(update_table_2, shift);
	
		// Swap those tables
		std::vector<int> flattened_table_2;
		flattenTable(update_table_2, flattened_table_2);
	
		std::vector<int> flattened_table_3;
		if(parties[0].isLocalParty()) {
			sm_ch->send(flattened_table_2);
			sm_ch->recv(flattened_table_3);
		} else if(parties[1].isLocalParty()) {
			sm_ch->recv(flattened_table_3);
			sm_ch->send(flattened_table_2);
		}
		std::vector<std::vector<int> > update_table_3;
		expandTable(flattened_table_3, num_pu, num_pr_per_pu, update_table_3);
	
		// Shift that table by -shift
		shiftTable(update_table_3, -shift);
	
		// Combine the tables
		std::vector<std::vector<int> > update_table_4(num_pu, std::vector<int>(num_pr_per_pu));
		for(unsigned int i = 0; i < update_table_1.size(); ++i) {
			for(unsigned int j = 0; j < update_table_1[i].size(); ++j) {
				update_table_4[i][j] = update_table_1[i][j] + update_table_3[i][j];
			}
		}
	
		// Swap the combined tables.
		std::vector<int> flattened_table_4;
		flattenTable(update_table_4, flattened_table_4);
	
		std::vector<int> flattened_table_5;
		if(parties[0].isLocalParty()) {
			sm_ch->send(flattened_table_4);
			sm_ch->recv(flattened_table_5);
		} else if(parties[1].isLocalParty()) {
			sm_ch->recv(flattened_table_5);
			sm_ch->send(flattened_table_4);
		}
		std::vector<std::vector<int> > update_table_5;
		expandTable(flattened_table_5, num_pu, num_pr_per_pu, update_table_5);
	
		// Shift the swapped and combined tables by -shift
		shiftTable(update_table_5, -shift);
	
		// Add to pu_table
		for(auto itr = pu_table->pus.begin(); itr != pu_table->pus.end(); ++itr) {
			int i = itr->first;
			for(unsigned int j = 0; j < itr->second.prs.size(); ++j) {
				itr->second.prs[j].threshold += update_table_5[i][j];
			}
		}
		// TODO check that the update is actually applied
	} else if(secure_write_algo == SpectrumManager::SecureWriteAlgo::SPC) {
		// TODO
		for(int i = 0; i < num_pu; ++i) {
			std::vector<sInt> this_updates;

			for(int j = 0; j < num_pr_per_pu; ++j) {
				this_updates.push_back(INPUT(parties, 0, 0, bit_count));
			}

			sInt this_index = INPUT(parties, 0, i, bit_count);
			for(unsigned int j = 0; j < updates.size(); ++j) {
				sInt eq_index = (this_index >= updates[j].first) & (this_index <= updates[j].first);
				for(unsigned int y = 0; y < this_updates.size(); ++y) {
					this_updates[y] = eq_index.ifelse(updates[j].second[y], this_updates[y]);
				}
			}

			for(unsigned int j = 0; j < this_updates.size(); ++j) {
				sInt u_a = INPUT(parties, 0, rand() % 100 - 50, bit_count);
				sInt u_b = this_updates[j] - u_a;

				parties[0].reveal(u_a);
				parties[1].reveal(u_b);

				int u = 0;
				if(parties[0].isLocalParty()) {
					u = u_a.getValue();
				} else if(parties[1].isLocalParty()) {
					u = u_b.getValue();
				}

				pu_table->pus[i].prs[j].threshold += u;
			}
			
		}
	}

	if(parties[0].isLocalParty() && secure_write_timer != nullptr) {
		secure_write_timer->end(Timer::secure_write);
	}
}

std::vector<float> SpectrumManager::plainTextRun(const std::vector<SU>& sus, const std::vector<PU>& pus, const std::vector<SS>& sss, Timer* timer, PathLossTable* path_loss_table) const {
	std::map<int, std::vector<const PU*> > pu_groups;
	std::map<int, std::vector<const SS*> > ss_groups;
	if(use_grid) {
		std::map<int, std::vector<int> > pu_int_groups;
		std::map<int, std::vector<int> > ss_int_groups;
		plainTextGrid(pus, sss, &pu_int_groups, &ss_int_groups);

		for(auto itr = pu_int_groups.begin(); itr != pu_int_groups.end(); ++itr) {
			for(unsigned int i = 0; i < itr->second.size(); ++i) {
				pu_groups[itr->first].push_back(&(pus[itr->second[i]]));
			}
		}

		for(auto itr = ss_int_groups.begin(); itr != ss_int_groups.end(); ++itr) {
			for(unsigned int i = 0; i < itr->second.size(); ++i) {
				ss_groups[itr->first].push_back(&(sss[itr->second[i]]));
			}
		}
	}

	std::vector<const PU*> sel_pus;
	std::vector<const SS*> sel_sss;
	if(!use_grid) {
		for(unsigned int j = 0; j < pus.size(); ++j) {
			sel_pus.push_back(&pus[j]);
		}

		for(unsigned int i = 0; i < sss.size(); ++i) {
			sel_sss.push_back(&sss[i]);
		}
	}
	std::vector<float> all_vals;
	for(unsigned int i = 0; i < sus.size(); ++i) {
		if(use_grid) {
			sel_pus.clear();
			sel_sss.clear();

			int group_id = int(sus[i].loc.x / grid_delta_x) + int(sus[i].loc.y / grid_delta_y) * grid_num_x;

			auto pu_itr = pu_groups.find(group_id);
			if(pu_itr == pu_groups.end()) {
				std::cerr << "Invalid group_id: " << group_id << std::endl;
				exit(1);
			}

			sel_pus = pu_itr->second;

			auto ss_itr = ss_groups.find(group_id);
			if(ss_itr == ss_groups.end()) {
				std::cerr << "Invalid group_id: " << group_id << std::endl;
				exit(1);
			}

			sel_sss = ss_itr->second;
		}

		float v = plainTextRadar(sus[i], sel_pus, sel_sss, path_loss_table);
		all_vals.push_back(v);
	}
	return all_vals;
}

void SpectrumManager::plainTextGrid(const std::vector<PU>& pus, const std::vector<SS>& sss,
		std::map<int, std::vector<int> >* pu_int_groups, std::map<int, std::vector<int> >* ss_int_groups) const {
	struct comp {
		bool operator()(const std::pair<float, int>& v1, const std::pair<float, int>& v2) const {
			return v1.first < v2.first;
		}
	};

	// PU Groups
	for(int x = 0; x < grid_num_x; ++x) {
		for(int y = 0; y < grid_num_y; ++y) {
			std::vector<std::pair<float, int> > pu_vals;

			for(unsigned int i = 0; i < pus.size(); ++i) {
				float this_dist = pus[i].loc.dist(Location((x + 0.5) * grid_delta_x, (y + 0.5) * grid_delta_y));

				if((int) pu_vals.size() < grid_min_num_pu || this_dist < pu_vals.front().first) {
					if((int) pu_vals.size() >= grid_min_num_pu) {
						std::pop_heap(pu_vals.begin(), pu_vals.end(), comp());
						pu_vals.pop_back();
					}

					pu_vals.push_back(std::make_pair(this_dist, i));
					std::push_heap(pu_vals.begin(), pu_vals.end(), comp());
				}
			}

			int group_id = x + y * grid_num_x;
			for(unsigned int i = 0; i < pu_vals.size(); ++i) {
				(*pu_int_groups)[group_id].push_back(pu_vals[i].second);
			}
		}
	}

	// SS Groups
	for(int x = 0; x < grid_num_x; ++x) {
		for(int y = 0; y < grid_num_y; ++y) {
			std::vector<std::pair<float, int>> ss_vals;

			for(unsigned int i = 0; i < sss.size(); ++i) {
				float this_dist = sss[i].loc.dist(Location((x + 0.5) * grid_delta_x, (y + 0.5) * grid_delta_y));

				if((int) ss_vals.size() < grid_min_num_ss || this_dist < ss_vals.front().first) {
					if((int) ss_vals.size() >= grid_min_num_ss) {
						std::pop_heap(ss_vals.begin(), ss_vals.end(), comp());
						ss_vals.pop_back();
					}

					ss_vals.push_back(std::make_pair(this_dist, i));
					std::push_heap(ss_vals.begin(), ss_vals.end(), comp());
				}
			}

			int group_id = x + y * grid_num_x;
			for(unsigned int i = 0; i < ss_vals.size(); ++i) {
				(*ss_int_groups)[group_id].push_back(ss_vals[i].second);
			}
		}
	}
}

float SpectrumManager::plainTextRadar(
		const SU& su, const std::vector<const PU*>& pus, const std::vector<const SS*>& sss, PathLossTable* path_loss_table) const {
	// PU selection.
	std::vector<int> pu_inds;
	unsigned int k_pu = numSelect(this->num_pu_selection, pus.size());

	// SS selection.
	std::vector<int> sss_inds;
	std::vector<float> sss_dists;
	unsigned int k_ss = numSelect(this->num_ss_selection, sss.size());
	if(selection_algo == SpectrumManager::SelectionAlgo::NONE) {
		for(unsigned int i = 0; i < k_pu; ++i) {
			pu_inds.push_back(i);
		}

		for(unsigned int i = 0; i < k_ss; ++i) {
			sss_inds.push_back(i);
			sss_dists.push_back(su.loc.dist(sss[i]->loc));
		}
	} else if(selection_algo == SpectrumManager::SelectionAlgo::SORT) {
		std::vector<float> pu_dists;
		for(unsigned int i = 0; i < pus.size(); ++i) {
			float dist = su.loc.dist(pus[i]->loc);
			int ind = i;

			if(i < k_pu) {
				pu_inds.push_back(ind);
				pu_dists.push_back(dist);
			} else {
				for(unsigned int x = 0; x < pu_inds.size(); ++x) {
					if(dist < pu_dists[x]) {
						float tmp_dist = dist;
						dist = pu_dists[x];
						pu_dists[x] = tmp_dist;

						float tmp_ind = ind;
						ind = pu_inds[x];
						pu_inds[x] = tmp_ind;
					}
				}
			}
		}

		for(unsigned int i = 0; i < sss.size(); ++i) {
			float dist = su.loc.dist(sss[i]->loc);
			int ind = i;

			if(i < k_ss) {
				sss_inds.push_back(ind);
				sss_dists.push_back(dist);
			} else {
				for(unsigned int x = 0; x < sss_inds.size(); ++x) {
					if(dist < sss_dists[x]) {
						float tmp_dist = dist;
						dist = sss_dists[x];
						sss_dists[x] = tmp_dist;

						int tmp_ind = ind;
						ind = sss_inds[x];
						sss_inds[x] = tmp_ind;
					}
				}
			}
		}
	} else if(selection_algo == SpectrumManager::SelectionAlgo::RANDOM) {
		std::vector<int> all_pu_inds;
		for(unsigned int i = 0; i < pus.size(); ++i) {
			all_pu_inds.push_back(i);
		}
		std::random_shuffle(all_pu_inds.begin(), all_pu_inds.end());

		for(unsigned int i = 0; i < k_pu; ++i) {
			sss_inds.push_back(all_pu_inds[i]);
		}

		std::vector<int> all_ss_inds;
		for(unsigned int i = 0; i < sss.size(); ++i) {
			all_ss_inds.push_back(i);
		}
		std::random_shuffle(all_ss_inds.begin(), all_ss_inds.end());

		for(unsigned int i = 0; i < k_ss; ++i) {
			sss_inds.push_back(all_ss_inds[i]);
			sss_dists.push_back(su.loc.dist(sss[all_ss_inds[i]]->loc));
		}
	} else {
		std::cerr << "Unsupported selection_algo" << std::endl;
		exit(1);
	}

	// Compute weight
	// w = 1 / d(SU, SS)
	std::vector<float> weights;
	float tmp_sum_weight = 0.0;
	for(unsigned int x = 0; x < sss_inds.size(); ++x) {
		float d = sss_dists[x];
		if (d < 0.0001) {
			d = 0.0001; 
		}

		float w = 1.0 / pow(d, pl_alpha);
		weights.push_back(w);
		tmp_sum_weight += w;
	}

	float transmit_power = std::numeric_limits<float>::infinity();
	// Compute the share of received power for each
	// Indexed by SS first then PU. received_powers[i][j] is the power received at 
	std::vector<std::vector<float> > received_powers;
	for(unsigned int i = 0; i < sss.size(); ++i) {
		std::vector<float> rp_weights;
		float sum_rp_weights = 0.0;
		for(unsigned int j = 0; j < pus.size(); ++j) {
			rp_weights.push_back(1.0 / pow(sss[i]->loc.dist(pus[j]->loc), rp_alpha));
			sum_rp_weights += rp_weights[j];
		}

		received_powers.push_back(std::vector<float>());

		for(unsigned int j = 0; j < pus.size(); ++j) {
			if(utils::unit_type == utils::UnitType::ABS) {
				received_powers[i].push_back(rp_weights[j] * sss[i]->received_power / sum_rp_weights);
			} else if(utils::unit_type == utils::UnitType::DB) {
				received_powers[i].push_back(utils::todBm(rp_weights[j] / sum_rp_weights) + sss[i]->received_power);
			} else {
				std::cerr << "Unsupported unit_type" << std::endl;
				exit(1);
			}
		}
	}
	
	// Compute SU transmit power
	// tp = thresh * sum(w(SS)) / sum(r(SS) / t(PU) * w(SS))
	for(unsigned int y = 0; y < pu_inds.size(); ++y) {
		unsigned int j = pu_inds[y];
		float sum_weight = 0.0;
		float sum_weighted_ratio = 0.0;
		for(unsigned int x = 0; x < sss_inds.size(); ++x) {
			int i = sss_inds[x];
			sum_weight = sum_weight + weights[x];
			if(utils::unit_type == utils::UnitType::ABS) {
				sum_weighted_ratio = sum_weighted_ratio + weights[x] * received_powers[i][j] / pus[j]->transmit_power;
			} else if(utils::unit_type == utils::UnitType::DB) {
				sum_weighted_ratio = sum_weighted_ratio + weights[x] * (received_powers[i][j] - pus[j]->transmit_power);
			}
		}
		
		float this_path_loss = sum_weighted_ratio / sum_weight;
		path_loss_table->addPlaintextPathLoss(su.index, j, 0, this_path_loss);

		float this_transmit_power;
		if(utils::unit_type == utils::UnitType::ABS) {
			this_transmit_power = pus[j]->prs[0].threshold * sum_weight / sum_weighted_ratio;
		} else if(utils::unit_type == utils::UnitType::DB) {
			this_transmit_power = pus[j]->prs[0].threshold - sum_weighted_ratio / sum_weight;
		} else {
			std::cerr << "Unsupported unit_type" << std::endl;
			exit(1);
		}

		if(this_transmit_power < transmit_power) {
			transmit_power = this_transmit_power;
		}
	}

	return transmit_power;
}

void SpectrumManager::initGroupDeviations() {

	std::map<float, std::vector<std::pair<int, int> > > deviations_by_dist_squared;
	for(int xd = 0; xd <= grid_num_x; ++xd) {
		for(int yd = 0; yd <= grid_num_y; ++yd) {
			float this_dist = xd * xd * grid_delta_x * grid_delta_x + yd * yd * grid_delta_y * grid_delta_y;

			deviations_by_dist_squared[this_dist].push_back(std::make_pair(xd, yd));
			
			if(xd > 0) {
				deviations_by_dist_squared[this_dist].push_back(std::make_pair(-xd, yd));
			}

			if(yd > 0) {
				deviations_by_dist_squared[this_dist].push_back(std::make_pair(xd, -yd));
			}

			if(xd > 0 && yd > 0) {
				deviations_by_dist_squared[this_dist].push_back(std::make_pair(-xd, -yd));
			}
		}
	}

	float prev_dist = -1.0;
	for(auto itr = deviations_by_dist_squared.begin(); itr != deviations_by_dist_squared.end(); ++itr) {
		if(fabs(prev_dist - itr->first) < 0.00001) {
			// If dist is close enough, merge the two lists.
			for(unsigned int i = 0; i < itr->second.size(); ++i) {
				all_deviations[all_deviations.size() - 1].push_back(itr->second[i]);
			}
		} else {
			// Else create a new entry.
			all_deviations.push_back(itr->second);
			prev_dist = itr->first;
		}
	}
}

const std::vector<std::pair<int, int> >& SpectrumManager::getDeviations(int iter) const {
	if(iter < 0 || (unsigned int)iter >= all_deviations.size()) {
		std::cerr << "Invalid iter: " << iter << std::endl;
		exit(1);
	}
	return all_deviations[iter];
}
