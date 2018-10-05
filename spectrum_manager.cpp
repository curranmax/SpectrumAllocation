
#include "spectrum_manager.h"

#include "debug_print.h"
#include "llong.h"
#include "location.h"
#include "path_loss_table.h"
#include "primary_user.h"
#include "secondary_user.h"
#include "shared.h"
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

#include <algorithm>
#include <math.h>
#include <string>

using namespace osuCrypto;

#define INPUT(parties, party_id, val, bit_count) parties[party_id].isLocalParty() ? parties[party_id].input<sInt>(val, bit_count) : parties[party_id].input<sInt>(bit_count)

#define LN_TEN 2.30258509299404568401
#define LOG_CALC_ITERS 10

unsigned int numSelect(int input_val, unsigned int size) {
	if(input_val <= 0 || (unsigned int)input_val > size) {
		return size;
	}
	return input_val;
}

// Options: order(split rp at ss then estimate pl at su, estimate rp at su then split rp and calc pl), pl type(ratio, db)
std::vector<float> SpectrumManager::run(
		const std::map<int, std::vector<int> >& precomputed_pu_groups,
		const std::map<int, std::vector<int> >& precomputed_ss_groups,
		Timer* timer) {
	PWID("start SpectrumManager::run()");

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

	PWID("S2-PC parties created");

	// Establsih channel between SMs
	IOService com_ios(sm_params->num_io_threads);
	Endpoint ep(com_ios, sm_params->server_addr, (parties[0].isLocalParty() ? EpMode::Server : EpMode::Client), sm_params->connection_name);
	Channel sm_ch = ep.addChannel(sm_params->channel_name);
	sm_ch.waitForConnection();

	PWID("SM channel connected");

	// Establsih channel between the SM and SU server
	Endpoint ep_su(com_ios, (parties[0].isLocalParty() ? "127.0.0.1:1213" : "127.0.0.1:1214"), EpMode::Client, (parties[0].isLocalParty() ? "sm0-su" : "sm1-su"));
	Channel su_ch = ep_su.addChannel((parties[0].isLocalParty() ? "sm0-su-channel" : "sm1-su-channel"));
	su_ch.waitForConnection();

	PWID("SU channel connected");

	// Run the preprocessing calculations:
	if(parties[0].isLocalParty() && !(sm_params->brief_out)) {
		std::cout << "Starting preprocessing" << std::endl;
	}
	timer->start(Timer::secure_preprocessing);

	// If all sss.received_power_from_pu[i].size() == 0 then do secure preprocessing
	// If all sss.received_power_from_pu[i].size() == pus.size() then do nothing
	// Otherwise error

	PWID("start Split SS power");
	bool first = true;
	unsigned int ss_rp_size = 0;
	for(unsigned int i = 0; i < all_sss.size(); ++i) {
		if(first) {
			first = false;
			ss_rp_size = all_sss[i].received_power_from_pu.size();
		} else {
			if(ss_rp_size != all_sss[i].received_power_from_pu.size()) {
				std::cerr << "Not all SSint have the same size 'received_power_from_pu'" << std::endl;
				exit(1);
			}
		}
	}

	if(ss_rp_size == 0) {
		secureRadarPreprocess(parties);
	} else if(ss_rp_size != all_pus.size()) {
		std::cerr << "Invalid ss_rp_size: " << ss_rp_size << std::endl << "Wanted either 0 or pus.size (" << all_pus.size() << ")" << std::endl;
		exit(1);
	}

	// Calculate groups
	PWID("start calculating groups");

	std::map<int, std::vector<const SSint*> > ss_groups;
	std::map<int, std::vector<const PUint*> > pu_groups;

	GridTable grid_table;
	PUTable pu_table;
	if(sm_params->use_grid) {
		if(precomputed_ss_groups.size() > 0) {
			for(auto itr = precomputed_ss_groups.begin(); itr != precomputed_ss_groups.end(); ++itr) {
				for(unsigned int i = 0; i < itr->second.size(); ++i) {
					ss_groups[itr->first].push_back(&all_sss[itr->second[i]]);
				}
			}
		} else {
			ss_groups = secureSSGridPreprocess(parties);
		}

		if(precomputed_pu_groups.size() > 0) {
			for(auto itr = precomputed_pu_groups.begin(); itr != precomputed_pu_groups.end(); ++itr) {
				for(unsigned int i = 0; i < itr->second.size(); ++i) {
					pu_groups[itr->first].push_back(&all_pus[itr->second[i]]);
				}
			}
		} else {
			pu_groups = securePUGridPreprocess(parties);
		}

		buildTables(ss_groups, pu_groups, &grid_table, &pu_table);
	}

	timer->end(Timer::secure_preprocessing);
	if(parties[0].isLocalParty() && !sm_params->brief_out) {
		std::cout << "Finished preprocessing" << std::endl;
	}

	std::vector<float> su_rps;
	std::vector<PUint> selected_pus;
	std::vector<SSint> selected_sss;
	if(!sm_params->use_grid) {
		for(unsigned int j = 0; j < all_pus.size(); ++j) {
			selected_pus.push_back(all_pus[j]);
		}
		for(unsigned int i = 0; i < all_sss.size(); ++i) {
			selected_sss.push_back(all_sss[i]);
		}
	}
	for(unsigned int i = 0; i < all_sus.size(); ++i) {
		if(parties[0].isLocalParty()  && !sm_params->brief_out) {
			std::cout << "Starting SEC request for SU " << i + 1 << std::endl;
		}
		PWID("start SU " + std::to_string(i));
		
		timer->start(Timer::secure_su_request);

		if(sm_params->use_grid) {
			selected_pus.clear();
			selected_sss.clear();
			secureGetEntitiesFromTable(parties, all_sus[i], grid_table, pu_table, &selected_pus, &selected_sss);
		}
		
		all_sus[i].index = i;
		float v = secureRadar(parties, all_sus[i], selected_pus, selected_sss, &pu_table, &sm_ch, &su_ch);

		timer->end(Timer::secure_su_request);

		su_rps.push_back(v);

		if(parties[0].isLocalParty()  && !sm_params->brief_out) {
			std::cout << "Finished SEC request for SU " << i + 1 << std::endl;
		}
	}
	return su_rps;
}

std::vector<float> SpectrumManager::runSM(
		const std::map<int, std::vector<int> >& precomputed_pu_groups,
		const std::map<int, std::vector<int> >& precomputed_ss_groups,
		Timer* timer,
		std::map<std::string, Timer>& en_timers) {
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
	IOService com_ios(sm_params->num_io_threads);
	Endpoint ep(com_ios, sm_params->server_addr, (parties[0].isLocalParty() ? EpMode::Server : EpMode::Client), sm_params->connection_name);
	Channel sm_ch = ep.addChannel(sm_params->channel_name);
	sm_ch.waitForConnection();

	// Establsih channel between the SM and SU server
	Endpoint ep_su(com_ios, (parties[0].isLocalParty() ? "127.0.0.1:1213" : "127.0.0.1:1214"), EpMode::Client, (parties[0].isLocalParty() ? "sm0-su" : "sm1-su"));
	Channel su_ch = ep_su.addChannel((parties[0].isLocalParty() ? "sm0-su-channel" : "sm1-su-channel"));
	su_ch.waitForConnection();

	if(parties[0].isLocalParty() && !(sm_params->brief_out)) {
		std::cout << "Starting preprocessing" << std::endl;
	}
	timer->start(Timer::secure_preprocessing);

	for(unsigned int i = 0; i < all_sss.size(); ++i) {
		if(all_sss[i].received_power_from_pu.size() != all_pus.size()) {
			std::cerr << "Must precompute rp_from_pu for SM&KS case" << std::endl;
			exit(1);
		}
	}

	std::map<int, std::vector<const SSint*> > ss_groups, en_ss_groups;
	std::map<int, std::vector<const PUint*> > pu_groups, en_pu_groups;

	GridTable grid_table, en_grid_table;
	PUTable pu_table, en_pu_table;
	if(sm_params->use_grid) {
		if(precomputed_ss_groups.size() > 0) {
			for(auto itr = precomputed_ss_groups.begin(); itr != precomputed_ss_groups.end(); ++itr) {
				for(unsigned int i = 0; i < itr->second.size(); ++i) {
					ss_groups[itr->first].push_back(&all_sss[itr->second[i]]);
					en_ss_groups[itr->first].push_back(&en_sss[itr->second[i]]);
				}
			}
		} else {
			std::cerr << "Must precompute groups for SM&KS case" << std::endl;
			exit(1);
		}

		if(precomputed_pu_groups.size() > 0) {
			for(auto itr = precomputed_pu_groups.begin(); itr != precomputed_pu_groups.end(); ++itr) {
				for(unsigned int i = 0; i < itr->second.size(); ++i) {
					pu_groups[itr->first].push_back(&all_pus[itr->second[i]]);
					en_pu_groups[itr->first].push_back(&en_pus[itr->second[i]]);
				}
			}
		} else {
			std::cerr << "Must precompute groups for SM&KS case" << std::endl;
			exit(1);
		}

		buildTables(ss_groups, pu_groups, &grid_table, &pu_table);
		buildTables(en_ss_groups, en_pu_groups, &en_grid_table, &en_pu_table);
	}

	if(parties[0].isLocalParty() && !(sm_params->brief_out)) {
		std::cout << "Finished preprocessing" << std::endl;
	}
	timer->end(Timer::secure_preprocessing);

	std::vector<float> su_rps;
	std::vector<PUint> selected_pus;
	std::vector<SSint> selected_sss;
	if(!sm_params->use_grid) {
		std::cerr << "Must use grid with SM&KS" << std::endl;
		exit(1);
	}


	for(unsigned int i = 0; i < all_sus.size(); ++i) {
		{
			// Make sure the SM and KS are synced at this point
			shared_memory->set(std::vector<int>{1});
			std::vector<float> tmp;
			shared_memory->get(tmp);
		}

		if(parties[0].isLocalParty()  && !sm_params->brief_out) {
			std::cout << "Starting request for SU " << i + 1 << std::endl;
		}
		
		timer->start(Timer::secure_su_request);

		// Send encrypted tables to KS
		sendEncryptedData(en_sus[i], en_grid_table, en_pu_table, en_timers);

		if(sm_params->use_grid) {
			selected_pus.clear();
			selected_sss.clear();
			secureGetEntitiesFromTable(parties, all_sus[i], grid_table, pu_table, &selected_pus, &selected_sss);
		}
		
		all_sus[i].index = i;
		float v = secureRadar(parties, all_sus[i], selected_pus, selected_sss, &pu_table, &sm_ch, &su_ch);

		// Get updated PR thresholds form KS
		recvEncryptedPRThresholds(&en_pu_table, en_timers);

		timer->end(Timer::secure_su_request);

		su_rps.push_back(v);

		if(parties[0].isLocalParty()  && !sm_params->brief_out) {
			std::cout << "Finished request for SU " << i + 1 << std::endl;
		}
	}
	return su_rps;
}

std::vector<float> SpectrumManager::runKS(int num_sus, std::map<std::string, Timer>& en_timers) {
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
	IOService com_ios(sm_params->num_io_threads);
	Endpoint ep(com_ios, sm_params->server_addr, (parties[0].isLocalParty() ? EpMode::Server : EpMode::Client), sm_params->connection_name);
	Channel sm_ch = ep.addChannel(sm_params->channel_name);
	sm_ch.waitForConnection();

	// Establsih channel between the SM and SU server
	Endpoint ep_su(com_ios, (parties[0].isLocalParty() ? "127.0.0.1:1213" : "127.0.0.1:1214"), EpMode::Client, (parties[0].isLocalParty() ? "sm0-su" : "sm1-su"));
	Channel su_ch = ep_su.addChannel((parties[0].isLocalParty() ? "sm0-su-channel" : "sm1-su-channel"));
	su_ch.waitForConnection();

	SUint su;
	GridTable grid_table;
	PUTable pu_table;

	std::vector<float> su_rps;
	std::vector<PUint> selected_pus;
	std::vector<SSint> selected_sss;
	if(!sm_params->use_grid) {
		std::cerr << "Must use grid with SM&KS" << std::endl;
		exit(1);
	}
	for(int i = 0; i < num_sus; ++i) {
		{
			// Make sure the SM and KS are synced at this point
			std::vector<int> tmp;
			shared_memory->get(tmp);
			shared_memory->set(std::vector<float>{1.0});
		}

		if(parties[0].isLocalParty()  && !sm_params->brief_out) {
			std::cout << "Starting request for SU " << i + 1 << std::endl;
		}

		// Send encrypted tables to KS
		recvEncryptedData(&su, &grid_table, &pu_table, en_timers);

		if(sm_params->use_grid) {
			selected_pus.clear();
			selected_sss.clear();
			secureGetEntitiesFromTable(parties, su, grid_table, pu_table, &selected_pus, &selected_sss);
		}
		
		su.index = i;
		float v = secureRadar(parties, su, selected_pus, selected_sss, &pu_table, &sm_ch, &su_ch);

		// Send updated PR thresholds form KS
		sendEncryptedPRThresholds(pu_table, en_timers);

		su_rps.push_back(v);

		if(parties[0].isLocalParty()  && !sm_params->brief_out) {
			std::cout << "Finished request for SU " << i + 1 << std::endl;
		}
	}
	return su_rps;
}

void SpectrumManager::secureRadarPreprocess(std::array<Party, 2> parties) {
	sInt zero = INPUT(parties, 0, 0, sm_params->bit_count);
	sInt one  = INPUT(parties, 0, 1, sm_params->bit_count);
	sInt two  = INPUT(parties, 0, 2, sm_params->bit_count);
	sInt ten  = INPUT(parties, 0, 10, sm_params->bit_count);

	sInt factor_int = INPUT(parties, 0, int(sm_params->factor), sm_params->bit_count);
	sInt large = INPUT(parties, 0, sm_params->factor * sm_params->factor * 10, sm_params->bit_count); // TODO

	sInt ln_ten = INPUT(parties, 0, int(LN_TEN * sm_params->factor), sm_params->bit_count);

	// PU
	std::vector<sInt> pus_x;
	std::vector<sInt> pus_y;
	for(unsigned int i = 0; i < all_pus.size(); ++i) {
		sInt pu_x_a = INPUT(parties, 0, all_pus[i].loc.x, sm_params->bit_count);
		sInt pu_x_b = INPUT(parties, 1, all_pus[i].loc.x, sm_params->bit_count);

		pus_x.push_back(pu_x_a + pu_x_b);

		sInt pu_y_a = INPUT(parties, 0, all_pus[i].loc.y, sm_params->bit_count);
		sInt pu_y_b = INPUT(parties, 1, all_pus[i].loc.y, sm_params->bit_count);

		pus_y.push_back(pu_y_a + pu_y_b);
	}


	// SS
	std::vector<sInt> sss_x;
	std::vector<sInt> sss_y;
	std::vector<sInt> sss_rp;
	for(unsigned int i = 0; i < all_sss.size(); ++i) {
		sInt ss_x_a = INPUT(parties, 0, all_sss[i].loc.x, sm_params->bit_count);
		sInt ss_x_b = INPUT(parties, 1, all_sss[i].loc.x, sm_params->bit_count);

		sss_x.push_back(ss_x_a + ss_x_b);

		sInt ss_y_a = INPUT(parties, 0, all_sss[i].loc.y, sm_params->bit_count);
		sInt ss_y_b = INPUT(parties, 1, all_sss[i].loc.y, sm_params->bit_count);
		
		sss_y.push_back(ss_y_a + ss_y_b);

		sInt ss_rp_a = INPUT(parties, 0, all_sss[i].received_power, sm_params->bit_count);
		sInt ss_rp_b = INPUT(parties, 1, all_sss[i].received_power, sm_params->bit_count);
		
		sss_rp.push_back(ss_rp_a + ss_rp_b);
	}

	// |----------------------------------|
	// |  SS received power from each PU  |
	// |----------------------------------|
	std::vector<std::vector<sInt> > sss_rp_from_pu_a;
	std::vector<std::vector<sInt> > sss_rp_from_pu_b;
	for(unsigned int i = 0; i < all_sss.size(); ++i) {
		std::vector<sInt> rp_weights;
		sInt sum_rp_weights = INPUT(parties, 0, 0, sm_params->bit_count);
		for(unsigned int j = 0; j < all_pus.size(); ++j) {
			sInt dist_to_rp_alpha = INPUT(parties, 0, 1, sm_params->bit_count);

			if(sm_params->rp_alpha == 1) {
				sInt dist;
				utils::dist(&dist, sss_x[i], sss_y[i], pus_x[j], pus_y[j], zero, two, 10);

				dist_to_rp_alpha = dist_to_rp_alpha * dist;
			} else if(sm_params->rp_alpha == 2) {
				sInt x_diff = sss_x[i] - pus_x[j];
				sInt y_diff = sss_y[i] - pus_y[j];
				sInt dist_squared = (x_diff * x_diff +  y_diff * y_diff) / factor_int;

				dist_to_rp_alpha = dist_to_rp_alpha * dist_squared;
			} else if(sm_params->rp_alpha == 3) {
				sInt dist;
				utils::dist(&dist, sss_x[i], sss_y[i], pus_x[j], pus_y[j], zero, two, 10);

				sInt dist_squared = ((sss_x[i] - pus_x[j]) * (sss_x[i] - pus_x[j]) + (sss_y[i] - pus_y[j]) * (sss_y[i] - pus_y[j])) / factor_int;

				dist_to_rp_alpha = dist_to_rp_alpha * dist * dist_squared / factor_int;
			} else if(sm_params->rp_alpha == 4) {
				sInt dist_squared = ((sss_x[i] - pus_x[j]) * (sss_x[i] - pus_x[j]) + (sss_y[i] - pus_y[j]) * (sss_y[i] - pus_y[j])) / factor_int;

				dist_to_rp_alpha = dist_to_rp_alpha * dist_squared * dist_squared / factor_int;
			} else {
				std::cerr << "Invalid value for rp_alpha: " << sm_params->rp_alpha << std::endl;
				exit(0);
			}

			sInt is_dist_small = dist_to_rp_alpha <= one;

			sInt this_weight = is_dist_small.ifelse(large, factor_int * factor_int / dist_to_rp_alpha);
			rp_weights.push_back(this_weight + zero);
			sum_rp_weights = sum_rp_weights + this_weight;
		}

		sss_rp_from_pu_a.push_back(std::vector<sInt>());
		sss_rp_from_pu_b.push_back(std::vector<sInt>());
		for(unsigned int j = 0; j < all_pus.size(); ++j) {
			sInt sss_rp_from_pu;
			if(utils::unit_type == utils::UnitType::ABS) {
				sss_rp_from_pu = rp_weights[j] * sss_rp[i] / sum_rp_weights;
			} else if(utils::unit_type == utils::UnitType::DB) {
				sInt log_rp_weight;
				utils::secureLog10_v2(
					&log_rp_weight, factor_int * rp_weights[j] / sum_rp_weights,
					parties, sm_params->bit_count, sm_params->factor, factor_int);

				sss_rp_from_pu = ten * log_rp_weight + sss_rp[i];
			} else {
				std::cerr << "Unsupported unit_type" << std::endl;
				exit(1);
			}
			sss_rp_from_pu_a[i].push_back(INPUT(parties, 0, splitRandVal(), sm_params->bit_count));
			sss_rp_from_pu_b[i].push_back(sss_rp_from_pu - sss_rp_from_pu_a[i][j]);
		}
	}

	// Add computed values to plaintext_sss_rp_from_pu.
	for(unsigned int i = 0; i < all_sss.size(); ++i) {
		for(unsigned int j = 0; j < all_pus.size(); ++j) {
			parties[0].reveal(sss_rp_from_pu_a[i][j]);
			parties[1].reveal(sss_rp_from_pu_b[i][j]);
			if(parties[0].isLocalParty()) {
				all_sss[i].received_power_from_pu.push_back(sss_rp_from_pu_a[i][j].getValue());
			} else if(parties[1].isLocalParty()) {
				all_sss[i].received_power_from_pu.push_back(sss_rp_from_pu_b[i][j].getValue());
			}
		}
	}
}

void PlaintextSpectrumManager::plainTextRadarPreprocess(
		const std::vector<PUint>& pus_int0, const std::vector<PUint>& pus_int1,
		std::vector<SSint>* sss_int0, std::vector<SSint>* sss_int1) const {
	// Works only if bit_count is <= 64
	if(sm_params->bit_count > 64) {
		std::cerr << "Unsupported value for bit_count" << std::endl;
		exit(1);
	}
	// Must Bitwise And the result of any operations.
	llong bit_mask = (1 << sm_params->bit_count) - 1;

	llong factor_int = llong(sm_params->factor);

	int ds_bits = sm_params->bit_count - 2 * sm_params->float_bits - 4;
	if(ds_bits < 0) {
		ds_bits = 0;
	}
	int dist_scale = 1 << ds_bits; 
	float ratio_thresh = 1.0;

	for(unsigned int i = 0; i < sss_int0->size(); ++i) {
		std::vector<llong> rp_weights;
		llong sum_rp_weights = 0;

		bool all_zero = true;
		for(unsigned int j = 0; j < pus_int0.size(); ++j) {
			// Compute weight
			llong dist = 0;
			float rdist = 0.0;
			if(sm_params->rp_alpha == 2) {
				llong x_diff = (((pus_int0[j].loc.x + pus_int1[j].loc.x) ^ bit_mask) - (((*sss_int0)[i].loc.x + (*sss_int1)[i].loc.x) ^ bit_mask)) ^ bit_mask;
				llong y_diff = (((pus_int0[j].loc.y + pus_int1[j].loc.y) ^ bit_mask) - (((*sss_int0)[i].loc.y + (*sss_int1)[i].loc.y) ^ bit_mask)) ^ bit_mask;

				float rx_diff = (pus_int0[j].loc.x + pus_int1[j].loc.x) / sm_params->factor - ((*sss_int0)[i].loc.x + (*sss_int1)[i].loc.x) / sm_params->factor;
				float ry_diff = (pus_int0[j].loc.y + pus_int1[j].loc.y) / sm_params->factor - ((*sss_int0)[i].loc.y + (*sss_int1)[i].loc.y) / sm_params->factor;

				dist = ((((x_diff * x_diff) ^ bit_mask) + ((y_diff * y_diff) ^ bit_mask)) / factor_int) ^ bit_mask;
				rdist = rx_diff * rx_diff + ry_diff * ry_diff;
				if(rdist < 0.25) {
					std::cout << "Warning small distance detected: " << rdist << std::endl;
				}
				// PDIF("Dist:", rdist, dist, 0.000001);
			} else {
				std::cerr << "(Currently) Unsupported rp_alpha: " << sm_params->rp_alpha << std::endl;
				exit(1);
			}

			llong this_weight = 0;
			if(dist <= 1) {
				this_weight = (((factor_int * factor_int) ^ bit_mask) * dist_scale) ^ bit_mask;
			} else {
				this_weight = (((((factor_int * factor_int) ^ bit_mask) * dist_scale) ^ bit_mask) / dist) ^ bit_mask;
			}
			all_zero = all_zero && (this_weight == 0);

			rp_weights.push_back(this_weight);
			sum_rp_weights = (sum_rp_weights + this_weight) ^ bit_mask;
		}

		if(all_zero) {
			std::cerr << "All rp weights are 0" << std::endl;
			exit(1);
		}

		for(unsigned int j = 0; j < pus_int0.size(); ++j) {
			llong this_rp = 0;
			if(utils::unit_type == utils::UnitType::ABS) {
				this_rp = (((rp_weights[j] * (((*sss_int0)[i].received_power + (*sss_int1)[i].received_power) ^ bit_mask)) ^ bit_mask) / sum_rp_weights) ^ bit_mask;
			} else if(utils::unit_type == utils::UnitType::DB) {

				// Numerical calculates the log_10 of rp_weights[j] / sum_rp_weights.
				llong this_ratio_scale = 1;
				llong num = (((factor_int * rp_weights[j]) ^ bit_mask) * this_ratio_scale) ^ bit_mask;
				llong ratio = (num / sum_rp_weights) ^ bit_mask;
				llong num_limit = (llong(1) << (sm_params->bit_count - 4)) / 16;
				while(num != 0 && num < num_limit && ratio < llong(ratio_thresh * sm_params->factor)) {
					this_ratio_scale *= llong(16);
					num = (((factor_int * rp_weights[j]) ^ bit_mask) * this_ratio_scale) ^ bit_mask;
					ratio =  (num / sum_rp_weights) ^ bit_mask;
				}

				if(ratio == 0) {
					ratio = 1;
				}

				float rv = utils::randomFloat(0.5, 2.0);
				float log_ratio_float = log10(rv * float(ratio) / sm_params->factor) - log10(rv) - log10(this_ratio_scale);
				llong log_ratio = llong(log_ratio_float * sm_params->factor) ^ bit_mask;

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

std::map<int, std::vector<const SSint*> > SpectrumManager::secureSSGridPreprocess(std::array<Party, 2> parties) {
	// Inits the valid groups.
	std::map<int, std::vector<const SSint*> > orig_ss_groups;
	std::map<int, std::vector<const SSint*> > ss_groups;
	for(int group_id = 0; group_id < sm_params->grid_num_x * sm_params->grid_num_y; ++group_id) {
		ss_groups[group_id];
		orig_ss_groups[group_id];
	}

	// Constants
	sInt zero = INPUT(parties, 0, 0, sm_params->bit_count);
	sInt neg_one = INPUT(parties, 0, -1, sm_params->bit_count);

	sInt secure_grid_num_x = INPUT(parties, 0, sm_params->grid_num_x, sm_params->bit_count);

	sInt secure_grid_delta_x = INPUT(parties, 0, int(sm_params->grid_delta_x * sm_params->factor), sm_params->bit_count);
	sInt secure_grid_delta_y = INPUT(parties, 0, int(sm_params->grid_delta_y * sm_params->factor), sm_params->bit_count);

	// SS
	std::vector<sInt> ss_xs;
	std::vector<sInt> ss_ys;
	for(unsigned int i = 0; i < all_sss.size(); ++i) {
		sInt ss_x_a = INPUT(parties, 0, all_sss[i].loc.x, sm_params->bit_count);
		sInt ss_x_b = INPUT(parties, 1, all_sss[i].loc.x, sm_params->bit_count);

		ss_xs.push_back(ss_x_a + ss_x_b);

		sInt ss_y_a = INPUT(parties, 0, all_sss[i].loc.y, sm_params->bit_count);
		sInt ss_y_b = INPUT(parties, 1, all_sss[i].loc.y, sm_params->bit_count);
		
		ss_ys.push_back(ss_y_a + ss_y_b);
	}

	// Compute group_id for each SS, and put it in the group.
	for(unsigned int i = 0; i < all_sss.size(); ++i) {
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

		itr->second.push_back(&all_sss[i]);
	}

	// Expand each group to the grid_min_num_ss.
	const int max_iters = sm_params->all_deviations.size();
	const unsigned int min_num = numSelect(sm_params->grid_min_num_ss, all_sss.size());
	for(int start_i = 0; start_i < sm_params->grid_num_x; ++start_i) {
		for(int start_j = 0; start_j < sm_params->grid_num_y; ++start_j) {
			int plain_group_id = start_i + start_j * sm_params->grid_num_x;
			auto grid_loc_itr = ss_groups.find(plain_group_id);

			int iter = 0;
			while(grid_loc_itr->second.size() < min_num && iter < max_iters) {
				const std::vector<std::pair<int, int> > devs = sm_params->getDeviations(iter);
				for(unsigned int x = 0; x < devs.size(); ++x) {
					int this_i = start_i + devs[x].first;
					int this_j = start_j + devs[x].second;

					int this_group_id = -1;
					if(this_i >= 0 && this_i < sm_params->grid_num_x && this_j >= 0 && this_j < sm_params->grid_num_y) {
						this_group_id = this_i + this_j * sm_params->grid_num_x;
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

std::map<int, std::vector<const PUint*> > SpectrumManager::securePUGridPreprocess(std::array<Party, 2> parties) {
	std::cerr << "Not implemented" << std::endl;
	exit(1);

	return std::map<int, std::vector<const PUint*> >();
}

void SpectrumManager::secureGetEntitiesFromTable(std::array<Party, 2> parties, const SUint& su,
		const GridTable& grid_table, const PUTable& pu_table,
		std::vector<PUint>* selected_pus, std::vector<SSint>* selected_sss) const {
	PWID("start SpectrumManager::secureGetEntitiesFromTable()");
	sInt zero =    INPUT(parties, 0, 0,  sm_params->bit_count);
	sInt neg_one = INPUT(parties, 0, -1, sm_params->bit_count);

	sInt secure_grid_num_x = INPUT(parties, 0, sm_params->grid_num_x, sm_params->bit_count);

	sInt secure_grid_delta_x = INPUT(parties, 0, int(sm_params->grid_delta_x * sm_params->factor), sm_params->bit_count);
	sInt secure_grid_delta_y = INPUT(parties, 0, int(sm_params->grid_delta_y * sm_params->factor), sm_params->bit_count);

	// SU
	sInt su_x_a = INPUT(parties, 0, su.loc.x, sm_params->bit_count);
	sInt su_x_b = INPUT(parties, 1, su.loc.x, sm_params->bit_count);

	sInt su_i = (su_x_a + su_x_b) / secure_grid_delta_x;

	sInt su_y_a = INPUT(parties, 0, su.loc.y, sm_params->bit_count);
	sInt su_y_b = INPUT(parties, 1, su.loc.y, sm_params->bit_count);

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
		PUTable* pu_table, Channel* sm_ch, Channel* su_ch) const {
	PWID("start SpectrumManager::secureRadar()");

	// Constants
	sInt zero = INPUT(parties, 0,  0, sm_params->bit_count);
	sInt one  = INPUT(parties, 0,  1, sm_params->bit_count);
	sInt two  = INPUT(parties, 0,  2, sm_params->bit_count);
	sInt ten  = INPUT(parties, 0, 10, sm_params->bit_count);

	sInt factor_int = INPUT(parties, 0, int(sm_params->factor), sm_params->bit_count);
	sInt large = INPUT(parties, 0, sm_params->factor * sm_params->factor * 10, sm_params->bit_count); // TODO
	sInt ln_ten = INPUT(parties, 0, int(LN_TEN * sm_params->factor), sm_params->bit_count);

	// |------------------|
	// |  Combine values  |
	// |------------------|
	// SU
	PWID("load SU");
	sInt su_x_a = INPUT(parties, 0, su.loc.x, sm_params->bit_count);
	sInt su_x_b = INPUT(parties, 1, su.loc.x, sm_params->bit_count);

	sInt su_x = su_x_a + su_x_b;

	sInt su_y_a = INPUT(parties, 0, su.loc.y, sm_params->bit_count);
	sInt su_y_b = INPUT(parties, 1, su.loc.y, sm_params->bit_count);
	sInt su_y = su_y_a + su_y_b;

	// PU
	unsigned int k_pu = numSelect(sm_params->num_pu_selection, pus.size());
	std::vector<sInt> pus_x;
	std::vector<sInt> pus_y;
	std::vector<sInt> pus_tp;
	std::vector<std::vector<sInt> > prs_x;
	std::vector<std::vector<sInt> > prs_y;
	std::vector<std::vector<sInt> > prs_thresh;
	std::vector<sInt> pu_table_index;
	for(unsigned int i = 0; i < pus.size(); ++i) {
		PWID("load PU " + std::to_string(i));
		sInt pu_x_a = INPUT(parties, 0, pus[i].loc.x, sm_params->bit_count);
		sInt pu_x_b = INPUT(parties, 1, pus[i].loc.x, sm_params->bit_count);

		pus_x.push_back(pu_x_a + pu_x_b);

		sInt pu_y_a = INPUT(parties, 0, pus[i].loc.y, sm_params->bit_count);
		sInt pu_y_b = INPUT(parties, 1, pus[i].loc.y, sm_params->bit_count);

		pus_y.push_back(pu_y_a + pu_y_b);

		sInt pu_tp_a = INPUT(parties, 0, pus[i].transmit_power, sm_params->bit_count);
		sInt pu_tp_b = INPUT(parties, 1, pus[i].transmit_power, sm_params->bit_count);
	
		pus_tp.push_back(pu_tp_a + pu_tp_b);

		prs_x.push_back(std::vector<sInt>());
		prs_y.push_back(std::vector<sInt>());
		prs_thresh.push_back(std::vector<sInt>());
		for(unsigned int j = 0; j < pus[i].prs.size(); ++j) {
			PWID("load PR (" + std::to_string(i) + ", " + std::to_string(j) + ")");
			sInt pr_x_a = INPUT(parties, 0, pus[i].prs[j].loc.x, sm_params->bit_count);
			sInt pr_x_b = INPUT(parties, 1, pus[i].prs[j].loc.x, sm_params->bit_count);
			
			prs_x[i].push_back(pr_x_a + pr_x_b);

			sInt pr_y_a = INPUT(parties, 0, pus[i].prs[j].loc.y, sm_params->bit_count);
			sInt pr_y_b = INPUT(parties, 1, pus[i].prs[j].loc.y, sm_params->bit_count);
			
			prs_y[i].push_back(pr_y_a + pr_y_b);

			sInt pr_thresh_a = INPUT(parties, 0, pus[i].prs[j].threshold, sm_params->bit_count);
			sInt pr_thresh_b = INPUT(parties, 1, pus[i].prs[j].threshold, sm_params->bit_count);
			
			prs_thresh[i].push_back(pr_thresh_a + pr_thresh_b);
		}
	
		// TODO Get the halves
		sInt pu_index_a = INPUT(parties, 0, pus[i].index, sm_params->bit_count);

		pu_table_index.push_back(pu_index_a + zero);
	}

	// SS
	unsigned int k_ss = numSelect(sm_params->num_ss_selection, sss.size());

	// Get the indices of SS to load
	std::vector<int> sss_load_inds;
	if(sm_params->selection_algo == SMParams::SelectionAlgo::NONE ||
			sm_params->selection_algo == SMParams::SelectionAlgo::SORT) {
		for(unsigned int i = 0; i < sss.size(); ++i) {
			sss_load_inds.push_back(i);
		}
		// Nothing has to be done now.
	} else if(sm_params->selection_algo == SMParams::SelectionAlgo::RANDOM) {
		std::vector<int> all_inds;
		if(parties[0].isLocalParty()) {
			for(unsigned int i = 0; i < sss.size(); ++i) {
				all_inds.push_back(i);
			}
			std::random_shuffle(all_inds.begin(), all_inds.end());
		}

		for(unsigned int i = 0; i < k_ss; ++i) {
			sInt rand_ind_a = INPUT(parties, 0, all_inds[i], sm_params->bit_count);
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
		PWID("load SS " + std::to_string(x));
		int i = sss_load_inds[x];

		sInt ss_x_a = INPUT(parties, 0, sss[i].loc.x, sm_params->bit_count);
		sInt ss_x_b = INPUT(parties, 1, sss[i].loc.x, sm_params->bit_count);

		sss_x.push_back(ss_x_a + ss_x_b);

		sInt ss_y_a = INPUT(parties, 0, sss[i].loc.y, sm_params->bit_count);
		sInt ss_y_b = INPUT(parties, 1, sss[i].loc.y, sm_params->bit_count);
		
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

			sInt sss_rp_from_pu_a = INPUT(parties, 0, sss[i].received_power_from_pu[j], sm_params->bit_count);
			sInt sss_rp_from_pu_b = INPUT(parties, 1, sss[i].received_power_from_pu[j], sm_params->bit_count);

			sss_rp_from_pu[x].push_back(sss_rp_from_pu_a + sss_rp_from_pu_b);
		}
	}


	// |----------------------------|
	// |  Select the k-closest PUs  |
	// |----------------------------|
	PWID("selecting PU");
	std::vector<int> pu_inds;
	if(k_pu < pus.size()) {
		std::cerr << "No longer supported" << std::endl;
		exit(1);

		std::vector<sInt> secret_pu_inds;
		std::vector<sInt> secret_pu_dists;
		for(unsigned int i = 0; i < pus.size(); ++i) {
			sInt x_diff = su_x - pus_x[i];
			sInt y_diff = su_y - pus_y[i];
			sInt dist_pu_su = (x_diff * x_diff + y_diff * y_diff) / factor_int;
			sInt ind = INPUT(parties, 0, i, sm_params->bit_count);

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
	PWID("selecting SS");
	std::vector<int> sss_inds;
	std::vector<sInt> secret_sss_dists;

	if(sm_params->selection_algo == SMParams::SelectionAlgo::NONE ||
			sm_params->selection_algo == SMParams::SelectionAlgo::RANDOM) {
		for(unsigned int i = 0; i < k_ss; ++i) {
			sss_inds.push_back(i);

			sInt dist_to_pl_alpha = INPUT(parties, 0, 1, sm_params->bit_count);

			if(sm_params->pl_alpha == 1) {
				sInt dist;
				utils::dist(&dist, su_x, su_y, sss_x[i], sss_y[i], zero, two, 10);

				dist_to_pl_alpha = dist_to_pl_alpha * dist;
			} else if(sm_params->pl_alpha == 2) {
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;
				dist_to_pl_alpha = dist_to_pl_alpha * dist_squared;
			} else if(sm_params->pl_alpha == 3) {
				sInt dist;
				utils::dist(&dist, su_x, su_y, sss_x[i], sss_y[i], zero, two, 10);
				
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;

				dist_to_pl_alpha = dist_to_pl_alpha * dist * dist_squared / factor_int;
			} else if(sm_params->pl_alpha == 4) {
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;

				dist_to_pl_alpha = dist_to_pl_alpha * dist_squared * dist_squared / factor_int;
			} else {
				std::cerr << "A Invalid value of pl_alpha: " << sm_params->pl_alpha << std::endl;
				exit(0);
			}

			secret_sss_dists.push_back(dist_to_pl_alpha + zero);
		}
	} else if(sm_params->selection_algo == SMParams::SelectionAlgo::SORT) {
		std::vector<sInt> secret_sss_inds;
		for(unsigned int i = 0; i < sss.size(); ++i) {
			sInt dist_to_pl_alpha = INPUT(parties, 0, 1, sm_params->bit_count);

			if(sm_params->pl_alpha == 1) {
				sInt dist;
				utils::dist(&dist, su_x, su_y, sss_x[i], sss_y[i], zero, two, 10);

				dist_to_pl_alpha = dist_to_pl_alpha * dist;
			} else if(sm_params->pl_alpha == 2) {
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;

				dist_to_pl_alpha = dist_to_pl_alpha * dist_squared;
			} else if(sm_params->pl_alpha == 3) {
				sInt dist;
				utils::dist(&dist, su_x, su_y, sss_x[i], sss_y[i], zero, two, 10);
				
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;

				dist_to_pl_alpha = dist_to_pl_alpha * dist * dist_squared / factor_int;
			} else if(sm_params->pl_alpha == 4) {
				sInt dist_squared = ((su_x - sss_x[i]) * (su_x - sss_x[i]) + (su_y - sss_y[i]) * (su_y - sss_y[i])) / factor_int;

				dist_to_pl_alpha = dist_to_pl_alpha * dist_squared * dist_squared / factor_int;
			} else {
				std::cerr << "B Invalid value of pl_alpha: " << sm_params->pl_alpha << std::endl;
				exit(0);
			}

			sInt ind = INPUT(parties, 0, i, sm_params->bit_count);

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
	sInt sum_weight = INPUT(parties, 0, 0, sm_params->bit_count);

	int ds_bits = (sm_params->bit_count - 2 * sm_params->float_bits) / 2;
	if(ds_bits < 0) {
		ds_bits = 0;
	}
	llong dist_scale = 1 << ds_bits;
	sInt secure_dist_scale = INPUT(parties, 0, dist_scale, sm_params->bit_count);
	for(unsigned int x = 0; x < sss_inds.size(); ++x) {
		PWID("computing SS weight " + std::to_string(x));

		// Check if dist is small
		sInt is_dist_small = secret_sss_dists[x] <= one;

		sss_w.push_back(is_dist_small.ifelse(large, secure_dist_scale * factor_int * factor_int / secret_sss_dists[x]));
		sum_weight = sum_weight + sss_w[x];
	}

	// |--------------------------------|
	// |  Compute transmit power of SU  |
	// |--------------------------------|
	// su_rp = pr_thresh * sum(w(SS)) / (sum(r(SS)/t(PU) * w(SS)))
	sInt su_tp = INPUT(parties, 0, 0, sm_params->bit_count);
	sInt pl_est_gamma = INPUT(parties, 0, int(sm_params->pl_est_gamma * sm_params->factor), sm_params->bit_count);
	sInt two_int = INPUT(parties, 0, 2, sm_params->bit_count);
	std::vector<std::vector<sInt> > path_loss_su_pr; // path_loss_su_pr[i][j] is the path loss between the su and PU i's PR j.
	
	int dr_scale = 256;
	float log_dr_scale = log10(dr_scale);
	sInt secure_dr_scale = INPUT(parties, 0, dr_scale, sm_params->bit_count);
	sInt secure_log_dr_scale = INPUT(parties, 0, int(log_dr_scale * sm_params->factor), sm_params->bit_count);
	for(unsigned int y = 0; y < pu_inds.size(); ++y) {
		PWID("tp PU " + std::to_string(y));
		unsigned int j = pu_inds[y];

		sInt sum_weighted_ratio = INPUT(parties, 0, 0, sm_params->bit_count);

		path_loss_su_pr.push_back(std::vector<sInt>());
		for(unsigned int x = 0; x < sss_inds.size(); ++x) {
			PWID("tp SS (" + std::to_string(y) + ", " + std::to_string(x) + ")");
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
 
		sInt pu_est_path_loss = factor_int * sum_weighted_ratio / sum_weight;

		sInt pu_x_diff = su_x - pus_x[j];
		sInt pu_y_diff = su_y - pus_y[j];
		sInt su_pu_dist_squared = (pu_x_diff * pu_x_diff + pu_y_diff * pu_y_diff) / factor_int;

		for(unsigned int x = 0; x < pus[j].prs.size(); ++x) {
			PWID("tp PU (" + std::to_string(y) + ", " + std::to_string(x) + ")");
			sInt pr_x_diff = su_x - prs_x[j][x];
			sInt pr_y_diff = su_y - prs_y[j][x];
			sInt su_pr_dist_squared = (pr_x_diff * pr_x_diff + pr_y_diff * pr_y_diff) / factor_int;

			sInt dist_ratio = secure_dr_scale * factor_int * su_pr_dist_squared / su_pu_dist_squared;

			sInt log_dist_ratio;
			utils::secureLog10_v2(&log_dist_ratio, dist_ratio, parties, sm_params->bit_count, sm_params->factor, factor_int);
			sInt pr_est_path_loss = pu_est_path_loss + ten * pl_est_gamma * (log_dist_ratio - secure_log_dr_scale) / two / factor_int;

			sInt this_su_tp;
			if(utils::unit_type == utils::UnitType::ABS) {
				sInt pl_zero = pr_est_path_loss < one;
				this_su_tp = pl_zero.ifelse(zero, factor_int * prs_thresh[j][x] / pr_est_path_loss);
			} else if(utils::unit_type == utils::UnitType::DB) {
				this_su_tp = prs_thresh[j][x] - pr_est_path_loss;
			} else {
				std::cerr << "Unsupported unit_type" << std::endl;
				exit(1);
			}
			path_loss_su_pr[j].push_back(pr_est_path_loss + zero);

			if(y == 0 && x == 0) {
				su_tp = this_su_tp + zero;
			} else {
				sInt min_su_tp = this_su_tp < su_tp;
				su_tp = min_su_tp.ifelse(this_su_tp, su_tp);
			}
		}
	}

	// |-------------------------------|
	// |  Reveal transmit power of SU  |
	// |-------------------------------|
	PWID("send SU max_tp");
	sInt su_tp_a = INPUT(parties, 0, splitRandVal(), sm_params->bit_count);
	sInt su_tp_b = su_tp - su_tp_a;

	// reveal this output to party 0.
	parties[0].reveal(su_tp_a);
	parties[1].reveal(su_tp_b);

	int max_su_tp_int;
	float max_su_tp;

	if (parties[0].isLocalParty()) {
		max_su_tp_int = su_tp_a.getValue();
	} else if(parties[1].isLocalParty()) {
		max_su_tp_int = su_tp_b.getValue();
	}
	max_su_tp = float(max_su_tp_int) / sm_params->factor;

	// Send data to SU, and get back actual transmit power
	su_ch->send(std::array<int, 3>{su.index, max_su_tp_int, int(sm_params->factor)});

	std::array<int, 2> su_vals;
	su_ch->recv(su_vals);

	int actual_su_tp_pt = su_vals[0];
	int is_su_transmitting = su_vals[1];

	// |------------------------|
	// |  Update PR thresholds  |
	// |------------------------|
	PWID("update PR thresholds");
	sInt actual_su_tp_a = INPUT(parties, 0, actual_su_tp_pt, sm_params->bit_count);
	sInt actual_su_tp_b = INPUT(parties, 1, actual_su_tp_pt, sm_params->bit_count);
	sInt actual_su_tp = actual_su_tp_a + actual_su_tp_b;

	// is_su_transmitting
	sInt is_su_transmitting_a = INPUT(parties, 0, is_su_transmitting, sm_params->bit_count);
	sInt is_su_transmitting_b = INPUT(parties, 1, is_su_transmitting, sm_params->bit_count);
	sInt is_su_transmitting_sec = (is_su_transmitting_a + is_su_transmitting_b) > zero;

	// Calculate updates
	// updates[i].first is the index of the PU to update, and updates[i].second is the list of updates to the PR threshholds to update.
	std::vector<std::pair<sInt, std::vector<sInt> > > updates;
	for(unsigned int y = 0; y < pu_inds.size(); ++y) {
		unsigned int j = pu_inds[y];
		updates.push_back(std::make_pair(pu_table_index[j] + zero, std::vector<sInt>()));

		for(unsigned int x = 0; x < prs_thresh[j].size(); ++x) {
			PWID("update PR (" + std::to_string(y) + ", " + std::to_string(x) + ")");

			// Calculate the update to the PR thresh: u = 10 * log_10(1.0 - 10 ^ ((rp - thresh) / 10))
			sInt diff_dbm = (actual_su_tp + path_loss_su_pr[j][x]) - prs_thresh[j][x];
			sInt diff_div_ten = diff_dbm / ten;

			sInt ten_to_diff;
			PWID("start pow10");
			utils::securePow10(&ten_to_diff, diff_div_ten, parties, sm_params->bit_count, sm_params->factor, factor_int);
			PWID("end pow10");

			ten_to_diff = factor_int - ten_to_diff;

			sInt log_diff;
			PWID("start log10");
			utils::secureLog10_v2(&log_diff, ten_to_diff, parties, sm_params->bit_count, sm_params->factor, factor_int);
			PWID("end log10");

			sInt update = ten * log_diff;

			// If SU isn't actually transmitting, make the update 0.
			update = is_su_transmitting_sec.ifelse(update, zero);

			if(sm_params->no_pr_thresh_update) {
				updates[y].second.push_back(zero + zero);
			} else {
				updates[y].second.push_back(update + zero);
			}
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

	PWID("start SpectrumManager::secureTableWrite()");

	if(parties[0].isLocalParty() && secure_write_timer != nullptr) {
		secure_write_timer->start(Timer::secure_write);
	}

	int num_pu = pu_table->pus.size();
	int num_pr_per_pu = pu_table->num_pr_per_pu;

	if(sm_params->secure_write_algo == SMParams::SecureWriteAlgo::PROPOSED){
		// TODO secure Table write.
		// Each side randomly shifts the "threshold table"
		int shift = rand() % pu_table->pus.size();
	
		sInt shift_a = INPUT(parties, 0, shift, sm_params->bit_count);
		sInt shift_b = INPUT(parties, 1, shift, sm_params->bit_count);
	
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
				sInt this_update_a = INPUT(parties, 0, rand() % 100 - 50, sm_params->bit_count);
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
	} else if(sm_params->secure_write_algo == SMParams::SecureWriteAlgo::SPC) {
		// TODO
		for(int i = 0; i < num_pu; ++i) {
			std::vector<sInt> this_updates;

			for(int j = 0; j < num_pr_per_pu; ++j) {
				this_updates.push_back(INPUT(parties, 0, 0, sm_params->bit_count));
			}

			sInt this_index = INPUT(parties, 0, i, sm_params->bit_count);
			for(unsigned int j = 0; j < updates.size(); ++j) {
				sInt eq_index = (this_index >= updates[j].first) & (this_index <= updates[j].first);
				for(unsigned int y = 0; y < this_updates.size(); ++y) {
					this_updates[y] = eq_index.ifelse(updates[j].second[y], this_updates[y]);
				}
			}

			for(unsigned int j = 0; j < this_updates.size(); ++j) {
				sInt u_a = INPUT(parties, 0, rand() % 100 - 50, sm_params->bit_count);
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

void SpectrumManager::sendEncryptedData(const SUint& su, const GridTable& grid_table, const PUTable& pu_table, std::map<std::string, Timer>& en_timers) {
	en_timers["total"].start("sendEncryptedData");

	// SU
	en_timers["entities"].start("SU");
	shared_memory->set(su.getValues());
	en_timers["entities"].end("SU");

	en_timers["entities"].start("GT");
	shared_memory->set(&grid_table);
	en_timers["entities"].end("GT");

	// // Grid Table - SS
	// std::vector<int> gt_ss_size = {int(grid_table.sss.size())};
	// shared_memory->set(gt_ss_size);

	// // Assume grid_table.sss keys go from 1 to n
	// for(unsigned int i = 0; i < grid_table.sss.size(); ++i) {
	// 	auto ss_itr = grid_table.sss.find(i);
	// 	if(ss_itr == grid_table.sss.end()) {
	// 		std::cerr << "Error in SS table" << std::endl;
	// 		exit(1);
	// 	}

	// 	shared_memory->set(std::vector<int>{int(ss_itr->second.size())});

	// 	for(unsigned int j = 0; j < ss_itr->second.size(); ++j) {
	// 		en_timers["entities"].start("SS");
	// 		shared_memory->set(ss_itr->second[j].getValues());
	// 		en_timers["entities"].end("SS");
	// 	}
	// }

	// // Grid Table - PU
	// shared_memory->set(std::vector<int>{int(grid_table.pu_refs.size())});

	// // Assume grid_table.pu_ref keys go from 1 to n
	// for(unsigned int i = 0; i < grid_table.pu_refs.size(); ++i) {
	// 	auto pu_ref_itr = grid_table.pu_refs.find(i);
	// 	if(pu_ref_itr == grid_table.pu_refs.end()) {
	// 		std::cerr << "Error in PU_REF table" << std::endl;
	// 		exit(1);
	// 	}

	// 	en_timers["entities"].start("PU inds");
	// 	shared_memory->set(pu_ref_itr->second);
	// 	en_timers["entities"].end("PU inds");
	// }

	// PU Table
	shared_memory->set(std::vector<int>{int(pu_table.pus.size()), pu_table.num_pr_per_pu});

	for(unsigned int i = 0; i < pu_table.pus.size(); ++i) {
		auto pu_itr = pu_table.pus.find(i);
		if(pu_itr == pu_table.pus.end()) {
			std::cerr << "Error in PU table" << std::endl;
			exit(1);
		}

		en_timers["entities"].start("PU");
		shared_memory->set(pu_itr->second.getValues());
		en_timers["entities"].end("PU");

		if(int(pu_itr->second.prs.size()) != pu_table.num_pr_per_pu) {
			std::cerr << "Invalid num pr: " <<  pu_itr->second.prs.size() << ", " << pu_table.num_pr_per_pu << std::endl;
			exit(1);
		}

		for(unsigned int j = 0; j < pu_itr->second.prs.size(); ++j) {
			en_timers["entities"].start("PR");
			shared_memory->set(pu_itr->second.prs[j].getValues());
			en_timers["entities"].end("PR");
		}
	}

	en_timers["total"].end("sendEncryptedData");
}

void SpectrumManager::recvEncryptedData(SUint* su, GridTable* grid_table, PUTable* pu_table, std::map<std::string, Timer>& en_timers) {
	std::vector<int> vals;

	en_timers["total"].start("recvEncryptedData-recv");

	// SU
	en_timers["entities"].start("SU");
	shared_memory->get(vals);
	su->setValues(vals);
	en_timers["entities"].end("SU");

	en_timers["entities"].start("GT");
	shared_memory->get(*grid_table);
	en_timers["entities"].end("GT");

	// // Grid Table - SS
	// std::vector<int> gt_ss_c;
	// shared_memory->get(gt_ss_c);
	// int gt_ss_entries = gt_ss_c[0];

	// for(int i = 0; i < gt_ss_entries; ++i) {
	// 	std::vector<int> gt_ss_entry_c;
	// 	shared_memory->get(gt_ss_entry_c);
	// 	int gt_ss_entry_size = gt_ss_entry_c[0];

	// 	std::vector<SSint> this_entry(gt_ss_entry_size);
	// 	for(int j = 0; j < gt_ss_entry_size; ++j) {
	// 		vals.clear();
	// 		en_timers["entities"].start("SS");
	// 		shared_memory->get(vals);
	// 		this_entry[j].setValues(vals);
	// 		en_timers["entities"].end("SS");
	// 	}

	// 	grid_table->sss[i] = this_entry;
	// }

	// // Grid Table - PU
	// std::vector<int> gt_pu_ref_c;
	// shared_memory->get(gt_pu_ref_c);
	// int gt_pu_ref_entries = gt_pu_ref_c[0];

	// for(int i = 0; i < gt_pu_ref_entries; ++i) {
	// 	vals.clear();

	// 	en_timers["entities"].start("PU inds");
	// 	shared_memory->get(vals);
	// 	grid_table->pu_refs[i] = vals;
	// 	en_timers["entities"].end("PU inds");
	// }

	// PU Table
	std::vector<int> put_c;
	shared_memory->get(put_c);
	int put_entries = put_c[0];
	int put_npr = put_c[1];

	pu_table->num_pr_per_pu = put_npr;

	for(int i = 0; i < put_entries; ++i) {
		vals.clear();

		en_timers["entities"].start("PU");
		shared_memory->get(vals);
		pu_table->pus[i].setValues(vals);
		en_timers["entities"].end("PU");

		std::vector<PRint> this_prs(put_npr);
		for(int j = 0; j < put_npr; ++j) {
			vals.clear();

			en_timers["entities"].start("PR");
			shared_memory->get(vals);
			this_prs[j].setValues(vals);
			en_timers["entities"].end("PR");
		}
		pu_table->pus[i].prs = this_prs;
	}

	en_timers["total"].end("recvEncryptedData-recv");

	// Decrypt Data
	en_timers["total"].start("recvEncryptedData-decrypt");
	
	// SU
	en_timers["entities"].start("SU-decrypt");
	key_server->decrypt(su);
	en_timers["entities"].end("SU-decrypt");

	// SS
	for(auto ss_itr = grid_table->sss.begin(); ss_itr != grid_table->sss.end(); ++ss_itr) {
		for(unsigned int i = 0; i < ss_itr->second.size(); ++i) {
			en_timers["entities"].start("SS-decrypt");
			key_server->decrypt(&(ss_itr->second[i]));
			en_timers["entities"].end("SS-decrypt");
		}
	}

	// PU
	for(auto pu_itr = pu_table->pus.begin(); pu_itr != pu_table->pus.end(); ++pu_itr) {
		en_timers["entities"].start("PU-decrypt");
		key_server->decrypt(&(pu_itr->second));
		en_timers["entities"].end("PU-decrypt");

		// PR
		for(unsigned int i = 0; i < pu_itr->second.prs.size(); ++i) {
			en_timers["entities"].start("PR-decrypt");
			key_server->decrypt(&(pu_itr->second.prs[i]));
			en_timers["entities"].end("PR-decrypt");
		}
	}
	en_timers["total"].end("recvEncryptedData-decrypt");
}

void SpectrumManager::sendEncryptedPRThresholds(const PUTable& pu_table, std::map<std::string, Timer>& en_timers) {
	en_timers["total"].start("sendEncryptedPRThresholds");
	for(auto pu_itr = pu_table.pus.begin(); pu_itr != pu_table.pus.end(); ++pu_itr) {
		for(unsigned int i = 0; i < pu_itr->second.prs.size(); ++i) {
			// Encrypt just the PR thresholds
			en_timers["entities"].start("PR thresh-encrypt");
			int v = key_server->encryptPRThreshold(pu_itr->second.prs[i]);
			en_timers["entities"].end("PR thresh-encrypt");

			// Send the encrypted threshold back to the SM
			std::vector<int> thresh_val{v};

			en_timers["entities"].start("PR thresh");
			shared_memory->set(thresh_val);
			en_timers["entities"].end("PR thresh");
		}
	}
	en_timers["total"].end("sendEncryptedPRThresholds");
}

void SpectrumManager::recvEncryptedPRThresholds(PUTable* pu_table, std::map<std::string, Timer>& en_timers) {
	en_timers["total"].start("recvEncryptedPRThresholds");
	for(auto pu_itr = pu_table->pus.begin(); pu_itr != pu_table->pus.end(); ++pu_itr) {
		for(unsigned int i = 0; i < pu_itr->second.prs.size(); ++i) {
			// Recv the new PR threshold
			std::vector<int> thresh_val;

			en_timers["entities"].start("PR thresh");
			shared_memory->get(thresh_val);

			// Update the pu_table
			pu_itr->second.prs[i].threshold = thresh_val[0];
			en_timers["entities"].end("PR thresh");
		}
	}
	en_timers["total"].end("recvEncryptedPRThresholds");
}

std::vector<float> PlaintextSpectrumManager::plainTextRun(const std::vector<SU>& sus, const std::vector<PU>& input_pus, const std::vector<SS>& sss, Timer* timer, PathLossTable* path_loss_table, std::vector<std::vector<float> >* rp_at_ss_from_pu_pt, std::vector<std::vector<float> >* su_pu_pl) const {
	P("start PlaintextSpectrumManager::plainTextRun");
	std::vector<PU> pus = input_pus;

	// Compute the share of received power for each
	// Indexed by SS first then PU. received_powers[i][j] is the power received at 
	std::vector<std::vector<float> > received_powers;
	if(rp_at_ss_from_pu_pt->size() == 0) {
		for(unsigned int i = 0; i < sss.size(); ++i) {
			std::vector<float> rp_weights;
			float sum_rp_weights = 0.0;
			for(unsigned int j = 0; j < pus.size(); ++j) {
				rp_weights.push_back(1.0 / pow(sss[i].loc.dist(pus[j].loc), sm_params->rp_alpha_f));
				sum_rp_weights += rp_weights[j];
			}

			received_powers.push_back(std::vector<float>());

			for(unsigned int j = 0; j < pus.size(); ++j) {
				if(utils::unit_type == utils::UnitType::ABS) {
					received_powers[i].push_back(rp_weights[j] * sss[i].received_power / sum_rp_weights);
				} else if(utils::unit_type == utils::UnitType::DB) {
					received_powers[i].push_back(utils::todBm(rp_weights[j] / sum_rp_weights) + sss[i].received_power);
				} else {
					std::cerr << "Unsupported unit_type" << std::endl;
					exit(1);
				}
			}
		}

		rp_at_ss_from_pu_pt->clear(); // rp_at_ss_from_pu_pt[j][i] is the estimated received power at SS i from PU j
		for(unsigned int j = 0; j < pus.size(); ++j) {
			rp_at_ss_from_pu_pt->push_back(std::vector<float>());
			for(unsigned int i = 0; i < sss.size(); ++i) {
				(*rp_at_ss_from_pu_pt)[j].push_back(received_powers[i][j]);
			}
		}
	} else {
		std::cout << "Using GT rp_at_ss_from_pu_pt" << std::endl;
		if(rp_at_ss_from_pu_pt->size() != pus.size()) {
			std::cerr << "GT rps don't match the number of pus" << std::endl;
			exit(1);
		}

		for(unsigned int j = 0; j < rp_at_ss_from_pu_pt->size(); ++j) {
			if((*rp_at_ss_from_pu_pt)[j].size() != sss.size()) {
				std::cerr << "GT rps don't match the number of sss" << std::endl;
				exit(1);
			}
		}

		for(unsigned int i = 0; i < sss.size(); ++i) {
			received_powers.push_back(std::vector<float>());
			for(unsigned int j = 0; j < pus.size(); ++j) {
				received_powers[i].push_back((*rp_at_ss_from_pu_pt)[j][i]);
			}
		}
	}

	std::map<int, std::vector<PU*> > pu_groups;
	std::map<int, std::vector<const SS*> > ss_groups;
	std::map<int, std::vector<std::vector<float> > > rp_powers_groups;
	if(sm_params->use_grid) {
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

				rp_powers_groups[itr->first].push_back(std::vector<float>());

				auto pu_itr = pu_int_groups.find(itr->first);
				if(pu_itr == pu_int_groups.end()) {
					std::cerr << "Error with rp and groups in PT" << std::endl;
					exit(1);
				}

				for(unsigned int j = 0; j < pu_itr->second.size(); ++j) {
					rp_powers_groups[itr->first][i].push_back(received_powers[itr->second[i]][pu_itr->second[j]]);
				}
			}
		}
	} else {
		std::cerr << "Not using grid is not supported" << std::endl;
		exit(1);
	}

	if(su_pu_pl == nullptr) {
		std::cerr << "No su_pu_pl supplied" << std::endl;
		exit(1);
	}

	if(su_pu_pl->size() == 0) {
		*su_pu_pl = std::vector<std::vector<float> >(sus.size());
	}

	std::vector<PU*> sel_pus;
	std::vector<const SS*> sel_sss;
	std::vector<std::vector<float> > this_rps;
	if(!sm_params->use_grid) {
		for(unsigned int j = 0; j < pus.size(); ++j) {
			sel_pus.push_back(&pus[j]);
		}

		for(unsigned int i = 0; i < sss.size(); ++i) {
			sel_sss.push_back(&sss[i]);
		}
	}
	std::vector<float> all_vals;
	for(unsigned int i = 0; i < sus.size(); ++i) {
		std::cout << "Starting PT request for SU " << i << std::endl;
		if(sm_params->use_grid) {
			sel_pus.clear();
			sel_sss.clear();

			int group_id = int(sus[i].loc.x / sm_params->grid_delta_x) + int(sus[i].loc.y / sm_params->grid_delta_y) * sm_params->grid_num_x;

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

			auto rp_itr = rp_powers_groups.find(group_id);
			if(rp_itr == rp_powers_groups.end()) {
				std::cerr << "Invalid group_id: " << group_id << std::endl;
				exit(1);
			}

			this_rps = rp_itr->second;
		}
		float v = plainTextRadar(sus[i], sel_pus, sel_sss, this_rps, path_loss_table, &((*su_pu_pl)[i]));
		all_vals.push_back(v);
		std::cout << "Finished PT request for SU " << i << std::endl;
	}

	P("end PlaintextSpectrumManager::plainTextRun");
	return all_vals;
}

void PlaintextSpectrumManager::plainTextGrid(const std::vector<PU>& pus, const std::vector<SS>& sss,
		std::map<int, std::vector<int> >* pu_int_groups, std::map<int, std::vector<int> >* ss_int_groups) const {
	struct comp {
		bool operator()(const std::pair<float, int>& v1, const std::pair<float, int>& v2) const {
			return v1.first < v2.first;
		}
	};

	// PU Groups
	for(int x = 0; x < sm_params->grid_num_x; ++x) {
		for(int y = 0; y < sm_params->grid_num_y; ++y) {
			std::vector<std::pair<float, int> > pu_vals;

			for(unsigned int i = 0; i < pus.size(); ++i) {
				float this_dist = pus[i].loc.dist(Location((x + 0.5) * sm_params->grid_delta_x, (y + 0.5) * sm_params->grid_delta_y));

				if((int) pu_vals.size() < sm_params->grid_min_num_pu || this_dist < pu_vals.front().first) {
					if((int) pu_vals.size() >= sm_params->grid_min_num_pu) {
						std::pop_heap(pu_vals.begin(), pu_vals.end(), comp());
						pu_vals.pop_back();
					}

					pu_vals.push_back(std::make_pair(this_dist, i));
					std::push_heap(pu_vals.begin(), pu_vals.end(), comp());
				}
			}

			int group_id = x + y * sm_params->grid_num_x;
			for(unsigned int i = 0; i < pu_vals.size(); ++i) {
				(*pu_int_groups)[group_id].push_back(pu_vals[i].second);
			}
		}
	}

	// SS Groups
	for(int x = 0; x < sm_params->grid_num_x; ++x) {
		for(int y = 0; y < sm_params->grid_num_y; ++y) {
			std::vector<std::pair<float, int>> ss_vals;

			for(unsigned int i = 0; i < sss.size(); ++i) {
				float this_dist = sss[i].loc.dist(Location((x + 0.5) * sm_params->grid_delta_x, (y + 0.5) * sm_params->grid_delta_y));

				if((int) ss_vals.size() < sm_params->grid_min_num_ss || this_dist < ss_vals.front().first) {
					if((int) ss_vals.size() >= sm_params->grid_min_num_ss) {
						std::pop_heap(ss_vals.begin(), ss_vals.end(), comp());
						ss_vals.pop_back();
					}

					ss_vals.push_back(std::make_pair(this_dist, i));
					std::push_heap(ss_vals.begin(), ss_vals.end(), comp());
				}
			}

			int group_id = x + y * sm_params->grid_num_x;
			for(unsigned int i = 0; i < ss_vals.size(); ++i) {
				(*ss_int_groups)[group_id].push_back(ss_vals[i].second);
			}
		}
	}
}

float PlaintextSpectrumManager::plainTextRadar(
		const SU& su, std::vector<PU*>& pus, const std::vector<const SS*>& sss, const std::vector<std::vector<float> >& received_powers, PathLossTable* path_loss_table, std::vector<float>* this_su_pu_pl) const {
	// PU selection.
	std::vector<int> pu_inds;
	unsigned int k_pu = numSelect(sm_params->num_pu_selection, pus.size());

	// SS selection.
	std::vector<int> sss_inds;
	std::vector<float> sss_dists;
	unsigned int k_ss = numSelect(sm_params->num_ss_selection, sss.size());
	if(sm_params->selection_algo == SMParams::SelectionAlgo::NONE) {
		for(unsigned int i = 0; i < k_pu; ++i) {
			pu_inds.push_back(i);
		}

		for(unsigned int i = 0; i < k_ss; ++i) {
			sss_inds.push_back(i);
			sss_dists.push_back(su.loc.dist(sss[i]->loc));
		}
	} else if(sm_params->selection_algo == SMParams::SelectionAlgo::SORT) {
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
	} else if(sm_params->selection_algo == SMParams::SelectionAlgo::RANDOM) {
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

		float w = 1.0 / pow(d, sm_params->pl_alpha);
		weights.push_back(w);
		tmp_sum_weight += w;
	}

	bool read_su_pu_pl = (this_su_pu_pl->size() != 0);
	if(!read_su_pu_pl) {
		*this_su_pu_pl = std::vector<float>(pus.size(), 0.0);
	}
	
	// Compute SU transmit power
	// tp = thresh * sum(w(SS)) / sum(r(SS) / t(PU) * w(SS))
	float max_transmit_power = std::numeric_limits<float>::infinity();
	std::vector<std::vector<float> > estimated_path_loss; // estimated_path_loss[i][j] is the path loss between the su and PU i's PR j.
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
		
		float this_pu_path_loss = sum_weighted_ratio / sum_weight;

		if(read_su_pu_pl) {
			this_pu_path_loss = (*this_su_pu_pl)[j];
		} else {
			(*this_su_pu_pl)[j] = this_pu_path_loss;
		}


		estimated_path_loss.push_back(std::vector<float>());
		for(unsigned int x = 0; x < pus[j]->prs.size(); ++x) {
			float this_pr_path_loss = this_pu_path_loss + 10.0 * sm_params->pl_est_gamma * log10(su.loc.dist(pus[j]->prs[x].loc) / su.loc.dist(pus[j]->loc));

			path_loss_table->addPlaintextPathLoss(su.index, j, x, this_pr_path_loss);
			estimated_path_loss[y].push_back(this_pr_path_loss);
			
			float this_transmit_power;
			if(utils::unit_type == utils::UnitType::ABS) {
				std::cerr << "Not implemented yet" << std::endl;
				exit(1);

				this_transmit_power = pus[j]->prs[x].threshold / this_pr_path_loss;
			} else if(utils::unit_type == utils::UnitType::DB) {
				this_transmit_power = pus[j]->prs[x].threshold - this_pr_path_loss;
			} else {
				std::cerr << "Unsupported unit_type" << std::endl;
				exit(1);
			}

			if(this_transmit_power < max_transmit_power) {
				max_transmit_power = this_transmit_power;
			}
		}
	}

	if(!(sm_params->no_pr_thresh_update)) {
		float actual_su_tp = max_transmit_power + su.less_max_tp;
		bool is_su_transmitting = actual_su_tp > su.min_tp;

		if(is_su_transmitting) {
			for(unsigned int y = 0; y < pu_inds.size(); ++y) {
				unsigned int j = pu_inds[y];
				if(j != y) {
					std::cerr << "False assumption made" << std::endl;
					exit(1);
				}
				for(unsigned int x = 0; x < pus[j]->prs.size(); ++x) {
					float rp = actual_su_tp + estimated_path_loss[y][x];

					float new_value = utils::todBm(utils::fromdBm(pus[j]->prs[x].threshold) - utils::fromdBm(rp));
					
					pus[j]->prs[x].threshold = new_value;
				}
			}
		}
	}

	return max_transmit_power;
}

void SMParams::initGroupDeviations() {

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

const std::vector<std::pair<int, int> >& SMParams::getDeviations(int iter) const {
	if(iter < 0 || (unsigned int)iter >= all_deviations.size()) {
		std::cerr << "Invalid iter: " << iter << std::endl;
		exit(1);
	}
	return all_deviations[iter];
}
