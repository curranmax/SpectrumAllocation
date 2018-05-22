
#include "spectrum_manager.h"

#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_sensor.h"
#include "split.h"
#include "utils.h"

#include "cryptoTools/Network/IOService.h"
#include "cryptoTools/Network/Session.h"
#include "cryptoTools/Crypto/PRNG.h"

#include "ivory/Runtime/ShGc/ShGcRuntime.h"
#include "ivory/Runtime/sInt.h"
#include "ivory/Runtime/Party.h"

#include <math.h>

using namespace osuCrypto;

#define INPUT(parties, party_id, val, bit_count) parties[party_id].isLocalParty() ? parties[party_id].input<sInt>(val, bit_count) : parties[party_id].input<sInt>(bit_count)

// Options: order(split rp at ss then estimate pl at su, estimate rp at su then split rp and calc pl), pl type(ratio, db)
std::vector<float> SpectrumManager::run(int party_id,
										const std::vector<PUint>& pus,
										const std::vector<SSint>& sss,
										const std::vector<SUint>& sus,
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

	// Run the preprocessing calculations:
	std::vector<std::vector<int> > plaintext_sss_rp_from_pu;

	if(parties[0].isLocalParty() && !brief_out) {
		std::cout << "Starting preprocessing" << std::endl;
	}
	timer->start(Timer::secure_preprocessing);
	secureRadarPreprocess(parties, pus, sss, &plaintext_sss_rp_from_pu);
	timer->end(Timer::secure_preprocessing);
	if(parties[0].isLocalParty()  && !brief_out) {
		std::cout << "Finished preprocessing" << std::endl;
	}

	std::vector<float> su_rps;
	for(unsigned int i = 0; i < sus.size(); ++i) {
		if(parties[0].isLocalParty()  && !brief_out) {
			std::cout << "Starting request for SU " << i + 1 << std::endl;
		}
		
		timer->start(Timer::secure_su_request);
		
		float v = secureRadar(parties, sus[i], pus, sss, plaintext_sss_rp_from_pu);

		timer->end(Timer::secure_su_request);

		su_rps.push_back(v);

		if(parties[0].isLocalParty()  && !brief_out) {
			std::cout << "Finished request for SU " << i + 1 << std::endl;
		}
	}
	return su_rps;
}

void SpectrumManager::secureRadarPreprocess(
		std::array<Party, 2> parties, const std::vector<PUint>& pus, const std::vector<SSint>& sss,
		std::vector<std::vector<int> >* plaintext_sss_rp_from_pu) const {
	sInt zero = INPUT(parties, 0, 0, bit_count);
	sInt one  = INPUT(parties, 0, 1, bit_count);
	sInt two  = INPUT(parties, 0, 2, bit_count);

	sInt factor_int = INPUT(parties, 0, int(factor), bit_count);
	sInt large = INPUT(parties, 0, factor * factor * 10, bit_count); // TODO

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
	for(unsigned int i = 0; i < sss.size(); ++i) {
		sInt ss_x_a = INPUT(parties, 0, sss[i].loc.x, bit_count);
		sInt ss_x_b = INPUT(parties, 1, sss[i].loc.x, bit_count);

		sss_x.push_back(ss_x_a + ss_x_b);

		sInt ss_y_a = INPUT(parties, 0, sss[i].loc.y, bit_count);
		sInt ss_y_b = INPUT(parties, 1, sss[i].loc.y, bit_count);
		
		sss_y.push_back(ss_y_a + ss_y_b);

		sInt ss_rp_a = INPUT(parties, 0, sss[i].received_power, bit_count);
		sInt ss_rp_b = INPUT(parties, 1, sss[i].received_power, bit_count);
		
		sss_rp.push_back(ss_rp_a + ss_rp_b);
	}

	// |----------------------------------|
	// |  SS received power from each PU  |
	// |----------------------------------|
	// TODO precompute this since none of it is related to the SU.
	std::vector<std::vector<sInt> > sss_rp_from_pu_a;
	std::vector<std::vector<sInt> > sss_rp_from_pu_b;
	for(unsigned int i = 0; i < sss.size(); ++i) {
		std::vector<sInt> rp_weights;
		sInt sum_rp_weights = INPUT(parties, 0, 0, bit_count);
		for(unsigned int j = 0; j < pus.size(); ++j) {
			sInt dist_to_rp_alpha = INPUT(parties, 0, 1, bit_count);

			if(rp_alpha == 1) {
				sInt dist;
				utils::dist(&dist, sss_x[i], sss_y[i], pus_x[j], pus_y[j], zero, two, 10);

				dist_to_rp_alpha = dist_to_rp_alpha * dist;
			} else if(rp_alpha == 2) {
				sInt dist_squared = ((sss_x[i] - pus_x[j]) * (sss_x[i] - pus_x[j]) + (sss_y[i] - pus_y[j]) * (sss_y[i] - pus_y[j])) / factor_int;

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

			rp_weights.push_back(is_dist_small.ifelse(large, factor_int * factor_int / dist_to_rp_alpha));
			sum_rp_weights = sum_rp_weights + rp_weights[j];
		}

		sss_rp_from_pu_a.push_back(std::vector<sInt>());
		sss_rp_from_pu_b.push_back(std::vector<sInt>());
		for(unsigned int j = 0; j < pus.size(); ++j) {
			sInt sss_rp_from_pu = rp_weights[j] * sss_rp[i] / sum_rp_weights;
			sss_rp_from_pu_a[i].push_back(INPUT(parties, 0, splitRandVal(), bit_count));
			sss_rp_from_pu_b[i].push_back(sss_rp_from_pu - sss_rp_from_pu_a[i][j]);
		}
	}

	// Add computed values to plaintext_sss_rp_from_pu.
	for(unsigned int i = 0; i < sss.size(); ++i) {
		plaintext_sss_rp_from_pu->push_back(std::vector<int>());
		for(unsigned int j = 0; j < pus.size(); ++j) {
			parties[0].reveal(sss_rp_from_pu_a[i][j]);
			parties[1].reveal(sss_rp_from_pu_b[i][j]);
			if(parties[0].isLocalParty()) {
				(*plaintext_sss_rp_from_pu)[i].push_back(sss_rp_from_pu_a[i][j].getValue());
			} else if(parties[1].isLocalParty()) {
				(*plaintext_sss_rp_from_pu)[i].push_back(sss_rp_from_pu_b[i][j].getValue());
			}
		}
	}
}

float SpectrumManager::secureRadar(
		std::array<Party, 2> parties, const SUint& su, const std::vector<PUint>& pus, const std::vector<SSint>& sss,
		const std::vector<std::vector<int> >& plaintext_sss_rp_from_pu) const {
	// Constants
	sInt zero = INPUT(parties, 0, 0, bit_count);
	sInt one  = INPUT(parties, 0, 1, bit_count);
	sInt two  = INPUT(parties, 0, 2, bit_count);

	sInt factor_int = INPUT(parties, 0, int(factor), bit_count);
	sInt large = INPUT(parties, 0, factor * factor * 10, bit_count); // TODO

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
	std::vector<sInt> pus_tp;
	std::vector<sInt> prs_thresh;
	for(unsigned int i = 0; i < pus.size(); ++i) {
		sInt pu_tp_a = INPUT(parties, 0, pus[i].transmit_power, bit_count);
		sInt pu_tp_b = INPUT(parties, 1, pus[i].transmit_power, bit_count);
	
		pus_tp.push_back(pu_tp_a + pu_tp_b);
	
		sInt pr_thresh_a = INPUT(parties, 0, pus[i].prs[0].threshold, bit_count);
		sInt pr_thresh_b = INPUT(parties, 1, pus[i].prs[0].threshold, bit_count);
	
		prs_thresh.push_back(pr_thresh_a + pr_thresh_b);
	}

	// SS
	std::vector<sInt> sss_x;
	std::vector<sInt> sss_y;
	for(unsigned int i = 0; i < sss.size(); ++i) {
		sInt ss_x_a = INPUT(parties, 0, sss[i].loc.x, bit_count);
		sInt ss_x_b = INPUT(parties, 1, sss[i].loc.x, bit_count);

		sss_x.push_back(ss_x_a + ss_x_b);

		sInt ss_y_a = INPUT(parties, 0, sss[i].loc.y, bit_count);
		sInt ss_y_b = INPUT(parties, 1, sss[i].loc.y, bit_count);
		
		sss_y.push_back(ss_y_a + ss_y_b);
	}

	// Get the preprocessing values.
	std::vector<std::vector<sInt> > sss_rp_from_pu;
	for(unsigned int i = 0; i < sss.size(); ++i) {
		sss_rp_from_pu.push_back(std::vector<sInt>());
		for(unsigned int j = 0; j < pus.size(); ++j) {
			sInt sss_rp_from_pu_a = INPUT(parties, 0, plaintext_sss_rp_from_pu[i][j], bit_count);
			sInt sss_rp_from_pu_b = INPUT(parties, 1, plaintext_sss_rp_from_pu[i][j], bit_count);

			sss_rp_from_pu[i].push_back(sss_rp_from_pu_a + sss_rp_from_pu_b);
		}
	}

	// |----------------------------|
	// |  Select the k-closest SSs  |
	// |----------------------------|
	// Calculate the k-closest SS
	std::vector<sInt> secret_sss_inds;
	std::vector<sInt> secret_sss_dists;
	// std::vector<sInt> other_dists;
	unsigned int k = sss.size();
	if(this->num_ss_selection > 0) {
		k = this->num_ss_selection;
	}	
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

		if(i < k) {
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
	std::vector<int> sss_inds;
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

	// |----------------------|
	// |  Compute PL weights  |
	// |----------------------|
	// w = (1 / d(SU, SS))
	std::vector<sInt> sss_w;
	for(unsigned int x = 0; x < sss_inds.size(); ++x) {
		sInt dist_to_pl_alpha = secret_sss_dists[x] + zero;

		// Check if dist is small
		sInt is_dist_small = dist_to_pl_alpha <= one;

		sss_w.push_back(is_dist_small.ifelse(large, factor_int * factor_int / dist_to_pl_alpha));
	}

	// |--------------------------------|
	// |  Compute transmit power of SU  |
	// |--------------------------------|
	// su_rp = pr_thresh * sum(w(SS)) / (sum(r(SS)/t(PU) * w(SS)))
	sInt su_tp = large + zero;

	for(unsigned int j = 0; j < pus.size(); ++j) {
		sInt sum_weight = INPUT(parties, 0, 0, bit_count);
		sInt sum_weighted_ratio = INPUT(parties, 0, 0, bit_count);
		for(unsigned int x = 0; x < sss_inds.size(); ++x) {
			unsigned int i = sss_inds[x];

			sum_weight = sum_weight + sss_w[x];
			if(pl_type == SpectrumManager::PathLossType::RATIO) {
				sum_weighted_ratio = sum_weighted_ratio + sss_w[x] * sss_rp_from_pu[i][j] / pus_tp[j];
			} else if(pl_type == SpectrumManager::PathLossType::DB) {
				std::cerr << "Not implemented yet" << std::endl;
				exit(0);
			}
		}
	
		sInt this_su_tp;
		if(pl_type == SpectrumManager::PathLossType::RATIO) {
			this_su_tp = prs_thresh[j] * sum_weight / sum_weighted_ratio;
		} else if(pl_type == SpectrumManager::PathLossType::DB) {
			std::cerr << "Not implemented yet" << std::endl;
			exit(0);
		}

		sInt min_su_tp = this_su_tp < su_tp;
		su_tp = min_su_tp.ifelse(this_su_tp, su_tp);
	}

	sInt su_tp_a = INPUT(parties, 0, splitRandVal(), bit_count);
	sInt su_tp_b = su_tp - su_tp_a;

	// reveal this output to party 0.
	parties[0].reveal(su_tp_a);
	parties[1].reveal(su_tp_b);

	// operations can get queued up in the background. Eventually this call should not
	// be required but in the mean time, if one party does not call getValue(), then
	// processesQueue() should be called.
	// parties[1].getRuntime().processesQueue();

	if (parties[0].isLocalParty()) {
		return float(su_tp_a.getValue()) / factor;
	} else if(parties[1].isLocalParty()) {
		return float(su_tp_b.getValue()) / factor;;
	}
	return 0.0;
}

float SpectrumManager::plainTextRadar(const SU& su,
										const std::vector<PU>& pus,
										const std::vector<SS>& sss) const {
	// SS selection.
	std::vector<int> sss_inds;
	std::vector<float> sss_dists;
	unsigned int k = sss.size();
	if(this->num_ss_selection > 0) {
		k = this->num_ss_selection;
	}	
	for(unsigned int i = 0; i < sss.size(); ++i) {
		float dist = su.loc.dist(sss[i].loc);
		int ind = i;

		if(i < k) {
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

	// Compute weight
	// w = 1 / d(SU, SS)
	std::vector<float> weights;
	for(unsigned int x = 0; x < sss_inds.size(); ++x) {
		float d = sss_dists[x];
		if (d < 0.0001) {
			d = 0.0001; 
		}
		weights.push_back(1.0 / pow(d, pl_alpha));
	}

	float transmit_power = std::numeric_limits<float>::infinity();
	if(order == SpectrumManager::AlgoOrder::SPLIT_THEN_IDW) {
		// Compute the share of received power for each
		// Indexed by SS first then PU. received_powers[i][j] is the power received at 
		std::vector<std::vector<float> > received_powers;
		for(unsigned int i = 0; i < sss.size(); ++i) {
			std::vector<float> rp_weights;
			float sum_rp_weights = 0.0;
			for(unsigned int j = 0; j < pus.size(); ++j) {
				rp_weights.push_back(1.0 / pow(sss[i].loc.dist(pus[j].loc), rp_alpha));
				sum_rp_weights += rp_weights[j];
			}

			received_powers.push_back(std::vector<float>());
			for(unsigned int j = 0; j < pus.size(); ++j) {
				received_powers[i].push_back(rp_weights[j] * sss[i].received_power / sum_rp_weights);
			}
		}
		
		// Compute SU transmit power
		// tp = thresh * sum(w(SS)) / sum(r(SS) / t(PU) * w(SS))
		for(unsigned int j = 0; j < pus.size(); ++j) {
			float sum_weight = 0.0;
			float sum_weighted_ratio = 0.0;
			float sum_weighted_ratio2 = 0.0;
			for(unsigned int x = 0; x < sss_inds.size(); ++x) {
				int i = sss_inds[x];
				sum_weight = sum_weight + weights[x];
				if(pl_type == SpectrumManager::PathLossType::RATIO) {
					sum_weighted_ratio = sum_weighted_ratio + weights[x] * received_powers[i][j] / pus[j].transmit_power;
				} else if(pl_type == SpectrumManager::PathLossType::DB) {
					sum_weighted_ratio2 = sum_weighted_ratio + weights[x] * received_powers[i][j] / pus[j].transmit_power;
					sum_weighted_ratio = sum_weighted_ratio + weights[x] * 10.0 * (log10(received_powers[i][j]) - log10(pus[j].transmit_power));
				}
			}
			


			float this_transmit_power;
			if(pl_type == SpectrumManager::PathLossType::RATIO) {
				this_transmit_power = pus[j].prs[0].threshold * sum_weight / sum_weighted_ratio;
			} else if(pl_type == SpectrumManager::PathLossType::DB) {
				this_transmit_power = utils::fromdBm(utils::todBm(pus[j].prs[0].threshold) - sum_weighted_ratio / sum_weight);
			}

			if(this_transmit_power < transmit_power) {
				transmit_power = this_transmit_power;
			}
		}
	} else if(order == SpectrumManager::AlgoOrder::IDW_THEN_SPLIT) {
		float sum_weight = 0.0;
		float sum_weighted_ratio = 0.0;
		for(unsigned int x = 0; x < sss_inds.size(); ++x) {
			int i = sss_inds[x];
			sum_weight = sum_weight + weights[x];
			sum_weighted_ratio = sum_weighted_ratio + sss[i].received_power * weights[x];
		}
		
		float received_power = sum_weighted_ratio / sum_weight;

		sum_weight = 0.0;
		std::vector<float> pu_weights;
		for(unsigned int j = 0; j < pus.size(); ++j) {
			float d = su.loc.dist(pus[j].loc);
			if(d < 0.0001) {
				d = 0.0001;
			}
			sum_weight = sum_weight + 1 / pow(d, rp_alpha);
			pu_weights.push_back(1 / pow(d, rp_alpha));
		}

		for(unsigned int j = 0; j < pus.size(); ++j) {
			float this_path_loss = pu_weights[j] * received_power / sum_weight / pus[j].transmit_power;
			float this_transmit_power = pus[j].prs[0].threshold / this_path_loss;

			if(this_transmit_power < transmit_power) {
				transmit_power = this_transmit_power;
			}
		}
	}

	return transmit_power;
}
