
#include "args.h"
#include "generator.h"
#include "location.h"
#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_manager.h"
#include "spectrum_sensor.h"
#include "split.h"
#include "timer.h"
#include "utils.h"

#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <math.h>
#include <sstream>

int main(int argc, char const *argv[]) {
	Args args(argc, argv);

	auto rs = time(nullptr);
	if(args.rand_seed != 0) {
		rs = args.rand_seed;
	}
	srand(rs);
	if(args.brief_out) {
		std::cout << "rand_seed|int|" << rs << std::endl;
	} else {
		std::cout << "Using random seed: " << rs << std::endl;
	}
	
	utils::setUnitType(args.unit_type);

	PropagationModel* pm;
	if(args.propagation_model == "log_distance") {
		pm = new LogDistancePM(args.ld_path_loss0, args.ld_dist0, args.ld_gamma);
	} else if (args.propagation_model == "longley_rice") {
		pm = new LongleyRicePM(args.ref_lat, args.ref_long,
								args.splat_dir, args.sdf_dir, args.return_dir);
	} else {
		std::cerr << "Unknown propagation model: " << args.propagation_model << std::endl;
		exit(0);
	}

	Generator gen(args.location_range, pm);

	std::vector<PU> pus;
	std::vector<SS> sss;
	std::vector<SU> sus;
	gen.generateEntities(args.num_pu, args.num_ss, args.num_su, &pus, &sss, &sus);
	
	float factor = 1 << args.num_float_bits;
	SM sm(factor, args.num_ss_selection, args.num_pu_selection, args.ss_receive_power_alpha, args.path_loss_alpha, args.s2_pc_bit_count, args.brief_out);
	sm.setAlgoOrder(args.algo_order);

	std::vector<float> secure_vs_a;
	std::vector<float> secure_vs_b;
	Timer t1;
	if(!args.skip_s2pc) {	
		// "Split values"
		std::vector<PUint> pus_int0;
		std::vector<PUint> pus_int1;
		for(unsigned int i = 0; i < pus.size(); ++i) {
			auto pu_split = splitPrimaryUserInt(pus[i], factor);
			pus_int0.push_back(pu_split.first);
			pus_int1.push_back(pu_split.second);
		}

		std::vector<SSint> sss_int0;
		std::vector<SSint> sss_int1;
		for(unsigned int i = 0; i < sss.size(); ++i) {
			auto ss_split = splitSpectrumSensorInt(sss[i], factor);
			sss_int0.push_back(ss_split.first);
			sss_int1.push_back(ss_split.second);
		}

		std::vector<SUint> sus_int0;
		std::vector<SUint> sus_int1;
		for(unsigned int i = 0; i < sus.size(); ++i) {
			auto su_split = splitSecondaryUserInt(sus[i], factor);
			sus_int0.push_back(su_split.first);
			sus_int1.push_back(su_split.second);
		}

		NullTimer null_timer;

		std::thread thrd0([&]() {
			int party_id = 0;
			secure_vs_a = sm.run(party_id, pus_int0, sss_int0, sus_int0, &t1);
		});

		std::thread thrd1([&]() {
			int party_id = 1;
			secure_vs_b = sm.run(party_id, pus_int1, sss_int1, sus_int1, &null_timer);
		});

		thrd0.join();
		thrd1.join();
	}

	std::vector<float> secure_vs;
	std::vector<float> plaintext_vs;
	std::vector<float> ground_truth_vs;
	for(unsigned int i = 0; i < sus.size(); ++i) {
		if(!args.skip_s2pc) {
			float secure_v = secure_vs_a[i] + secure_vs_b[i];
			secure_vs.push_back(secure_v);
		}

		float plaintext_v = sm.plainTextRadar(sus[i], pus, sss);
		float ground_truth_v = gen.computeGroundTruth(sus[i], pus);

		plaintext_vs.push_back(plaintext_v);
		ground_truth_vs.push_back(ground_truth_v);
	}

	if(args.brief_out) {
		if(!args.skip_s2pc) {
			std::cout << "su_transmit_power(secure,plain,ground)|list(float,float,float)|";
		} else {
			std::cout << "su_transmit_power(plain,ground)|list(float,float)|";
		}
		for(unsigned int i = 0; i < plaintext_vs.size(); ++i) {
			if(!args.skip_s2pc) {
				std::cout << secure_vs[i] << ":";
			}
			std::cout << plaintext_vs[i] << ":" << ground_truth_vs[i];
			if(i < plaintext_vs.size() - 1) {
				std::cout << ",";
			}
		}
		std::cout << std::endl;

		if(!args.skip_s2pc) {
			std::cout << "preprocess_time|float|" << t1.getAverageDuration(Timer::secure_preprocessing) << std::endl;
			std::cout << "time_per_request|float|" << t1.getAverageDuration(Timer::secure_su_request) << std::endl;
		}
	} else {
		for(unsigned int i = 0; i < plaintext_vs.size(); ++i) {
			std::cout << "--------------------" << std::endl;
			if(!args.skip_s2pc) {
				std::cout << "Secure:     " << secure_vs[i] << std::endl;
			}
			std::cout << "Plain:      " << plaintext_vs[i] << std::endl;
			std::cout << "Ground:     " << ground_truth_vs[i] << std::endl;
			if(!args.skip_s2pc) {
				std::cout << "Diff PvS:   " << fabs(secure_vs[i] - plaintext_vs[i]) << std::endl;
			}
			std::cout << "Diff GvP:   " << fabs(plaintext_vs[i] - ground_truth_vs[i]) << std::endl;
			if(!args.skip_s2pc) {
				std::cout << "% Diff PvS: " << fabs(secure_vs[i] - plaintext_vs[i]) / plaintext_vs[i] * 100 << "%" << std::endl;
			}
			std::cout << "% Diff GvP: " << fabs(plaintext_vs[i] - ground_truth_vs[i]) / ground_truth_vs[i] * 100 << "%" << std::endl;
			std::cout << "--------------------" << std::endl;
		}

		if(!args.skip_s2pc) {
			std::cout << "Average duration of " << Timer::secure_preprocessing << ": " << t1.getAverageDuration(Timer::secure_preprocessing) << std::endl;
			std::cout << "Average duration of " << Timer::secure_su_request << ": " << t1.getAverageDuration(Timer::secure_su_request) << std::endl;
		}
	}

	return 0;
}

