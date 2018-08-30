
#include "args.h"
#include "generator.h"
#include "key_server.h"
#include "location.h"
#include "path_loss_table.h"
#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_manager.h"
#include "spectrum_sensor.h"
#include "split.h"
#include "su_server.h"
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
#include <algorithm>

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

	// Set up ground truth propagation model
	PropagationModel* pm;
	if(args.propagation_model == "log_distance") {
		pm = new LogDistancePM(args.ld_path_loss0, args.ld_dist0, args.ld_gamma);
	} else if (args.propagation_model == "longley_rice") {
		pm = new LongleyRicePM(args.splat_cmd, args.ref_lat, args.ref_long,
								args.splat_dir, args.sdf_dir, args.return_dir,
								args.use_itwom_pl);
	} else if(args.propagation_model == "single_lr") {
		pm = new TransmitterOnlySplatPM(args.splat_cmd, args.ref_lat, args.ref_long,
										args.splat_dir, args.sdf_dir, args.return_dir);
	} else {
		std::cerr << "Unknown propagation model: " << args.propagation_model << std::endl;
		exit(0);
	}

	// Generate enetities
	Generator gen(args.location_range, pm);

	std::vector<PU> pus;
	std::vector<SS> sss;
	std::vector<SU> sus;
	gen.generateEntities(args.num_pu, args.num_ss, args.num_su, args.num_pr_per_pu, args.pr_range, args.out_filename, &pus, &sss, &sus);
	
	// Set up Spectrum Manager params
	float factor = 1 << args.num_float_bits;
	SMParams sm_params(factor, args.s2_pc_bit_count, args.num_pu_selection, args.num_ss_selection, args.ss_receive_power_alpha, args.path_loss_alpha, args.brief_out);
	sm_params.setSelectionAlgo(args.selection_algo);

	sm_params.setSecureWriteAlgo(args.secure_write_algo);
	sm_params.no_pr_thresh_update = args.no_pr_thresh_update;
	sm_params.pl_est_gamma = args.pl_est_gamma;

	if(args.grid_num_x > 0 && args.grid_num_y > 0) {
		int k_pu = args.num_pu_selection;
		int k_ss = args.num_ss_selection;

		if(k_pu <= 0) {
			k_pu = args.num_pu;
		}

		if(k_ss <= 0) {
			k_ss = args.num_ss;
		}

		sm_params.setGridParams(
				k_pu, k_ss,
				args.grid_num_x, args.grid_num_y,
				args.location_range / args.grid_num_x, args.location_range / args.grid_num_y);
	}

	PtSM pt_sm(&sm_params);

	std::vector<float> secure_vs_a;
	std::vector<float> secure_vs_b;
	Timer t1, secure_write_timer;
	if(!args.skip_s2pc) {
		// "Split values"
		// PU
		std::vector<PUint> pus_int0;
		std::vector<PUint> pus_int1;
		for(unsigned int i = 0; i < pus.size(); ++i) {
			auto pu_split = splitPrimaryUserInt(pus[i], factor);
			pus_int0.push_back(pu_split.first);
			pus_int1.push_back(pu_split.second);
		}

		// Set all PUint's indexes
		for(unsigned int i = 0; i < pus.size(); ++i) {
			pus_int0[i].setInd(i);
			pus_int1[i].setInd(i);
		}

		// SS
		std::vector<SSint> sss_int0;
		std::vector<SSint> sss_int1;
		for(unsigned int i = 0; i < sss.size(); ++i) {
			auto ss_split = splitSpectrumSensorInt(sss[i], factor);
			sss_int0.push_back(ss_split.first);
			sss_int1.push_back(ss_split.second);
		}

		// SU
		std::vector<SUint> sus_int0;
		std::vector<SUint> sus_int1;
		for(unsigned int i = 0; i < sus.size(); ++i) {
			auto su_split = splitSecondaryUserInt(sus[i], factor);
			sus_int0.push_back(su_split.first);
			sus_int1.push_back(su_split.second);
		}

		if(args.do_plaintext_split) {
			t1.start(Timer::plaintext_split_preprocessing);

			pt_sm.plainTextRadarPreprocess(pus_int0, pus_int1, &sss_int0, &sss_int1);

			t1.end(Timer::plaintext_split_preprocessing);
		}

		std::map<int, std::vector<int> > precomputed_pu_int_groups;
		std::map<int, std::vector<int> > precomputed_ss_int_groups;
		if(sm_params.use_grid) {
			t1.start(Timer::plaintext_grid_preprocessing);
			pt_sm.plainTextGrid(pus, sss, &precomputed_pu_int_groups, &precomputed_ss_int_groups);
			t1.end(Timer::plaintext_grid_preprocessing);
		}

		KeyServer ks;
		if(args.central_entities == "two_sms") {
			// Do nothing
		} else if(args.central_entities == "sm_ks") {
			// Encrypt the '1's
			for(unsigned int i = 0; i < pus_int1.size(); ++i) {
				ks.encryptAndInit(&(pus_int1[i]));

				for(unsigned int j = 0; j < pus_int1[i].prs.size(); ++j) {
					ks.encryptAndInit(&(pus_int1[i].prs[j]));
				}
			}

			for(unsigned int i = 0; i < sss_int1.size(); ++i) {
				ks.encryptAndInit(&(sss_int1[i]));
			}

			for(unsigned int i = 0; i < sus_int1.size(); ++i) {
				ks.encryptAndInit(&(sus_int1[i]));
			}
		} else {
			std::cerr << "Unknown central_entities: " << args.central_entities << std::endl;
			exit(1);
		}

		NullTimer null_timer;

		// Used for communication between the two SMs.
		int num_io_threads = 1;
		const std::string server_addr = "127.0.0.1:1212";
		const std::string connection_name = "sm_connection";
		const std::string channel_name = "sm-to-sm_channel";
		sm_params.setCommunicationValues(num_io_threads, server_addr, connection_name, channel_name);

		std::thread thrd0([&]() {
			int party_id = 0;
			if(args.central_entities == "two_sms") {
				SM sm(party_id, &sm_params, pus_int0, sss_int0, sus_int0);
				sm.setSecureWriteTimer(&secure_write_timer);
	
				secure_vs_a = sm.run(precomputed_pu_int_groups, precomputed_ss_int_groups, &t1);
			} else if(args.central_entities == "sm_ks") {
				SM sm(party_id, &sm_params, pus_int0, sss_int0, sus_int0, pus_int1, sss_int1, sus_int1);
				sm.setSecureWriteTimer(&secure_write_timer);
				
				secure_vs_a = sm.runSM(precomputed_pu_int_groups, precomputed_ss_int_groups, &t1);
			}
		});

		std::thread thrd1([&]() {
			int party_id = 1;
			if(args.central_entities == "two_sms") {
				SM sm(party_id, &sm_params, pus_int1, sss_int1, sus_int1);
	
				secure_vs_b = sm.run(precomputed_pu_int_groups, precomputed_ss_int_groups, &null_timer);
			} else if(args.central_entities == "sm_ks") {
				SM sm(party_id, &sm_params, &ks); // Probably size information

				secure_vs_b = sm.runKS(int(sus_int1.size()));
			}
		});

		std::thread su_thrd([&]() {
			SUServer su_server(sus);
			su_server.run();
		});

		thrd0.join();
		thrd1.join();
		su_thrd.join();
	}

	PathLossTable path_loss_table;

	std::vector<float> plaintext_vs = pt_sm.plainTextRun(sus, pus, sss, &t1, &path_loss_table);
	std::vector<float> ground_truth_vs = gen.computeGroundTruth(sus, pus, &path_loss_table, args.no_pr_thresh_update);

	std::vector<float> secure_vs;
	for(unsigned int i = 0; i < sus.size(); ++i) {
		if(!args.skip_s2pc) {
			float secure_v = secure_vs_a[i] + secure_vs_b[i];
			secure_vs.push_back(secure_v);
		}
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

		std::cout << "path_loss(plain,ground)|list(float,float)|";
		unsigned int x = 0;
		for(auto itr = path_loss_table.table.begin(); itr != path_loss_table.table.end(); ++itr) {
			if(!itr->second.pt_set || !itr->second.gt_set) {
				std::cerr << "Path loss not set" << std::endl;
				exit(1);
			}

			std::cout << itr->second.pt_pl << ":" << itr->second.gt_pl;
			if(x < path_loss_table.table.size() - 1) {
				std::cout << ",";
			}
			++x;
		}
		std::cout << std::endl;

		if(!args.skip_s2pc) {
			std::cout << "preprocess_time|float|" << t1.getAverageDuration(Timer::secure_preprocessing) + t1.getAverageDuration(Timer::plaintext_split_preprocessing) + t1.getAverageDuration(Timer::plaintext_grid_preprocessing) << std::endl;
			std::cout << "time_per_request|float|" << t1.getAverageDuration(Timer::secure_su_request) << std::endl;
			std::cout << "secure_write_time|float|" << secure_write_timer.getAverageDuration(Timer::secure_write) << std::endl;
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
				std::cout << "% Diff PvS: " << fabs(secure_vs[i] - plaintext_vs[i]) / fabs(plaintext_vs[i]) * 100 << "%" << std::endl;
			}
			std::cout << "% Diff GvP: " << fabs(plaintext_vs[i] - ground_truth_vs[i]) / fabs(ground_truth_vs[i]) * 100 << "%" << std::endl;
			std::cout << "--------------------" << std::endl;
		}

		if(!args.skip_s2pc) {
			std::cout << "Average duration of " << Timer::secure_preprocessing << ": " << t1.getAverageDuration(Timer::secure_preprocessing) + t1.getAverageDuration(Timer::plaintext_split_preprocessing) + t1.getAverageDuration(Timer::plaintext_grid_preprocessing) << std::endl;
			std::cout << "Average duration of " << Timer::secure_su_request << ": " << t1.getAverageDuration(Timer::secure_su_request) << std::endl;
			std::cout << "Average duration of " << Timer::secure_write << ": " << secure_write_timer.getAverageDuration(Timer::secure_write) << std::endl;
		}
	}

	return 0;
}

