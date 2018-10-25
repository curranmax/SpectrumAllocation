
#include "args.h"
#include "debug_print.h"
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
#include <set>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <math.h>
#include <sstream>
#include <algorithm>

#define PDIF(prefix, v1, v2, pdif_limit) {if(fabs(v1 - v2) / fabs(v1) >= pdif_limit) { std::cout << prefix << " " << v1 << ", " << v2 << ", " << fabs(v1 - v2) / fabs(v1) << std::endl; }}

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
	PropagationModel* pm = nullptr;
	if(args.propagation_model == "log_distance") {
		pm = new LogDistancePM(args.ld_path_loss0, args.ld_dist0, args.ld_gamma);
	} else if (args.propagation_model == "longley_rice") {
		pm = new LongleyRicePM(args.splat_cmd, args.ref_lat, args.ref_long,
								args.splat_dir, args.sdf_dir, args.return_dir,
								args.use_itwom_pl);
	} else if(args.propagation_model == "single_lr") {
		pm = new TransmitterOnlySplatPM(args.splat_cmd, args.ref_lat, args.ref_long,
										args.splat_dir, args.sdf_dir, args.return_dir);
	} else if(args.propagation_model == "input_file") {
		// Nothing for now
		pm = new NullPropagationModel();
	} else {
		std::cerr << "Unknown propagation model: " << args.propagation_model << std::endl;
		exit(0);
	}

	// Generate enetities
	Generator gen(args.location_range, args.su_buffer, pm);

	std::vector<PU> pus;
	std::vector<SS> sss;
	std::vector<SU> sus;

	PathLossTable path_loss_table;
	std::vector<std::vector<float> > rp_at_ss_from_pu, rp_at_ss_from_pu_pt, rp_at_ss_from_pu_uo_pt;
	std::vector<std::vector<float> > gt_su_pu_pl, pt_su_pu_pl, uo_pt_su_pu_pl; // gt_su_pu_pl[i][j] is the GT path loss between SU i and PU j
	if(args.in_filename == "") {
		gen.generateEntities(args.num_pu, args.num_ss, args.num_su, args.num_pr_per_pu, args.pr_range, args.out_filename, &pus, &sss, &sus, &rp_at_ss_from_pu, &gt_su_pu_pl);
	} else {
		gen.getEntitiesFromFile(args.num_pu, args.num_ss, args.num_su, args.num_pr_per_pu, args.in_filename, &pus, &sss, &sus, &rp_at_ss_from_pu, &gt_su_pu_pl);
	}

	if(args.use_gt_rp_at_ss_from_pu) {
		rp_at_ss_from_pu_uo_pt = rp_at_ss_from_pu;

	}

	if(args.use_gt_su_pu_pl) {
		pt_su_pu_pl = gt_su_pu_pl;
		uo_pt_su_pu_pl = gt_su_pu_pl;
	}
	
	// Set up Spectrum Manager params
	float factor = 1 << args.num_float_bits;
	SMParams sm_params(factor, args.num_float_bits, args.s2_pc_bit_count, args.num_pu_selection, args.num_ss_selection, args.ss_receive_power_alpha, args.path_loss_alpha, args.brief_out);
	sm_params.setSelectionAlgo(args.selection_algo);

	sm_params.setSecureWriteAlgo(args.secure_write_algo);
	sm_params.no_pr_thresh_update = args.no_pr_thresh_update;
	sm_params.pl_est_gamma = args.pl_est_gamma;
	sm_params.pt_record_split_power = args.pt_record_split_power;

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

	SMParams uo_params(factor, args.num_float_bits, args.s2_pc_bit_count, 0, 0, args.ss_receive_power_alpha, args.path_loss_alpha, args.brief_out);
	uo_params.setSelectionAlgo("none");

	uo_params.no_pr_thresh_update = args.no_pr_thresh_update;
	uo_params.pl_est_gamma = args.pl_est_gamma;
	uo_params.pt_record_split_power = args.pt_record_split_power;

	PtSM uo_pt_sm(&uo_params);

	std::vector<float> secure_vs_a;
	std::vector<float> secure_vs_b;
	Timer t1, secure_write_timer;
	std::map<std::string, Timer> sm_timers, ks_timers;
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
			pt_sm.plainTextGrid(sus, pus, sss, &precomputed_pu_int_groups, &precomputed_ss_int_groups);
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

		Shared shared_memory;

		std::thread thrd0([&]() {
			int party_id = 0;
			if(args.central_entities == "two_sms") {
				SM sm(party_id, &sm_params, pus_int0, sss_int0, sus_int0);
				sm.setSecureWriteTimer(&secure_write_timer);
	
				secure_vs_a = sm.run(precomputed_pu_int_groups, precomputed_ss_int_groups, &t1);
			} else if(args.central_entities == "sm_ks") {
				SM sm(party_id, &sm_params, pus_int0, sss_int0, sus_int0, pus_int1, sss_int1, sus_int1, &shared_memory);
				sm.setSecureWriteTimer(&secure_write_timer);
				
				secure_vs_a = sm.runSM(precomputed_pu_int_groups, precomputed_ss_int_groups, &t1, sm_timers);
			}
		});

		std::thread thrd1([&]() {
			int party_id = 1;
			if(args.central_entities == "two_sms") {
				SM sm(party_id, &sm_params, pus_int1, sss_int1, sus_int1);
	
				secure_vs_b = sm.run(precomputed_pu_int_groups, precomputed_ss_int_groups, &null_timer);
			} else if(args.central_entities == "sm_ks") {
				SM sm(party_id, &sm_params, &ks, &shared_memory);

				secure_vs_b = sm.runKS(int(sus_int1.size()), ks_timers);
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

	std::vector<float> plaintext_vs = pt_sm.plainTextRun(sus, pus, sss, &t1, &path_loss_table, &rp_at_ss_from_pu_pt, &pt_su_pu_pl);
	std::vector<float> ground_truth_vs = gen.computeGroundTruth(sus, pus, &path_loss_table, args.no_pr_thresh_update);

	std::vector<float> uo_plaintext_vs;
	if(args.run_unoptimized_plaintext) {
		uo_plaintext_vs = uo_pt_sm.unoptimizedPlaintextRun(sus, pus, sss, &t1, &path_loss_table, &rp_at_ss_from_pu_uo_pt, &uo_pt_su_pu_pl);
	}

	{
		// Add ground truth SU-PU path losses to the path loss table
		for(unsigned int i = 0; i < gt_su_pu_pl.size(); ++i) {
			for(unsigned int j = 0; j < gt_su_pu_pl[i].size(); ++j) {
				path_loss_table.addGroundTruthPathLoss(i, j, gt_su_pu_pl[i][j]);
			}
		}
	}

	std::vector<float> secure_vs;
	for(unsigned int i = 0; i < sus.size(); ++i) {
		if(!args.skip_s2pc) {
			float secure_v = secure_vs_a[i] + secure_vs_b[i];
			secure_vs.push_back(secure_v);
		}
	}

	if(args.brief_out) {
		std::stringstream names, types;
		if(!args.skip_s2pc) {
			names << "secure,";
			types << "float,";
		}

		names << "plain,";
		types << "float,";

		if(args.run_unoptimized_plaintext) {
			names << "uo,";
			types << "float,";
		}

		names << "ground";
		types << "float";

		std::cout << "su_transmit_power(" << names.str() << ")|list(" << types.str() << ")|";
		for(unsigned int i = 0; i < plaintext_vs.size(); ++i) {
			if(!args.skip_s2pc) {
				std::cout << secure_vs[i] << ":";
			}

			std::cout << plaintext_vs[i] << ":";

			if(args.run_unoptimized_plaintext) {
				std::cout << uo_plaintext_vs[i] << ":";
			}

			std::cout << ground_truth_vs[i];
			if(i < plaintext_vs.size() - 1) {
				std::cout << ",";
			}
		}
		std::cout << std::endl;

		if(path_loss_table.table.size() > 0) {
			if(args.run_unoptimized_plaintext) {
				std::cout << "su_pr_path_loss(plain,uo,ground)|list(float,float,float)|";
			} else {
				std::cout << "su_pr_path_loss(plain,ground)|list(float,float)|";
			}
			unsigned int x = 0;
			for(auto itr = path_loss_table.table.begin(); itr != path_loss_table.table.end(); ++itr) {
				if(!itr->second.pt_set || !itr->second.gt_set || (args.run_unoptimized_plaintext && !itr->second.uo_pt_set)) {
					std::cerr << std::endl << "Path loss not set: " << !itr->second.pt_set << ", " << !itr->second.gt_set << ", " << args.run_unoptimized_plaintext << ", " << !itr->second.uo_pt_set << std::endl;
					exit(1);
				}

				std::cout << itr->second.pt_pl << ":";

				if(args.run_unoptimized_plaintext) {
					std::cout << itr->second.uo_pt_pl << ":";
				}

				std::cout << itr->second.gt_pl;
				if(x < path_loss_table.table.size() - 1) {
					std::cout << ",";
				}
				++x;
			}
			std::cout << std::endl;
		}

		if(path_loss_table.pu_table.size() > 0) {
			std::cout << "su_pu_path_loss(plain,ground)|list(float,float)|";
			unsigned int x = 0;
			for(auto itr = path_loss_table.pu_table.begin(); itr != path_loss_table.pu_table.end(); ++itr) {
				if(!itr->second.pt_set || !itr->second.gt_set) {
					std::cerr << "Path loss not set" << std::endl;
					exit(1);
				}

				std::cout << itr->second.pt_pl << ":" << itr->second.gt_pl;
				if(x < path_loss_table.pu_table.size() - 1) {
					std::cout << ",";
				}
				++x;
			}
			std::cout << std::endl;
		}

		if(!args.skip_s2pc) {
			std::cout << "preprocess_time|float|" << t1.getAverageDuration(Timer::secure_preprocessing) + t1.getAverageDuration(Timer::plaintext_split_preprocessing) + t1.getAverageDuration(Timer::plaintext_grid_preprocessing) << std::endl;
			std::cout << "time_per_request|float|" << t1.getAverageDuration(Timer::secure_su_request) << std::endl;
			std::cout << "secure_write_time|float|" << secure_write_timer.getAverageDuration(Timer::secure_write) << std::endl;
		}

		std::cout << Timer::opt_pt_request << "|float|" << t1.getAverageDuration(Timer::opt_pt_request) << std::endl;
		if(args.run_unoptimized_plaintext) {
			std::cout << Timer::unopt_pt_request << "|float|" << t1.getAverageDuration(Timer::unopt_pt_request) << std::endl;
		}

		if(!args.use_gt_rp_at_ss_from_pu && args.pt_record_split_power) {
			std::map<int, std::vector<int> > precomputed_pu_int_groups;
			std::map<int, std::vector<int> > precomputed_ss_int_groups;

			P("start grid comp for rp reporting");
			pt_sm.plainTextGrid(sus, pus, sss, &precomputed_pu_int_groups, &precomputed_ss_int_groups);
			P("end grid comp for rp reporting");
			

			P("start revelant SS-PU pairs");
			std::set<std::pair<int, int> > relevant_ss_pu;
			for(auto pu_itr = precomputed_pu_int_groups.begin(); pu_itr != precomputed_pu_int_groups.end(); ++pu_itr) {
				auto ss_itr = precomputed_ss_int_groups.find(pu_itr->first);
				if(ss_itr == precomputed_ss_int_groups.end()) {
					std::cerr << "Mismatch groups: " << pu_itr->first << std::endl;
					exit(1);
				}

				for(unsigned int x = 0; x < pu_itr->second.size(); ++x) {
					for(unsigned int y = 0; y < ss_itr->second.size(); ++y) {
						relevant_ss_pu.insert(std::make_pair(pu_itr->second[x], ss_itr->second[y]));
					}
				}
			}
			P("end revelant SS-PU pairs");

			std::cout << "rp_at_ss_from_pu|list(float)|";
			unsigned int x = 0;
			for(auto itr = relevant_ss_pu.begin(); itr != relevant_ss_pu.end(); ++itr) {
				std::cout << rp_at_ss_from_pu[itr->first][itr->second];
				if(x < relevant_ss_pu.size() - 1) {
					std::cout << ",";
				}
				++x;
			}
			std::cout << std::endl;

			std::cout << "rp_at_ss_from_pu_pt|list(float)|";
			x = 0;
			for(auto itr = relevant_ss_pu.begin(); itr != relevant_ss_pu.end(); ++itr) {
				std::cout << rp_at_ss_from_pu_pt[itr->first][itr->second];
				if(x < relevant_ss_pu.size() - 1) {
					std::cout << ",";
				}
				++x;
			}
			std::cout << std::endl;

			if(args.run_unoptimized_plaintext) {
				std::cout << "rp_at_ss_from_pu_uo_pt|list(float)|";
				x = 0;
				for(auto itr = relevant_ss_pu.begin(); itr != relevant_ss_pu.end(); ++itr) {
					std::cout << rp_at_ss_from_pu_uo_pt[itr->first][itr->second];
					if(x < relevant_ss_pu.size() - 1) {
						std::cout << ",";
					}
					++x;
				}
				std::cout << std::endl;
			}
		}
	} else {
		for(unsigned int i = 0; i < plaintext_vs.size(); ++i) {
			std::cout << "--------------------" << std::endl;
			if(!args.skip_s2pc) {
				std::cout << "Secure:     " << secure_vs[i] << std::endl;
			}
			std::cout << "Plain:      " << plaintext_vs[i] << std::endl;
			if(args.run_unoptimized_plaintext) {
				std::cout << "UnOpt PT:   " << uo_plaintext_vs[i] << std::endl;
			}
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

		std::cout << "Average duration of " << Timer::opt_pt_request << ": " << t1.getAverageDuration(Timer::opt_pt_request) << std::endl;
		if(args.run_unoptimized_plaintext) {
			std::cout << "Average duration of " << Timer::unopt_pt_request << ": " << t1.getAverageDuration(Timer::unopt_pt_request) << std::endl;
		}
	}

	if(!args.skip_s2pc && args.central_entities == "sm_ks") {
		// SM - "total"
		auto sm_total = sm_timers["total"];
		std::cout << "sendEncryptedData|float|" << sm_total.getAverageDuration("sendEncryptedData") << std::endl;
		std::cout << "recvEncryptedPRThresholds|float|" << sm_total.getAverageDuration("recvEncryptedPRThresholds") << std::endl;
	
		// KS - "total"
		auto ks_total = ks_timers["total"];
		std::cout << "recvEncryptedData_recv|float|" << ks_total.getAverageDuration("recvEncryptedData-recv") << std::endl;
		std::cout << "recvEncryptedData_decrypt|float|" << ks_total.getAverageDuration("recvEncryptedData-decrypt") << std::endl;
		std::cout << "sendEncryptedPRThresholds|float|" << ks_total.getAverageDuration("sendEncryptedPRThresholds") << std::endl;

		// SM - "entities"
		auto sm_entity = sm_timers["entities"];
		std::cout << "send_SU|float|" << sm_entity.getAverageDuration("SU") << std::endl;
		std::cout << "send_SU_num|float|" << sm_entity.numDurations("SU") << std::endl;

		std::cout << "send_GT|float|" << sm_entity.getAverageDuration("GT") << std::endl;
		std::cout << "send_GT_num|float|" << sm_entity.numDurations("GT") << std::endl;

		// std::cout << "send_SS|float|" << sm_entity.getAverageDuration("SS") << std::endl;
		// std::cout << "send_SS_num|float|" << sm_entity.numDurations("SS") << std::endl;

		// std::cout << "send_PU_inds|float|" << sm_entity.getAverageDuration("PU inds") << std::endl;
		// std::cout << "send_PU_inds_num|float|" << sm_entity.numDurations("PU inds") << std::endl;

		std::cout << "send_PU|float|" << sm_entity.getAverageDuration("PU") << std::endl;
		std::cout << "send_PU_num|float|" << sm_entity.numDurations("PU") << std::endl;

		std::cout << "send_PR|float|" << sm_entity.getAverageDuration("PR") << std::endl;
		std::cout << "send_PR_num|float|" << sm_entity.numDurations("PR") << std::endl;

		std::cout << "recv_PR_thresh|float|" << sm_entity.getAverageDuration("PR thresh") << std::endl;
		std::cout << "recv_PR_thresh_num|float|" << sm_entity.numDurations("PR thresh") << std::endl;

		// KS - "entities"
		auto ks_entity = ks_timers["entities"];
		std::cout << "recv_SU|float|" << ks_entity.getAverageDuration("SU") << std::endl;
		std::cout << "recv_SU_num|float|" << ks_entity.numDurations("SU") << std::endl;

		std::cout << "recv_GT|float|" << ks_entity.getAverageDuration("GT") << std::endl;
		std::cout << "recv_GT_num|float|" << ks_entity.numDurations("GT") << std::endl;

		// std::cout << "recv_SS|float|" << ks_entity.getAverageDuration("SS") << std::endl;
		// std::cout << "recv_SS_num|float|" << ks_entity.numDurations("SS") << std::endl;

		// std::cout << "recv_PU_inds|float|" << ks_entity.getAverageDuration("PU inds") << std::endl;
		// std::cout << "recv_PU_inds_num|float|" << ks_entity.numDurations("PU inds") << std::endl;

		std::cout << "recv_PU|float|" << ks_entity.getAverageDuration("PU") << std::endl;
		std::cout << "recv_PU_num|float|" << ks_entity.numDurations("PU") << std::endl;

		std::cout << "recv_PR|float|" << ks_entity.getAverageDuration("PR") << std::endl;
		std::cout << "recv_PR_num|float|" << ks_entity.numDurations("PR") << std::endl;

		std::cout << "decrypt_SU|float|" << ks_entity.getAverageDuration("SU-decrypt") << std::endl;
		std::cout << "decrypt_SU_num|float|" << ks_entity.numDurations("SU-decrypt") << std::endl;

		std::cout << "decrypt_SS|float|" << ks_entity.getAverageDuration("SS-decrypt") << std::endl;
		std::cout << "decrypt_SS_num|float|" << ks_entity.numDurations("SS-decrypt") << std::endl;

		std::cout << "decrypt_PU|float|" << ks_entity.getAverageDuration("PU-decrypt") << std::endl;
		std::cout << "decrypt_PU_num|float|" << ks_entity.numDurations("PU-decrypt") << std::endl;

		std::cout << "decrypt_PR|float|" << ks_entity.getAverageDuration("PR-decrypt") << std::endl;
		std::cout << "decrypt_PR_num|float|" << ks_entity.numDurations("PR-decrypt") << std::endl;

		std::cout << "encrypt_PR_thresh|float|" << ks_entity.getAverageDuration("PR thresh-encrypt") << std::endl;
		std::cout << "encrypt_PR_thresh_num|float|" << ks_entity.numDurations("PR thresh-encrypt") << std::endl;

		std::cout << "send_PR_thresh|float|" << ks_entity.getAverageDuration("PR thresh") << std::endl;
		std::cout << "send_PR_thresh_num|float|" << ks_entity.numDurations("PR thresh") << std::endl;
	}

	return 0;
}

