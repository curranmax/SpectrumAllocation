
#ifndef _SPECTRUM_MANAGER_H_
#define _SPECTRUM_MANAGER_H_

#include "ivory/Runtime/Party.h"

#include "key_server.h"
#include "path_loss_table.h"
#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_sensor.h"
#include "tables.h"
#include "timer.h"

#include <cryptoTools/Network/Channel.h>

#include <vector>

using namespace osuCrypto;

class SMParams {
public:
	enum SelectionAlgo {
		NONE,
		SORT,
		RANDOM
	};

	enum SecureWriteAlgo {
		PROPOSED,
		SPC
	};

	SMParams() = delete;
	SMParams(float _factor, int _bit_count, int _num_pu_selection, int _num_ss_selection, int _rp_alpha, int _pl_alpha, bool _brief_out) :
			brief_out(_brief_out), factor(_factor), bit_count(_bit_count), num_ss_selection(_num_ss_selection), num_pu_selection(_num_pu_selection),
			rp_alpha(_rp_alpha), pl_alpha(_pl_alpha), pl_est_gamma(0.0), no_pr_thresh_update(false), selection_algo(NONE), secure_write_algo(PROPOSED),
			use_grid(false), grid_min_num_pu(0), grid_min_num_ss(0), grid_num_x(0), grid_num_y(0), grid_delta_x(0), grid_delta_y(0),
			num_io_threads(0), server_addr(""), connection_name(""), channel_name("") {}

	SMParams(const SMParams& sm_params) = delete;
	const SMParams& operator=(const SMParams& sm_params) = delete;

	void setSelectionAlgo(const std::string& algo_string) {
		if(algo_string == "none") { selection_algo = SelectionAlgo::NONE; }
		else if(algo_string == "sort") { selection_algo = SelectionAlgo::SORT; }
		else if(algo_string == "random") { selection_algo = SelectionAlgo::RANDOM; }
		else { std::cerr << "Unknown selection_algo: " << algo_string << std::endl; exit(1); }
	}

	void setSecureWriteAlgo(const std::string algo_string) {
		if(algo_string == "proposed") { secure_write_algo = SecureWriteAlgo::PROPOSED; }
		else if(algo_string == "spc") { secure_write_algo = SecureWriteAlgo::SPC; }
		else { std::cerr << "Unknown selection_algo: " << algo_string << std::endl; exit(1); }
	}

	void setGridParams(int _grid_min_num_pu, int _grid_min_num_ss, int _grid_num_x, int _grid_num_y, float _grid_delta_x, float _grid_delta_y) {
		use_grid = true;
		grid_min_num_pu = _grid_min_num_pu;
		grid_min_num_ss = _grid_min_num_ss;
		grid_num_x = _grid_num_x;
		grid_num_y = _grid_num_y;
		grid_delta_x = _grid_delta_x;
		grid_delta_y = _grid_delta_y;

		initGroupDeviations();
	}

	void setCommunicationValues(int _num_io_threads, const std::string& _server_addr, const std::string& _connection_name, const std::string& _channel_name) {
		num_io_threads = _num_io_threads;
		server_addr = _server_addr;
		connection_name = _connection_name;
		channel_name = _channel_name;
	}

	void initGroupDeviations();
	const std::vector<std::pair<int, int> >& getDeviations(int iter) const;

	std::vector<std::vector<std::pair<int, int> > > all_deviations;

	bool brief_out;

	float factor;
	int bit_count;

	int num_ss_selection, num_pu_selection;
	int rp_alpha, pl_alpha;

	float pl_est_gamma;

	bool no_pr_thresh_update;

	SelectionAlgo selection_algo;
	SecureWriteAlgo secure_write_algo;

	bool use_grid;
	int grid_min_num_pu, grid_min_num_ss;
	int grid_num_x, grid_num_y;
	float grid_delta_x, grid_delta_y;

	int num_io_threads;
	std::string server_addr, connection_name, channel_name;
};

class SpectrumManager {
public:
	SpectrumManager() = delete;
	SpectrumManager(const SpectrumManager& sm) = delete;
	
	SpectrumManager(int _party_id, const SMParams* _sm_params, const std::vector<PUint>& _pus, const std::vector<SSint>& _sss, const std::vector<SUint>& _sus)
			: party_id(_party_id), sm_params(_sm_params), all_pus(_pus), all_sss(_sss), all_sus(_sus), secure_write_timer(nullptr) {}
	
	SpectrumManager(int _party_id, const SMParams* _sm_params, const std::vector<PUint>& _pus, const std::vector<SSint>& _sss, const std::vector<SUint>& _sus,
					const std::vector<PUint>& _en_pus, const std::vector<SSint>& _en_sss, const std::vector<SUint>& _en_sus)
			: party_id(_party_id), sm_params(_sm_params), all_pus(_pus), en_pus(_en_pus), all_sss(_sss), en_sss(_en_sss), all_sus(_sus), en_sus(_en_sus), secure_write_timer(nullptr) {}
	
	SpectrumManager(int _party_id, const SMParams* _sm_params, KeyServer* _key_server)
			: party_id(_party_id), sm_params(_sm_params), key_server(_key_server) {}

	void setSecureWriteTimer(Timer* _secure_write_timer) {
		secure_write_timer = _secure_write_timer;
	}

	std::vector<float> run(
			const std::map<int, std::vector<int> >& precomputed_pu_groups, const std::map<int, std::vector<int> >& precomputed_ss_groups,
			Timer* timer);

	std::vector<float> runSM(
			const std::map<int, std::vector<int> >& precomputed_pu_groups, const std::map<int, std::vector<int> >& precomputed_ss_groups,
			Timer* timer);

	std::vector<float> runKS(int num_sus);

	// Assumes that each PU has a single PR and they are at the same location.
	void secureRadarPreprocess(std::array<Party, 2> parties);

	std::map<int, std::vector<const SSint*> > secureSSGridPreprocess(std::array<Party, 2> parties);

	std::map<int, std::vector<const PUint*> > securePUGridPreprocess(std::array<Party, 2> parties);

	void secureGetEntitiesFromTable(
			std::array<Party, 2> parties, const SUint& su,
			const GridTable& grid_table, const PUTable& pu_table,
			std::vector<PUint>* selected_pus, std::vector<SSint>* selected_sss) const;

	float secureRadar(std::array<Party, 2> parties,
						const SUint& su,
						const std::vector<PUint>& pus,
						const std::vector<SSint>& sss,
						PUTable* pu_table,
						Channel* sm_ch,
						Channel* su_ch) const;
	void secureTableWrite(
			std::array<Party, 2> parties,
			PUTable* pu_table,
			const std::vector<std::pair<sInt, std::vector<sInt> > >& updates,
			Channel* sm_ch) const;

	void sendEncryptedData(Channel* sm_ch, const SUint& su, const GridTable& grid_table, const PUTable& pu_table);
	void recvEncryptedData(Channel* sm_ch, SUint* su, GridTable* grid_table, PUTable* pu_table);

	void sendEncryptedPRThresholds(Channel* sm_ch, const PUTable& pu_table);
	void recvEncryptedPRThresholds(Channel* sm_ch, PUTable* pu_table);

	bool useGrid() const { return sm_params->use_grid; }
private:
	int party_id;

	const SMParams* sm_params;

	std::vector<PUint> all_pus, en_pus;
	std::vector<SSint> all_sss, en_sss;
	std::vector<SUint> all_sus, en_sus;

	KeyServer* key_server;

	Timer* secure_write_timer;
};

class PlaintextSpectrumManager {
public:
	PlaintextSpectrumManager() = delete;

	PlaintextSpectrumManager(const SMParams* _sm_params) : sm_params(_sm_params) {}

	PlaintextSpectrumManager(const PlaintextSpectrumManager& pt_sm) = delete;
	const PlaintextSpectrumManager& operator=(const PlaintextSpectrumManager& pt_sm) = delete;

	void plainTextRadarPreprocess(
			const std::vector<PUint>& pus_int0, const std::vector<PUint>& pus_int1,
			std::vector<SSint>* sss_int0, std::vector<SSint>* sss_int1) const;

	std::vector<float> plainTextRun(
			const std::vector<SU>& sus,
			const std::vector<PU>& input_pus, const std::vector<SS>& sss,
			Timer* timer,
			PathLossTable* path_loss_table) const;

	void plainTextGrid(const std::vector<PU>& pus, const std::vector<SS>& sss,
			std::map<int, std::vector<int> >* pu_int_groups, std::map<int, std::vector<int> >* ss_int_groups) const;

	float plainTextRadar(const SU& su,
							std::vector<PU*>& pus,
							const std::vector<const SS*>& sss,
							PathLossTable* path_loss_table) const;

	bool useGrid() const { return sm_params->use_grid; }

private:
	const SMParams* sm_params;
};

typedef SpectrumManager SM;
typedef PlaintextSpectrumManager PtSM;

#endif
