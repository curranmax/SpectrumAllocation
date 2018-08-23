
#ifndef _SPECTRUM_MANAGER_H_
#define _SPECTRUM_MANAGER_H_

#include "ivory/Runtime/Party.h"

#include "path_loss_table.h"
#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_sensor.h"
#include "tables.h"
#include "timer.h"

#include <cryptoTools/Network/Channel.h>

#include <vector>

using namespace osuCrypto;

class SpectrumManager {
public:
	enum AlgoOrder{
		SPLIT_THEN_IDW,
		IDW_THEN_SPLIT
	};

	enum SelectionAlgo {
		NONE,
		SORT,
		RANDOM
	};

	enum SecureWriteAlgo {
		PROPOSED,
		SPC
	};

	SpectrumManager() = delete;
	SpectrumManager(float _factor, int _num_ss_selection, int _num_pu_selection, float _rp_alpha, float _pl_alpha, int _bit_count, bool _brief_out)
			: factor(_factor), num_ss_selection(_num_ss_selection), num_pu_selection(_num_pu_selection), rp_alpha(_rp_alpha), pl_alpha(_pl_alpha), bit_count(_bit_count), brief_out(_brief_out),
			order(AlgoOrder::SPLIT_THEN_IDW), selection_algo(SelectionAlgo::SORT), secure_write_algo(SecureWriteAlgo::PROPOSED), secure_write_timer(nullptr), use_grid(false), grid_min_num_pu(0), grid_min_num_ss(0), grid_num_x(0), grid_num_y(0), grid_delta_x(0.0), grid_delta_y(0.0),
			num_io_threads(0), server_addr(""), connection_name(""), channel_name("") {}
	SpectrumManager(const SpectrumManager& sm) = delete;

	void setSelectionAlgo(const std::string& algo_string) {
		if(algo_string == "none") { selection_algo = SelectionAlgo::NONE; }
		else if(algo_string == "sort") { selection_algo = SelectionAlgo::SORT; }
		else if(algo_string == "random") { selection_algo = SelectionAlgo::RANDOM; }
		else { std::cerr << "Unknown selection_algo: " << algo_string << std::endl; exit(1); }
	}

	void setSecureWriteAlgo(const std::string algo_string, Timer* _secure_write_timer) {
		secure_write_timer = _secure_write_timer;
		if(algo_string == "proposed") { secure_write_algo = SecureWriteAlgo::PROPOSED; }
		else if(algo_string == "spc") { secure_write_algo = SecureWriteAlgo::SPC; }
		else { std::cerr << "Unknown selection_algo: " << algo_string << std::endl; exit(1); }
	}

	void setGridParams(int _grid_min_num_pu, int _grid_min_num_ss, int _grid_num_x, int _grid_num_y, float _grid_delta_x, float _grid_delta_y);

	void setCommunicationValues(int _num_io_threads, const std::string& _server_addr, const std::string& _connection_name, const std::string& _channel_name);

	std::vector<float> run(int party_id,
							const std::vector<PUint>& pus,
							std::vector<SSint>& sss,
							const std::vector<SUint>& sus,
							const std::map<int, std::vector<int> >& precomputed_pu_groups,
							const std::map<int, std::vector<int> >& precomputed_ss_groups,
							Timer* timer) const;

	// Assumes that each PU has a single PR and they are at the same location.
	void secureRadarPreprocess(std::array<Party, 2> parties,
								const std::vector<PUint>& pus,
								std::vector<SSint>* sss) const;

	void plainTextRadarPreprocess(
			const std::vector<PUint>& pus_int0, const std::vector<PUint>& pus_int1,
			std::vector<SSint>* sss_int0, std::vector<SSint>* sss_int1) const;

	std::map<int, std::vector<const SSint*> > secureSSGridPreprocess(
			std::array<Party, 2> parties, const std::vector<SSint>& sss) const;

	std::map<int, std::vector<const PUint*> > securePUGridPreprocess(
			std::array<Party, 2> parties, const std::vector<PUint>& sss) const;

	void secureGetEntitiesFromTable(
			std::array<Party, 2> parties, const SUint& su,
			const GridTable& grid_table, const PUTable& pu_table,
			std::vector<PUint>* selected_pus, std::vector<SSint>* selected_sss) const;

	float secureRadar(std::array<Party, 2> parties,
						const SUint& su,
						const std::vector<PUint>& pus,
						const std::vector<SSint>& sss,
						PUTable* pu_table,
						Channel* sm_ch) const;
	void secureTableWrite(
			std::array<Party, 2> parties,
			PUTable* pu_table,
			const std::vector<std::pair<sInt, std::vector<sInt> > >& updates,
			Channel* sm_ch) const;

	std::vector<float> plainTextRun(
			const std::vector<SU>& sus,
			const std::vector<PU>& pus, const std::vector<SS>& sss,
			Timer* timer,
			PathLossTable* path_loss_table) const;

	void plainTextGrid(const std::vector<PU>& pus, const std::vector<SS>& sss,
			std::map<int, std::vector<int> >* pu_int_groups, std::map<int, std::vector<int> >* ss_int_groups) const;

	float plainTextRadar(const SU& su,
							const std::vector<const PU*>& pus,
							const std::vector<const SS*>& sss,
							PathLossTable* path_loss_table) const;

	bool useGrid() const { return use_grid; }
private:
	unsigned int numSelect(int input_val, unsigned int size) const {
		if(input_val <= 0 || (unsigned int)input_val > size) {
			return size;
		}
		return input_val;
	}

	float factor;

	int num_ss_selection;
	int num_pu_selection;

	int rp_alpha;
	int pl_alpha;

	int bit_count;

	bool brief_out;

	AlgoOrder order;
	SelectionAlgo selection_algo;
	SecureWriteAlgo secure_write_algo;

	Timer* secure_write_timer;

	bool use_grid;
	int grid_min_num_pu, grid_min_num_ss;
	int grid_num_x, grid_num_y;
	float grid_delta_x, grid_delta_y;

	// Deviations for group expansion.
	void initGroupDeviations();
	const std::vector<std::pair<int, int> >& getDeviations(int iter) const;

	std::vector<std::vector<std::pair<int, int> > > all_deviations;

	int num_io_threads;
	std::string server_addr, connection_name, channel_name;
};

typedef SpectrumManager SM;

#endif
