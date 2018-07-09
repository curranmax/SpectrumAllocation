
#ifndef _SPECTRUM_MANAGER_H_
#define _SPECTRUM_MANAGER_H_

#include "ivory/Runtime/Party.h"

#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_sensor.h"
#include "timer.h"

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

	SpectrumManager() = delete;
	SpectrumManager(float _factor, int _num_ss_selection, int _num_pu_selection, float _rp_alpha, float _pl_alpha, int _bit_count, bool _brief_out)
			: factor(_factor), num_ss_selection(_num_ss_selection), num_pu_selection(_num_pu_selection), rp_alpha(_rp_alpha), pl_alpha(_pl_alpha), bit_count(_bit_count), brief_out(_brief_out),
			order(AlgoOrder::SPLIT_THEN_IDW), selection_algo(SORT), use_grid(false), grid_min_num(0), grid_num_x(0), grid_num_y(0), grid_delta_x(0.0), grid_delta_y(0.0) {}
	SpectrumManager(const SpectrumManager& sm) = delete;

	void setAlgoOrder(const std::string& order_string) {
		if(order_string == "split_then_idw") { order = AlgoOrder::SPLIT_THEN_IDW; }
		else if(order_string == "idw_then_split") { order = AlgoOrder::IDW_THEN_SPLIT; }
		else { std::cerr << "Unknown algo_order: " << order_string << std::endl; exit(1); }
	}

	void setSelectionAlgo(const std::string& algo_string) {
		if(algo_string == "none") { selection_algo = SelectionAlgo::NONE; }
		else if(algo_string == "sort") { selection_algo = SelectionAlgo::SORT; }
		else if(algo_string == "random") { selection_algo = SelectionAlgo::RANDOM; }
		else { std::cerr << "Unknown selection_algo: " << algo_string << std::endl; exit(1); }
	}

	void setGridParams(int _grid_min_num, int _grid_num_x, int _grid_num_y, float _grid_delta_x, float _grid_delta_y);

	std::vector<float> run(int party_id,
							const std::vector<PUint>& pus,
							std::vector<SSint>& sss,
							const std::vector<SUint>& sus,
							const std::map<int, std::vector<int> >& precomputed_ss_groups,
							Timer* timer) const;

	// Assumes that each PU has a single PR and they are at the same location.
	void secureRadarPreprocess(std::array<Party, 2> parties,
								const std::vector<PUint>& pus,
								std::vector<SSint>* sss) const;

	void plainTextRadarPreprocess(
			const std::vector<PU>& pus, const std::vector<SS>& sss,
			std::vector<SSint>* sss_int0, std::vector<SSint>* sss_int1) const;

	std::map<int, std::vector<const SSint*> > secureGridPreprocess(
			std::array<Party, 2> parties, const std::vector<SSint>& sss) const;

	void secureGetSSForSUFromGrid(std::array<Party, 2> parties,
									const SUint& su,
									const std::map<int, std::vector<const SSint*> >& ss_groups,
									std::vector<const SSint*>* selected_sss) const;

	float secureRadar(std::array<Party, 2> parties,
						const SUint& su,
						const std::vector<PUint>& pus,
						const std::vector<const SSint*>& sss) const;

	float secureRadarAlternateOrder(
						std::array<Party, 2> parties,
						const SUint& su,
						const std::vector<PUint>& pus,
						const std::vector<SSint>& sss) const;

	std::vector<float> plainTextRun(
			const std::vector<SU>& sus,
			const std::vector<PU>& pus, const std::vector<SS>& sss,
			Timer* timer) const;

	std::map<int, std::vector<int> > plainTextGrid(const std::vector<SS>& sss) const;

	float plainTextRadar(const SU& su,
							const std::vector<PU>& pus,
							const std::vector<const SS*>& sss) const;

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

	bool use_grid;
	int grid_min_num;
	int grid_num_x, grid_num_y;
	float grid_delta_x, grid_delta_y;

	// Deviations for group expansion.
	void initGroupDeviations();
	const std::vector<std::pair<int, int> >& getDeviations(int iter) const;

	std::vector<std::vector<std::pair<int, int> > > all_deviations;
};

typedef SpectrumManager SM;

#endif
