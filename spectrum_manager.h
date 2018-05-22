
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

	enum PathLossType {
		DB,
		RATIO
	};

	SpectrumManager() = delete;
	SpectrumManager(float _factor, int _num_ss_selection, float _rp_alpha, float _pl_alpha, int _bit_count, bool _brief_out)
			: factor(_factor), num_ss_selection(_num_ss_selection), rp_alpha(_rp_alpha), pl_alpha(_pl_alpha), bit_count(_bit_count), brief_out(_brief_out),
			order(AlgoOrder::SPLIT_THEN_IDW), pl_type(PathLossType::RATIO) {}
	SpectrumManager(const SpectrumManager& sm) = delete;

	void setAlgoOrder(const std::string order_string) {
		if(order_string == "split_then_idw") { order = AlgoOrder::SPLIT_THEN_IDW; }
		else if(order_string == "idw_then_split") { order = AlgoOrder::IDW_THEN_SPLIT; }
		else { std::cerr << "Unknown algo_order: " << order_string << std::endl; }
	}

	void setPathLossType(const std::string pl_type_string) {
		if(pl_type_string == "ratio") { pl_type = PathLossType::RATIO; }
		else if(pl_type_string == "db") { pl_type = PathLossType::DB; }
		else { std::cerr << "Unknown path_loss_type: " << pl_type_string << std::endl; }
	}


	std::vector<float> run(int party_id,
							const std::vector<PUint>& pus,
							const std::vector<SSint>& sss,
							const std::vector<SUint>& sus,
							Timer* timer) const;

	// Assumes that each PU has a single PR and they are at the same location.
	void secureRadarPreprocess(std::array<Party, 2> parties,
								const std::vector<PUint>& pus,
								const std::vector<SSint>& sss,
								std::vector<std::vector<int> >* plaintext_sss_rp_from_pu) const;
	float secureRadar(std::array<Party, 2> parties,
						const SUint& su,
						const std::vector<PUint>& pus,
						const std::vector<SSint>& sss,
						const std::vector<std::vector<int> >& plaintext_sss_rp_from_pu) const;
	float plainTextRadar(const SU& su,
							const std::vector<PU>& pus,
							const std::vector<SS>& sss) const;
private:
	float factor;

	int num_ss_selection;

	int rp_alpha;
	int pl_alpha;

	int bit_count;

	bool brief_out;

	AlgoOrder order;
	PathLossType pl_type;
};

typedef SpectrumManager SM;

#endif
