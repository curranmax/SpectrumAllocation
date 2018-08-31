
#ifndef _GENERATOR_H_
#define _GENERATOR_H_

#include "location.h"
#include "path_loss_table.h"
#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_sensor.h"
#include "utils.h"

#include <string>
#include <vector>

class PropagationModel {
public:
	PropagationModel() {}

	virtual void preprocessPathLoss(PU* pu, float ss_height, int pu_id) {}
	virtual void loadANOFile(const PU& pu) {}

	virtual void preprocessPathLoss(PR* pr, float ss_height, int pu_id, int pr_id) {}
	virtual void loadANOFile(const PR& pr) {}

	virtual float getPathLoss(const Location& loc1, const Location& loc2) const = 0;
	virtual ~PropagationModel() {}
};

class Generator {
public:
	Generator() = delete;
	Generator(float location_range_, PropagationModel* pm_) : location_range(location_range_), pm(pm_), pr_pls() {}
	Generator(const Generator& gen) = delete;

	~Generator() {
		delete pm;
	}

	void generateEntities(int num_pu, int num_ss, int num_su, int num_pr_per_pu, float pr_range, const std::string& out_filename,
						std::vector<PU>* pus, std::vector<SS>* sss, std::vector<SU>* sus) const;

	void getEntitiesFromFile(int num_pu, int num_ss, int num_su, int num_pr_per_pu, const std::string& in_filename,
						std::vector<PU>* pus, std::vector<SS>* sss, std::vector<SU>* sus);

	void outputEntities(const std::string& out_filename, std::vector<PU>& pus, const std::vector<SS>& sss, const std::vector<SU>& sus) const;

	std::vector<float> computeGroundTruth(const std::vector<SU>& su, const std::vector<PU>& input_pus, PathLossTable* path_loss_table, bool no_pr_thresh_update) const;

private:
	float location_range;

	PropagationModel* pm;

	std::map<std::pair<std::pair<int, int>, int>, float> pr_pls;
};

class NullPropagationModel : public PropagationModel {
public:
	NullPropagationModel() {}

	float getPathLoss(const Location& loc1, const Location& loc2) const { return 0.0; }
};

class LogDistancePM : public PropagationModel {
public:
	LogDistancePM() = delete;
	LogDistancePM(float pl0_, float d0_, float gamma_) : PropagationModel(), pl0(pl0_), d0(d0_), gamma(gamma_) {}
	LogDistancePM(const LogDistancePM& ld) = delete;
	~LogDistancePM() {}

	float getPathLoss(const Location& loc1, const Location& loc2) const;

private:
	float pl0, d0, gamma;
};

class LongleyRicePM : public PropagationModel {
public:
	LongleyRicePM() = delete;
	LongleyRicePM(const std::string _splat_cmd, float _ref_lat, float _ref_long,
					const std::string _splat_dir, const std::string _sdf_dir, const std::string _return_dir,
					bool _use_itwom_pl)
			: PropagationModel(), splat_cmd(_splat_cmd), ref_lat(_ref_lat), ref_long(_ref_long), splat_dir(_splat_dir), sdf_dir(_sdf_dir), return_dir(_return_dir), use_itwom_pl(_use_itwom_pl) {
		if(splat_cmd != "splat" && splat_cmd != "splat-hd") {
			std::cerr << "LR command must be either 'splat' or 'splat-hd', got: " << splat_cmd << std::endl;
			exit(1);
		}
	}
	LongleyRicePM(const LongleyRicePM& lr) = delete;
	~LongleyRicePM() {}

	float getPathLoss(const Location& loc1, const Location& loc2) const;

private:
	std::string splat_cmd;

	float ref_lat, ref_long;

	std::string splat_dir, sdf_dir, return_dir;

	bool use_itwom_pl;
};

class TransmitterOnlySplatPM :public PropagationModel {
public:
	TransmitterOnlySplatPM() = delete;
	TransmitterOnlySplatPM(const std::string _splat_cmd, float _ref_lat, float _ref_long,
					const std::string _splat_dir, const std::string _sdf_dir, const std::string _return_dir)
			: PropagationModel(), splat_cmd(_splat_cmd), ref_lat(_ref_lat), ref_long(_ref_long), splat_dir(_splat_dir), sdf_dir(_sdf_dir), return_dir(_return_dir),
			k(10), alpha(4.0), gamma(0.25), cur_loc(0, 0), cur_ano_data()
	{
		if(splat_cmd != "splat") {
			std::cerr << "LR command must be 'splat', got: " << splat_cmd << std::endl;
			exit(1);
		}	
	}

	TransmitterOnlySplatPM(const LongleyRicePM& lr) = delete;
	~TransmitterOnlySplatPM() {}

	void preprocessPathLoss(PU* pu, float ss_height, int pu_id);
	void loadANOFile(const PU& pu);

	void preprocessPathLoss(PR* pr, float su_height, int pu_id, int pr_id);
	void loadANOFile(const PR& pr);

	float getPathLoss(const Location& loc1, const Location& loc2) const;
private:
	std::string splat_cmd;

	float ref_lat, ref_long;

	std::string splat_dir, sdf_dir, return_dir;

	int k;
	float alpha, gamma;

	Location cur_loc;
	std::vector<std::pair<Location, float> > cur_ano_data;
};

#endif
