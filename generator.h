
#ifndef _GENERATOR_H_
#define _GENERATOR_H_

#include "location.h"
#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_sensor.h"

#include <string>
#include <vector>

class PropagationModel {
public:
	virtual float getPathLossRatio(const Location& loc1, const Location& loc2) const = 0;
	virtual ~PropagationModel() {}
};

class Generator {
public:
	Generator() = delete;
	Generator(float location_range_, PropagationModel* pm_) : location_range(location_range_), pm(pm_) {}
	Generator(const Generator& gen) = delete;

	~Generator() {
		delete pm;
	}

	void generateEntities(int num_pu, int num_ss, int num_su,
						std::vector<PU>* pus, std::vector<SS>* sss, std::vector<SU>* sus) const;

	float computeGroundTruth(const SU& su, const std::vector<PU>& pus) const;

private:
	float location_range;

	PropagationModel* pm;
};

class LogDistancePM : public PropagationModel {
public:
	LogDistancePM() = delete;
	LogDistancePM(float pl0_, float d0_, float gamma_) : PropagationModel(), pl0(pl0_), d0(d0_), gamma(gamma_) {}
	LogDistancePM(const LogDistancePM& ld) = delete;
	~LogDistancePM() {}

	float getPathLossRatio(const Location& loc1, const Location& loc2) const;

private:
	float pl0, d0, gamma;
};

class LongleyRicePM : public PropagationModel {
public:
	LongleyRicePM() = delete;
	LongleyRicePM(int x) : y(x) {}
	LongleyRicePM(const LongleyRicePM& lr) = delete;
	~LongleyRicePM() {}

	float getPathLossRatio(const Location& loc1, const Location& loc2) const;

private:
	int y;
};

#endif
