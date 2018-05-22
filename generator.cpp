
#include "generator.h"

#include "itm.hpp"

#include "location.h"
#include "primary_user.h"
#include "secondary_user.h"
#include "spectrum_sensor.h"
#include "utils.h"

#include <string>

void Generator::generateEntities(
		int num_pu, int num_ss, int num_su,
		std::vector<PU>* pus, std::vector<SS>* sss, std::vector<SU>* sus) const {
	// PU
	for(int j = 0; j < num_pu; ++j) {
		float pu_x = utils::randomFloat(0.0, location_range);
		float pu_y = utils::randomFloat(0.0, location_range);
	
		// Both in mW.
		float transmit_power = 2.0;
		float threshold = 1.0;
		PU pu(Location(pu_x, pu_y), transmit_power, threshold);
		(*pus).push_back(pu);
	}

	// SS
	for(int i = 0; i < num_ss; ++i) {
		float x = utils::randomFloat(0.0, location_range);
		float y = utils::randomFloat(0.0, location_range);
		Location loc(x, y);

		float rp = 0.0;
		for(int j = 0; j < num_pu; ++j) {
			float pl_ratio = pm->getPathLossRatio((*pus)[j].loc, loc);
			rp += pl_ratio * (*pus)[j].transmit_power;
		}

		SS ss(loc, rp);
		(*sss).push_back(ss);
	}

	// SU
	for(int i = 0; i < num_su; ++i) {
		float x = utils::randomFloat(location_range / 5.0, 4.0 * location_range / 5.0);
		float y = utils::randomFloat(0.0, location_range);
		SU su(Location(x, y));
		(*sus).push_back(su);
	}
}

float Generator::computeGroundTruth(const SU& su, const std::vector<PU>& pus) const {
	float tp = 0.0;
	bool first = true;
	for(unsigned int i = 0; i < pus.size(); ++i) {
		for(unsigned int j = 0; j < pus[i].prs.size(); ++j) {
			float this_tp = pus[i].prs[j].threshold / pm->getPathLossRatio(su.loc, pus[i].prs[j].loc);
			if(first) {
				tp = this_tp;
				first = false;
			} else if(this_tp < tp) {
				tp = this_tp;
			}
		}
	}
	return tp;
}

float LogDistancePM::getPathLossRatio(const Location& loc1, const Location& loc2) const {
	float pl_dBm = pl0 + 10 * gamma * log10(loc1.dist(loc2) / d0);
	return utils::fromdBm(-pl_dBm);
}

float LongleyRicePM::getPathLossRatio(const Location& loc1, const Location& loc2) const {

	double elev[11] = {8, 10, 10, 20, 30, 40, 50, 40, 30, 20, 10};
	double height1 = 100;
	double height2 = 10;
	double eps_dielect = 10.0;
	double sgm_conductivity = 1.0;
	double eno_ns_surfref = 0.0;
	double frq_mhz = 2000.0;
	int radio_climate = 1;
	int pol = 0;
	double conf = 0.5;
	double rel = 0.5;
	double dbloss = 0.0;
	char strmode[256];
	int errnum = 0;

	point_to_point(elev, height1, height2, eps_dielect, sgm_conductivity,
					eno_ns_surfref, frq_mhz, radio_climate, pol, conf, rel,
					dbloss, strmode, errnum);

	std::cout << dbloss << " | " << strmode << " | " << errnum << std::endl;
	exit(0);
	return utils::fromdBm(-dbloss);
}

