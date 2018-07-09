
#include "split.h"

#include <stdlib.h>
#include <utility>

std::pair<float, float> splitFloat(float v) {
	// TODO
	float a = v;
	float b = 0.0;
	return std::make_pair(a, b);
}

std::pair<int, int> splitInt(int v) {
	// TODO might be a better way to do this
	int a = splitRandVal();
	int b = v - a;
	return std::make_pair(a, b);
}

int splitRandVal() {
	return rand() % (1 << 16) - (1 << 15);
}

std::pair<Location, Location> splitLocation(const Location& loc) {
	auto x_split = splitFloat(loc.x);
	auto y_split = splitFloat(loc.y);

	return std::make_pair(Location(x_split.first, y_split.first),
							Location(x_split.second, y_split.second));
}

std::pair<LocInt, LocInt> splitLocInt(const Location& loc, float factor) {
	auto x_split = splitInt(int(loc.x * factor));
	auto y_split = splitInt(int(loc.y * factor));

	return std::make_pair(LocInt(x_split.first, y_split.first),
							LocInt(x_split.second, y_split.second));
}

std::pair<PRint, PRint> splitPrimaryReceiverInt(const PR& pr, float factor) {
	auto loc_split = splitLocInt(pr.loc, factor);
	auto thresh_split = splitInt(int(pr.threshold * factor));

	return std::make_pair(PRint(loc_split.first, thresh_split.first),
							PRint(loc_split.second, thresh_split.second));
}

std::pair<PUint, PUint> splitPrimaryUserInt(const PU& pu, float factor) {
	auto loc_split = splitLocInt(pu.loc, factor);
	auto tp_split = splitInt(int(pu.transmit_power * factor));

	auto pu_split = std::make_pair(PUint(loc_split.first, tp_split.first),
									PUint(loc_split.second, tp_split.second));

	for(unsigned int i = 0; i < pu.prs.size(); ++i) {
		auto pr_split = splitPrimaryReceiverInt(pu.prs[i], factor);
		pu_split.first.prs.push_back(pr_split.first);
		pu_split.second.prs.push_back(pr_split.second);
	}

	return pu_split;
}

std::pair<SSint, SSint> splitSpectrumSensorInt(const SS& ss, float factor) {
	auto loc_split = splitLocInt(ss.loc, factor);
	auto rp_split = splitInt(int(ss.received_power * factor));

	return std::make_pair(SSint(loc_split.first, rp_split.first),
							SSint(loc_split.second, rp_split.second));
}

std::pair<SUint, SUint> splitSecondaryUserInt(const SU& su, float factor) {
	auto loc_split = splitLocInt(su.loc, factor);
	
	return std::make_pair(SUint(loc_split.first),
							SUint(loc_split.second));
}
