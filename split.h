
#ifndef _SPLIT_H_
#define _SPLIT_H_

#include "location.h"
#include "primary_user.h"
#include "spectrum_sensor.h"
#include "secondary_user.h"

#include <utility>

std::pair<float, float> splitFloat(float v);
std::pair<int, int> splitInt(int v);
int splitRandVal();

std::pair<Location, Location> splitLocation(const Location& loc);
std::pair<LocInt, LocInt> splitLocInt(const Location& loc, float factor);

std::pair<PRint, PRint> splitPrimaryReceiverInt(const PR& pr, float factor);
std::pair<PUint, PUint> splitPrimaryUserInt(const PU& pu, float factor);

std::pair<SSint, SSint> splitSpectrumSensorInt(const SS& ss, float factor);

std::pair<SUint, SUint> splitSecondaryUserInt(const SU& su, float factor);

#endif
