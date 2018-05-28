
#ifndef _CUSTOM_IVORY_UTILS_H_
#define _CUSTOM_IVORY_UTILS_H_

#include "ivory/Runtime/sInt.h"

using namespace osuCrypto;


namespace utils {

enum UnitType {
	ABS,
	DB
};

void setUnitType(const std::string& unit_type_str);

extern UnitType unit_type;

// abs
// v is the value to calculate the abs on, the argument is modified in place.
void abs(sInt* v, const sInt& zero);

// sqrt
// guess should be set to an appropriate initial guess, and will be modified in place.
void sqrt(sInt* guess, const sInt& v, int num_iters, const sInt& two);

void dist(sInt* dist, const sInt& x1, const sInt& y1, const sInt& x2, const sInt& y2,
			const sInt& zero, const sInt& two, int num_iters);

float todBm(float rp_in_mW);
float fromdBm(float rp_in_dBm);

float randomFloat(float min, float max);

}

#endif
