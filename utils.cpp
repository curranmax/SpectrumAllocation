
#include "utils.h"

#include "ivory/Runtime/sInt.h"

#include <math.h>
#include <random>
#include <stdlib.h>

using namespace osuCrypto;

utils::UnitType utils::unit_type = utils::UnitType::ABS;

void utils::setUnitType(const std::string& unit_type_str) {
	if(unit_type_str == "abs" || unit_type_str == "ABS") {
		utils::unit_type = UnitType::ABS;
	} else if(unit_type_str == "db" || unit_type_str == "dB" || unit_type_str == "DB") {
		utils::unit_type = UnitType::DB;
	} else {
		std::cerr << "Invalid unit_type_str: " << unit_type_str << std::endl;
		exit(1);
	}
}

void utils::abs(sInt* v, const sInt& zero) {
	sInt is_positive = *v >= zero;
	*v = is_positive.ifelse(*v, zero - *v);
}

void utils::sqrt(sInt* guess, const sInt& v, int num_iters, const sInt& two) {
	for(int i = 0; i < num_iters; ++i) {
		*guess = (*guess + v / *guess) / two;
	}
}

void utils::dist(sInt* dist, const sInt& x1, const sInt& y1, const sInt& x2, const sInt& y2,
			const sInt& zero, const sInt& two, int num_iters) {
	// perform some computation
	sInt xv = (x1 - x2);
	sInt yv = (y1 - y2);
	utils::abs(&xv, zero);
	utils::abs(&yv, zero);

	// Initialize dist to a good initial guess, then calculate the sqrt.
	*dist = xv + yv;
	utils::sqrt(dist, xv * xv + yv * yv,
				num_iters, two);
}

void utils::secureLog10(
		sInt* ans, const sInt& v,
		const sInt& zero, const sInt& factor_int, const sInt& ln_10,
		int num_iters) {
	*ans = zero + zero;

	sInt num   = v - factor_int;
	sInt denom = factor_int + zero;
	sInt sign  = factor_int + zero;

	for(int i = 0; i < num_iters; ++i) {
		*ans = (*ans) + sign * num / denom;

		if(i < num_iters - 1) {
			num   = num * (v - factor_int) / factor_int;
			denom = denom + factor_int;
			sign  = zero - sign;
		}
	}

	*ans = factor_int * (*ans) / ln_10;
}

float utils::todBm(float rp_in_mW) {
	return 10.0 * log10(rp_in_mW);
}

float utils::fromdBm(float rp_in_dBm) {
	return pow(10, rp_in_dBm / 10.0);
}

float utils::randomFloat(float min, float max) {
	return ((float) rand()) / ((float) RAND_MAX) * (max - min) + min;
}
