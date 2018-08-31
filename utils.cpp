
#include "utils.h"

#include "ivory/Runtime/Party.h"
#include "ivory/Runtime/sInt.h"

#include <math.h>
#include <random>
#include <stdlib.h>

using namespace osuCrypto;

#define INPUT(parties, party_id, val, bit_count) parties[party_id].isLocalParty() ? parties[party_id].input<sInt>(val, bit_count) : parties[party_id].input<sInt>(bit_count)

#define POW_RAND_MIN -1.0
#define POW_RAND_MAX -1.0

#define LOG_RAND_MIN 0.5
#define LOG_RAND_MAX 2.0

utils::UnitType utils::unit_type = utils::UnitType::ABS;

void utils::setUnitType(const std::string& unit_type_str) {
	if(unit_type_str == "abs" || unit_type_str == "ABS") {
		std::cerr << "ABS unit_type is not fully implemented" << std::endl;
		exit(1);

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

void utils::secureLog10_v1(
		sInt* ans, const sInt& v,
		const sInt& zero, const sInt& factor_int, const sInt& ln_10,
		int num_iters) {
	*ans = zero + zero;

	sInt one = factor_int / factor_int;

	sInt ratio = factor_int * (v - factor_int) / (v + factor_int);
	sInt ratio_sq = ratio * ratio / factor_int;

	sInt k = ratio + zero;
	sInt d = one + zero;

	for(int i = 0; i < num_iters; ++i) {
		*ans = (*ans) + k / d;

		if(i < num_iters - 1) {
			k = k * ratio_sq / factor_int;
			d = d + one + one;
		}
	}

	*ans = (factor_int + factor_int) * (*ans) / ln_10;
}

void utils::secureLog10_v2(
		sInt* ans, const sInt& v,
		std::array<Party, 2>& parties, int bit_count,
		float factor, const sInt& factor_int) {
	float rv = randomFloat(LOG_RAND_MIN, LOG_RAND_MAX);
	sInt rv_secure = INPUT(parties, 0, int(rv * factor), bit_count);

	sInt v_times_rv_secure = v * rv_secure / factor_int;
	
	parties[1].reveal(v_times_rv_secure);
	float v_times_rv = 0.0;
	if(parties[1].isLocalParty()) {
		v_times_rv = float(v_times_rv_secure.getValue()) / factor;
	}
	parties[0].getRuntime().processesQueue();

	float log_rv = 0.0;
	float log_v_times_rv = 0.0;
	if(parties[0].isLocalParty()) {
		log_rv = log10(rv);
	}
	if(parties[1].isLocalParty()) {
		log_v_times_rv = log10(v_times_rv);
	}

	sInt log_rv_secure = INPUT(parties, 0, int(log_rv * factor), bit_count);
	sInt log_v_times_rv_secure = INPUT(parties, 0, int(log_v_times_rv * factor), bit_count);

	*ans = log_v_times_rv_secure - log_rv_secure;
}

void utils::securePow10(
		sInt* ans, const sInt& v,
		std::array<Party, 2>& parties, int bit_count,
		float factor, const sInt& factor_int) {
	// Party 0 generates a random value rv
	float rv = randomFloat(POW_RAND_MIN, POW_RAND_MAX);
	sInt rv_secure = INPUT(parties, 0, int(rv * factor), bit_count);

	// Party 1 gets v + rv
	sInt v_plus_rv_secure = v + rv_secure;

	parties[1].reveal(v_plus_rv_secure);
	float v_plus_rv = 0.0;
	if(parties[1].isLocalParty()) {
		v_plus_rv = float(v_plus_rv_secure.getValue()) / factor;
	}
	parties[0].getRuntime().processesQueue();

	// Party 0 calculates 10^rv, and Party 1 calculates 1o^(v + rv)
	float ten_to_rv = 0.0;
	float ten_to_v_plus_rv = 0.0;
	if(parties[0].isLocalParty()) {
		ten_to_rv = pow(10.0, rv);
	}
	if(parties[1].isLocalParty()) {
		ten_to_v_plus_rv = pow(10.0, v_plus_rv);
	}

	// Input values to S2-PC and calculate desired result
	sInt ten_to_rv_secure = INPUT(parties, 0, ten_to_rv, bit_count);
	sInt ten_to_v_plus_rv_secure = INPUT(parties, 1, ten_to_v_plus_rv, bit_count);

	*ans = factor_int * ten_to_rv_secure / ten_to_v_plus_rv_secure;
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

bool utils::deleteFile(const std::string& filename) {
	return remove(filename.c_str()) == 0;
}
