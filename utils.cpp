
#include "utils.h"

#include "ivory/Runtime/sInt.h"

#include <math.h>
#include <random>

using namespace osuCrypto;

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

float utils::todBm(float rp_in_mW) {
	return 10.0 * log10(rp_in_mW);
}

float utils::fromdBm(float rp_in_dBm) {
	return pow(10, rp_in_dBm / 10.0);
}

float utils::randomFloat(float min, float max) {
	return ((float) rand()) / ((float) RAND_MAX) * (max - min) + min;
}
