
#include "location.h"

#include <math.h>

float Location::dist(const Location& other) const {
	return sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
}
