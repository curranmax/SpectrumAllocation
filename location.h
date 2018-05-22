
#ifndef _LOCATION_H_
#define _LOCATION_H_

class Location {
public:
	Location() : x(0.0), y(0.0) {}
	Location(float _x, float _y) : x(_x), y(_y) {}
	Location(const Location& loc) : x(loc.x), y(loc.y) {}

	const Location& operator=(const Location& loc) {
		this->x = loc.x;
		this->y = loc.y;
		return *this;
	}

	float dist(const Location& other) const;

	float x, y;
};

class LocInt {
public:
	LocInt() : x(0), y(0) {}
	LocInt(int _x, int _y) : x(_x), y(_y) {}
	LocInt(const LocInt& loc) : x(loc.x), y(loc.y) {}

	const LocInt& operator=(const LocInt& loc) {
		this->x = loc.x;
		this->y = loc.y;
		return *this;
	}

	int x, y;
};

#endif
