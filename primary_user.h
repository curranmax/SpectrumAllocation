
#ifndef _PRIMARY_USER_H_
#define _PRIMARY_USER_H_

#include "location.h"

#include <vector>

class PrimaryReceiver {
public:
	PrimaryReceiver() : loc(), threshold(0.0) {}
	PrimaryReceiver(const Location& _loc, float _threshold) : loc(_loc), threshold(_threshold) {}
	PrimaryReceiver(const PrimaryReceiver& pr) : loc(pr.loc), threshold(pr.threshold) {}

	const PrimaryReceiver& operator=(const PrimaryReceiver& pr) {
		this->loc = pr.loc;
		this->threshold = pr.threshold;
		return *this;
	}

	Location loc;

	// In mW.
	float threshold;
};

class PrimaryReceiverInt {
public:
	PrimaryReceiverInt() : loc(), threshold(0) {}
	PrimaryReceiverInt(const LocInt& _loc, int _threshold) : loc(_loc), threshold(_threshold) {}
	PrimaryReceiverInt(const PrimaryReceiverInt& pr) : loc(pr.loc), threshold(pr.threshold) {}

	const PrimaryReceiverInt& operator=(const PrimaryReceiverInt& pr) {
		this->loc = pr.loc;
		this->threshold = pr.threshold;
		return *this;
	}

	int id;

	LocInt loc;

	// In mW.
	int threshold;
};

class PrimaryUser {
public:
	PrimaryUser() : loc(), transmit_power(0.0), prs() {}
	PrimaryUser(const Location& _loc, float _transmit_power) : loc(_loc), transmit_power(_transmit_power), prs() {}
	PrimaryUser(const Location& _loc, float _transmit_power, float _threshold) :
			loc(_loc), transmit_power(_transmit_power),
			prs({PrimaryReceiver(_loc, _threshold)}) {}
	PrimaryUser(const PrimaryUser& pu) : loc(pu.loc), transmit_power(pu.transmit_power), prs(pu.prs) {}

	const PrimaryUser& operator=(const PrimaryUser& pu) {
		this->loc = pu.loc;
		this->transmit_power = pu.transmit_power;
		this->prs = pu.prs;
		return *this;
	}

	Location loc;

	// In mW.
	float transmit_power;

	std::vector<PrimaryReceiver> prs;
};

class PrimaryUserInt {
public:
	PrimaryUserInt() : loc(), transmit_power(0 ), prs() {}
	PrimaryUserInt(const LocInt& _loc, int _transmit_power) : loc(_loc), transmit_power(_transmit_power), prs() {}
	PrimaryUserInt(const LocInt& _loc, int _transmit_power, int _threshold) :
			loc(_loc), transmit_power(_transmit_power),
			prs({PrimaryReceiverInt(_loc, _threshold)}) {}
	PrimaryUserInt(const PrimaryUserInt& pu) : loc(pu.loc), transmit_power(pu.transmit_power), prs(pu.prs) {}

	const PrimaryUserInt& operator=(const PrimaryUserInt& pu) {
		this->loc = pu.loc;
		this->transmit_power = pu.transmit_power;
		this->prs = pu.prs;
		return *this;
	}

	LocInt loc;

	// In mW.
	int transmit_power;

	std::vector<PrimaryReceiverInt> prs;
};

typedef PrimaryReceiver PR;
typedef PrimaryReceiverInt PRint;

typedef PrimaryUser PU;
typedef PrimaryUserInt PUint;

#endif
