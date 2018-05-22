
#ifndef _SECONDARY_USER_H_
#define _SECONDARY_USER_H_

#include "location.h"

// TODO this model doesn't work with current implementation of SM, since it is just two threads.
// Will eventually need to have this support network updates from the SM
class SecondaryUser {
public:
	SecondaryUser() : loc() {}
	SecondaryUser(const Location& _loc) : loc(_loc) {}
	SecondaryUser(const SecondaryUser& su) : loc(su.loc) {}

	const SecondaryUser& operator=(const SecondaryUser& su) {
		this->loc = su.loc;
		return *this;
	}

	Location loc;
};

class SecondaryUserInt {
public:
	SecondaryUserInt() : loc() {}
	SecondaryUserInt(const LocInt& _loc) : loc(_loc) {}
	SecondaryUserInt(const SecondaryUserInt& su) : loc(su.loc) {}

	const SecondaryUserInt& operator=(const SecondaryUserInt& su) {
		this->loc = su.loc;
		return *this;
	}

	LocInt loc;
};

typedef SecondaryUser SU;
typedef SecondaryUserInt SUint;

#endif
