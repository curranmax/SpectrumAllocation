
#ifndef _SECONDARY_USER_H_
#define _SECONDARY_USER_H_

#include "location.h"

// TODO this model doesn't work with current implementation of SM, since it is just two threads.
// Will eventually need to have this support network updates from the SM
class SecondaryUser {
public:
	SecondaryUser() : index(0), loc() {}
	SecondaryUser(int _index, const Location& _loc) : index(_index), loc(_loc) {}
	SecondaryUser(const SecondaryUser& su) : index(su.index), loc(su.loc) {}

	const SecondaryUser& operator=(const SecondaryUser& su) {
		this->index = su.index;
		this->loc = su.loc;
		return *this;
	}

	int index;

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
