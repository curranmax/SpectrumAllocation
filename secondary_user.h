
#ifndef _SECONDARY_USER_H_
#define _SECONDARY_USER_H_

#include "entity.h"
#include "location.h"

#include <iostream>
#include <vector>

// TODO this model doesn't work with current implementation of SM, since it is just two threads.
// Will eventually need to have this support network updates from the SM
class SecondaryUser {
public:
	SecondaryUser() : index(0), loc() {}
	SecondaryUser(int _index, const Location& _loc) : index(_index), loc(_loc), less_max_tp(0.0), min_tp(0.0) {}
	SecondaryUser(const SecondaryUser& su) : index(su.index), loc(su.loc), less_max_tp(su.less_max_tp), min_tp(su.min_tp), pl_id(su.pl_id) {}

	const SecondaryUser& operator=(const SecondaryUser& su) {
		this->index = su.index;
		this->loc = su.loc;
		this->less_max_tp = su.less_max_tp;
		this->min_tp = su.min_tp;
		this->pl_id = su.pl_id;
		return *this;
	}

	int index;

	Location loc;

	float less_max_tp, min_tp;

	int pl_id;
};

class SecondaryUserInt : public Entity {
public:
	SecondaryUserInt() : index(-1), loc(), global_id(-1) {}
	~SecondaryUserInt() {}
	SecondaryUserInt(const LocInt& _loc) : index(-1), loc(_loc), global_id(-1) {}
	SecondaryUserInt(const SecondaryUserInt& su) : index(su.index), loc(su.loc), global_id(su.global_id) {}

	const SecondaryUserInt& operator=(const SecondaryUserInt& su) {
		this->index = su.index;
		this->loc = su.loc;
		this->global_id = su.global_id;
		return *this;
	}

	std::string getType() const { return "SU"; }
	int getID() const { return global_id; }
	void setID(int new_id) { global_id = new_id; }

	std::vector<int> getValues() const {
		std::vector<int> values;
		values.push_back(index);

		values.push_back(loc.x);
		values.push_back(loc.y);

		values.push_back(global_id);
		return values;
	}

	void setValues(const std::vector<int>& values) {
		if(values.size() != 4) {
			std::cerr << "SU setValues error" << std::endl;
			exit(1);
		}

		index = values[0];

		loc.x = values[1];
		loc.y = values[2];

		global_id = values[3];
	}

	std::vector<int> getKsValues() const {
		std::vector<int> values;
		values.push_back(loc.x);
		values.push_back(loc.y);
		return values;
	}

	void setKsValues(const std::vector<int>& values) {
		if(values.size() != 2) {
			std::cerr << "SU setKsValues error" << std::endl;
			exit(1);
		}

		loc.x = values[0];
		loc.y = values[1];
	}
 
	int index;

	LocInt loc;

	int global_id;
};

typedef SecondaryUser SU;
typedef SecondaryUserInt SUint;

#endif
