
#ifndef _PRIMARY_USER_H_
#define _PRIMARY_USER_H_

#include "entity.h"
#include "location.h"

#include <iostream>
#include <string>
#include <vector>

class PrimaryReceiver {
public:
	PrimaryReceiver() : loc(), threshold(0.0), splat_ano_filename("") {}
	PrimaryReceiver(const Location& _loc, float _threshold) : loc(_loc), threshold(_threshold), splat_ano_filename("") {}
	PrimaryReceiver(const PrimaryReceiver& pr) : loc(pr.loc), threshold(pr.threshold), splat_ano_filename(pr.splat_ano_filename) {}

	const PrimaryReceiver& operator=(const PrimaryReceiver& pr) {
		this->loc = pr.loc;
		this->threshold = pr.threshold;
		this->splat_ano_filename = pr.splat_ano_filename;
		return *this;
	}

	Location loc;

	// In mW.
	float threshold;

	std::string splat_ano_filename;
};

class PrimaryReceiverInt : public Entity {
public:
	PrimaryReceiverInt() : id(-1), loc(), threshold(0), global_id(-1) {}
	~PrimaryReceiverInt() {}
	PrimaryReceiverInt(const LocInt& _loc, int _threshold)
			: id(-1), loc(_loc), threshold(_threshold), global_id(-1) {}
	PrimaryReceiverInt(const PrimaryReceiverInt& pr)
			: id(pr.id), loc(pr.loc), threshold(pr.threshold), global_id(pr.global_id) {}

	const PrimaryReceiverInt& operator=(const PrimaryReceiverInt& pr) {
		this->id = pr.id;
		this->loc = pr.loc;
		this->threshold = pr.threshold;
		this->global_id = pr.global_id;
		return *this;
	}

	std::string getType() const { return "PR"; }
	int getID() const { return global_id; }
	void setID(int new_id) { global_id = new_id;}

	std::vector<int> getValues() const {
		std::vector<int> values;

		// ID
		values.push_back(id);

		// Loc
		values.push_back(loc.x);
		values.push_back(loc.y);

		// threshold
		values.push_back(threshold);

		// global_id
		values.push_back(global_id);
		return values;
	}

	void setValues(const std::vector<int>& values) {
		if(values.size() != 5) {
			std::cerr << "PR setValues() error" << std::endl;
			exit(1);
		}

		id = values[0];

		loc.x = values[1];
		loc.y = values[2];

		threshold = values[3];

		global_id = values[4];
	}

	std::vector<int> getKsValues() const {
		std::vector<int> values;
		// Loc
		values.push_back(loc.x);
		values.push_back(loc.y);

		// threshold
		values.push_back(threshold);
		return values;
	}

	void setKsValues(const std::vector<int>& values) {
		if(values.size() != 3) {
			std::cerr << "PR setValues() error" << std::endl;
			exit(1);
		}

		loc.x = values[0];
		loc.y = values[1];

		threshold = values[2];
	}

	int id;

	LocInt loc;

	// In mW.
	int threshold;

	int global_id;
};

class PrimaryUser {
public:
	PrimaryUser() : loc(), transmit_power(0.0), prs(), splat_ano_filename("") {}
	PrimaryUser(const Location& _loc, float _transmit_power) : loc(_loc), transmit_power(_transmit_power), prs(), splat_ano_filename("") {}
	PrimaryUser(const Location& _loc, float _transmit_power, float _threshold) :
			loc(_loc), transmit_power(_transmit_power),
			prs({PrimaryReceiver(_loc, _threshold)}),
			splat_ano_filename("") {}
	PrimaryUser(const PrimaryUser& pu) : loc(pu.loc), transmit_power(pu.transmit_power), prs(pu.prs), splat_ano_filename(pu.splat_ano_filename) {}

	const PrimaryUser& operator=(const PrimaryUser& pu) {
		this->loc = pu.loc;
		this->transmit_power = pu.transmit_power;
		this->prs = pu.prs;
		this->splat_ano_filename = pu.splat_ano_filename;
		return *this;
	}

	Location loc;

	// In mW.
	float transmit_power;

	std::vector<PrimaryReceiver> prs;

	std::string splat_ano_filename;
};

class PrimaryUserInt : public Entity {
public:
	PrimaryUserInt() : index(-1), loc(), transmit_power(0 ), prs(), global_id(-1) {}
	PrimaryUserInt(const LocInt& _loc, int _transmit_power)
			: index(-1), loc(_loc), transmit_power(_transmit_power), prs(), global_id(-1) {}
	PrimaryUserInt(const LocInt& _loc, int _transmit_power, int _threshold) :
			index(-1), loc(_loc), transmit_power(_transmit_power),
			prs({PrimaryReceiverInt(_loc, _threshold)}) {}
	PrimaryUserInt(const PrimaryUserInt& pu)
		: index(pu.index), loc(pu.loc), transmit_power(pu.transmit_power), prs(pu.prs), global_id(pu.global_id) {}

	const PrimaryUserInt& operator=(const PrimaryUserInt& pu) {
		this->index = pu.index;
		this->loc = pu.loc;
		this->transmit_power = pu.transmit_power;
		this->prs = pu.prs;
		this->global_id = pu.global_id;
		return *this;
	}

	std::string getType() const { return "PU"; }
	int getInd() const { return index; }
	void setInd(int _index) { index = _index; }

	int getID() const { return global_id; }
	void setID(int new_id) { global_id = new_id;}

	std::vector<int> getValues() const {
		std::vector<int> values;
		
		// Index
		values.push_back(index);

		// Loc
		values.push_back(loc.x);
		values.push_back(loc.y);

		// tp
		values.push_back(transmit_power);

		values.push_back(global_id);
		return values;
	}

	void setValues(const std::vector<int>& values) {
		if(values.size() != 5) {
			std::cerr << "PU setValues() error" << std::endl;
			exit(1);
		}

		// Index
		index = values[0];

		// Loc
		loc.x = values[1];
		loc.y = values[2];

		// tp
		transmit_power = values[3];

		// global_id
		global_id = values[4];
	}

	std::vector<int> getKsValues() const {
		std::vector<int> values;

		// Loc
		values.push_back(loc.x);
		values.push_back(loc.y);

		// tp
		values.push_back(transmit_power);
		return values;
	}

	void setKsValues(const std::vector<int>& values) {
		if(values.size() != 3) {
			std::cerr << "PU setValues() error" << std::endl;
			exit(1);
		}

		// Loc
		loc.x = values[0];
		loc.y = values[1];

		// tp
		transmit_power = values[2];
	}

	int index;

	LocInt loc;

	// In mW.
	int transmit_power;

	std::vector<PrimaryReceiverInt> prs;

	int global_id;
};

typedef PrimaryReceiver PR;
typedef PrimaryReceiverInt PRint;

typedef PrimaryUser PU;
typedef PrimaryUserInt PUint;

#endif
