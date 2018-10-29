
#ifndef _SPECTRUM_SENSOR_H_
#define _SPECTRUM_SENSOR_H_

#include "entity.h"
#include "location.h"

#include <iostream>
#include <vector>

class SpectrumSensor {
public:
	SpectrumSensor() : loc(), received_power(0.0) {}
	SpectrumSensor(const Location& _loc, float _received_power) : loc(_loc), received_power(_received_power) {}
	SpectrumSensor(const SpectrumSensor& ss) : loc(ss.loc), received_power(ss.received_power), pl_id(ss.pl_id) {}

	const SpectrumSensor& operator=(const SpectrumSensor& ss) {
		this->loc = ss.loc;
		this->received_power = ss.received_power;
		this->pl_id = ss.pl_id;
		return *this;
	}

	Location loc;

	float received_power;

	int pl_id;
};

class SpectrumSensorInt : public Entity {
public:
	SpectrumSensorInt()
			: loc(), received_power(0.0), received_power_from_pu(), global_id(-1) {}
	~SpectrumSensorInt() {}
	SpectrumSensorInt(const LocInt& _loc, int _received_power)
			: loc(_loc), received_power(_received_power), received_power_from_pu(), global_id(-1) {}
	SpectrumSensorInt(const SpectrumSensorInt& ss)
			: loc(ss.loc), received_power(ss.received_power),
			received_power_from_pu(ss.received_power_from_pu), global_id(ss.global_id) {}

	const SpectrumSensorInt& operator=(const SpectrumSensorInt& ss) {
		this->loc = ss.loc;
		this->received_power = ss.received_power;
		this->received_power_from_pu = ss.received_power_from_pu;
		this->global_id = global_id;
		return *this;
	}

	std::string getType() const { return "SS"; }
	int getID() const { return global_id; }
	void setID(int new_id) { global_id = new_id; }

	bool fullEncrypt() const { return true; }
	
	std::vector<int> getValues() const {
		std::vector<int> values;

		// Location
		values.push_back(loc.x);
		values.push_back(loc.y);

		// RP
		values.push_back(received_power);

		// Split RP
		for(unsigned int i = 0; i < received_power_from_pu.size(); ++i) {
			values.push_back(received_power_from_pu[i]);
		}

		// global_id
		values.push_back(global_id);

		return values;
	}

	void setValues(const std::vector<int>& values) {
		loc.x = values[0];
		loc.y = values[1];

		received_power = values[2];

		received_power_from_pu = std::vector<int>(values.size() - 4);
		for(unsigned int i = 3; i < values.size() - 1; ++i) {
			received_power_from_pu[i - 3] = values[i];
		}
		global_id = values[values.size() - 1];
	}

	std::vector<int> getKsValues() const {
		std::vector<int> values;

		// Location
		values.push_back(loc.x);
		values.push_back(loc.y);

		// RP
		values.push_back(received_power);

		// Split RP
		for(unsigned int i = 0; i < received_power_from_pu.size(); ++i) {
			values.push_back(received_power_from_pu[i]);
		}

		return values;
	}

	void setKsValues(const std::vector<int>& values) {
		loc.x = values[0];
		loc.y = values[1];

		received_power = values[2];

		received_power_from_pu = std::vector<int>(values.size() - 3);
		for(unsigned int i = 3; i < values.size(); ++i) {
			received_power_from_pu[i - 3] = values[i];
		}
	}

	LocInt loc;

	int received_power;

	std::vector<int> received_power_from_pu;

	int global_id;
};

typedef SpectrumSensor SS;
typedef SpectrumSensorInt SSint;

#endif
