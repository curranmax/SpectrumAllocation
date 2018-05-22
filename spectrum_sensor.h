
#ifndef _SPECTRUM_SENSOR_H_
#define _SPECTRUM_SENSOR_H_

#include "location.h"

class SpectrumSensor {
public:
	SpectrumSensor() : loc(), received_power(0.0) {}
	SpectrumSensor(const Location& _loc, float _received_power) : loc(_loc), received_power(_received_power) {}
	SpectrumSensor(const SpectrumSensor& ss) : loc(ss.loc), received_power(ss.received_power) {}

	const SpectrumSensor& operator=(const SpectrumSensor& ss) {
		this->loc = ss.loc;
		this->received_power = ss.received_power;
		return *this;
	}

	Location loc;

	// In mW
	float received_power;
};

class SpectrumSensorInt {
public:
	SpectrumSensorInt() : loc(), received_power(0.0) {}
	SpectrumSensorInt(const LocInt& _loc, int _received_power) : loc(_loc), received_power(_received_power) {}
	SpectrumSensorInt(const SpectrumSensorInt& ss) : loc(ss.loc), received_power(ss.received_power) {}

	const SpectrumSensorInt& operator=(const SpectrumSensorInt& ss) {
		this->loc = ss.loc;
		this->received_power = ss.received_power;
		return *this;
	}

	LocInt loc;

	// In mW
	int received_power;
};

typedef SpectrumSensor SS;
typedef SpectrumSensorInt SSint;

#endif
