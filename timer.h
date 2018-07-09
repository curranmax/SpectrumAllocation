
#ifndef _TIMER_H_
#define _TIMER_H_

#include <chrono>
#include <iostream>
#include <map>
#include <string>
#include <vector>

class Timer {
public:
	Timer() : current_tag(""), current_start(), durations() {}

	virtual void start(const std::string& tag);
	virtual void end(const std::string& tag);

	virtual float getAverageDuration(const std::string& tag) const;

	static const std::string secure_preprocessing;
	static const std::string secure_su_request;

	static const std::string plaintext_split_preprocessing;
	static const std::string plaintext_grid_preprocessing;
private:
	std::string current_tag;
	std::chrono::time_point<std::chrono::high_resolution_clock> current_start;

	std::map<std::string, std::vector<float> > durations;
};

class NullTimer: public Timer {
public:
	NullTimer() : Timer() {}

	void start(const std::string& tag) {}
	void end(const std::string& tag) {}

	float getAverageDuration(const std::string& tag) const { return 0.0; }
};

#endif