
#include "timer.h"

#include <chrono>
#include <map>
#include <stdlib.h>
#include <string>

const std::string Timer::secure_preprocessing = "secure_preprocessing";
const std::string Timer::secure_su_request = "secure_su_request";

const std::string Timer::secure_write = "secure_write";

const std::string Timer::plaintext_split_preprocessing = "plaintext_split_preprocessing";
const std::string Timer::plaintext_grid_preprocessing = "plaintext_grid_preprocessing";

void Timer::start(const std::string& tag) {
	if(current_tag != "") {
		std::cerr << "Trying to start timer with mismatched tags: (" << tag << ", " << current_tag << ")" << std::endl;
		exit(0);
	}

	current_start = std::chrono::high_resolution_clock::now();
	current_tag = tag;
}

void Timer::end(const std::string& tag) {
	if(current_tag != tag) {
		std::cerr << "Trying to end timer with mismatched tags: (" << tag << ", " << current_tag << ")" << std::endl;
		exit(0);
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = end - current_start;

	durations[tag].push_back(elapsed.count());

	current_tag = "";
}

float Timer::getAverageDuration(const std::string& tag) const {
	auto itr = durations.find(tag);
	if(itr == durations.cend()) {
		return 0.0;
	}

	float sum_durs = 0.0;
	int num_durs = 0;
	for(unsigned int i = 0; i < itr->second.size(); ++i) {
		num_durs++;
		sum_durs += itr->second[i];
	}

	if(num_durs == 0) {
		return 0.0;
	}
	return sum_durs / num_durs;
}

int Timer::numDurations(const std::string& tag) const {
	auto itr = durations.find(tag);
	if(itr == durations.cend()) {
		return 0;
	}

	return itr->second.size();
}
