
#include "debug_print.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

const std::chrono::time_point<std::chrono::high_resolution_clock> start_test = std::chrono::high_resolution_clock::now();

std::string toString(std::chrono::duration<float> dur) {
	auto h =  std::chrono::duration_cast<std::chrono::hours>(dur);
	auto m =  std::chrono::duration_cast<std::chrono::minutes>(dur -= h);
	auto s =  std::chrono::duration_cast<std::chrono::seconds>(dur -= m);
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur -= s);

	std::stringstream sstr;
	sstr << h.count() << ":" << (m.count() < 10 ? "0" : "") << m.count() << ":" << (s.count() < 10 ? "0" : "") << s.count() << "." << (ms.count() < 10 ? "00" : (ms.count() < 100 ? "0" : "")) << ms.count();

	return sstr.str();
}

void debug::print(const std::string& filename, int line_number, const std::string& msg) {
	// Get time since last print
	auto now = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = now - start_test;

	std::cout << filename << ":" << line_number << " " << toString(elapsed) << " -> " << msg << std::endl;
}

void debug::print(const std::string& filename, int line_number, int party_id, const std::string& msg) {
	// Get time since last print
	auto now = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = now - start_test;

	std::cout << filename << ":" << line_number << " " << toString(elapsed) << " -> Party " << party_id << ": " << msg << std::endl;
}
