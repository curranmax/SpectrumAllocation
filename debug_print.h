
#ifndef _DEBUG_PRINT_H_
#define _DEBUG_PRINT_H_

#include <chrono>
#include <string>

#define P(msg) {} //debug::print(__FILE__, __LINE__, msg)
#define PWID(msg) {} //debug::print(__FILE__, __LINE__, party_id, msg) 

namespace debug {

void print(const std::string& filename, int line_number, const std::string& msg);
void print(const std::string& filename, int line_number, int party_id, const std::string &msg);

}

#endif
