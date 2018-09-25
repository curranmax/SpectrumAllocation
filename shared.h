
#ifndef _SHARED_H_
#define _SHARED_H_

#include <mutex>
#include <vector>

class Shared {
public:
	Shared() : mtx(), float_set(false), vs_float(), int_set(false), vs_int() {}
	~Shared() {}

	void set(const std::vector<int>& vs);
	void get(std::vector<int>& vs);

	void set(const std::vector<float>& vs);
	void get(std::vector<float>& vs);

private:
	std::mutex mtx;

	bool float_set;
	std::vector<float> vs_float;

	bool int_set;
	std::vector<int> vs_int;
};

#endif
