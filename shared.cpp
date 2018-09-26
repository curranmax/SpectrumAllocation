
#include "shared.h"

void Shared::set(const std::vector<int>& vs) {
	// TODO not a busy wait or switch to a queue of 
	while(int_set) {}

	mtx.lock();
	vs_int = vs;	

	int_set = true;
	mtx.unlock();
}

void Shared::get(std::vector<int>& vs) {
	while(!int_set) {}

	mtx.lock();
	vs = vs_int;

	int_set = false;
	mtx.unlock();
}

void Shared::set(const std::vector<float>& vs) {
	// TODO not a busy wait or switch to a queue of 
	while(float_set) {}

	mtx.lock();
	vs_float = vs;	

	float_set = true;
	mtx.unlock();
}

void Shared::get(std::vector<float>& vs) {
	while(!float_set) {}

	mtx.lock();
	vs = vs_float;

	float_set = false;
	mtx.unlock();
}

void Shared::set(const GridTable* gt) {
	while(gt_set) {}

	mtx.lock();
	vs_gt = gt;

	gt_set = true;
	mtx.unlock();
}

void Shared::get(GridTable& gt) {
	while(!gt_set) {}

	mtx.lock();
	gt = *vs_gt;
	vs_gt = nullptr;

	gt_set = false;
	mtx.unlock();
}
