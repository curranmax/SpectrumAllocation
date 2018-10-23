
#ifndef _GRID_TABLE_H_
#define _GRID_TABLE_H_

#include "spectrum_sensor.h"
#include "primary_user.h"

#include <map>
#include <vector>

class GridTable {
public:
	GridTable() : sss(), pu_refs() {}
	~GridTable() {}
	
	std::map<int, std::vector<SSint> > sss;
	std::map<int, std::vector<int> > pu_refs;
};

class PUTable {
public:
	PUTable() : pus(), num_pr_per_pu(1) {}
	~PUTable() {}
	
	std::map<int, PUint> pus;
	int num_pr_per_pu;
};

void buildTables(const std::map<int, std::vector<const SSint*> >& ss_groups,
					const std::map<int, std::vector<const PUint*> >& pu_groups,
					const std::vector<PUint>& pus, int grid_size,
					GridTable* grid_table,
					PUTable* pu_table);

#endif
