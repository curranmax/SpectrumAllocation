
#include "path_loss_table.h"

bool operator<(const PLTKey& a, const PLTKey& b) {
	if(a.su_ind != b.su_ind) {
		return a.su_ind < b.su_ind;
	}
	if(a.pu_ind != b.pu_ind) {
		return a.pu_ind < b.pu_ind;
	}
	return a.pr_ind < b.pr_ind;
}
