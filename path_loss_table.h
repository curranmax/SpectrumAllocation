
#ifndef _PATH_LOSS_TABLE_H_
#define _PATH_LOSS_TABLE_H_

#include <map>

class PLTKey {
public:
	PLTKey() = delete;
	~PLTKey() {}
	PLTKey(int _su_ind, int _pu_ind, int _pr_ind) :
			su_ind(_su_ind), pu_ind(_pu_ind), pr_ind(_pr_ind) {}
	PLTKey(const PLTKey& key) :
			su_ind(key.su_ind), pu_ind(key.pu_ind), pr_ind(key.pr_ind) {}

	const PLTKey& operator=(const PLTKey& key) {
		su_ind = key.su_ind; pu_ind = key.pu_ind; pr_ind = key.pr_ind;
		return *this;
	}

	int su_ind, pu_ind, pr_ind;
};

bool operator<(const PLTKey& a, const PLTKey& b);

class PLTValue {
public:
	PLTValue() :
			sec_set(false), pt_set(false), gt_set(false),
			sec_pl(0.0), pt_pl(0.0), gt_pl(0.0) {}
	~PLTValue() {}
	PLTValue(const PLTValue& val) :
			sec_set(val.sec_set), pt_set(val.pt_set), gt_set(val.gt_set),
			sec_pl(val.sec_pl), pt_pl(val.pt_pl), gt_pl(val.gt_pl) {}
	const PLTValue& operator=(const PLTValue& val) {
		sec_set = val.sec_set; pt_set = val.pt_set; gt_set = val.gt_set;
		sec_pl = val.sec_pl; pt_pl = val.pt_pl; gt_pl = val.gt_pl;
		return *this;
	}

	bool setSecurePathLoss(float v) {
		if(sec_set) { return false; }
		sec_pl = v;
		sec_set = true;
		return true;
	}

	bool setPlaintextPathLoss(float v) {
		if(pt_set) { return false; }
		pt_pl = v;
		pt_set = true;
		return true;
	}

	bool setUnoptimizedPlaintextPathLoss(float v) {
		if(uo_pt_set) { return false; }
		uo_pt_pl = v;
		uo_pt_set = true;
		return true;
	}

	bool setGroundTruthPathLoss(float v) {
		if(gt_set) { return false; }
		gt_pl = v;
		gt_set = true;
		return true;
	}

	bool sec_set, pt_set, uo_pt_set, gt_set;
	float sec_pl, pt_pl, uo_pt_pl, gt_pl;
};

class PathLossTable {
public:
	PathLossTable() : table(), pu_table(), uo_mode(false) {}
	~PathLossTable() {}

	PathLossTable(const PathLossTable& plt) = delete;
	const PathLossTable& operator=(const PathLossTable& plt) = delete;
	
	// All pt path losses get added
	bool addPlaintextPathLoss(int su_ind, int pu_ind, int pr_ind, float path_loss) {
		if(uo_mode) {
			return this->addUnoptimizedPlaintextPathLoss(su_ind, pu_ind, pr_ind, path_loss);
		} else {
			return table[PLTKey(su_ind, pu_ind, pr_ind)].setPlaintextPathLoss(path_loss);
		}
	}

	bool addUnoptimizedPlaintextPathLoss(int su_ind, int pu_ind, int pr_ind, float path_loss) {
		PLTKey key(su_ind, pu_ind, pr_ind);
		if(table.count(key) <= 0) {
			return false;
		}
		return table[key].setUnoptimizedPlaintextPathLoss(path_loss);
	}

	// Only add gt path losses, if there is already an entry in the table from pt.
	bool addGroundTruthPathLoss(int su_ind, int pu_ind, int pr_ind, float path_loss) {
		PLTKey key(su_ind, pu_ind, pr_ind);
		if(table.count(key) <= 0) {
			return false;
		}
		return table[key].setGroundTruthPathLoss(path_loss);
	}

	// All pt path losses get added
	bool addPlaintextPathLoss(int su_ind, int pu_ind, float path_loss) {
		if(uo_mode) {
			return this->addUnoptimizedPlaintextPathLoss(su_ind, pu_ind, path_loss);
		} else {
			return pu_table[PLTKey(su_ind, pu_ind, -1)].setPlaintextPathLoss(path_loss);
		}
	}

	bool addUnoptimizedPlaintextPathLoss(int su_ind, int pu_ind, float path_loss) {
		PLTKey key(su_ind, pu_ind, -1);
		if(pu_table.count(key) <= 0) {
			return false;
		}
		return pu_table[key].setUnoptimizedPlaintextPathLoss(path_loss);
	}

	// Only add gt path losses, if there is already an entry in the table from pt.
	bool addGroundTruthPathLoss(int su_ind, int pu_ind, float path_loss) {
		PLTKey key(su_ind, pu_ind, -1);
		if(pu_table.count(key) <= 0) {
			return false;
		}
		return pu_table[key].setGroundTruthPathLoss(path_loss);
	}

	// (SU, PR) -> (pt path loss, gt path loss)
	std::map<PLTKey, PLTValue> table;
	std::map<PLTKey, PLTValue> pu_table;

	bool uo_mode;
};

#endif
