#include "Ant.h"

void Ant::mask_visit(const size_t idx, const vector<IdxValPt>& sz_p_pt)
{
	_visit_mask[idx] = true;
	size_t sz_id = sz_p_pt[idx].first;
	for (size_t i = idx - 1; i > 1; i--) {
		if (sz_id == sz_p_pt[i].first) { _visit_mask[i] = true; }
		else { break; }
	}
	for (size_t i = idx + 1; i < sz_p_pt.size(); i++) {
		if (sz_id == sz_p_pt[i].first) { _visit_mask[i] = true; }
		else { break; }
	}
}

void Ant::unmask_visit(const size_t idx, const vector<IdxValPt>& sz_p_pt)
{
	_visit_mask[idx] = false;
	size_t sz_id = sz_p_pt[idx].first;
	for (size_t i = idx - 1; i > 1; i--) {
		if (sz_id == sz_p_pt[i].first) { _visit_mask[i] = false; }
		else { break; }
	}
	for (size_t i = idx + 1; i < sz_p_pt.size(); i++) {
		if (sz_id == sz_p_pt[i].first) { _visit_mask[i] = false; }
		else { break; }
	}
}

void Ant::reset()
{
	cur_node = 0;  // start node
	cost = prize = 0;
	feasible = false;
	if (!path.empty()) { path.clear(); } path.reserve(_num_path);
	path.push_back(0);
	// Reset visit mask
	_visit_mask = vector<bool>(_num_nodes, false);
	_visit_mask[0] = _visit_mask[1] = true;
}
