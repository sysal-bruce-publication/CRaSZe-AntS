#include "StnZoneSrch.h"

void StnZoneSrch::_calc_dist_mat()
{
	_dists = vector<vector<double>>(_num_circs, vector<double>(_num_circs));
	for (size_t i = 0; i < _num_circs - 1; ++i) {
		for (size_t j = i + 1; j < _num_circs; ++j) {
			_dists[i][j] = _dists[j][i] = (_cs[i] - _cs[j]).norm2();
		}
	}
}

void StnZoneSrch::init(const bool shuffle_circs)
{
	if (shuffle_circs) {
		unique_ptr<Rd> rd = make_unique<Rd>();
		rd->vec_shuffle_skip_first_n(_num_depots, _cs);
	}
	// Make circle id align with its index in _cs.
	for (size_t i = 0; i < _num_circs; i++) { _cs[i].set_id(i); }
	_calc_dist_mat();
	// Clear SZ vector, here start node is excluded
	if (!_szs.empty()) { _szs.clear(); } _szs.reserve(_num_circs - _num_depots);
	// Clear indices of available circles
	_feas_idx = vector<size_t>(_num_circs - _num_depots);
	std::iota(_feas_idx.begin(), _feas_idx.end(), _num_depots);
}

void StnZoneSrch::construct_stn_zones()
{
	vector<Point> intx_pts; vector<DualIdx> cc_ids;
	size_t sz_cnt = _num_depots; // Note that Steiner Zone vector won't save the start node.
	for (size_t i = _num_depots; i < _num_circs; ++i) {
		// Get information of first circle to be added in a SZ, and its index in _cs.
		Circle c1 = _cs[i];
		if (std::find(_feas_idx.begin(), _feas_idx.end(), i) != _feas_idx.end()) {
			// First remove this circle from available circles
			Common::vec_erase(_feas_idx, i);
			// Initialize the Steiner Zone
			StnZone cur_sz = StnZone(sz_cnt, _max_circs);
			cur_sz.push_circ(c1); sz_cnt++;
			// Clear intersection points and circle IDs from last iteration
			if (!intx_pts.empty()) { intx_pts.clear(); }
			// Any SZ with degree K has K intersection points and an additional
			// intersection point to be removed.
			intx_pts.reserve(_max_circs + 1);
			// Same for cc_idxs.
			if (!cc_ids.empty()) { cc_ids.clear(); }
			cc_ids.reserve(_max_circs + 1);
			// Then find all xs involved within the circle's interval, i.e.,
			// (x-, x+), NOTICE all xs within the interval should be left extreme
			// point, i.e., x-, of RHS circles.
			if (!_feas_idx.empty()) {
				// Any candidate circles successfully added, remove it from _feas_cs. 
				auto rm_cs = remove_if(_feas_idx.begin(), _feas_idx.end(),
					[this, &cur_sz, &intx_pts, &cc_ids](const size_t idx) -> bool
					{  // First check if two circles are close.
						if (!cur_sz.is_full()) {
							if (cur_sz.c1().is_close(_cs[idx])) {
								// If their rectangular overlaps
								if (cur_sz.check_necessity(_dists, _cs[idx])) {
									// Check necessity
									bool is_ovlp = false;
									if (cur_sz.check_sufficiency(_dists, _cs[idx],
										is_ovlp, intx_pts, cc_ids)) {
										// Check the sufficiency
										if (is_ovlp) { cur_sz.push_ovlp_circ(_cs[idx]); }
										else { cur_sz.push_circ(_cs[idx]); }
										return true;
									}
								}
							}
						}
						return false;
					});
				_feas_idx.erase(rm_cs, _feas_idx.end());
			}
			// If there exists intersection points, push them to the SZ.
			if (!intx_pts.empty()) {
				for (size_t i = 0; i < intx_pts.size(); ++i) {
					cur_sz.push_vertex(intx_pts[i]);
					cur_sz.push_cc_id(cc_ids[i]);
				}
			}
			// After updating the centroid of the SZ, push it to the vector.
			_szs.push_back(cur_sz);
		}
	}
	_szs.shrink_to_fit();
}

vector<StnZone> StnZoneSrch::randomized_search()
{
	std::multimap<size_t, vector<StnZone>> num_sz_map;
	for (size_t t = 0; t < _num_iters_szs; ++t) {
		if (t == 0) { init(false); } // First iteration no need to shuffle circles.
		else { init(true); }  // Later iterations we randomize the circle sequence.
		construct_stn_zones();
		num_sz_map.insert(std::make_pair(_szs.size(), _szs));
	}
	// Lastly set vertices for singletons
	vector<StnZone> final_szs = num_sz_map.begin()->second;
	for (StnZone& sz : final_szs) {
		if (sz.is_singleton()) { sz.set_disc_pts_for_singleton(_num_vtxs); }
	}
	return final_szs;
}
