#include "StnZone.h"

using std::min, std::max;
using std::make_pair;

bool StnZone::check_necessity(const vector<vector<double>>& dist_mat,
	const Circle& c2test) const
{
	for (const Circle& c : _cs) {
		// Here two circles should be strictly "close" to each other (i.e. overlap), 
		// otherwise it may happen that two circles only intersect at one point. 
		if (dist_mat[c.id()][c2test.id()] >= 0.95 * (c.r() + c2test.r())) {
			return false;
		}
	}
	return true;
}

bool StnZone::check_sufficiency(const vector<vector<double>>& dist_mat,
	const Circle& c2test, bool& is_ovlp, vector<Point>& intx_pts,
	vector<DualIdx>& cc_ids) const
{
	if (_cs.size() == 1) {  // If only 1 circle, necessity is enough.
		PtPair intx_pt_pair = c2test.intx_pts(
			_cs[0], dist_mat[_cs[0].id()][c2test.id()]);
		intx_pts.push_back(intx_pt_pair.first);
		intx_pts.push_back(intx_pt_pair.second);
		cc_ids.push_back({ _cs[0].id(), c2test.id() });
		cc_ids.push_back({ _cs[0].id(), c2test.id() });
		return true;
	}
	vector<size_t> idx2rm; // Indices of intersection points to be removed
	// Because some intersection points may not be feasible any more. Considering 
	// a SZ of degree 2 turned to degree 3, one intersection point is replaced with
	// other 2 new intersection points.
	if (!intx_pts.empty()) {  // Non-empty means some circles already in this SZ
		size_t num_intx_pts = intx_pts.size();
		idx2rm.reserve(num_intx_pts);
		// Here we need strictly larger because intersection points can  
		// be on the circle.
		double r_sq = c2test.r_sq() + EPS;
		for (size_t i = 0; i < num_intx_pts; ++i) {
			if ((intx_pts[i] - c2test.ctr()).sq_sum() > r_sq) { idx2rm.push_back(i); }
		}
		// If all intersection points are located in this circle to be added, 
		// the existing Steiner Zone is completely covered by the circle to test,
		// sufficiency is automatically satisfied, and no need to add any
		// intersection points.
		if (idx2rm.empty()) { is_ovlp = true; return true; }
	}
	// If some intersection points are outside this circle to be added,
	size_t num_cs = _cs.size();
	vector<Point> found_intx_pts; found_intx_pts.reserve(2 * num_cs);
	vector<DualIdx> cc_id_pairs; cc_id_pairs.reserve(2 * num_cs);
	for (size_t i = 0; i < num_cs; ++i) { // Test with all circles already in the SZ
		bool is_p1_in_all_circs = true, is_p2_in_all_circs = true;
		double d_p1c_sq = 0, d_p2c_sq = 0;
		// Find p1 and p2 between the circle to test and the circle i.
		PtPair intx_pt_pair = c2test.intx_pts(_cs[i], dist_mat[_cs[i].id()][c2test.id()]);
		// Then check if p1 or p2 is in ALL OTHER CIRCLES (except circle i)  
		// in this SZ. Note that both p1 and p2 can be in all other circles. 
		for (size_t j = 0; j < num_cs; ++j) {
			if (j == i) { continue; } // No need to check same circle.
			// Only when p1 or p2 is strictly outside the other circle, i.e.,
			// exclude boundary, set in_all_circs flag to false. 
			double r_sq = _cs[j].r_sq() + EPS;
			if (is_p1_in_all_circs) {
				d_p1c_sq = (intx_pt_pair.first - _cs[j].ctr()).sq_sum();
				if (d_p1c_sq > r_sq) { is_p1_in_all_circs = false; }
			}
			if (is_p2_in_all_circs) {
				d_p2c_sq = (intx_pt_pair.second - _cs[j].ctr()).sq_sum();
				if (d_p2c_sq > r_sq) { is_p2_in_all_circs = false; }
			}
			if (!is_p1_in_all_circs && !is_p2_in_all_circs) { break; }
		}
		if (is_p1_in_all_circs) {
			found_intx_pts.push_back(intx_pt_pair.first);
			cc_id_pairs.push_back(DualIdx(_cs[i].id(), c2test.id()));
		}
		if (is_p2_in_all_circs) {
			found_intx_pts.push_back(intx_pt_pair.second);
			cc_id_pairs.push_back(DualIdx(_cs[i].id(), c2test.id()));
		}
	}
	if (!found_intx_pts.empty()) {
		if (!intx_pts.empty()) { // Now remove former infeasible intersection points.
			for (auto it = idx2rm.rbegin(); it != idx2rm.rend(); ++it) {
				intx_pts.erase(intx_pts.begin() + *it);
				cc_ids.erase(cc_ids.begin() + *it);
			}
		}
		for (size_t i = 0; i < found_intx_pts.size(); ++i) {
			intx_pts.push_back(found_intx_pts[i]); cc_ids.push_back(cc_id_pairs[i]);
		}
		return true;
	}
	return false;
}

void StnZone::update_sz()
{
	if (_cs.size() == 1) {  // Singleton
		_ctr = _cs[0].ctr();
		_vx_max = _vy_max = _cs[0].r();
	}
	else {  // Overlapped circles
		size_t num_vs = _vs.size();
		// Find the max and min x and y coordinates of all vertices.
		double x_max = _vs[num_vs - 1].x, x_min = x_max;
		double y_max = _vs[num_vs - 1].y, y_min = y_max;
		_sq_dists = vector<vector<double>>(num_vs, vector<double>(num_vs));
		for (size_t i = 0; i < num_vs - 1; ++i) {
			_ctr += _vs[i];
			x_min = min(x_min, _vs[i].x); x_max = max(x_max, _vs[i].x);
			y_min = min(y_min, _vs[i].y); y_max = max(y_max, _vs[i].y);
			for (size_t j = i + 1; j < num_vs; ++j) {
				_sq_dists[i][j] = _sq_dists[j][i] = (_vs[j] - _vs[i]).sq_sum();
			}
		}
		// Still miss the last one
		_ctr += _vs.back();
		_ctr /= static_cast<double>(num_vs);
		const double k = 0.5;
		double x_diff = (x_max - x_min) * k, y_diff = (y_max - y_min) * k;
		_vx_max = abs(x_diff) <= EPS ? y_diff : x_diff;
		_vy_max = abs(y_diff) <= EPS ? x_diff : y_diff;
	}
}

void StnZone::_circ_feas_range(const Circle& c, double& theta1, double& theta2) const
{
	size_t cid = c.id();
	vector<size_t> idxs; idxs.reserve(2);
	for (size_t i = 0; i < _cc_ids.size(); ++i) {
		if (_cc_ids[i].first == cid || _cc_ids[i].second == cid) {
			idxs.push_back(i);
		}
	}
	if (!idxs.empty()) {
		theta1 = c.cartesian2polar(_vs[idxs[0]]);
		theta2 = c.cartesian2polar(_vs[idxs[1]]);
	}
}

void StnZone::restrict_pos(const Point& pt, const Point& vel, Point& new_pt) const
{
	Point p2 = pt + vel;
	if (_vs.empty()) {  // If singleton, restrict point to circumference.
	  if (!_cs[0].is_close_enough(p2)) { new_pt = _cs[0].ext_pt2bound(p2); }
	  else { new_pt = p2; }
	  return;
	}
	// If the SZ includes >= 2 circles, check each circle's feasible range.
	// There will only be one valid intersection point because the SZ is convex.
	// Thus, we only need to find the closest feasible intersection point.
	size_t num_cs = _cs.size();
	vector<Point::IdxPtPair> cid_intx_pts; cid_intx_pts.reserve(num_cs * 2);
	for (size_t i = 0; i < num_cs; i++) {
		// If the circle is not the completely overlap circle:
		if (!_cs[i].is_close_enough(p2)) { // If p2 is outside the circle.
			// Get the intersection points between the line and the circle.	
			PtPair intx1_2 = _cs[i].intx_pts(pt, p2);
			Point intx1 = intx1_2.first, intx2 = intx1_2.second;
			// Convert intersection points to polar coordinates.
			double intx1_angle = _cs[i].cartesian2polar(intx1);
			double intx2_angle = _cs[i].cartesian2polar(intx2);
			// Find the feasible range of the circle.
			double angle1 = 0, angle2 = 0;
			_circ_feas_range(_cs[i], angle1, angle2);
			double max_angle = max(angle1, angle2);
			double min_angle = min(angle1, angle2);
			// Check if the intersection points are within the feasible range.
			if (Common::is_same_sign(angle1, angle2)) { // If at same side 
				if (min_angle <= intx1_angle && intx1_angle <= max_angle) {
					cid_intx_pts.push_back(make_pair(i, intx1));
				}
				if (min_angle <= intx2_angle && intx2_angle <= max_angle) {
					cid_intx_pts.push_back(make_pair(i, intx2));
				}
			}
			else {  // Otherwise, min_angle in [-pi, 0], max_angle in [0, pi]
				double angle_diff = max_angle - min_angle;
				if (angle_diff > PI) {
					// If > pi, then feasible arc pass through pi or -pi.
					if (intx1_angle <= min_angle || max_angle <= intx1_angle) {
						// p1 in [-pi, min_angle] or in [max_angle, pi]
						cid_intx_pts.push_back(make_pair(i, intx1));
					}
					if (intx2_angle <= min_angle || max_angle <= intx2_angle) {
						// p2 in [-pi, min_angle] or in [max_angle, pi]
						cid_intx_pts.push_back(make_pair(i, intx2));
					}
				}
				else {
					if (min_angle <= intx1_angle && intx1_angle <= max_angle) {
						cid_intx_pts.push_back(make_pair(i, intx1));
					}
					if (min_angle <= intx2_angle && intx2_angle <= max_angle) {
						cid_intx_pts.push_back(make_pair(i, intx2));
					}
				}
			}
		}
	}
	// Then update new point and new velocity (possibly).
	if (cid_intx_pts.empty()) { new_pt = p2; }
	else {
		// If has more that one intersection point, sort them by distance to p2.
		if (cid_intx_pts.size() > 1) { p2.sort_pt_vec(cid_intx_pts); }
		// Then find out the closest feasible intersection point.
		new_pt = cid_intx_pts[0].second;
		// If the new point hits the boundary, reflect the velocity.
		//vel = _cs[intx_pts[0].first].reflect_vec(pt, vel, new_pt);
	}
}

void StnZone::restrict_vel(Point& vel) const
{
	Common::restrict_value_in_range(-_vx_max, _vx_max, vel.x);
	Common::restrict_value_in_range(-_vy_max, _vy_max, vel.y);
}

double StnZone::_dist_pt2line(const double len_sq, const Point& pt, 
	const Point& p1, const Point& p2) const
{
	double A = pt.x - p1.x, B = pt.y - p1.y;
	double C = p2.x - p1.x, D = p2.y - p1.y;
	double dot = A * C + B * D;
	double param = dot / len_sq;
	double xx = 0, yy = 0;
	if (param < 0) { xx = p1.x; yy = p1.y; }
	else if (param > 1) { xx = p2.x; yy = p2.y; }
	else { xx = p1.x + param * C; yy = p1.y + param * D; }
	double dx = pt.x - xx, dy = pt.y - yy;
	return sqrt(pow(dx, 2) + pow(dy, 2));
}

Point StnZone::rand_pt() const
{
	if (_vs.empty()) { return _cs[0].rand_pt(); }
	else if (_vs.size() == 2) {  // Only works for circle with same radius
		double d = sqrt(_sq_dists[0][1]) / 2.;
		double internal_r = _cs[0].r() - sqrt(_cs[0].r_sq() - pow(d, 2));
		return Circle(_ctr, internal_r).rand_pt();
	}
	else {
		vector<double> ds; ds.reserve(_vs.size());
		vector<bool> trav_mask = vector<bool>(_vs.size(), false);
		size_t cur_idx = 0;
		trav_mask[cur_idx] = true;
		size_t num_true = 1;
		do {
			std::multimap<double, size_t> dist_map;
			for (size_t j = 1; j < _vs.size(); j++) {
				if (!trav_mask[j]) { 
					dist_map.insert(make_pair(_sq_dists[cur_idx][j], j)); 
				}
			}
			size_t next_idx = dist_map.begin()->second;
			Point p1 = _vs[cur_idx], p2 = _vs[next_idx];
			ds.push_back(_dist_pt2line(_sq_dists[cur_idx][next_idx], _ctr, p1, p2));
			cur_idx = next_idx;
			trav_mask[cur_idx] = true;
			num_true++;
		} while (num_true < _vs.size());
		// Complete the last one
		ds.push_back(_dist_pt2line(_sq_dists[cur_idx][0], _ctr, _vs[cur_idx], _vs[0]));
		Point rd_pt = Circle(_ctr, Common::vec_min<double>(ds)).rand_pt();
		return rd_pt;
	}
}

double StnZone::nbhd_c(const Point& pt) const
{
	if (_vs.empty()) {
		double d_pt_c = (pt - _cs[0].ctr()).norm2();
		double v_uav_depart = _cs[0].w() * _v_uav;
		double this_cost = d_pt_c / v_uav_depart + _t_service + d_pt_c / _v_uav;
		return this_cost;
	}
	else {
		vector<double> costs; costs.reserve(_cs.size());
		for (const Circle& c : _cs) {
			double d_pt_c = (pt - c.ctr()).norm2();
			double v_depart = c.w() * _v_uav;
			double this_cost = d_pt_c / v_depart + _t_service + d_pt_c / _v_uav;
			costs.push_back(this_cost);
		}
		if (!_ovlp_cs.empty()) {
			for (const Circle& c : _ovlp_cs) {
				double d_pt_c = (pt - c.ctr()).norm2();
				double v_depart = c.w() * _v_uav;
				double this_cost = d_pt_c / v_depart + _t_service + d_pt_c / _v_uav;
				costs.push_back(this_cost);
			}
		}
		return Common::vec_max<double>(costs);
	}
}

void StnZone::free_memory()
{
	for (size_t i = 0; i < _sq_dists.size(); i++) { _sq_dists[i].clear(); }
	_sq_dists.clear(); _sq_dists.shrink_to_fit();
}
