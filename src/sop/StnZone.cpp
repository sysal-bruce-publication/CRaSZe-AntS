#include "StnZone.h"

using std::min, std::max;
using std::make_pair;

bool StnZone::check_necessity(const vector<vector<double>>& dist_mat,
	const Circle& c2test) const
{
	for (const Circle& c : _cs) {
		// Here two circles should be strictly "close" to each other (i.e. overlap), 
		// otherwise it may happen that two circles only intersect at one point. 
		if (dist_mat[c.id()][c2test.id()] >= c.r() + c2test.r() - EPS) {
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

void StnZone::set_disc_pts_for_singleton(const size_t num_disc_pts)
{
	// Get the circle's center and its radius
	Point ctr0 = _cs[0].ctr(); double r0 = _cs[0].r();
	vector<Point> disc_pts;
	if (num_disc_pts == 2) {       // 90, 270 degree
		_vs.push_back(Point(ctr0.x, ctr0.y + 0.5 * r0));
		_vs.push_back(Point(ctr0.x, ctr0.y - 0.5 * r0));
	}
	else if (num_disc_pts == 3) {  // 30, 150, 270 degree
		_vs.push_back(Point(ctr0.x + 0.4330127 * r0, ctr0.y + 0.25 * r0));
		_vs.push_back(Point(ctr0.x - 0.4330127 * r0, ctr0.y + 0.25 * r0));
		_vs.push_back(Point(ctr0.x, ctr0.y - 0.5 * r0));
	}
	else if (num_disc_pts == 4) {  // 45, 135, 225, 315 degree
		_vs.push_back(Point(ctr0.x + 0.5 * r0, ctr0.y + 0.5 * r0));
		_vs.push_back(Point(ctr0.x - 0.5 * r0, ctr0.y + 0.5 * r0));
		_vs.push_back(Point(ctr0.x - 0.5 * r0, ctr0.y - 0.5 * r0));
		_vs.push_back(Point(ctr0.x + 0.5 * r0, ctr0.y - 0.5 * r0));
	}
	else { _vs.push_back(ctr0); return; }
}

bool StnZone::_is_pt_in_arc(const double theta0,
	const double min_theta, const double max_theta) const
{
	if (Common::is_same_sign(min_theta, max_theta)) {  // If at same side 
		if (min_theta <= theta0 && theta0 <= max_theta) { return true; }
		else { return false; }
	}
	else {  // Otherwise, min_angle in [-pi, 0], max_angle in [0, pi]
		double angle_diff = max_theta - min_theta;
		if (angle_diff > PI) {  // If > pi, then feasible arc pass through pi or -pi.
			if (theta0 <= min_theta || max_theta <= theta0) { return true; }
			else { return false; }
		}
		else {
			if (min_theta <= theta0 && theta0 <= max_theta) { return true; }
			else { return false; }
		}
	}
}

bool StnZone::_common_circ(const size_t i, const size_t j, Circle& circ) const
{
	size_t common_cid = 0;
	if (_cc_ids[j].first == _cc_ids[i].first ||
		_cc_ids[j].first == _cc_ids[i].second) {
		common_cid = _cc_ids[j].first;
	}	
	if (_cc_ids[j].second == _cc_ids[i].first
		|| _cc_ids[j].second == _cc_ids[i].second) {
		common_cid = _cc_ids[j].second;
	}
	auto it = std::find_if(_cs.begin(), _cs.end(), 
		[common_cid](const Circle& c) -> bool { return common_cid == c.id(); });
	if (it != _cs.end()) { circ = *it; return true; }
	// For SZ > 3, it may happen 2 vertices on two pairs of different circles.
	else { return false; }  
}

bool StnZone::arc_search(const Point& prev_pt,
	const Point& next_pt, Point& opt_pt) const
{
	// For singleton
	if (is_singleton()) {  opt_pt = _cs[0].opt_pt(prev_pt, next_pt); return true; }
	// For SZ with degree >= 2
	if (_cs.size() == 2) {  // Degree 2
		for (const Circle& c : _cs) {
			unique_ptr<Point> circ_opt = make_unique<Point>();
			// Get the relation between Line segment (prev_pt, next_pt) and circle.
			vector<Point> close_pts;
			c.closest_pt_line_segment(prev_pt, next_pt, *circ_opt, close_pts);
			// Get the feasible range of the circle.
			unique_ptr<Point> vec1 = make_unique<Point>(_vs[0] - c.ctr());
			unique_ptr<Point> vec2 = make_unique<Point>(_vs[1] - c.ctr());
			double theta1 = atan2(vec1->y, vec1->x), theta2 = atan2(vec2->y, vec2->x);
			double min_theta = std::min(theta1, theta2);
			double max_theta = std::max(theta1, theta2);
			if (close_pts.empty()) {  // Line segment misses circle.
				Point ctr2opt = *circ_opt - c.ctr();
				double opt_theta = atan2(ctr2opt.y , ctr2opt.x);
				if (_is_pt_in_arc(opt_theta, min_theta, max_theta)) {
					opt_pt = *circ_opt; return true;
				}
			}
			else {  // Line segment tangent/intersect circle
				for (const Point& pt : close_pts) {
					Point ctr2pt = pt - c.ctr();
					double pt_theta = atan2(ctr2pt.y, ctr2pt.x);
					if (_is_pt_in_arc(pt_theta, min_theta, max_theta)) {
						// If the found closest point is in range of two vertices,
						// then update the opt_pt and stop this program
						opt_pt = pt; return true;
					}
					// Unlike above code, here else scenario will not stop
					// the program (in case more than one close points)
				}
			}
		}
	}
	else {  // Degree >2
		for (size_t i = 0; i < _vs.size() - 1; i++) {
			for (size_t j = i + 1; j < _vs.size(); j++) {
				Circle c;
				if (_common_circ(i, j, c)) {
					unique_ptr<Point> circ_opt = make_unique<Point>();
					// Get relation between Line segment (prev_pt, next_pt) and circle.
					vector<Point> close_pts;
					c.closest_pt_line_segment(prev_pt, next_pt, *circ_opt, close_pts);
					// Get the feasible range of the circle.
					unique_ptr<Point> vec1 = make_unique<Point>(_vs[i] - c.ctr());
					unique_ptr<Point> vec2 = make_unique<Point>(_vs[j] - c.ctr());
					double theta1          = atan2(vec1->y, vec1->x);
					double theta2          = atan2(vec2->y, vec2->x);
					double min_theta       = std::min(theta1, theta2);
					double max_theta       = std::max(theta1, theta2);
					if (close_pts.empty()) {  // Line segment misses circle.
						Point ctr2opt = *circ_opt - c.ctr();
						double opt_theta = atan2(ctr2opt.y, ctr2opt.x);
						if (_is_pt_in_arc(opt_theta, min_theta, max_theta)) {
							opt_pt = *circ_opt; return true;
						}
					}
					else {  // Line segment tangent/intersect circle
						for (const Point& pt : close_pts) {
							Point ctr2pt = pt - c.ctr();
							double pt_theta = atan2(ctr2pt.y, ctr2pt.x);
							if (_is_pt_in_arc(pt_theta, min_theta, max_theta)) {
								// If the found closest point is in range of arc,
								// then update the opt_pt and stop this program.
								opt_pt = pt; return true;
							}
							// Unlike above code, here else scenario will not stop
							// the program (in case more than one close points)
						}
					}
				}
			}
		}
	}
	return false;
}

bool StnZone::find_closest_vertex(const Point& p1, const Point& p2,
	const double init_dist, Point& closest_vert) const
{
	double min_d = 1e6; Point temp_pt = Point();
	for (const Point& v : _vs) {
		double d = (v - p1).norm2() + (v - p2).norm2();
		if (d < min_d) { min_d = d; temp_pt = v; }
	}
	if (abs(init_dist - min_d) < EPS) { return false; }
	else { closest_vert = temp_pt; return true; }
}
