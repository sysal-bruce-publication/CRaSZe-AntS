#include "Circle.h"

PtPair Circle::intx_pts(const Circle& c2, double cc_dist) const
{
	double c2_x = c2._ctr.x, c2_y = c2._ctr.y, c2_r = c2._r;
	// Get SQUARED circle-circle distance if not given, i.e.,
	// D = \sqrt{ (x_2 - x_1)^2 + (y_2 - y_1)^2 }
	double cc_dist_sq = 0;
	if (cc_dist < 0) { cc_dist_sq = (*this - c2).sq_sum(); }
	else { cc_dist_sq = pow(cc_dist, 2); }
	if (cc_dist_sq < EPS) { return PtPair(); }
	double c1_r_sq = pow(_r, 2), c2_r_sq = pow(c2_r, 2);
	double r_sq_sum = c1_r_sq + c2_r_sq, r_sq_sub = c1_r_sq - c2_r_sq;
	// 2nd coefficient (r_1^2 - r_2^2) / (D^2)
	double coeff2 = r_sq_sub / cc_dist_sq;
	// 3rd coefficient 
	// \sqrt{ 2 * (r_1^2 + r_2^2) / D^2 - (r_1^2 - r_2^2)^2 / D^4 - 1 }
	double coeff3 = sqrt(2 * r_sq_sum / cc_dist_sq
		- pow(r_sq_sub, 2) / pow(cc_dist_sq, 2) - 1);
	// The fixed x term, i.e., x_1 + x_2 + coeff2 * (x_2 - x_1)
	double fx = _ctr.x + c2_x + coeff2 * (c2_x - _ctr.x);
	// The last x term, i.e., coeff3 * (y_2 - y_1)
	double gx = coeff3 * (c2_y - _ctr.y);
	// The fixed y term, i.e., y_1 + y_2 + coeff2 * (y_2 - y_1)
	double fy = _ctr.y + c2_y + coeff2 * (c2_y - _ctr.y);
	// The last y term, i.e., coeff3 * (x_1 - x_2)
	double gy = coeff3 * (_ctr.x - c2_x);

	return PtPair(Point(0.5 * (fx + gx), 0.5 * (fy + gy)),
		Point(0.5 * (fx - gx), 0.5 * (fy - gy)));
}

PtPair Circle::intx_pts(const Point& p1, const Point& p2) const
{
	double A = p2.y - p1.y, B = p1.x - p2.x;
	double C = p2.x * p1.y - p1.x * p2.y;
	double A_sq = pow(A, 2), B_sq = pow(B, 2), C_sq = pow(C, 2);
	double x0_sq = pow(_ctr.x, 2), y0_sq = pow(_ctr.y, 2);
	double a = A_sq + B_sq, b = 0, c = 0;
	double twice_a = 2 * a;
	if (std::abs(B) > EPS) {
		double B_y0_prod = B * _ctr.y;
		b = 2 * A * C + 2 * A * B_y0_prod - 2 * B_sq * _ctr.x;
		c = C_sq + 2 * C * B_y0_prod - B_sq * (_r_sq - x0_sq - y0_sq);
		// x has solution (-b +/- sqrt(b^2 - 4ac)) / 2a
		double delta = sqrt(pow(b, 2) - 4 * a * c) / twice_a;
		double fixed = -b / twice_a;
		double p1_x = fixed + delta, p2_x = fixed - delta;
		return PtPair(Point(p1_x, (A * p1_x + C) / -B),
			Point(p2_x, (A * p2_x + C) / -B));
	}
	else { // In case B is close to 0
		double A_x0_prod = A * _ctr.x;
		b = 2 * B * C + 2 * B * A_x0_prod - 2 * A_sq * _ctr.y;
		c = C_sq + 2 * C * A_x0_prod - A_sq * (_r_sq - x0_sq - y0_sq);
		// y has solution (-b +/- sqrt(b^2 - 4ac)) / 2a
		double delta = sqrt(pow(b, 2) - 4 * a * c) / twice_a;
		double fixed = -b / twice_a;
		double p1_y = fixed + delta, p2_y = fixed - delta;
		return PtPair(Point((B * p1_y + C) / -A, p1_y),
			Point((B * p2_y + C) / -A, p2_y));
	}
}

Point Circle::opt_pt(const Point& p1, const Point& p2) const
{
	Point p1p2 = p2 - p1;
	double d_p1p2 = p1p2.norm2();
	Point proj_pt;
	if (d_p1p2 <= EPS) { proj_pt = p1; }
	else {
		// Get p1->p2 's unit vector (direction)
		Point unit_vec = p1p2 / d_p1p2;
		// Get the projection of unit_vec on p1->c, 0 <= |proj| <= d_p1p2 
		double d = std::min(d_p1p2, std::max(0., unit_vec * (_ctr - p1)));
		// Then add this distance to unit vector. 
		proj_pt = p1 + unit_vec * d;
	}
	// Get c->pt_on_p1p2
	Point diff = proj_pt - _ctr;
	double theta = atan2(diff.y, diff.x);
	double d_diff_sq = diff.sq_sum();
	// if proj_pt is inside (or on) the circle
	if (d_diff_sq <= _r_sq + EPS) { return proj_pt; }
	// Otherwise, the optimal point is set on the circle 
	// (along the direction of c->opt_pt).
	return _ctr + (Point(cos(theta), sin(theta)) * _r);
}

void Circle::closest_pt_line_segment(const Point& p1, const Point& p2,
	Point& opt_pt, vector<Point>& close_pts) const
{
	Point p1p2 = p2 - p1;
	double d_p1p2 = p1p2.norm2();
	Point proj_pt;
	double d_p1_proj = 0;
	Point unit_vec;
	if (d_p1p2 <= EPS) { proj_pt = p1; }
	else {
		// Get p1->p2 's unit vector (direction)
		unit_vec = p1p2 / d_p1p2;
		// Get the projection of unit_vec on p1->c, 0 <= |proj| <= d_p1p2 
		d_p1_proj = std::min(d_p1p2, std::max(0., unit_vec * (_ctr - p1)));
		// Then add this distance to unit vector. 
		proj_pt = p1 + unit_vec * d_p1_proj;
	}
	Point ctr2proj = proj_pt - _ctr;
	double d_diff_sq = ctr2proj.sq_sum();
	if (!close_pts.empty()) { close_pts.clear(); } close_pts.reserve(2);
	if (d_diff_sq > _r_sq + EPS) {  	// The circle miss the line segment
		// The closest point is on the circle
		double theta = atan2(ctr2proj.y, ctr2proj.x);
		opt_pt = _ctr + Point(cos(theta) * _r, sin(theta) * _r);
	}
	else if (d_diff_sq < _r_sq - EPS) {  // Inside
		// Projection point the intersection of vertical line and line segment,
		// but here we want the intersection point(s) of p1p2 and the circle
		// because proj point's polar coordinate is "incorrect". See further
		// arc search algorithm for more info.
		double d_proj_intx = sqrt(_r_sq - d_diff_sq);
		double d_p1_intx = d_p1_proj - d_proj_intx;
		close_pts.push_back(p1 + unit_vec * d_p1_intx);
		// There exists two intersect points if p1p2 go through the circle.
		double d_next_intx = d_p1_proj + d_proj_intx;
		if (d_next_intx <= d_p1p2 + EPS) {
			close_pts.push_back(p1 + unit_vec * d_next_intx);
		}
	}
	else { close_pts.push_back(proj_pt); }  // tangent
	close_pts.shrink_to_fit(); 
}
