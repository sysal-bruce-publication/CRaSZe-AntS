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
