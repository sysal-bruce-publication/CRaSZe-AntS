/*****************************************************************//**
 * \file   Point.h
 * \brief  Basic class for Point(x, y) object and its operations.
 * 
 * \author Qiuchen Qian
 * \date   October 2023
 *********************************************************************/
#ifndef POINT_H
#define POINT_H

#include "Utils.h"

using std::abs, std::size_t;
typedef std::pair<size_t, size_t> DualIdx;

struct Point
{
	Point() = default;
	Point(double x_, double y_) : x(x_), y(y_) {}
	Point(const Point& pt) : x(pt.x), y(pt.y) {}
	~Point() {}

	double x = 0;
	double y = 0;

	typedef std::pair<size_t, Point> IdxPtPair;

	double sq_sum() const { return pow(x, 2) + pow(y, 2); }
	double norm2() const { return sqrt(sq_sum()); }

	Point& operator=(const Point& pt) { x = pt.x; y = pt.y; return *this; }
	bool operator==(const Point& pt) const
	{
		if (abs(pt.x - x) < EPS && abs(pt.y - y) < EPS) { return true; }
		else { return false; }
	}
	Point operator+(const Point& pt) const { return Point(x + pt.x, y + pt.y); }
	Point operator-(const Point& pt) const { return Point(x - pt.x, y - pt.y); }
	/**
	 * Dot product (i.e. Hadamard product) of points.
	 *
	 * \param pt
	 * \return
	 */
	double operator*(const Point& pt) const { return x * pt.x + y * pt.y; }
	Point operator*(const double val) const { return Point(x * val, y * val); }
	Point operator/(const double val) const { return Point(x / val, y / val); }
	Point& operator+=(const Point& pt) { x += pt.x; y += pt.y; return *this; }
	Point& operator-=(const Point& pt) { x -= pt.x; y -= pt.y; return *this; }
	Point& operator*=(const double val) { x *= val; y *= val; return *this; }
	Point& operator/=(const double val) { x /= val; y /= val; return *this; }
	/**
	 * Cross product of two point (vector w.r.t. the origin).
	 *
	 * \param other
	 * \return
	 */
	double cross_prod(const Point& pt) const { return x * pt.y - pt.y * x; }
	/**
	 * Given a vector of points, sort it according the distance.
	 *
	 * \param idx_pts
	 */
	void sort_pt_vec(vector<IdxPtPair>& idx_pts) const
	{
		return std::sort(idx_pts.begin(), idx_pts.end(),
			[this](const IdxPtPair& l_pt, const IdxPtPair& r_pt) -> bool
			{
				return (l_pt.second - *this).sq_sum() < (r_pt.second - *this).sq_sum();
			});
	}

	double vec_angle(const Point& pt) const
	{
		return acos((*this * pt) / (norm2() * pt.norm2()));
	}
};

#endif // !POINT_H

