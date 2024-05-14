/*****************************************************************//**
 * \file   Circle.h
 * \brief  Basic class for Circle(center, radius) object and its operations.
 * 
 * \author Qiuchen Qian
 * \date   October 2023
 *********************************************************************/
#ifndef CIRCLE_H
#define CIRCLE_H

#include <memory>
#include "Point.h"

using std::unique_ptr, std::make_unique;
typedef std::pair<Point, Point> PtPair;

class Circle
{
public:
	Circle() = default;
	/**
	 * Copy constructor.
	 *
	 * \param c
	 */
	Circle(const Circle& c) : _id(c._id), _ctr(c._ctr), _w(c._w), _p(c._p),
		_r(c._r), _r_sq(c._r_sq), _x1(c._x1), _x2(c._x2), _y1(c._y1), _y2(c._y2) {}

	/**
	 * For start and end depot.
	 * 
	 * \param id
	 * \param center
	 */
	Circle(const size_t id, const Point& center) : _id(id), _ctr(center) {}

	/**
	 * Create a virtual circle. For circle geometry or calculation only.
	 * 
	 * \param center
	 * \param radius
	 */
	Circle(const Point& center, const double radius) : _ctr(center), _r(radius) {}

	/**
	 * For target circles.
	 * 
	 * \param id
	 * \param radius
	 * \param center
	 */
	Circle(const size_t id, const Point& center, const double radius, 
		const double prize, const double weight) : _id(id), _ctr(center), 
		_w(weight), _p(prize), _r(radius), _r_sq(pow(radius, 2)),
		_x1(center.x - radius), _x2(center.x + radius), 
		_y1(center.y - radius), _y2(center.y + radius) {}

	~Circle() {}

	size_t id() const { return _id; }
	Point ctr() const { return _ctr; }
	double r() const { return _r; }
	double r_sq() const { return _r_sq; }
	double w() const { return _w; }
	double p() const { return _p; }

	void set_id(const size_t cid) { _id = cid; }

	bool operator==(const Circle& c2) const { return _id == c2._id; }
	Point operator-(const Circle& c2) const { return _ctr - c2._ctr; }
	Point operator+(const Circle& c2) const { return _ctr + c2._ctr; }

	/**
	 * Convert a Point on this circle with polar coordinate (here angle only)
	 * to Cartesian coordinate. Angle range [-pi, pi].
	 *
	 * \param angle
	 * \return
	 */
	Point polar2cartesian(const double& angle) const
	{
		return Point(_ctr.x + _r * cos(angle), _ctr.y + _r * sin(angle));
	}

	/**
	 * Convert a Point in euclidean coordinate to polar coordinate of this circle
	 * with atan2 function, range [-pi, pi].
	 * 
	 * \param pt
	 * \return 
	 */
	double cartesian2polar(const Point& pt) const
	{
		return atan2(pt.y - _ctr.y, pt.x - _ctr.x);
	}

	/**
	 * Abstract two circles to two square. Compare if two squares have overlap.
	 *
	 * \param c2 Another Circle object.
	 * \return true if overlapped. Otherwise false.
	 */
	bool is_close(const Circle& c2) const
	{
		if (_x1 < c2._x2 - EPS && _x2 > c2._x1 + EPS
			&& _y2 > c2._y1 + EPS && _y1 < c2._y2 - EPS) {
			return true;
		}
		else { return false; }
	}

	/**
	 * Check if the point is within the circle or on the circle. Here we use
	 * strictly larger to determine the boundary equal.
	 *
	 * \param pt
	 * \return false if outside the circle. Otherwise true.
	 */
	bool is_close_enough(const Point& pt) const
	{
		return (pt - _ctr).sq_sum() > _r_sq + EPS ? false : true;
	}

	/**
	 * Set an external point to the boundary of this circle.
	 *
	 * \param pt
	 * \return
	 */
	inline Point ext_pt2bound(const Point& pt) const
	{
		Point diff = pt - _ctr;
		double theta = atan2(diff.y, diff.x);
		return Point(_ctr.x + _r * cos(theta), _ctr.y + _r * sin(theta));
	}

	/**
	 * Generate a random point within (or on) the circle.
	 *
	 * \return
	 */
	Point rand_pt() const
	{
		std::unique_ptr<Rd> rand = std::make_unique<Rd>();
		double theta = rand->uniform_rand(0., 2 * PI);
		double r = rand->uniform_rand(0., _r + EPS); // Allow equal
		return Point(_ctr.x + r * cos(theta), _ctr.y + r * sin(theta));
	}

	/**
	 * Get the intersect points of another circle. See below link for more details:
	 * https://math.stackexchange.com/a/1367732/735790
	 *
	 * \param c2: Another Circle object
	 * \param cc_dist: Optional parameter, the distance between two circle centers.
	 *	If negative, it will be calculated automatically.
	 * \return
	 */
	PtPair intx_pts(const Circle& c2, double cc_dist) const;

	/**
	 * Find the intersection points between the line (characterized by
	 * p1 & p2) and this circle.
	 * 
	 * \param p1
	 * \param p2
	 * \return 
	 */
	PtPair intx_pts(const Point& p1, const Point& p2) const;

private:
	size_t _id   = 0;			// ID
	Point _ctr   = Point();		// Center
	double _w    = 0;			// Parcel weight
	double _p    = 0;			// Prize (priority) 
	double _r    = 0;			// Radius
	double _r_sq = 0;			// Radius squared
	double _x1   = 0;			// Left extreme point of the circle
	double _x2   = 0;			// Right extreme point of the circle
	double _y1   = 0;			// Bottom extreme point of the circle
	double _y2   = 0;			// Top extreme point of the circle
};

#endif // !CIRCLE_H
