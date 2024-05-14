/*****************************************************************//**
 * \file   StnZone.h
 * \brief  Steiner Zone
 *
 * \author Bruce
 * \date   August 2023
 *********************************************************************/
#ifndef STNZONE_H
#define STNZONE_H

#include "Circle.h"

class StnZone
{
public:
	StnZone() = default;
	StnZone(const size_t sz_id, const size_t max_circs) : 
		_id(sz_id), _max_cs(max_circs)
	{
		_vs.reserve(max_circs);
		_cs.reserve(max_circs);
		_cc_ids.reserve(max_circs);
	}
	~StnZone() {}

	/**
	 * Get Steiner Zone ID.
	 *
	 * \return
	 */
	size_t id() const { return _id; }

	/**
	 * Get the first added circle.
	 *
	 * \return
	 */
	Circle c1() const { return _cs[0]; }

	/**
	 * Get the number of vertices.
	 *
	 * \return
	 */
	size_t num_vetxs() const { return _vs.size(); }

	/**
	 * Get the number of circles.
	 *
	 * \return
	 */
	size_t num_circs() const { return _cs.size(); }

	/**
	 * Get the number of maximum circles can be involved in a SZ.
	 * 
	 * \return
	 */
	size_t num_max_circs() const { return _max_cs; }

	/**
	 * Get the prize sum.
	 * 
	 * \return 
	 */
	double p() const { return _p; }

	/**
	 * Get the vector of all vertices.
	 *
	 * \return
	 */
	vector<Point> vs() const { return _vs; }

	/**
	 * Check if this SZ is full.
	 *
	 * \return
	 */
	bool is_full() const { return _cs.size() + _ovlp_cs.size() == _max_cs; }

	/**
	 * Check if this SZ is a singleton, i.e. a circle.
	 *
	 * \return
	 */
	bool is_singleton() const { return _cs.size() == 1; }

	/**
	 * Necessity of a circle to be added into a Steiner Zone. It checks if
	 * all circles in this SZ intersect with the Circle to test.
	 *
	 * \param dist_mat
	 * \param c
	 * \return
	 */
	bool check_necessity(const vector<vector<double>>& dist_mat,
		const Circle& c2test) const;

	/**
	 * Check if the circle to add meets the sufficiency condition, i.e., for its
	 * every intersection point with one circle in the existing list, there must
	 * be at least one intersection point located within (or on the boundary) each
	 * circle in the existing list. Note that this check is unnecessary when only
	 * one circle in the existing list.
	 *
	 * \param dist_mat
	 * \param c2test
	 * \param intx_pts
	 * \param circ_ids One intersection point corresponds to two circles. Thus, each
	 * element in intx_pts has a corresponding DualIdx in cc_ids.
	 * \return
	 */
	bool check_sufficiency(const vector<vector<double>>& dist_mat, const Circle& c2test,
		bool& is_ovlp, vector<Point>& intx_pts, vector<DualIdx>& cc_ids) const;

	/**
	 * Push one Point into vector of vertices.
	 *
	 * \param pt
	 */
	void push_vertex(const Point& pt) { _vs.push_back(pt); }

	/**
	 * One vertex corresponds to two circle ids. Push the circle-circle id pair.
	 *
	 * \param cc_id
	 */
	void push_cc_id(const DualIdx& cc_id) { _cc_ids.push_back(cc_id); }

	/**
	 * Push one ID of the covered circle into vector of Circles and sum up prizes.
	 *
	 * \param c The circle ID and prize will be kept in this Steiner Zone.
	 */
	void push_circ(const Circle& c) { _cs.push_back(c); _p += c.p(); }

	/**
	 * There exists one scenario that the circle can be added to this SZ as well,
	 * i.e. the circle completely overlaps this SZ. But in this case, the circle
	 * doesn't have any valid intersection points with this SZ.
	 *
	 * \param c
	 */
	void push_ovlp_circ(const Circle& c) { _ovlp_cs.push_back(c); _p += c.p(); }

	/**
	 * Set discrete points for singletons.
	 * 
	 * \param num_disc_pts
	 */
	void set_disc_pts_for_singleton(const size_t num_disc_pts);

	/**
	 * Arc search based on the line characterized by prev_pt and next_pt. Find the 
	 * optimal point on the boundary of this SZ.
	 * 
	 * \param prev_pt
	 * \param next_pt
	 * \param opt_pt
	 * \return 
	 */
	bool arc_search(const Point& prev_pt, const Point& next_pt, Point& opt_pt) const;

	/**
	 * .
	 * 
	 * \param p1
	 * \param p2
	 * \param init_dist
	 * \param closest_vert
	 * \return 
	 */
	bool find_closest_vertex(const Point& p1, const Point& p2,
		const double init_dist, Point& closest_vert) const;

private:
	double _p = 0;
	vector<Point> _vs;			// Vertices of this SZ
	vector<Circle> _cs;			// Circles covered by this SZ
	vector<Circle> _ovlp_cs;	// Circles that completely overlapped this SZ
	vector<DualIdx> _cc_ids;	// Circle ID pairs matched for all intersection points

	const size_t _id     = 0;	 // ID of Steiner Zone
	const size_t _max_cs = 0;    // Maximum number of circles that can be covered

	bool _is_pt_in_arc(const double theta0,
		const double min_theta, const double max_theta) const;

	/**
	 * Get the common circle of two intersection points.
	 * 
	 * \param i: Index i in _vs
	 * \param j: Index j in _vs
	 * \param circ: The common circle
	 * \return true if found, false otherwise.
	 */
	bool _common_circ(const size_t i, const size_t j, Circle& circ) const;
};

#endif // !STNZONE_H
