/*****************************************************************//**
 * \file   StnZone.h
 * \brief  basic class for Steiner Zone (multiple circles) object and its operations.
 * 
 * \author Bruce
 * \date   August 2023
 *********************************************************************/
#ifndef STNZONE_H
#define STNZONE_H

#include <map>
#include "Circle.h"


class StnZone
{
public:
	StnZone() = default;
	StnZone(const size_t sz_id, const size_t max_circs, const double v_uav) :
		_id(sz_id), _max_cs(max_circs), _v_uav(v_uav)
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
	 * Get Steiner Zone prize.
	 * 
	 * \return 
	 */
	double p() const { return _p; }

	/**
	 * Get the centroid of this SZ.
	 *
	 * \return
	 */
	Point ctr() const { return _ctr; }

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
	 * Get the vector of all vertices.
	 *
	 * \return
	 */
	vector<Point> vs() const { return _vs; }

	vector<Circle> cs() const { return _cs; }

	/**
	 * Check if this SZ is full.
	 *
	 * \return
	 */
	bool is_full() const { return _cs.size() + _ovlp_cs.size() >= _max_cs; }

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
	 * Update the centroid of this SZ, corresponding bound and distance matrix
	 * of all vertices.
	 *
	 * \return
	 */
	void update_sz();

	/**
	 * Restrict a random point within the constraint
	 *
	 * \param pt: Current pos (not moved yet).
	 * \param vel: Velocity of the pos.
	 * \param new_pt: New pos after restriction.
	 * \return
	 */
	void restrict_pos(const Point& pt, const Point& vel, Point& new_pt) const;

	/**
	 * Restrict vel of a point
	 *
	 * \param v
	 * \return
	 */
	void restrict_vel(Point& vel) const;

	/**
	 * Generate a random point based on the centroid.
	 *
	 * \return
	 */
	Point rand_pt() const;

	/**
	 * With given Point, get the neighborhood cost.
	 *
	 * \param pt
	 * \return
	 */
	double nbhd_c(const Point& pt) const;

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
	 * Push one ID of the covered circle into vector of circ_ids.
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
	 * Free memory of distance matrix.
	 * 
	 * \return 
	 */
	void free_memory();

private:
	double _p      = 0;
	double _vx_max = 1e6;		// Minimum x-coordinate bound of this Steiner Zone
	double _vy_max = 1e6;			// Maximum x-coordinate bound of this Steiner Zone
	Point _ctr     = Point();
	vector<Point> _vs;
	vector<Circle> _cs;
	vector<Circle> _ovlp_cs;	// IDs of circles completely overlapped this SZ
	vector<DualIdx> _cc_ids;	// Circle ID pairs matched for all intersection points
	vector<vector<double>> _sq_dists;	// squared distance matrix of all vertices

	const size_t _id        = 0;		  // ID of Steiner Zone
	const size_t _max_cs    = 0;		  // Maximum number of circles that can be covered
	const double _v_uav     = 0;		  // UAV velocity
	const double _t_service = 0.08333333; // UAV service time

	/**
	 * Get the feasible range of the given Circle, range [-pi, pi].
	 *
	 * \param c: The given Circle.
	 * \param theta1
	 * \param theta2
	 */
	void _circ_feas_range(const Circle& c, double& theta1, double& theta2) const;

	double _dist_pt2line(const double len_sq, const Point& pt,
		const Point& p1, const Point& p2) const;
};

#endif // !STNZONE_H
