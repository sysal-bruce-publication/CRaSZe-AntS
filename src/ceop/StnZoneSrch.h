/*****************************************************************//**
 * \file   StnZoneSrch.h
 * \brief  Randomized Steiner Zone Search operator
 *
 * \author Bruce
 * \date   August 2023
 *********************************************************************/
#ifndef STNZONESRCH_H
#define STNZONESRCH_H

#include <map>
#include "StnZone.h"

typedef std::pair<size_t, Circle> IdxCirc;

class StnZoneSrch
{
public:
	StnZoneSrch() = default;
	StnZoneSrch(const size_t num_iters_szs, const size_t max_circs, 
		const size_t num_vtxs_singleton, vector<Circle>& all_circs) : 
		_num_iters_szs(num_iters_szs),  _num_circs(all_circs.size()), 
		_max_circs(max_circs), _num_vtxs(num_vtxs_singleton)
	{
		_cs = std::move(all_circs);
	}
	~StnZoneSrch() {}

	/**
	 * Get the start point.
	 * 
	 * \return 
	 */
	Point start_pt() const { return _cs[0].ctr(); }

	/**
	 * Get the end point.
	 * 
	 * \return 
	 */
	Point end_pt() const { return _cs[1].ctr(); }

	/**
	 * Get the Steiner Zones.
	 * 
	 * \return 
	 */
	vector<StnZone> stn_zones() const { return _szs; }

	/**
	 * Initialize this operator, calculate distance matrix of (shuffled) circles,
	 * clear stn zone vector and reinitialize _feas_idx.
	 *
	 * \param shuffle_circs: whether to shuffle circles.
	 * \return
	 */
	void init(const bool shuffle_circs);

	/**
	 * Construct the vector of Steiner Zones.
	 *
	 * \return
	 */
	void construct_stn_zones();

	/**
	 * The main body of the Randomized Steiner Zone Search operator.
	 *
	 * \return The Steiner Zone vector with minimum number of SZs.
	 */
	vector<StnZone> randomized_search();

private:
	vector<size_t> _feas_idx; // Index of feasible circles.
	vector<Circle> _cs;	      // Vector of shuffled circles.
	vector<StnZone> _szs;	  // Vector of constructed Steiner Zones.
	vector<vector<double>> _dists;		// Distance matrix of all circle centers
	const size_t _num_iters_szs = 0;
	const size_t _num_vtxs      = 0;  // Number of vertices in a singleton.
	const size_t _num_circs     = 0;
	const size_t _max_circs     = 0;  // Maximum circle number of a SZ.
	const size_t _num_depots    = 2;

	void _calc_dist_mat();
};


#endif // !STNZONESRCH_H
