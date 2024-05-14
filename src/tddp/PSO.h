/*****************************************************************//**
 * \file   PSO.h
 * \brief  Particle Swarm Optimization with Steiner Zone boundary constraints
 * 
 * \author Qiuchen Qian
 * \date   October 2023
 *********************************************************************/
#ifndef PSO_H
#define PSO_H

#include "Particle.h"
#include "DataIO.h"
#include "Timer.h"

class PSO
{
public:
	PSO() = default;
	PSO(const size_t num_ptcls, const size_t max_no_impr_pso, 
		const size_t num_iters_pso, const double c1,
		const double c2, const double inertia_min, const double inertia_max,
		const double tol_pso) : 
		_num_ptcls(num_ptcls), 
		_max_no_impr(max_no_impr_pso),	
		_num_iters_pso(num_iters_pso),
		_c1(c1), 
		_c2(c2), 
		_w_max(inertia_max), 
		_w_const((inertia_max - inertia_min) / num_iters_pso),
		_tol_pso(tol_pso)
	{}
	~PSO() {}

	vector<StnZone> get_stn_zones() const { return _stn_zones; }

	/**
	 * Get the final output of PSO.
	 *
	 * \param fitness
	 * \param best_pos
	 */
	void get_output(Fitness& fit, vector<Point>& best_pos) const;

	/**
	 * Create uninitialized swarm. This function will only push empty Particle
	 * with ACS parameters into the swarm.
	 * 
	 * \param num_ants
	 * \param num_iters_acs
	 * \param num_iters_2opt
	 * \param v_truck
	 * \param q0
	 * \param alpha
	 * \param beta
	 * \param rho
	 * \param tol_acs
	 */
	void create_swarm_with_acs(const size_t num_ants, const size_t num_iters_acs,
		const size_t num_iters_2opt, const size_t max_no_impr_acs, const double v_truck,
		const double q0, const double alpha, const double beta, const double rho,
		const double tol_acs, const double budget);

	/**
	 * Initialize the swarm with the start point and vector of Steiner Zones.
	 * For each particle, this function will initialize its position, velocity,
	 * neighborhood fitness, and update history best fitness and position. Finally,
	 * this function will update globally best particle.
	 *
	 * \param start_pt
	 * \param
	 * \param stn_zones
	 */
	void init_swarm_with_stn_zones(Point& start_pt, Point& end_pt,
		vector<StnZone>& stn_zones);

	/**
	 * Evolve until meet stop criterion.
	 * 
	 * \return
	 */
	void evolve_until_stop_criterion();

private:
	vector<StnZone> _stn_zones;
	vector<Particle> _ptcls;
	unique_ptr<GBPtcl> _gb_ptcl;
	unique_ptr<Timer> _timer;
	size_t _lb_ptcl_idx         = 0;	// Index of local best particle
	size_t _no_impr_cnt         = 0;	// Number of iterations without improvement
	size_t _num_waypts          = 0;	// Number of way-points
	Point _start_pt             = Point();	// Start depot coordinate
	Point _end_pt               = Point();  // End depot coordinate
	const size_t _num_ptcls     = 0;	// Number of particles
	const size_t _num_iters_pso = 0;	// Number of PSO iterations
	const size_t _max_no_impr   = 0;	// Maximum number of iterations without improvement
	const size_t _num_depots    = 2;
	const size_t _start_idx     = 0;
	const size_t _end_idx       = 1;
	const double _c1            = 0;	// Acceleration coefficient for individual
	const double _c2            = 0;	// Acceleration coefficient for swarm
	const double _w_max         = 0;	// Initial inertia weight
	const double _w_const       = 0;    // Constant for computing inertia
	const double _tol_pso       = 0;	// Tolerance for improvement
	const double _time_limit    = 595000;	 // 10 minutes main loop time limit [ms]

	/**
	 * Update global best particle's fitness, position, and path.
	 * 
	 * \return
	 */
	void _update_gb_ptcl();
};

#endif // !PSO_H
