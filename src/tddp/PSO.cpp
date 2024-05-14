#include "PSO.h"

void PSO::get_output(Fitness& fit, vector<Point>& best_pos) const
{
	fit = _gb_ptcl->hb_fit;
	size_t path_len = _gb_ptcl->hb_path.size();
	if (!best_pos.empty()) { best_pos.clear(); } best_pos.reserve(_num_waypts);
	best_pos.push_back(_start_pt);
	for (size_t j = 1; j < path_len - 1; ++j) {
		// hb_pos doesn't contain start and end nodes. But hb_path contains
		// start and end nodes. So here we need to subtract _num_depots.
		best_pos.push_back(_gb_ptcl->hb_pos[_gb_ptcl->hb_path[j] - _num_depots]);
	}
	best_pos.push_back(_end_pt);
}

void PSO::create_swarm_with_acs(const size_t num_ants, const size_t num_iters_acs,
	const size_t num_iters_2opt, const size_t max_no_impr_acs, const double v_truck,
	const double q0, const double alpha, const double beta, const double rho, 
	const double tol_acs, const double budget)
{
	_ptcls.reserve(_num_ptcls);
	for (size_t k = 0; k < _num_ptcls; ++k) {
		_ptcls.push_back(Particle(num_ants, num_iters_acs, num_iters_2opt, 
			max_no_impr_acs, q0, alpha, beta, rho, tol_acs, budget, v_truck));
	}
	_gb_ptcl = make_unique<GBPtcl>();
}

void PSO::_update_gb_ptcl()
{
	double prize_diff = _ptcls[_lb_ptcl_idx].hb_fit.prize - _gb_ptcl->hb_fit.prize;
	if (prize_diff > 0) { // General ptcl has higher prize
		_gb_ptcl->hb_fit = _ptcls[_lb_ptcl_idx].hb_fit;
		_gb_ptcl->hb_path = _ptcls[_lb_ptcl_idx].hb_path;
		_gb_ptcl->hb_pos = _ptcls[_lb_ptcl_idx].hb_pos;
		if (prize_diff > _tol_pso) { _no_impr_cnt = 0; return; }
	}
	else if (abs(prize_diff) <= EPS) {
		double cost_diff = _gb_ptcl->hb_fit.cost - _ptcls[_lb_ptcl_idx].hb_fit.cost;
		if (cost_diff > 0) {  // General ptcl has less cost
			_gb_ptcl->hb_fit = _ptcls[_lb_ptcl_idx].hb_fit;
			_gb_ptcl->hb_path = _ptcls[_lb_ptcl_idx].hb_path;
			_gb_ptcl->hb_pos = _ptcls[_lb_ptcl_idx].hb_pos;
			if (cost_diff > _tol_pso) { _no_impr_cnt = 0; return; }
		}
	}
	_no_impr_cnt++;
}

void PSO::init_swarm_with_stn_zones(Point& start_pt, Point& end_pt,
	vector<StnZone>& stn_zones)
{
  _timer = make_unique<Timer>(); _timer->tick();
	_start_pt = std::move(start_pt); _end_pt = std::move(end_pt);
	_stn_zones = std::move(stn_zones);
	_num_waypts = _stn_zones.size() + _num_depots;
	Fitness lb_hb_fit = Fitness(1e8, -1);  // Locally best history-best fitness
	for (size_t k = 0; k < _num_ptcls; ++k) {
		_ptcls[k].pos.reserve(_num_waypts - _num_depots);
		// Create neighborhood way-points
		vector<NbhdWayPt> nbhd_waypts; nbhd_waypts.reserve(_num_waypts);
		nbhd_waypts.push_back(NbhdWayPt(0, 0, _start_pt));
		nbhd_waypts.push_back(NbhdWayPt(0, 0, _end_pt));
		for (size_t j = 0; j < _stn_zones.size(); j++) {
			Point next_pos = _stn_zones[j].rand_pt();
 			_ptcls[k].pos.push_back(next_pos);
			nbhd_waypts.push_back(NbhdWayPt(
				_stn_zones[j].nbhd_c(next_pos), _stn_zones[j].p(), next_pos));
		}
		_ptcls[k].vel = vector<Point>(_num_waypts - _num_depots);  
		// Pass the vector of NhbdWayPt s to update sequence fitness, 
		// note this is the 1st iteration. Thus, ACS use nearest neighbor 
		// to determine the initial pheromone matrix.
		_ptcls[k].update_seq_fit(k, true, nbhd_waypts);
		// Because it's 1st iteration, here directly assign value.
		_ptcls[k].hb_fit = _ptcls[k].fit;
		_ptcls[k].hb_pos = _ptcls[k].pos;
		_ptcls[k].hb_path = _ptcls[k].path;

		// Lastly update the index of particle whose history best fitness is 
		// best in this iteration, i.e. locally best hb_fit.
		if (_ptcls[k].hb_fit > lb_hb_fit) {
			lb_hb_fit = _ptcls[k].hb_fit; _lb_ptcl_idx = k;
		}
	}
	// Clear internal distance matrix
	for (StnZone& sz : _stn_zones) { sz.free_memory(); }
	_update_gb_ptcl();
}

void PSO::evolve_until_stop_criterion()
{
	unique_ptr<Rd> rd = make_unique<Rd>();
	for (size_t t = 0; t < _num_iters_pso; ++t) {
	  
		if (_no_impr_cnt > _max_no_impr) { break; }
		
		Fitness lb_hb_fit = Fitness(1e8, -1);
		double w = _w_max - t * _w_const;
		for (size_t k = 0; k < _num_ptcls; ++k) {
			//! If want to set computing time limit, uncomment the following lines.
			//_timer->tock();
			//if (_timer->duration().count() > _time_limit) { break; }

			double lr1 = _c1 * rd->uniform_rand(), lr2 = _c2 * rd->uniform_rand();
			vector<NbhdWayPt> nbhd_waypts; nbhd_waypts.reserve(_num_waypts);
			nbhd_waypts.push_back(NbhdWayPt(0, 0, _start_pt));
			nbhd_waypts.push_back(NbhdWayPt(0, 0, _end_pt));
			for (size_t j = 0; j < _ptcls[k].pos.size(); ++j) {
				// The attraction velocity caused by individual best particle
				Point v_ib = (_ptcls[k].hb_pos[j] - _ptcls[k].pos[j]) * lr1;
				// The attraction velocity caused by global best particle
				Point v_gb = (_gb_ptcl->hb_pos[j] - _ptcls[k].pos[j]) * lr2;
				// The velocity to push the particle to the position at the next iteration.
				Point push_vel = _ptcls[k].vel[j] * w + v_ib + v_gb;
				// Restrict the current velocity first
				_stn_zones[j].restrict_vel(push_vel);
				// Then restrict next iteration's position in case it's out of the SZ.
				Point next_pos = Point();
				_stn_zones[j].restrict_pos(_ptcls[k].pos[j], push_vel, next_pos);
				// Then update pos and vel.
				_ptcls[k].pos[j] = next_pos; _ptcls[k].vel[j] = push_vel;
				nbhd_waypts.push_back(NbhdWayPt(
					_stn_zones[j].nbhd_c(next_pos), _stn_zones[j].p(), next_pos));
			}

			// Pass the vector of NhbdWayPt s to update sequence fitness, now it's
			// not the initialization phase, so we inherit ACS to update pheromone.
			_ptcls[k].update_seq_fit(k, false, nbhd_waypts);
			// Then update history best information.
			if (_ptcls[k].fit > _ptcls[k].hb_fit) {
				_ptcls[k].hb_fit = _ptcls[k].fit;
				_ptcls[k].hb_pos = _ptcls[k].pos;
				_ptcls[k].hb_path = _ptcls[k].path;
			}
			// Lastly update the index of particle whose history best fitness is 
			// best in this iteration, i.e. local-best hist_best_fit.
			if (_ptcls[k].hb_fit > lb_hb_fit) {
				lb_hb_fit = _ptcls[k].hb_fit; _lb_ptcl_idx = k;
			}
		}
		// At the end of this iteration, update global best particle.
		_update_gb_ptcl();
	}
}
