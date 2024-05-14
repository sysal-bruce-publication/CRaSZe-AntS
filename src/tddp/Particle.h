#ifndef PARTICLE_H
#define PARTICLE_H

#include "DataIO.h"
#include "ACS.h"

struct GBPtcl
{
	GBPtcl() = default;

	Fitness hb_fit = Fitness(1e8, -1);
	vector<size_t> hb_path;
	vector<Point> hb_pos;
};

using std::shared_ptr, std::make_shared;

class Particle
{
public:
	Particle() = default;
	Particle(const size_t num_ants, const size_t num_iters_acs, 
		const size_t num_iters_2opt, const size_t max_no_impr_acs, 
		const double q0, const double alpha, const double beta, const double rho, 
		const double tol_acs, const double budget, const double v_truck)
	{
		_acs = make_shared<ACS>(num_ants, max_no_impr_acs, num_iters_acs, 
			num_iters_2opt, q0, alpha, beta, rho, tol_acs, budget, v_truck);
	}
	~Particle() {}

	Fitness fit    = Fitness();			// Fitness of current iteration
	Fitness hb_fit = Fitness(1e8, -1);	// History best fitness
	vector<size_t> path;				// Path sequence 
	vector<size_t> hb_path;				// History best path sequence
	vector<Point> pos;					// Position
	vector<Point> hb_pos;				// History best position
	vector<Point> vel;					// Velocity

	/**
	 * Apply ACS solver to evaluate the fitness of this particle. Given particle's
	 * way-points, the IACS solver will find the best path sequence and its fitness.
	 * Note that this function must be called after way-point fitness update.
	 * 
	 * \param is_init_phase: Whether in initialization phase.
	 */
	void update_seq_fit(const size_t idx, const bool is_init_phase, 
		const vector<NbhdWayPt>& nbhd_waypts)
	{
		_acs->init_colony(is_init_phase, nbhd_waypts);
		_acs->evolve_until_stop_criterion();
		
		//vector<double> gb_fit_hist = _acs->gb_cost_history();
		//DataIO::write_ant_fit(fs::current_path() / "data" / "cetsp",
		//	"team5_499.cetsp", idx, gb_fit_hist);
		
		_acs->inherit_acs();
		fit = _acs->best_fit();
		path = _acs->best_path();
	}

private:
	shared_ptr<ACS> _acs;
};

#endif // !PARTICLE_H

