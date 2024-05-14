/*****************************************************************//**
 * \file   ACS.h
 * \brief  Ant Colony System with inheritance mechanism.
 * 
 * \author Qiuchen Qian
 * \date   October 2023
 *********************************************************************/
#ifndef ACS_H
#define ACS_H

#include <memory>
#include <map>
#include "Ant.h"

typedef std::pair<size_t, double> IdxVal;
typedef std::pair<double, double> DualVal;
typedef std::pair<DualIdx, double> DualidxVal;
using std::unique_ptr, std::make_unique;
using std::multimap;

struct NbhdWayPt
{
	NbhdWayPt(const double cost, const double prize, const Point& coord) :
		c(cost), p(prize), pt(coord) {}
	double c = 0;			// Cost of overlapped neighborhood
	double p = 0;			// Prize of overlapped neighborhood
	Point pt = Point();		// Coordinate of this way point
};

class ACS
{
public:
	ACS() = default;
	ACS(const size_t num_ants, const size_t max_no_impr, const size_t num_iters_acs, 
		const size_t num_iters_2opt, const double q0,  const double alpha, 
		const double beta, const double rho, const double tol_acs, 
		const double budget, const double v_truck) :
		_num_ants(num_ants), 
		_num_iters_acs(num_iters_acs), 
		_num_iters_2opt(num_iters_2opt),
		_max_no_impr(max_no_impr),
		_q0(q0), 
		_alpha(alpha), 
		_beta(beta),
		_rho(rho),
		_tol_acs(tol_acs),
		_bgt(budget), 
		_v_truck(v_truck) {}
	~ACS() {}

	/**
	 * Get best cost based on Euclidean distance.
	 * 
	 * \return 
	 */
	Fitness best_fit() const { return _gb_ant->fit; }

	/**
	 * Get best path sequence.
	 * 
	 * \return 
	 */
	vector<size_t> best_path() const { return _gb_ant->path; }

	/**
	 * Based on whether this is the first iteration, (re-)initialize the entire
	 * colony. If 1st, use nearest neighbor cost as guidance, otherwise, use the 
	 * global best ant's path as guidance.
	 * 
	 * \param is_init_phase
	 * \param way_pts: Vector of NbhdWayPt from PSO.
	 */
	void init_colony(const bool is_init_phase, const vector<NbhdWayPt>& way_pts);

	/**
	 * Evolve the ACS until meets the stopping criterion. 
	 * 
	 * \return
	 */
	void evolve_until_stop_criterion();

	/**
	 * Inherit the global best path as the guidance for next iteration. And
	 * reset the pheromone trail, distance matrix, etc.
	 * 2. No need to reset tau0 because any very rough approximation of the 
	 * optimal tour length would suffice.
	 * 3. No need to reset _num_waypts because it won't change.
	 * 5. Do not reset _gb_ant because we'll use it to initialize the pheromone 
	 * matrix at next iteration. And in the next iteration's init_colony(), it
	 * will be reset.
	 * 
	 *
	 * \return
	 */
	void inherit_acs();

private:
	size_t _lb_ant_id   = 0;		// The index of local best Ant
	size_t _num_nodes   = 0;		// Number of nodes, including the depot
	size_t _no_impr_cnt = 0;		// Number of iterations without improvement
	size_t _num_ants    = 0;		// Number of ants, may change based on number of target nodes
	double _tau0        = 0;		// The initial pheromone 
	unique_ptr<GBAnt> _gb_ant;		// Global best ant
	vector<Ant> _ants;				// Vector of all Ants
	// Vector of way-points from PSO, 0 is start, 1 is end, 2... are targets.
	vector<NbhdWayPt> _way_pts;
	vector<vector<double>> _costs;  // Cost matrix
	vector<vector<double>> _taus;	// Pheromone matrix

	const size_t _num_iters_acs  = 0;	// Number of iterations of ACS
	const size_t _num_iters_2opt = 0;	// Number of iterations of 2-opt
	const size_t _max_no_impr    = 0;   // Maximum number of iterations without improvement
	const size_t _num_depots     = 2;	// Number of depot nodes, OP is 2
	const size_t _start_idx      = 0;	// Start node index in _sz_p_pt
	const size_t _end_idx        = 1;	// End node index in _sz_p_pt
	const double _q0             = 0;	// Probability of choosing the best next node
	const double _alpha          = 0;	// Pheromone importance
	const double _beta           = 0;	// Heuristic importance
	const double _rho            = 0;	// Pheromone evaporation rate
	const double _tol_acs        = 0;	// Tolerance of improvement
	const double _bgt            = 0;	// Budget
	const double _v_truck        = 0;	// Truck speed

	/**
	 * Update the cost matrix with given way-points. 
	 * cost := dist[i][j] / v_truck;
	 * 
	 * \param way_pts
	 */
	void _update_costs(const vector<NbhdWayPt>& way_pts);

	/**
	 * Get the prize-cost ratio by greedy nearest neighbor strategy.
	 * 
	 * \return 
	 */
	double _nn_nbr_ratio() const;

	/**
	 * Get the path cost.
	 * 
	 * \param path
	 * \return 
	 */
	double _path_cost(const vector<size_t>& path) const;

	/**
	 * ACS State Transition Rule.
	 * If q <= q0, exploitation
	 * Else, biased exploration
	 *
	 * \param ant_idx
	 * \return The picked node to visit
	 */
	size_t _state_transition_rule(const Ant& ant_k) const;

	/**
	 * Check if the ant's path is feasible.
	 *
	 * \param ant_idx
	 * \param is_path_complete
	 * \return
	 */
	bool _is_ant_path_feasible(const size_t ant_idx,
		const bool is_path_complete) const;

	/**
	 * Visit a specific node. Update Ant's cost, cur_node, path, and visit_mask.
	 * 
	 * \param node2visit
	 * \param ant_idx
	 */
	void _set_node_for_ant(const size_t node2visit, const size_t ant_idx);

	/**
	 * Construct a feasible path for the ant.
	 *
	 * \param ant_idx
	 */
	void _construct_path_for_ant(const size_t ant_idx);

	/**
	 * Apply 2-opt search operator for the ant.
	 *
	 * \param ant_idx
	 */
	void _2opt_for_ant(const size_t ant_idx);

	/**
	 * Iterative drop process after obtained drop value map. Note
	 * The given dropval_node_cost has ALREADY DELETED the 1ST dropped node.
	 *
	 * \param ant_idx: Ant index
	 * \param idx2drop: At which index in path should be dropped
	 * \param drop_node: Which node exactly should be dropped
	 * \param dropval_node_cost multimap of drop value, {node, cost_diff}
	 */
	void _iterative_drop(const size_t k, size_t& idx2drop, size_t& drop_node,
		multimap<double, IdxVal>& dropval_node_cost);

	/**
	 * For ant k's solution, keep dropping low-profit nodes until the solution
	 * is feasible, the metric of each node in the solution is defined as:
	 *
	 * drop[v_i] = p[i] / (d[i-1][i] + d[i][i+1] + nbhd_cost[i] - d[i-1][i+1])
	 *
	 * Reference:
	 * An efficient evolutionary algorithm for the orienteering problem
	 *
	 * Note that:
	 * 1. This function WILL DO an in-place modification for passed
	 * Ant object;
	 * 2. This function WILL update ant_k's solution feasibility;
	 * 3. The dropped node will be added back to ant's visit mask.
	 *
	 * \param ant_idx
	 */
	void _drop_nodes_until_feasible(const size_t ant_idx);

	/**
	 * Get the add cost based on node index.
	 * 
	 * \param adj0
	 * \param adj1
	 * \param node2add
	 * \return 
	 */
	double _add_cost(const size_t adj0, 
		const size_t adj1, const size_t node2add) const;

	/**
	 * Compute the first addvalue_nodepaircost map for the ant.
	 *
	 * \param ant_idx
	 */
	multimap<double, DualidxVal, std::greater<double>> _calc_init_addval_map(
		const size_t ant_idx) const;

	/**
	 * For ant k's solution, keep adding high-profit nodes until the path cannot
	 * add any nodes any more. The metric of each node in the solution is defined:
	 *
	 * add[v_i] = p[i] / (d[i-1][i] + d[i][i+1] + nbhd_cost[i] - d[i-1][i+1])
	 *
	 * Reference:
	 * An efficient evolutionary algorithm for the orienteering problem
	 *
	 * \param k
	 * \param idx2add
	 * \param add_node
	 * \param addval_map: {Addval, {{add_node, idx2add}, addcost}}
	 */
	void _iterative_add(const size_t ant_idx, size_t& idx2add, size_t& add_node,
		multimap<double, DualidxVal, std::greater<double>>& addval_map);

	/**
	 * Decide where to insert a node in the path,
	 * and which node to be inserted.
	 *
	 * Note that this function WILL do an in-place modification
	 * for the Ant object.
	 *
	 * \param ant_idx
	 */
	void _add_nodes_until_budget(const size_t ant_idx);

	/**
	 * Local updating rule:
	 * tau(i, j) = (1 - rho) * tau(i, j) + rho * tau0
	 * for each ant.
	 *
	 * \param i
	 * \param j
	 */
	void _local_updating_rule(const size_t i, const size_t j);

	/**
	 * Update globally best ant's cost, path, and record
	 * pheromone trail for later global updating rule.
	 *
	 * \return
	 */
	void _update_gb_ant();

	/**
	 * Global updating rule:
	 * tau(i, j) = (1 - alpha) * tau(i, j) + alpha / gb_path_cost
	 *
	 * \return
	 */
	void _global_updating_rule();

	/**
	 * Check whether the global best solution is valid.
	 *
	 * \return: 0:correct; -1:exceed budget; -2:cost not match; -3:prize not match.
	 */
	int _is_gb_sol_valid() const;
};

#endif // !ACS_H
