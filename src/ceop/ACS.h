#ifndef ACS_H
#define ACS_H

#include <map>
#include <string>
#include "Ant.h"
//#include "StnZone.h"
#include "DataIO.h"

typedef std::pair<size_t, double> IdxValPair;
typedef std::pair<DualIdx, double> DualidxVal;

using std::multimap;
using std::greater;

class ACS
{
public:
  ACS() = default;
  ACS(const size_t num_ants, const size_t max_no_impr, const size_t num_iters_acs, 
      const size_t num_iters_2opt, const size_t num_iters_arc, const double budget, 
      const double q0, const double alpha, const double beta, const double rho, 
      const double tol_acs) : 
    _num_ants(num_ants), 
    _max_no_impr(max_no_impr),
    _num_iters_acs(num_iters_acs), 
    _num_iters_2opt(num_iters_2opt),
    _num_iters_arc(num_iters_arc), 
    _bgt(budget), 
    _q0(q0), 
    _alpha(alpha), 
    _beta(beta), 
    _rho(rho), 
    _tol_acs(tol_acs) {}
  ~ACS() {}

	/**
	 * Get best cost based on Euclidean distance.
	 *
	 * \return
	 */
	double best_cost() const { return _gb_ant->cost; }

	/**
	 * Get best path sequence.
	 *
	 * \return
	 */
	vector<size_t> best_path() const { return _gb_ant->path; }

	vector<IdxValPt> sz_p_pt() const { return _sz_p_pt; }

	/**
	 * Initialize the entire colony.
	 *
	 * \param start: The start Point
	 * \param end: The end Point
	 * \param stn_zones: Vector of Steiner Zones, excluding the start and end node.
	 */
	void init_colony(const Point& start, const Point& end, vector<StnZone>& stn_zones);

	/**
	 * Evolve the ACS until meets the stopping criterion.
	 *
	 * \return
	 */
	void evolve_until_stop_criterion();


	/**
	 * Get final information.
	 *
	 * \param best_prize
	 * \param best_cost
	 * \param best_pts
	 */
	void get_gb_info(double& best_prize,
		double& best_cost, vector<Point>& best_pts) const;

	/**
	 * Apply arc search to SOP solution.
	 * 
	 * \return
	 */
	void arc_search(double& best_prize, double& best_cost,
		vector<Point>& best_pts) const;

private:
	size_t _lb_ant_id   = 0;		// The index of local best Ant
	size_t _num_nodes	= 0;		// Number of way points, including the depot
	size_t _no_impr_cnt = 0;		// Number of iterations without improvement
	size_t _num_ants    = 0;		// Number of ants, may change based on number of target nodes
	size_t _num_groups  = 0;
	size_t _num_targets = 0;
	double _tau0        = 0;		// The initial pheromone 
	unique_ptr<GBAnt> _gb_ant;		// Global best ant
	vector<Ant> _ants;				// Vector of all Ants
	vector<StnZone> _stn_zones;		// Vector of StnZones
	vector<IdxValPt> _sz_p_pt;		// Vector of StnZone pairs and their probabilities
	vector<vector<double>> _dists;  // Distance matrix
	vector<vector<double>> _taus;	// Pheromone matrix

	const size_t _num_iters_acs  = 0;	// Number of iterations of ACS
	const size_t _num_iters_2opt = 0;	// Number of iterations of 2-opt
	const size_t _num_iters_arc  = 0;	// Number of iterations of arc search
	const size_t _max_no_impr    = 0;  // Maximum number of iterations without improvement
	const size_t _num_depots     = 2;	// Number of depot nodes, OP is 2
	const size_t _start_idx      = 0;	// Start node index in _sz_p_pt
	const size_t _end_idx        = 1;	// End node index in _sz_p_pt
	const double _bgt            = 0;	// Budget
	const double _q0             = 0;	// Probability of choosing the best next node
	const double _alpha          = 0;	// Pheromone importance
	const double _beta           = 0;	// Heuristic importance
	const double _rho            = 0;	// Pheromone evaporation rate
	const double _tol_acs        = 0;	// Tolerance of improvement

	/**
	 * Update _dists of _sz_p_pt.
	 * 
	 * \return
	 */
	void _update_dists();

	/**
	 * Get the prize / cost by greedy nearest neighbor strategy.
	 * 
	 * \return
	 */
	double _nn_nbr_rate() const;

	/**
	 * With given SZ ID, random sample its corresponding index in _sz_p_pt.
	 * 
	 * \param sz_id
	 * \return 
	 */
	int _rand_sample_start_node(const size_t sz_id) const;

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
	 * ACS State Transition Rule.
	 * If q <= q0, exploitation
	 * Else, biased exploration
	 *
	 * \param ant_idx
	 * \return The picked node to visit
	 */
	size_t _state_transition_rule(const Ant& ant_k) const;

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
		multimap<double, IdxValPair>& dropval_node_cost);

	/**
	 * For ant k's solution, keep dropping low-profit nodes until the solution
	 * is feasible, the metric of each node in the solution is defined as:
	 * 
	 * drop[v_i] = p[i] / (d[i-1][i] + d[i][i+1] - d[i-1][i+1])
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

	double _addcost(const size_t adj0, const size_t adj1, 
		const size_t node2add) const;

	/**
	 * Compute the first addvalue_nodepaircost map for the ant.
	 *
	 * \param ant_idx
	 */
	multimap<double, DualidxVal, greater<double>> _calc_init_addval_map(
		const size_t ant_idx) const;

	/**
	 * For ant k's solution, keep adding high-profit nodes until the path cannot
	 * add any nodes any more. The metric of each node in the solution is defined:
	 *
	 * add[v_i] = p[i] / (d[i-1][i] + d[i][i+1] - d[i-1][i+1])
	 *
	 * Reference:
	 * An efficient evolutionary algorithm for the orienteering problem
	 * 
	 * \param k
	 * \param idx2add
	 * \param add_node
	 * \param addval_map: {Addval, {{add_node, idx2add}, addcost}}  
	 */
	void _iterative_add(const size_t k, size_t& idx2add, size_t& add_node,
		multimap<double, DualidxVal, greater<double>>& addval_map);

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
	 * \return: 0:correct; -1:exceed budget; -2:cost not match; -3:prize not match;
	 * -4:exist same SZ vertices in path
	 */
	int _is_gb_sol_valid() const;

	/**
	 * One iteration of arc search. Note that arc search won't change prize.
	 * 
	 * \param is_1st_iter
	 * \param waypts: Vector of SOP solution way-points
	 * \param cost
	 * \param is_waypt_changed
	 */
	void _arc_srch_one_iter(const bool is_1st_iter, const vector<size_t>& path, 
		vector<Point>& waypts, double& cost, bool& is_waypt_changed) const;

	/**
	 * Continuous version of addcost.
	 * 
	 * \param adj0
	 * \param adj1
	 * \param pt2add
	 * \return 
	 */
	double _addcost(const Point& adj0, const Point& adj1,
		const Point& pt2add) const;

	void _mask_visit(const size_t node_idx, vector<bool>& visit_mask) const;

	/**
	 * Continuous version of _calc_first_addval_map. Here the point to add
	 * has NOT BEEN added yet, so in this function we don't need path_pts,
	 * visit_mask, etc.
	 *
	 * \param path_pts
	 * \param
	 */
	void _calc_addval_map(const vector<Point>& waypts, 
		const vector<bool>& visit_mask, const double init_cost, 
		multimap<double, DualidxVal, greater<double>>& addval_map) const;

	void _add_nodes_until_budget(double& path_cost, double& path_prize,
		vector<size_t>& path, vector<bool>& visit_mask, vector<Point>& waypts) const;
};

#endif // !ACS_H
