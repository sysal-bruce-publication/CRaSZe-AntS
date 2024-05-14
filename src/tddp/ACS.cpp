#include "ACS.h"

using std::make_pair;

void ACS::_update_costs(const vector<NbhdWayPt>& way_pts)
{
	_costs = vector<vector<double>>(_num_nodes, vector<double>(_num_nodes));
	for (size_t i = 0; i < _num_nodes - 1; ++i) {
		for (size_t j = i + 1; j < _num_nodes; ++j) {
			_costs[i][j] = _costs[j][i] 
				= std::max((way_pts[j].pt - way_pts[i].pt).norm2() / _v_truck, 1e-5);
		}
	}
}

double ACS::_nn_nbr_ratio() const
{
	unique_ptr<Ant> nn_ant = make_unique<Ant>(_num_nodes);
	multimap<double, size_t> next_dists;
	while (!nn_ant->is_all_visit()) {
		next_dists.clear();
		for (size_t i = _num_depots; i < _num_nodes; i++) {
			if (nn_ant->is_node_visited(i)) { continue; }
			next_dists.insert(
				make_pair(_costs[nn_ant->cur_node][i] + _way_pts[i].c, i));
		}
		size_t node2visit = next_dists.begin()->second;
		double needed_cost = next_dists.begin()->first + _costs[node2visit][_end_idx];
		if (nn_ant->fit.cost + needed_cost > _bgt) { break; }
		nn_ant->cur_node = next_dists.begin()->second;
		nn_ant->fit.prize += _way_pts[nn_ant->cur_node].p;
		nn_ant->fit.cost += next_dists.begin()->first + _way_pts[nn_ant->cur_node].c;
		// Remove min_idx in unvisits
		nn_ant->mask_visit(nn_ant->cur_node);
	}
	nn_ant->fit.cost += _costs[nn_ant->cur_node][0];  // No extra cost for end node.
	return nn_ant->fit.prize / nn_ant->fit.cost;
}

double ACS::_path_cost(const vector<size_t>& path) const
{
	double total_cost = 0.;
	// extra cost at start and end is 0
	for (size_t i = 0; i < path.size() - 1; i++) {  
		total_cost += _costs[path[i]][path[i + 1]] + _way_pts[i].c;
	}
	return total_cost;
}

void ACS::init_colony(const bool is_init_phase, const vector<NbhdWayPt>& way_pts)
{
	_way_pts = way_pts;
	_num_nodes = way_pts.size();
	_num_ants = std::min(_num_ants, _num_nodes - _num_depots);
	// In ACS, we don't need original way-points, heuristic information is enough.
	_update_costs(way_pts);
	if (is_init_phase) { // If 1st iteration, use nearest neighbor distance as tau0.
		// Ref: Ant Colony System Paper: tau0 = 1 / (num_city * nn_nbr_cost)
		_tau0 = _nn_nbr_ratio() / _num_nodes;
		_taus = vector<vector<double>>(_num_nodes, vector<double>(_num_nodes, _tau0));
	}
	else {
		_gb_ant->fit.cost = _path_cost(_gb_ant->path);  // Update cost first.
		// Prize won't change, so we can directly use it to calculate tau0.
		_tau0 = _gb_ant->fit.prize / (_num_nodes * _gb_ant->fit.cost);
		_taus = vector<vector<double>>(_num_nodes, vector<double>(_num_nodes, _tau0));
		double tau_delta = _tau0 * (1 - _alpha) 
			+ _alpha * _gb_ant->fit.prize / _gb_ant->fit.cost;
		for (size_t j = 0; j < _gb_ant->path.size() - 1; ++j) {
			_taus[_gb_ant->path[j]][_gb_ant->path[j + 1]] = tau_delta;
		}
	}
	// Then create all ants
	_ants.reserve(_num_ants);
	for (size_t k = 0; k < _num_ants; ++k) { _ants.push_back(Ant(_num_nodes)); }
	_gb_ant = make_unique<GBAnt>();
}

size_t ACS::_state_transition_rule(const Ant& ant_k) const
{
	size_t node2visit = _start_idx; 
	multimap<double, size_t, std::greater<double>> prod_nodes_map;
	for (size_t j = _num_depots; j < _num_nodes; j++) {
		if (!ant_k.is_node_visited(j)) {
			double eta_rcpl = pow(_way_pts[j].p / _costs[ant_k.cur_node][j], _beta);
			prod_nodes_map.insert(make_pair(_taus[ant_k.cur_node][j] * eta_rcpl, j));
		}
	}
	unique_ptr<Rd> rand = make_unique<Rd>();
	if (rand->uniform_rand() <= _q0) { node2visit = prod_nodes_map.begin()->second; }
	else {
		vector<double> prob_js; prob_js.resize(_num_nodes);
		for (auto it = prod_nodes_map.begin(); it != prod_nodes_map.end(); ++it) {
			prob_js[it->second] = it->first;
		}
		double neighb_sum = Common::vec_sum(prob_js);
		std::transform(prob_js.begin(), prob_js.end(), prob_js.begin(),
			[neighb_sum](double& c) -> double { return c / neighb_sum; });
		rand->vec_biased_sample_one_idx(prob_js, prob_js, node2visit);
	}
	return node2visit;
}

bool ACS::_is_ant_path_feasible(const size_t ant_idx,
	const bool is_path_complete) const
{
	size_t last_node = _ants[ant_idx].path.back();
	double cost = _ants[ant_idx].fit.cost
		+ (is_path_complete ? 0 : _costs[last_node][_end_idx]);
	if (cost > _bgt + EPS) { return false; } // Allow equal
	else { return true; }
}

void ACS::_set_node_for_ant(const size_t node2visit, const size_t ant_idx)
{
	_ants[ant_idx].fit.cost += 
		_costs[_ants[ant_idx].cur_node][node2visit] + _way_pts[node2visit].c;
	_ants[ant_idx].fit.prize += _way_pts[node2visit].p;
	_ants[ant_idx].cur_node = node2visit;
	_ants[ant_idx].path.push_back(node2visit);
	_ants[ant_idx].mask_visit(node2visit);
	_ants[ant_idx].feasible = _is_ant_path_feasible(ant_idx, false);
}

void ACS::_construct_path_for_ant(const size_t ant_idx)
{
	while ((!_ants[ant_idx].is_all_visit()) && _ants[ant_idx].feasible) {
		size_t node2visit = _state_transition_rule(_ants[ant_idx]);
		_set_node_for_ant(node2visit, ant_idx);
	}
	// Add the last distance back to the end node
	_ants[ant_idx].fit.cost += _costs[_ants[ant_idx].cur_node][0];
	_ants[ant_idx].path.push_back(_end_idx);
	_ants[ant_idx].feasible = _is_ant_path_feasible(ant_idx, true);
}

void ACS::_2opt_for_ant(const size_t ant_idx)
{
	auto delta_cost = [this](const vector<size_t>& path,
		const size_t v1, const size_t v2) -> double
	{
		size_t path_len = path.size();
		size_t v1_next = (v1 + 1) % path_len, v2_next = (v2 + 1) % path_len;
		return _costs[path[v1]][path[v2]] - _costs[path[v1]][path[v1_next]]
			+ _costs[path[v1_next]][path[v2_next]] - _costs[path[v2]][path[v2_next]];
	};

	bool improved = true;
	size_t cnt = 0;
	size_t path_len = _ants[ant_idx].path.size();
	while (improved && cnt < _num_iters_2opt) {
		improved = false;
		for (size_t i = 0; i < path_len - 2; ++i) {
			for (size_t j = i + 1; j < path_len - 1; ++j) {
				double cost_diff = delta_cost(_ants[ant_idx].path, i, j);
				if (cost_diff < 0) {
					Common::vec_reverse(_ants[ant_idx].path, i, j);
					_ants[ant_idx].fit.cost += cost_diff;
					improved = true;
					cnt++;
				}
			}
		}
	}
}

void ACS::_iterative_drop(const size_t ant_idx, size_t& idx2drop, size_t& drop_node,
	multimap<double, IdxVal>& dropval_node_cost)
{
	// This is for the first drop value map update
	_ants[ant_idx].fit.prize -= _way_pts[drop_node].p;
	_ants[ant_idx].fit.cost -= dropval_node_cost.begin()->second.second;
	// Before erase the index, record its previous and next neighbor node.
	size_t idx0 = idx2drop - 1;
	_ants[ant_idx].path.erase(_ants[ant_idx].path.begin() + idx2drop);
	// ALTHOUGH the dropped node performs POOR in current path, we still
	// need to add it back to unvisited node list because it may be a good 
	// candidate after adding some other nodes. 
	_ants[ant_idx].unmask_visit(drop_node);
	_ants[ant_idx].feasible = _is_ant_path_feasible(ant_idx, true);
	if (_ants[ant_idx].feasible) { return; }
	// Otherwise we need to drop more nodes. First update the dropped nodes.
	dropval_node_cost.erase(dropval_node_cost.begin());
	while (!_ants[ant_idx].feasible) {
		// Traverse dropval_node_cost and update previous and next node.
		for (size_t v = 0; v < dropval_node_cost.size(); ++v) {
			auto it = dropval_node_cost.begin(); std::advance(it, v);
			double cost_diff = 0, drop_val = 0;
			// if the previous node to be updated
			size_t cur_node = _ants[ant_idx].path[idx0];
			if (it->second.first == cur_node) {  
				size_t prev_node = _ants[ant_idx].path[idx0 - 1];
				size_t next_node = _ants[ant_idx].path[idx0 + 1];
				cost_diff = _costs[prev_node][cur_node] + _costs[cur_node][next_node]
					+ _way_pts[cur_node].c - _costs[prev_node][next_node];
				drop_val = _way_pts[cur_node].p / cost_diff;
				dropval_node_cost.insert(
					make_pair(drop_val, make_pair(cur_node, cost_diff)));
				dropval_node_cost.erase(it);
				continue;  // The iterator has been deleted
			}
			// if the next node to be updated
			cur_node = _ants[ant_idx].path[idx2drop];
			if (it->second.first == cur_node) {
				size_t prev_node = _ants[ant_idx].path[idx0];
				size_t next_node = _ants[ant_idx].path[idx2drop + 1];
				cost_diff = _costs[prev_node][cur_node] + _costs[cur_node][next_node]
					+ _way_pts[cur_node].c - _costs[prev_node][next_node];
				drop_val = _way_pts[cur_node].p / cost_diff;
				dropval_node_cost.insert(
					make_pair(drop_val, make_pair(cur_node, cost_diff)));
				dropval_node_cost.erase(it);
			}
		}
		// Update where2drop and prev_idx for next usage.
		drop_node = dropval_node_cost.begin()->second.first;
		_ants[ant_idx].fit.prize -= _way_pts[drop_node].p;
		_ants[ant_idx].fit.cost -= dropval_node_cost.begin()->second.second;
		idx2drop = Common::find_index_by_member(_ants[ant_idx].path, drop_node);
		_ants[ant_idx].path.erase(_ants[ant_idx].path.begin() + idx2drop);
		idx0 = idx2drop - 1; // Update previous index
		_ants[ant_idx].unmask_visit(drop_node);
		_ants[ant_idx].feasible = _is_ant_path_feasible(ant_idx, true);
		// This node has been dropped, remove it
		dropval_node_cost.erase(dropval_node_cost.begin());
	}
}

void ACS::_drop_nodes_until_feasible(const size_t ant_idx)
{
	if (!_ants[ant_idx].feasible) {
		// drop value, {node, cost difference}
		multimap<double, IdxVal> dropval_node_cost;
		// Calculate the initial dropval_node_cost map
		for (size_t j = 1; j < _ants[ant_idx].path.size() - 1; j++) {
			size_t cur_node = _ants[ant_idx].path[j];
			size_t prev_node = _ants[ant_idx].path[j - 1];
			size_t next_node = _ants[ant_idx].path[j + 1];
			double cost_diff = _costs[cur_node][next_node] + _costs[prev_node][cur_node]
				+ _way_pts[cur_node].c - _costs[prev_node][next_node];
			cost_diff = std::max(EPS, cost_diff);  // Avoid i is on line (i-1)(i+1)
			dropval_node_cost.insert(make_pair(_way_pts[cur_node].p / cost_diff,
					make_pair(cur_node, cost_diff)));
		}
		size_t node2drop = dropval_node_cost.begin()->second.first;
		size_t idx2drop = Common::find_index_by_member(_ants[ant_idx].path, node2drop);
		// Here we HAVEN'T DONE ANY MODIFICATION FOR ant_k
		_iterative_drop(ant_idx, idx2drop, node2drop, dropval_node_cost);
	}
}

double ACS::_add_cost(const size_t adj0, 
	const size_t adj1, const size_t node2add) const
{
	return _costs[adj0][node2add] + _costs[node2add][adj1] + _way_pts[node2add].c 
		- _costs[adj0][adj1];
}

multimap<double, DualidxVal, std::greater<double>> ACS::_calc_init_addval_map(
	const size_t ant_idx) const
{
	typedef std::pair<DualIdx, DualIdx> DualPair;
	multimap<double, DualidxVal, std::greater<double>> addval_map;
	for (size_t v = _num_depots; v < _num_nodes; ++v) {
		if (!_ants[ant_idx].is_node_visited(v)) {
			// distance, {node, node index in path}
			multimap<double, DualIdx> cost_nbrs;
			// Find neighbor nodes of node, here we don't consider the 
			// neighborhood cost because here is to find the closest neighbor.
			for (size_t j = 0; j < _ants[ant_idx].path.size(); j++) {
				cost_nbrs.insert(make_pair(_costs[_ants[ant_idx].path[j]][v],
					make_pair(_ants[ant_idx].path[j], j)));
			}
			// {node, node index in path}
			vector<DualIdx> nbr_pairs; nbr_pairs.reserve(3);
			for (auto it = cost_nbrs.begin(); it != cost_nbrs.end(); ++it) {
				nbr_pairs.push_back(make_pair(it->second.first, it->second.second));
				if (nbr_pairs.size() > 2) { cost_nbrs.clear();  break; }
			}
			// add cost, {{nbr, adj}, {nbr index in path, adj index in path}}
			multimap<double, DualPair> cost_nbr_pairs;
			// Now fill neighbor node pairs, their indices in path and addcost
			// Note that we DO NOT check whether picked neighbors are adjacent
			// because addval_map will be ITERATIVELY computed ANYWAY, if we 
			// directly return adjacent neighbor pair, it will raise a problem later.
			for (size_t i = 0; i < nbr_pairs.size(); i++) {
				size_t nbr = nbr_pairs[i].first;
				size_t nbr_idx = nbr_pairs[i].second;
				// Here we DO NOT add pair of start and end node, i.e. ABANDON
				// start_node -> end_node AND end_node -> start_node
				// because the first and last node in path 
				// MUST be at start and end node.
				if (nbr_idx < _ants[ant_idx].path.size() - 1) {
					size_t next_idx = nbr_idx + 1;
					cost_nbr_pairs.insert(make_pair(
						_add_cost(nbr, _ants[ant_idx].path[next_idx], v), make_pair(
							make_pair(nbr, _ants[ant_idx].path[next_idx]),
							make_pair(nbr_idx, next_idx))));
				}
				if (nbr_idx > 0) {
					size_t prev_idx = nbr_idx - 1;
					cost_nbr_pairs.insert(make_pair(
						_add_cost(_ants[ant_idx].path[prev_idx], nbr, v), make_pair(
							make_pair(_ants[ant_idx].path[prev_idx], nbr),
							make_pair(prev_idx, nbr_idx))));
				}
			}
			nbr_pairs.clear();  // We don't need it any more.
			// Get the smallest addcost of this node v.
			double addcost_v = 1e6; size_t idx2add_v = 0;
			// If any two nodes adjacent in path, then find the pair w/
			// the SMALLEST addcost and find index to insert node v.
			addcost_v = cost_nbr_pairs.begin()->first;
			idx2add_v = std::max(  // left, node to add, right
				cost_nbr_pairs.begin()->second.second.first,
				cost_nbr_pairs.begin()->second.second.second);
			// If the addcost within budget constraint, add it into the addval_map
			if (addcost_v + _ants[ant_idx].fit.cost <= _bgt + EPS) {
				addval_map.insert(make_pair(_way_pts[v].p / addcost_v,
					make_pair(make_pair(v, idx2add_v), addcost_v)));
			}
		}
	}
	return addval_map;
}

void ACS::_iterative_add(const size_t ant_idx, size_t& idx2add, size_t& add_node,
	multimap<double, DualidxVal, std::greater<double>>& addval_map)
{
	// add_value, {{node2add, add_idx}, add_cost}
	vector<std::pair<double, DualidxVal>> map2update; 
	map2update.reserve(addval_map.size());
	while (add_node != _start_idx) { // true means there exists feasible nodes to add
		// ONLY update prize and cost FOR NOW
		_ants[ant_idx].fit.prize += _way_pts[add_node].p;
		_ants[ant_idx].fit.cost += addval_map.begin()->second.second;
		addval_map.erase(addval_map.begin());  // Erase added node
		// Update map2update	
		Common::erase_if(addval_map,
			[add_node, idx2add, this, ant_idx, &map2update](const auto& pair) -> bool
			{
				size_t unvis_i = pair.second.first.first;  // Get unvisited node
				// Update this unvis_i's addvalue
				multimap<double, DualidxVal, std::greater<double>> temp_map;
				size_t unvis_idx = pair.second.first.second;  // Get index to insert
				double unvis_p = _way_pts[unvis_i].p;
				// Only when the previous add node is NOT LOCATED ADJACENT to the 
				// current add node, then insert previous calculated addval_map.
				if (unvis_idx != idx2add && unvis_idx != idx2add + 1) {
					if (idx2add < unvis_idx) {
						temp_map.insert(make_pair(pair.first, make_pair(
							make_pair(unvis_i, unvis_idx + 1), pair.second.second)));
					}
					else { temp_map.insert(pair); }
				}
				// Otherwise, there exists two scenarios:
				// 1. Previous idx2add is located before or after current unvis_idx;
				// 2. Insertion of previous add node will lead to shorter addcost.
				// BOTH scenario will need to compute the addcost of inserting 
				// unvis_i before and after add_node. 
				double addcost_before_node2add = _add_cost(
					_ants[ant_idx].path[idx2add - 1], add_node, unvis_i);
				temp_map.insert(make_pair(unvis_p / addcost_before_node2add,
					make_pair(make_pair(unvis_i, idx2add), addcost_before_node2add)));
				double addcost_aftr_node2add = _add_cost(
					add_node, _ants[ant_idx].path[idx2add], unvis_i);
				temp_map.insert(make_pair(unvis_p / addcost_aftr_node2add,
					make_pair(make_pair(unvis_i, idx2add + 1), addcost_aftr_node2add)));
				// If meets the constraint, add it into addval_node_pair_cost and 
				// erase the iterator (because the value MAY not be updated)
				if (temp_map.begin()->second.second + _ants[ant_idx].fit.cost <= _bgt + EPS) {
					map2update.push_back(*temp_map.begin());
				}
				// If the smallest addcost cannot meet the budget 
				// OR it has been updated, erase the iterator.
				return true;
			});
		// Insert back updated unvisited node
		for (const auto& pair : map2update) { addval_map.insert(pair); }
		map2update.clear();  // Has done its job.
		map2update.reserve(addval_map.size());
		// This is updating the PREVIOUS node2add !
		Common::vec_insert(_ants[ant_idx].path, idx2add, add_node);
		_ants[ant_idx].mask_visit(add_node);
		_ants[ant_idx].feasible = _is_ant_path_feasible(ant_idx, true);
		// Then get the CURRENT node2add
		if ((!addval_map.empty()) && addval_map.begin()->first > 0) {
			add_node = addval_map.begin()->second.first.first;
			idx2add = addval_map.begin()->second.first.second;
		}
		else { add_node = _start_idx; }
	}
}

void ACS::_add_nodes_until_budget(const size_t ant_idx)
{
	if ((!_ants[ant_idx].is_all_visit()) && _ants[ant_idx].feasible) {
		size_t node2add = _start_idx, idx2add = _start_idx;
		// add value, {{node to add, path index to add}, add cost}
		multimap<double, DualidxVal, std::greater<double>> addval_map = 
			_calc_init_addval_map(ant_idx);
		if ((!addval_map.empty()) && addval_map.begin()->first > 0) {
			// Get the node v and index to add in current path (whose 
			// addval is positive and the largest).
			node2add = addval_map.begin()->second.first.first;
			idx2add = addval_map.begin()->second.first.second;
		}
		else { node2add = _start_idx; }
		// Then iteratively add nodes until no room to add.
		_iterative_add(ant_idx, idx2add, node2add, addval_map);
	}
}

void ACS::_local_updating_rule(const size_t i, const size_t j)
{
	_taus[i][j] = (1 - _rho) * _taus[i][j] + _rho * _tau0;
}

void ACS::_update_gb_ant()
{
	// If local best ant not feasible, we will see error results at the end.
	// But luckily I didn't find any :)
	double prize_diff = _ants[_lb_ant_id].fit.prize - _gb_ant->fit.prize;
	if (prize_diff > 0) {  // General ant higher prize
		_gb_ant->fit.cost = _ants[_lb_ant_id].fit.cost;
		_gb_ant->fit.prize = _ants[_lb_ant_id].fit.prize;
		_gb_ant->path = _ants[_lb_ant_id].path;
		size_t path_len = _gb_ant->path.size() - 1;
		_gb_ant->delta_tau = _alpha * _gb_ant->fit.prize / _gb_ant->fit.cost;
		if (prize_diff > _tol_acs) { _no_impr_cnt = 0; return; }
	}
	else if (abs(prize_diff) <= EPS) {
		double cost_diff = _gb_ant->fit.cost - _ants[_lb_ant_id].fit.cost;
		if (cost_diff > 0) {  // General ant less cost
			_gb_ant->fit.cost = _ants[_lb_ant_id].fit.cost;
			_gb_ant->fit.prize = _ants[_lb_ant_id].fit.prize;
			_gb_ant->path = _ants[_lb_ant_id].path;
			_gb_ant->delta_tau = _alpha * _gb_ant->fit.prize / _gb_ant->fit.cost;
			if (cost_diff > _tol_acs) { _no_impr_cnt = 0; return; }
		}
	}
	_no_impr_cnt++;
}

void ACS::_global_updating_rule()
{
	double evap_ratio = 1 - _alpha;
	for (size_t i = 0; i < _num_nodes; i++) {
		if (i == _end_idx) { continue; }  // Don't care tau[_end_idx][x]!
		for (size_t j = 1; j < _num_nodes; j++) {
			if (i != j) { _taus[i][j] *= evap_ratio; }
		}
	}
	for (size_t j = 0; j < _gb_ant->path.size() - 1; ++j) {
		_taus[_gb_ant->path[j]][_gb_ant->path[j+1]] += _gb_ant->delta_tau;
	}
}

int ACS::_is_gb_sol_valid() const
{
	double prize = 0, cost = 0;
	for (size_t i = 0; i < _gb_ant->path.size() - 1; i++) {
		size_t node = _gb_ant->path[i], next_node = _gb_ant->path[i + 1];
		prize += _way_pts[node].p;
		cost += _costs[node][next_node] + _way_pts[next_node].c;
	}
	if (cost > _bgt) {
		return -1;
	}
	if (abs(cost - _gb_ant->fit.cost) > EPS) {
		return -2;
	}
	if (abs(prize - _gb_ant->fit.prize) > EPS) {
		return -3;
	}
	return 0;
}

void ACS::evolve_until_stop_criterion()
{
	unique_ptr<Rd> rd = make_unique<Rd>();
	vector<int> range_nodes = vector<int>(_num_nodes - _num_depots);
	std::iota(range_nodes.begin(), range_nodes.end(), static_cast<int>(_num_depots));
	vector<int> start_nodes;
	for (size_t t = 0; t < _num_iters_acs; ++t) {
		if (_no_impr_cnt > _max_no_impr) { break; }
		Fitness lb_fit = Fitness(1e8, -1);
		// If same length, shuffle node indices
		if (_num_ants == _num_nodes - _num_depots) { 
			start_nodes = range_nodes; rd->vec_shuffle_skip_first_n(0, start_nodes);
		}
		// Otherwise, sample _num_ants start nodes from range_nodes.
		else { rd->vec_uniform_rand_sample(range_nodes, _num_ants, start_nodes); }
		for (size_t k = 0; k < _num_ants; k++) {
			// Set the first random start node.
			_set_node_for_ant(static_cast<size_t>(start_nodes[k]), k);
			// Construct the solution path.
			_construct_path_for_ant(k);
			if (_num_iters_2opt > 0) {  // Apply 2-opt to the solution path.
				_2opt_for_ant(k); 
				_ants[k].feasible = _is_ant_path_feasible(k, true);
			}
			// [DROP NODE PHASE] Keep dropping node until solution is feasible.
			_drop_nodes_until_feasible(k);
			// [ADD NODE PHASE] Keep adding node until no unvisited nodes,
			// reach budget constraint, or reach number of searches.
			_add_nodes_until_budget(k);
			// Determine the locally best ant for later updating globally best ant.
			if (_ants[k].fit > lb_fit) { lb_fit = _ants[k].fit; _lb_ant_id = k; }
			// Local updating rule, here the path is start -> path -> start.
			for (size_t j = 0; j < _ants[k].path.size() - 1; ++j) {
				_local_updating_rule(_ants[k].path[j], _ants[k].path[j + 1]);
			}
		}
		_update_gb_ant();			// Update globally best ant.
		_global_updating_rule();    // Global updating rule.
		if (t < _num_iters_acs - 1) {  // Lastly reset every individual ant.
			for (size_t k = 0; k < _num_ants; k++) { _ants[k].reset(); }
		}
	}
}

void ACS::inherit_acs()
{
	_lb_ant_id = _no_impr_cnt = 0;
	// Then clear the pheromone and distance matrix
	for (size_t i = 0; i < _num_nodes; i++) { _taus[i].clear(); _costs[i].clear(); }
	_taus.clear(); _taus.shrink_to_fit();
	_costs.clear(); _costs.shrink_to_fit();
	_ants.clear(); _ants.shrink_to_fit();
}
