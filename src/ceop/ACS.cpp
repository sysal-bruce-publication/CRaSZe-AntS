#include "ACS.h"

using std::make_pair;

void ACS::_update_dists()
{
	_dists = vector<vector<double>>(_num_nodes, vector<double>(_num_nodes));
	for (size_t i = 0; i < _num_nodes - 1; i++) {
		for (size_t j = i + 1; j < _num_nodes; j++) {
			// No need to calculate the distance between nodes in the same SZ.
			if (_sz_p_pt[i].first == _sz_p_pt[j].first) continue;
			double d = (_sz_p_pt[i].second.second - _sz_p_pt[j].second.second).norm2();
			_dists[i][j] = _dists[j][i] = d;
		}
	}
}

double ACS::_nn_nbr_rate() const
{
	unique_ptr<Ant> nn_ant = make_unique<Ant>(_num_nodes, _num_groups + _num_depots);
	nn_ant->feasible = true;
	multimap<double, size_t> next_dists;
	while (!nn_ant->is_all_visit()) {
		next_dists.clear();
		for (size_t i = _num_depots; i < _num_nodes; i++) {
			if (nn_ant->is_node_visit(i)) { continue; }
			next_dists.insert(make_pair(_dists[nn_ant->cur_node][i], i));
		}
		size_t node2visit = next_dists.begin()->second;
		double visit_cost = next_dists.begin()->first;
		if (nn_ant->cost + visit_cost + _dists[node2visit][_end_idx] > _bgt) { 
			break; 
		}
		nn_ant->cur_node = next_dists.begin()->second;
		nn_ant->prize += _sz_p_pt[nn_ant->cur_node].second.first;
		nn_ant->cost += next_dists.begin()->first;
		// Remove min_idx in unvisits
		nn_ant->mask_visit(nn_ant->cur_node, _sz_p_pt);
	}
	nn_ant->cost += _dists[nn_ant->cur_node][0];
	return nn_ant->prize / nn_ant->cost;
}

void ACS::init_colony(const Point& start, const Point& end, vector<StnZone>& stn_zones)
{
	_stn_zones = std::move(stn_zones);
	_num_groups = _stn_zones.size();
	_sz_p_pt.reserve(_stn_zones.front().num_max_circs() * _num_groups + _num_depots);
	// First push two depots.
	_sz_p_pt.push_back(make_pair(_start_idx, make_pair(0, start)));
	_sz_p_pt.push_back(make_pair(_end_idx,make_pair(0, end)));
	for (const StnZone& sz : _stn_zones) {
		for (const Point& v : sz.vs()) {
			_sz_p_pt.push_back(make_pair(sz.id(), make_pair(sz.p(), v)));
		}
	}
	_sz_p_pt.shrink_to_fit();
	// Distance matrix and pheromone matrix
	_num_nodes = _sz_p_pt.size();  // All vertices and depots
	_update_dists();
	_tau0 = _nn_nbr_rate() / (_num_groups + _num_depots);
	_taus = vector<vector<double>>(
		_num_nodes, vector<double>(_num_nodes, _tau0));
	// Initialize ants
	_num_ants = std::min(_num_ants, _num_groups);
	_ants.reserve(_num_ants);
	for (size_t k = 0; k < _num_ants; k++) { 
		_ants.push_back(Ant(_num_nodes, _num_groups + _num_depots)); 
	}
	_gb_ant = make_unique<GBAnt>();
}

int ACS::_rand_sample_start_node(const size_t sz_id) const
{
	size_t left = _num_depots;  // Exclude start node and end node
	// group_id stands for steiner zone _id, which starts from 2
	for (size_t i = _num_depots; i < sz_id; i++) {
		left += _stn_zones[i - _num_depots].num_vetxs();
	}
	unique_ptr<Rd> rd = make_unique<Rd>();
	return rd->uniform_rand(static_cast<int>(left), 
		static_cast<int>(left + _stn_zones[sz_id - _num_depots].num_vetxs() - 1));
}

size_t ACS::_state_transition_rule(const Ant& ant_k) const
{
	size_t node2visit = _start_idx;  
	std::multimap<double, size_t, greater<double>> prod_nodes_map;
	for (size_t i = _num_depots; i < _num_nodes; i++) {
		if (!ant_k.is_node_visit(i)) {
			double eta = pow(_sz_p_pt[i].second.first / _dists[ant_k.cur_node][i], _beta); 
			prod_nodes_map.insert(make_pair(_taus[ant_k.cur_node][i] * eta, i));
		}
	}
	unique_ptr<Rd> rd = make_unique<Rd>();
	if (rd->uniform_rand() <= _q0) { node2visit = prod_nodes_map.begin()->second; }
	else {
		vector<double> prob_js; prob_js.resize(_num_nodes);
		for (auto it = prod_nodes_map.begin(); it != prod_nodes_map.end(); ++it) {
			prob_js[it->second] = it->first;
		}
		double neighb_sum = Common::vec_sum(prob_js);
		std::transform(prob_js.begin(), prob_js.end(), prob_js.begin(),
			[neighb_sum](double& c) -> double { return c / neighb_sum; });
		rd->vec_biased_sample_one_idx(prob_js, prob_js, node2visit);
	}
	return node2visit;
}

bool ACS::_is_ant_path_feasible(const size_t ant_idx,
	const bool is_path_complete) const
{
	size_t last_node = _ants[ant_idx].path.back();
	double cost = _ants[ant_idx].cost 
		+ (is_path_complete ? 0 : _dists[last_node][_end_idx]);
	if (cost > _bgt + EPS) { return false; } // Allow equal
	else { return true; }
}

void ACS::_set_node_for_ant(const size_t node2visit, const size_t ant_idx)
{
	_ants[ant_idx].cost += _dists[_ants[ant_idx].cur_node][node2visit];
	_ants[ant_idx].prize += _sz_p_pt[node2visit].second.first;
	_ants[ant_idx].cur_node = node2visit;
	_ants[ant_idx].path.push_back(node2visit);
	_ants[ant_idx].mask_visit(node2visit, _sz_p_pt);
	_ants[ant_idx].feasible = _is_ant_path_feasible(ant_idx, false);
}

void ACS::_construct_path_for_ant(const size_t ant_idx)
{
	while (_ants[ant_idx].feasible && (!_ants[ant_idx].is_all_visit())) {
		size_t node2visit = _state_transition_rule(_ants[ant_idx]);
		_set_node_for_ant(node2visit, ant_idx);
	}
	// Add the last distance back to the end node
	_ants[ant_idx].cost += _dists[_ants[ant_idx].cur_node][_end_idx];
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
		return _dists[path[v1]][path[v2]] - _dists[path[v1]][path[v1_next]]
			+ _dists[path[v1_next]][path[v2_next]] - _dists[path[v2]][path[v2_next]];
	};

	bool improved = true;
	size_t cnt = 0;
	double cost_diff = 0;
	size_t path_len = _ants[ant_idx].path.size();
	while (improved && cnt < _num_iters_2opt) {
		improved = false;
		for (size_t i = 0; i < path_len - 2; ++i) {
			for (size_t j = i + 1; j < path_len - 1; ++j) {
				cost_diff = delta_cost(_ants[ant_idx].path, i, j);
				if (cost_diff < 0) {
					Common::vec_reverse(_ants[ant_idx].path, i, j);
					_ants[ant_idx].cost += cost_diff;
					improved = true;
					cnt++;
				}
			}
		}
	}
}

void ACS::_iterative_drop(const size_t ant_idx, size_t& idx2drop, size_t& drop_node,
	multimap<double, IdxValPair>& dropval_node_cost)
{
	// This is for the first dropval map update
	_ants[ant_idx].prize -= _sz_p_pt[drop_node].second.first;
	_ants[ant_idx].cost -= dropval_node_cost.begin()->second.second;
	// Before erase the index, record its previous and next neighbor node.
	size_t idx0 = idx2drop - 1;
	_ants[ant_idx].path.erase(_ants[ant_idx].path.begin() + idx2drop);
	// ALTHOUGH the dropped node performs POOR in current path, we still
	// need to add it back to unvisited node list because it may be a good 
	// candidate after adding some other nodes. 
	_ants[ant_idx].unmask_visit(drop_node, _sz_p_pt);
	_ants[ant_idx].feasible = _is_ant_path_feasible(ant_idx, true);
	if (_ants[ant_idx].feasible) { return; }
	// This node has been dropped, remove it
	dropval_node_cost.erase(dropval_node_cost.begin());
	while (!_ants[ant_idx].feasible) {  // Keep dropping nodes until feasible
		// Traverse dropval_node_cost and modify previous and next node.
		for (size_t v = 0; v < dropval_node_cost.size(); ++v) {
			auto it = dropval_node_cost.begin(); std::advance(it, v);
			double cost_diff = 0, drop_val = 0;
			size_t cur_node = _ants[ant_idx].path[idx0];
			if (it->second.first == cur_node) {  // Previous node
				size_t prev_node = _ants[ant_idx].path[idx0 - 1];
				size_t next_node = _ants[ant_idx].path[idx0 + 1];
				cost_diff = _dists[prev_node][cur_node] + _dists[cur_node][next_node]
					- _dists[prev_node][next_node];
				drop_val = _sz_p_pt[cur_node].second.first / cost_diff;
				dropval_node_cost.insert(
					make_pair(drop_val, make_pair(cur_node, cost_diff)));
				dropval_node_cost.erase(it);
				continue;  // The iterator has been deleted
			}
			cur_node = _ants[ant_idx].path[idx2drop];
			if (it->second.first == cur_node) {  // Next node
				size_t prev_node = _ants[ant_idx].path[idx0];
				size_t next_node = _ants[ant_idx].path[idx2drop + 1];
				cost_diff = _dists[prev_node][cur_node] + _dists[cur_node][next_node]
					- _dists[prev_node][next_node];
				drop_val = _sz_p_pt[cur_node].second.first / cost_diff;
				dropval_node_cost.insert(
					make_pair(drop_val, make_pair(cur_node, cost_diff)));
				dropval_node_cost.erase(it);
			}
		}
		// Update where2drop and prev_idx for next usage.
		drop_node = dropval_node_cost.begin()->second.first;
		_ants[ant_idx].prize -= _sz_p_pt[drop_node].second.first;
		_ants[ant_idx].cost -= dropval_node_cost.begin()->second.second;
		idx2drop = Common::find_index_by_member(_ants[ant_idx].path, drop_node);
		_ants[ant_idx].path.erase(_ants[ant_idx].path.begin() + idx2drop);
		idx0 = idx2drop - 1; // Update previous index
		_ants[ant_idx].unmask_visit(drop_node, _sz_p_pt);
		_ants[ant_idx].feasible = _is_ant_path_feasible(ant_idx, true);
		// This node has been dropped, remove it
		dropval_node_cost.erase(dropval_node_cost.begin());
	}
}

void ACS::_drop_nodes_until_feasible(const size_t ant_idx)
{
	if (!_ants[ant_idx].feasible) {
		// drop value, {node, cost difference}
		multimap<double, IdxValPair> dropval_node_cost;
		// Calculate the initial dropval_node_cost map
		for (size_t j = 1; j < _ants[ant_idx].path.size() - 1; j++) {
			size_t cur_node = _ants[ant_idx].path[j];
			size_t prev_node = _ants[ant_idx].path[j - 1];
			size_t next_node = _ants[ant_idx].path[j + 1];
			double cost_diff = _dists[prev_node][cur_node] 
				+ _dists[cur_node][next_node] - _dists[prev_node][next_node];
			cost_diff = std::max(EPS, cost_diff);  // Avoid i is on line (i-1)(i+1)
			dropval_node_cost.insert(
				make_pair(_sz_p_pt[cur_node].second.first / cost_diff,
				make_pair(cur_node, cost_diff)));
		}
		size_t node2drop = dropval_node_cost.begin()->second.first;
		size_t idx2drop = Common::find_index_by_member(_ants[ant_idx].path, node2drop);
		// Here we HAVEN'T DONE ANY MODIFICATION FOR ant_k
		_iterative_drop(ant_idx, idx2drop, node2drop, dropval_node_cost);
	}
}

double ACS::_addcost(const size_t adj0, const size_t adj1, const size_t node2add) const
{
	return _dists[adj0][node2add] + _dists[node2add][adj1] - _dists[adj0][adj1];
}

multimap<double, DualidxVal, greater<double>> ACS::_calc_init_addval_map(
	const size_t ant_idx) const
{
	multimap<double, DualidxVal, greater<double>> addval_map;
	for (size_t v = _num_depots; v < _num_nodes; ++v) {
		if (!_ants[ant_idx].is_node_visit(v)) {
			// distance, {node, node index in path}
			multimap<double, DualIdx> cost_nbrs;
			// Find neighbor nodes of node v.
			for (size_t j = 0; j < _ants[ant_idx].path.size(); j++) {
				cost_nbrs.insert(make_pair(_dists[v][_ants[ant_idx].path[j]],
					make_pair(_ants[ant_idx].path[j], j)));
			}
			// {node, node index in path}
			vector<DualIdx> nbr_pairs; nbr_pairs.reserve(3); 
			for (auto it = cost_nbrs.begin(); it != cost_nbrs.end(); ++it) {
				nbr_pairs.push_back(make_pair(it->second.first, it->second.second));
				if (nbr_pairs.size() > 2) { cost_nbrs.clear();  break; }
			}
			// {{nbr, adj}, {nbr index in path, adj index in path}}, add cost
			typedef std::pair<DualIdx, DualIdx> DualPair;
			std::multimap<double, DualPair> cost_nbr_pairs;
			// Now fill neighbor node pairs, their indices in path and addcost
			// Note that we DO NOT check whether picked neighbors are adjacent
			// because addval_map will be ITERATIVELY computed ANYWAY, if we 
			// directly return adjacent neighbor pair, it will raise a problem later.
			for (size_t i = 0; i < nbr_pairs.size(); i++) {
				size_t nbr = nbr_pairs[i].first; size_t nbr_idx = nbr_pairs[i].second;
				// Here we DO NOT add pair of start and end node, i.e. ABANDON
				// start_node -> end_node AND end_node -> start_node
				// because the first and last node MUST be start and end node
				if (nbr_idx < _ants[ant_idx].path.size() - 1) {
					size_t next_idx = nbr_idx + 1;
					cost_nbr_pairs.insert(make_pair(
						_addcost(nbr, _ants[ant_idx].path[next_idx], v), make_pair(
							make_pair(nbr, _ants[ant_idx].path[next_idx]), 
							make_pair(nbr_idx, next_idx))));
				}
				if (nbr_idx > 0) {
					size_t prev_idx = nbr_idx - 1;
					cost_nbr_pairs.insert(make_pair(
						_addcost(_ants[ant_idx].path[prev_idx], nbr, v), make_pair(
							make_pair(_ants[ant_idx].path[prev_idx], nbr), 
							make_pair(prev_idx, nbr_idx))));
				}
			}
			nbr_pairs.clear();
			// Get the smallest addcost of this node v.
			double addcost_v = 1e6; size_t idx2add_v = 0;
			// If any two nodes adjacent in path, then find the pair w/
			// the SMALLEST addcost and find index to insert node v.
			addcost_v = cost_nbr_pairs.begin()->first;
			idx2add_v = std::max(  // left, node to add, right
				cost_nbr_pairs.begin()->second.second.first,
				cost_nbr_pairs.begin()->second.second.second);
			// If the addcost within budget constraint, add it into the addval_map
			if (addcost_v + _ants[ant_idx].cost <= _bgt) {
				addval_map.insert(make_pair(_sz_p_pt[v].second.first / addcost_v,
					make_pair(make_pair(v, idx2add_v), addcost_v)));
			}
		}
	}
	return addval_map;
}

void ACS::_iterative_add(const size_t ant_idx, size_t& idx2add, size_t& add_node,
	multimap<double, DualidxVal, greater<double>>& addval_map)
{
	vector<std::pair<double, DualidxVal>> map2update; 
	map2update.reserve(addval_map.size());
	while (add_node != _start_idx) { // It means there exists feasible nodes to add
		// ONLY update prize and cost FOR NOW
		_ants[ant_idx].prize += _sz_p_pt[add_node].second.first;
		_ants[ant_idx].cost += addval_map.begin()->second.second;
		// Remove any node that belongs to the same SZ.
		Common::erase_if(addval_map, 
			[add_node, this](const auto& pair) -> bool
			{
				size_t this_node = pair.second.first.first;
				if (_sz_p_pt[add_node].first == _sz_p_pt[this_node].first) {
					return true;
				}
				else { return false; }
			});
		Common::erase_if(addval_map,
			[add_node, idx2add, this, ant_idx, &map2update](const auto& pair) -> bool
			{
				size_t unvis_i = pair.second.first.first;
				// Update this unvis_i's addvalue
				multimap<double, DualidxVal, greater<double>> temp_map;
				size_t unvis_idx = pair.second.first.second;
				double unvis_p = _sz_p_pt[unvis_i].second.first;
				// Only when the previous add node is NOT LOCATED ADJACENT to the 
				// current add node, then insert previous calculated addval_map.
				if (unvis_idx != idx2add && unvis_idx != idx2add + 1) {
					if (idx2add < unvis_idx) {
						temp_map.insert(make_pair(pair.first, make_pair(
							make_pair(unvis_i, unvis_idx + 1), pair.second.second)));
					}
					else temp_map.insert(pair);
				}
				// Otherwise, there exists two scenarios:
				// 1. Previous idx2add is located before or after current unvis_idx;
				// 2. Insertion of previous add node will lead to shorter addcost.
				// BOTH scenario will need to compute the addcost of inserting 
				// unvis_i before and after add_node. 
				double addcost_before_node2add = _addcost(
					_ants[ant_idx].path[idx2add - 1], add_node, unvis_i);
				temp_map.insert(make_pair(unvis_p / addcost_before_node2add,
					make_pair(make_pair(unvis_i, idx2add), addcost_before_node2add)));
				double addcost_aftr_node2add = _addcost(
					add_node, _ants[ant_idx].path[idx2add], unvis_i);
				temp_map.insert(make_pair(unvis_p / addcost_aftr_node2add,
					make_pair(make_pair(unvis_i, idx2add + 1), addcost_aftr_node2add)));
				// If meets the constraint, add it into addval_node_pair_cost and 
				// erase the iterator (because the value MAY not be updated)
				if (temp_map.begin()->second.second + _ants[ant_idx].cost <= _bgt) {
					map2update.push_back(*temp_map.begin());
				}
				// If the smallest addcost cannot meet the budget 
				// OR it has been updated, erase the iterator.
				return true;
			});
		// Insert back updated unvisited node
		for (const auto& pair : map2update) { addval_map.insert(pair); }
		map2update.clear();
		map2update.reserve(addval_map.size());
		// This is updating the PREVIOUS node2add !
		Common::vec_insert(_ants[ant_idx].path, idx2add, add_node);
		_ants[ant_idx].mask_visit(add_node, _sz_p_pt);
		_ants[ant_idx].feasible = _is_ant_path_feasible(ant_idx, true);
		// Then get the CURRENT node2add
		if (!addval_map.empty() && addval_map.begin()->first > 0) {
			add_node = addval_map.begin()->second.first.first;
			idx2add = addval_map.begin()->second.first.second;
		}
		else { add_node = _start_idx; }
	}
}

void ACS::_add_nodes_until_budget(const size_t ant_idx)
{
	if (!_ants[ant_idx].is_all_visit() && _ants[ant_idx].feasible) {
		size_t node2add = 0; size_t idx2add = 0;
		// add value, {{node to add, path index to add}, add cost}
		multimap<double, DualidxVal, greater<double>> addval_map = _calc_init_addval_map(ant_idx);
		if (!addval_map.empty() && addval_map.begin()->first > 0) {
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
	if (_ants[_lb_ant_id].feasible) {
		double prize_diff = _ants[_lb_ant_id].prize - _gb_ant->prize;
		if (prize_diff > 0) {  // General ant higher prize
			_gb_ant->cost = _ants[_lb_ant_id].cost;
			_gb_ant->prize = _ants[_lb_ant_id].prize;
			_gb_ant->path = _ants[_lb_ant_id].path;
			if (!_gb_ant->gb_delta_taus.empty()) { _gb_ant->gb_delta_taus.clear(); }
			size_t path_len = _gb_ant->path.size() - 1;
			_gb_ant->gb_delta_taus.reserve(path_len);
			double delta_tau = _alpha * _gb_ant->prize / _gb_ant->cost;
			for (size_t i = 0; i < path_len; i++) {
				// Record the rhs term of the global updating rule.
				_gb_ant->gb_delta_taus.push_back(make_pair(make_pair(
					_gb_ant->path[i], _gb_ant->path[i + 1]), delta_tau));
			}
			if (prize_diff > _tol_acs) { _no_impr_cnt = 0; return; }
		}
		else if (abs(prize_diff) <= EPS) {
			double cost_diff = _gb_ant->cost - _ants[_lb_ant_id].cost;
			if (cost_diff > 0) {  // General ant less cost
				_gb_ant->cost = _ants[_lb_ant_id].cost;
				_gb_ant->prize = _ants[_lb_ant_id].prize;
				_gb_ant->path = _ants[_lb_ant_id].path;
				if (!_gb_ant->gb_delta_taus.empty()) { _gb_ant->gb_delta_taus.clear(); }
				size_t path_len = _gb_ant->path.size() - 1;
				_gb_ant->gb_delta_taus.reserve(path_len);
				double delta_tau = _alpha * _gb_ant->prize / _gb_ant->cost;
				for (size_t i = 0; i < path_len; i++) {
					// Record the rhs term of the global updating rule.
					_gb_ant->gb_delta_taus.push_back(make_pair(make_pair(
						_gb_ant->path[i], _gb_ant->path[i + 1]), delta_tau));
				}
				if (cost_diff > _tol_acs) { _no_impr_cnt = 0; return; }
			}
		}
	}
	_no_impr_cnt++;
}

void ACS::_global_updating_rule()
{
	double evap_ratio = 1 - _alpha;
	for (size_t i = 0; i < _num_nodes; i++) {
		if (i == _end_idx) { continue; }  // Don't care tau[_end_idx][x]!
		for (size_t j = 1; j < _num_nodes; j++) {  // don't care tau[x][_start_idx]!
			if (i != j) { _taus[i][j] *= evap_ratio; }
		}
	}
	for (const GBAnt::DualidxVal& x : _gb_ant->gb_delta_taus) {
		_taus[x.first.first][x.first.second] += x.second;
	}
}

int ACS::_is_gb_sol_valid() const
{
	double prize = 0, cost = 0;
	for (size_t i = 0; i < _gb_ant->path.size() - 1; i++) {
		size_t node = _gb_ant->path[i], next_node = _gb_ant->path[i + 1];
		prize += _sz_p_pt[node].second.first;
		cost += _dists[node][next_node];
		for (size_t j = i + 1; j < _gb_ant->path.size(); j++) {
			size_t test_node = _gb_ant->path[j];
			if (_sz_p_pt[node].first == _sz_p_pt[test_node].first) {
				return -4; 
			}
		}
	}
	if (cost > _bgt) { 
		return -1; 
	}
	if (abs(cost - _gb_ant->cost) > EPS) {
		return -2; 
	}
	if (abs(prize - _gb_ant->prize) > EPS) {
		return -3; 
	}
	return 0;
}

void ACS::evolve_until_stop_criterion()
{
	unique_ptr<Rd> rd = make_unique<Rd>();
	vector<int> groups = vector<int>(_num_groups);
	std::iota(groups.begin(), groups.end(), static_cast<int>(_num_depots));
	for (size_t t = 0; t < _num_iters_acs; t++) {
		if (_no_impr_cnt > _max_no_impr) { break; }
		double lb_prize = 0, lb_cost = 1e6;
		vector<int> start_groups;
		if (_num_ants == _num_groups) { // If same length, shuffle group indices
			start_groups = groups;
			rd->vec_shuffle_skip_first_n(0, start_groups);
		}
		else { rd->vec_uniform_rand_sample(groups, _num_ants, start_groups); }
		for (size_t k = 0; k < _num_ants; k++) {
			// Set the first random start node.
			int start_node_k = _rand_sample_start_node(start_groups[k]);
			_set_node_for_ant(static_cast<size_t>(start_node_k), k);
			// Construct the solution path.
			_construct_path_for_ant(k);
			// Apply 2-opt to the solution path.
			if (_num_iters_2opt > 0) { _2opt_for_ant(k); }
			_ants[k].feasible = _is_ant_path_feasible(k, true);
			// [DROP NODE PHASE] Keep dropping node until solution is feasible.
			_drop_nodes_until_feasible(k);
			// [ADD NODE PHASE] Keep adding node until no unvisited nodes,
			// reach budget constraint, or reach number of searches.
			_add_nodes_until_budget(k);
			// Determine the locally best ant for later updating globally best ant.
			if (_ants[k].prize > lb_prize || 
				(_ants[k].prize == lb_prize && _ants[k].cost < lb_cost)) {
				lb_prize = _ants[k].prize;
				lb_cost = _ants[k].cost;
				_lb_ant_id = k;
			}
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
	// Uncomment below code to test solution validity.
	//int return_status = _is_gb_sol_valid();
	//if (return_status < 0) { 
	//	throw std::runtime_error(std::to_string(return_status)); 
	//}
}

void ACS::_arc_srch_one_iter(const bool is_1st_iter, const vector<size_t>& path, 
	vector<Point>& waypts, double& cost, bool& is_waypt_changed) const
{
	double cost_diff = 0; Point opt_pt; bool is_opt_found = false;
	for (size_t j = 1; j < waypts.size() - 1; j++) {
		// Get the raw cost first
		if (is_1st_iter) {  // Some indices in solution are still valid.
			if (!is_waypt_changed) {
				// When THERE IS NO CHANGE MADE FOR INDEX j-1, we can 
				// directly use d[p1][p2] and d[p2][p3] from previous cost_mat 
				cost_diff = _dists[path[j - 1]][path[j]] 
					+ _dists[path[j]][path[j + 1]];
			}
			else {
				// If way-point change happened, we can still use d[p2][p3], BUT
				// p1 has changed to an unregistered point in cost_mat.
				cost_diff = (waypts[j - 1] - waypts[j]).norm2()
					+ _dists[path[j]][path[j + 1]];
			}
		}
		else {  // Otherwise cost difference has to be manually computed.
			cost_diff = (waypts[j - 1] - waypts[j]).norm2()
				+ (waypts[j] - waypts[j + 1]).norm2();
		}
		size_t sz_id = _sz_p_pt[path[j]].first - _num_depots;
		is_opt_found = _stn_zones[sz_id].arc_search(
			waypts[j - 1], waypts[j + 1], opt_pt);
		if (!is_opt_found) {
			is_opt_found = _stn_zones[sz_id].find_closest_vertex(
				waypts[j - 1], waypts[j + 1], cost_diff, opt_pt);
			// Note opt_pt may be not changed, i.e. is_opt_pt_found = false
		}
		if (is_opt_found) {
			// The updated opt_pt's cost
			double opt_diff = (opt_pt - waypts[j - 1]).norm2()
				+ (opt_pt - waypts[j + 1]).norm2() - cost_diff;
			if (opt_diff < 0) {
				// If improved, update total cost and point position.
				cost += opt_diff;
				waypts[j] = opt_pt;
				if (is_1st_iter) { is_waypt_changed = true; }
			}
		}
		// Note only 1st iteration has is_changed true or false, 
		// later iterations will always update based on manual cost computation.
		else { if (is_1st_iter) { is_waypt_changed = false; } }
	}
}

double ACS::_addcost(const Point& adj0, const Point& adj1, const Point& pt2add) const
{
	return (adj0 - pt2add).norm2() + (pt2add - adj1).norm2() - (adj0 - adj1).norm2();
}

void ACS::_calc_addval_map(const vector<Point>& waypts,
	const vector<bool>& visit_mask, const double init_cost,
	multimap<double, DualidxVal, greater<double>>& addval_map) const
{
	typedef std::pair<Point, size_t> PtIdx;
	typedef std::pair<Point, Point> DualPt;
	typedef std::pair<DualPt, DualIdx> DualptDualidx;

	for (size_t v = _num_depots; v < _num_nodes; v++) {
		if (!visit_mask[v]) {
			Point pt_v = _sz_p_pt[v].second.second;
			// distance, {Point, node index in path}
			multimap<double, PtIdx> cost_nbrs;
			// Find neighbor nodes of node v.
			for (size_t i = 0; i < waypts.size(); i++) {
				cost_nbrs.insert(
					make_pair((pt_v - waypts[i]).sq_sum(), make_pair(waypts[i], i)));
			}
			vector<PtIdx> nbr_pairs;	// {Point, node index in path}
			for (auto it = cost_nbrs.begin(); it != cost_nbrs.end(); ++it) {
				nbr_pairs.push_back(make_pair(it->second.first, it->second.second));
				if (nbr_pairs.size() > 2) { cost_nbrs.clear(); break; }
			}
			// add cost, {{nbr, adj}, {nbr index in path, adj index in path}}
			multimap<double, DualptDualidx> cost_dualnbrs;
			// Now fill neighbor Point pairs, their indices in path and addcost
			// Note that we DO NOT check whether picked neighbors are adjacent
			// because addval_map will be ITERATIVELY computed ANYWAY, if we 
			// directly return adjacent neighbor pair, it will raise a problem later.
			for (size_t i = 0; i < nbr_pairs.size(); i++) {
				Point nbr = nbr_pairs[i].first;
				size_t nbr_idx = nbr_pairs[i].second;
				// Here we DO NOT add pair of start and end node, i.e. ABANDON
				// start_node -> end_node AND end_node -> start_node
				// because the first and last node MUST be start and end node.
				if (nbr_idx < waypts.size() - 1) {
					size_t next_idx = nbr_idx + 1;
					cost_dualnbrs.insert(make_pair(
						_addcost(nbr, waypts[next_idx], pt_v), make_pair(
							make_pair(nbr, waypts[next_idx]), 
							make_pair(nbr_idx, next_idx))));
				}
				if (nbr_idx > 0) {
					size_t prev_idx = nbr_idx - 1;
					cost_dualnbrs.insert(make_pair(
						_addcost(waypts[prev_idx], nbr, pt_v), make_pair(
							make_pair(waypts[prev_idx], nbr), 
							make_pair(prev_idx, nbr_idx))));
				}
			}
			// Get the smallest addcost of this node v.
			double addcost_v = cost_dualnbrs.begin()->first;
			// If the addcost within budget constraint, add it into the addval_map
			if (addcost_v + init_cost <= _bgt) {
				size_t idx2add_v = std::max(cost_dualnbrs.begin()->second.second.first,
				cost_dualnbrs.begin()->second.second.second);
				// Here we can still use _sz_p_pt because v is in _sz_p_pt
				addval_map.insert(make_pair(_sz_p_pt[v].second.first / addcost_v,
					make_pair(make_pair(v, idx2add_v), addcost_v)));
			}
		}
	}
}

void ACS::_mask_visit(const size_t node_idx, vector<bool>& visit_mask) const
{
	visit_mask[node_idx] = true;
	size_t sz_id = _sz_p_pt[node_idx].first;
	for (size_t i = node_idx - 1; i > 1; i--) {
		if (sz_id == _sz_p_pt[i].first) { visit_mask[i] = true; }
		else { break; }
	}
	for (size_t i = node_idx + 1; i < _sz_p_pt.size(); i++) {
		if (sz_id == _sz_p_pt[i].first) { visit_mask[i] = true; }
		else { break; }
	}
}

void ACS::_add_nodes_until_budget(double& path_cost, double& path_prize,
	vector<size_t>& path, vector<bool>& visit_mask, vector<Point>& waypts) const
{
	if (!Common::is_bool_vec_all_true(visit_mask) && path_cost < _bgt - EPS) {
		size_t idx2add = 0, node2add = 0;
		// add value, {{node to add, path index to add}, add cost}
		multimap<double, DualidxVal, greater<double>> addval_map;
		do {
			_calc_addval_map(waypts, visit_mask, path_cost, addval_map);
			if (!addval_map.empty()) {
				if (addval_map.begin()->first > 0) {
					// Get the node v and index to add in current path (whose 
					// addval is positive and the largest).
					node2add = addval_map.begin()->second.first.first;
					idx2add = addval_map.begin()->second.first.second;
					path_prize += _sz_p_pt[node2add].second.first;
					path_cost += addval_map.begin()->second.second;
					Common::vec_insert(waypts, idx2add, 
						_sz_p_pt[node2add].second.second);
					Common::vec_insert(path, idx2add, node2add);
					_mask_visit(node2add, visit_mask);
				}
				else { node2add = _start_idx; }
			}
			else { node2add = _start_idx; }
			addval_map.clear();
		} while (node2add != _start_idx && path_cost < _bgt - EPS);
	}
}

void ACS::get_gb_info(double& best_prize,
	double& best_cost, vector<Point>& best_pts) const
{
	best_cost = _gb_ant->cost;
	best_prize = _gb_ant->prize;
	best_pts.reserve(_gb_ant->path.size());
	for (size_t i = 0; i < _gb_ant->path.size(); i++) {
		best_pts.push_back(_sz_p_pt[_gb_ant->path[i]].second.second);
	}
}

void ACS::arc_search(double& best_prize, double& best_cost,
	vector<Point>& best_pts) const
{
	vector<size_t> path_seq = _gb_ant->path;
	vector<bool> gb_visit_mask(_num_nodes, false);
	gb_visit_mask[0] = gb_visit_mask[1] = true;
	best_cost = _gb_ant->cost, best_prize = _gb_ant->prize;
	if (!best_pts.empty()) { best_pts.clear(); } 
	best_pts.reserve(_num_groups + _num_depots);
	for (size_t i = 0; i < path_seq.size(); i++) {
		size_t node_id = path_seq[i];
		size_t sz_id = _sz_p_pt[node_id].first;
		best_pts.push_back(_sz_p_pt[node_id].second.second);
		if (node_id != _start_idx && node_id != _end_idx) {
			_mask_visit(node_id, gb_visit_mask);
		}
	}
	// DataIO::write_ceop_sol(fs::current_path(), "team5_499.ceop", 
	// 	".impr0", _gb_ant->cost, _gb_ant->prize, best_pts);

	bool is_waypt_changed = false;
	for (size_t t = 0; t < _num_iters_arc; t++) {
		_arc_srch_one_iter(t == 0, path_seq, best_pts, best_cost, is_waypt_changed);
		_add_nodes_until_budget(best_cost, best_prize, path_seq, gb_visit_mask, best_pts);
		// DataIO::write_ceop_sol(fs::current_path(), "team5_499.ceop", 
		// 	".impr" + to_string(t + 1), best_cost, best_prize, best_pts);
	}
}
