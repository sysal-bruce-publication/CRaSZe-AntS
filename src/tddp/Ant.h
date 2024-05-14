/*****************************************************************//**
 * \file   Ant.h
 * \brief  Class for single population in ACS
 * 
 * \author Qiuchen Qian
 * \date   October 2023
 *********************************************************************/
#ifndef ANT_H
#define ANT_H

#include "Point.h"

struct Fitness
{
	Fitness() = default;
	Fitness(const double cost, const double prize) : cost(cost), prize(prize) {}

	double cost = 0;
	double prize = 0;

	Fitness& operator=(const Fitness& other)
	{
		cost = other.cost; prize = other.prize; return *this;
	}

	bool operator>(const Fitness& other) const
	{
		if (prize > other.prize) { return true; }
		else if (abs(prize - other.prize) <= EPS && cost < other.cost) { 
			return true; 
		}
		else { return false; }
	}
};

struct GBAnt
{
	//typedef std::pair<DualIdx, double> DualidxVal;

	GBAnt() = default;

	double delta_tau = 0;
	Fitness fit      = Fitness(1e8, -1);
	vector<size_t> path;
	//vector<DualidxVal> gb_delta_taus;

	void reset()
	{
		fit = Fitness(1e8, -1);
		if (!path.empty()) { path.clear(); }
		//if (!gb_delta_taus.empty()) { gb_delta_taus.clear(); }
	}
};

class Ant
{
public:
	Ant() = default;
	Ant(const size_t num_nodes) : _num_nodes(num_nodes)
	{
		path.reserve(_num_nodes);
		path.push_back(0);  // Push start node index.
		_visit_mask = vector<bool>(_num_nodes, false);
		_visit_mask[0] = _visit_mask[1] = true;
	}
	~Ant() {}

	bool feasible   = false;		// Whether the path is feasible.
	size_t cur_node = 0;			// The current node index
	Fitness fit     = Fitness();	// The fitness of the path
	vector<size_t> path;			// The path of the ant
	
	/**
	 * Check whether node with index idx has been visited.
	 * 
	 * \param idx
	 * \return 
	 */
	bool is_node_visited(const size_t idx) const { return _visit_mask[idx]; }

	/**
	 * Check whether all nodes have been visited.
	 *
	 * \return
	 */
	bool is_all_visit() const { return Common::is_bool_vec_all_true(_visit_mask); }

	/**
	 * If node visited, mask its index to true.
	 * 
	 * \param node_idx
	 */
	void mask_visit(const size_t node_idx) { _visit_mask[node_idx] = true; }

	/**
	 * If node dropped, unmask its index to false.
	 * 
	 * \param node_idx
	 */
	void unmask_visit(const size_t node_idx) { _visit_mask[node_idx] = false; }

	void reset()
	{
		cur_node = 0;  // start node
		fit = Fitness();
		if (!path.empty()) { path.clear(); } path.reserve(_num_nodes);
		path.push_back(0);
		// Reset visit mask
		_visit_mask = vector<bool>(_num_nodes, false);
		_visit_mask[0] = _visit_mask[1] = true;
	}

private:
	const size_t _num_nodes = 0;	// Number of all nodes including start and end.
	vector<bool> _visit_mask;
};

#endif // !ANT_H
