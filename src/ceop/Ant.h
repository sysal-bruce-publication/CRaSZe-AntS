#ifndef ANT_H
#define ANT_H

#include "Point.h"

typedef std::pair<std::size_t, std::pair<double, Point>> IdxValPt;

struct GBAnt
{
	typedef std::pair<DualIdx, double> DualidxVal;

	GBAnt() = default;
	double prize = 0;
	double cost = 1e6;
	vector<size_t> path;
	vector<DualidxVal> gb_delta_taus;

	void reset()
	{
		cost = 1e6;
		if (!path.empty()) { path.clear(); }
		if (!gb_delta_taus.empty()) { gb_delta_taus.clear(); }
	}
};

class Ant
{
public:
	Ant() = default;
	/**
	 * .
	 * 
	 * \param num_waypts: All nodes including start node and end node.
	 */
	Ant(const size_t num_nodes, const size_t num_waypts) : 
		_num_nodes(num_nodes), _num_path(num_waypts)
	{
		path.reserve(num_waypts);
		path.push_back(0);
		_visit_mask = vector<bool>(num_nodes, false);
		_visit_mask[0] = _visit_mask[1] = true;
	}
	~Ant() {}

	size_t cur_node = 0;
	bool feasible = false;
	double prize = 0;
	double cost = 0;		// The distance cost of the path
	vector<size_t> path;	// The path of the ant

	/**
	 * Check whether node with index idx has been visited.
	 *
	 * \param idx
	 * \return
	 */
	bool is_node_visit(const size_t idx) const { return _visit_mask[idx]; }

	/**
	 * Check whether all nodes have been visited.
	 *
	 * \return
	 */
	bool is_all_visit() const
	{
		return Common::is_bool_vec_all_true(_visit_mask);
	}

	void mask_visit(const size_t idx) { _visit_mask[idx] = true; }
	void mask_visit(const size_t idx, const vector<IdxValPt>& sz_p_pt);
	void unmask_visit(const size_t idx, const vector<IdxValPt>& sz_p_pt);

	void reset();

private:
	const size_t _num_path  = 0;  // Number of elements in path first -> ... -> end
	const size_t _num_nodes = 0;  // Number of all nodes
	vector<bool> _visit_mask;
};

#endif // !ANT_H

