/*****************************************************************//**
 * \file   Utils.h
 * \brief  Utils functions
 *
 * \author Qiuchen Qian
 * \date   August 2023
 *********************************************************************/
#ifndef UTILS_H
#define UTILS_H

#include <set>
#include <cmath>
#include <cstddef>		// std::size_t
#include <algorithm>	// std::shuffle
#include <numeric>		// std::accumulate
#include <vector>
#include <random>		// std::default_random_engine

constexpr double EPS = 1e-7;
constexpr double PI = 3.14159265;

using std::vector;

namespace Common
{
	/**
	 * Erase if condition meets.
	 * 
	 * \param items
	 * \param predicate
	 */
	template<typename ContainerT, typename PredicateT>
	void erase_if(ContainerT& items, const PredicateT& predicate) {
		for (auto it = items.begin(); it != items.end(); ) {
			if (predicate(*it)) it = items.erase(it);
			else ++it;
		}
	}

	/**
	 * Check whether all elements in the vector are unique.
	 *
	 * \param vec
	 * \return If unique, return true; otherwise, return false.
	 */
	template <typename T, typename A>
	inline bool is_elems_in_vec_unique(const vector<T, A>& vec)
	{
		std::set<T> s(vec.begin(), vec.end());
		return s.size() == vec.size();
	}

	/**
	 * Check whether the bool vector is all true.
	 * 
	 * \param vec
	 * \return 
	 */
	inline bool is_bool_vec_all_true(const vector<bool>& vec)
	{
		return std::all_of(vec.begin(), vec.end(), [](bool b) { return b; });
	}

	/**
	 * Find the index of specific member in the vector.
	 *
	 * \param vec
	 * \param member
	 * \return
	 */
	template <typename T, typename A>
	inline int find_index_by_member(const vector<T, A>& vec, const T member)
	{
		auto it = std::find(vec.begin(), vec.end(), member);
		if (it != vec.cend()) {
			return static_cast<int>(std::distance(vec.begin(), it));
		}
		else { return -1; }
	}

	template <typename T, typename A>
	inline T vec_min(const vector<T, A>& vec)
	{
		T min_elem = *std::min_element(vec.begin(), vec.end());
		return min_elem;
	}

	template <typename T, typename A>
	inline T vec_max(const vector<T, A>& vec)
	{
		T max_elem = *std::max_element(vec.begin(), vec.end());
		return max_elem;
	}

	template <typename T, typename A>
	inline T vec_sum(const vector<T, A>& vec)
	{
		T sum_res = 0;
		for (const T& elem : vec) { sum_res += elem; }
		return sum_res;
	}

	/**
	 * Erase one element from the vector.
	 *
	 * \param vec
	 * \param member
	 */
	template <typename T, typename A>
	inline void vec_erase(vector<T, A>& vec, const T member)
	{
		vec.erase(std::remove(vec.begin(), vec.end(), member), vec.end());
	}

	template <typename T, typename A>
	inline void vec_insert(vector<T, A>& vec, size_t idx_in_vec, const T item2insert)
	{
		vec.insert(vec.begin() + idx_in_vec, item2insert);
	}

	/**
	 * Take vec[v1+1] to vec[v2] and add them in reverse order.
	 *
	 * \param vec
	 * \param v1
	 * \param v2
	 */
	template <typename T, typename A>
	inline void vec_reverse(std::vector<T, A>& vec,
		const std::size_t v1, const std::size_t v2)
	{
		std::reverse(vec.begin() + v1 + 1, vec.begin() + v2 + 1);
	}

	/**
	 * Re-range atan2 to range [0, 2pi].
	 *
	 * \return
	 */
	inline double rerange_atan2(const double atan2_res)
	{
		return atan2_res >= 0 ? atan2_res : (2 * PI + atan2_res);
	}

	/**
	 * Restrict a value in range [val_min, val_max].
	 *
	 * \param val_min
	 * \param val_max
	 * \param value
	 */
	template<typename T>
	inline void restrict_value_in_range(const T val_min, const T val_max, T& value)
	{
		value = value < val_min ? val_min : (value > val_max ? val_max : value);
	}
	/**
	 * Check if two variables have same sign. If one of them is 0, we
	 * assume they have same sign.
	 *
	 * \param a
	 * \param b
	 * \return
	 */
	template<typename T = double>
	inline bool is_same_sign(const T a, const T b)
	{
		return (a * b < 0) ? false : true;
	}
}

using std::random_device, std::mt19937;
using std::uniform_real_distribution;
using std::uniform_int_distribution;
using std::discrete_distribution;

class Rd
{
public:
	Rd()
	{
		_gen = mt19937((random_device())());
	}
	~Rd() {}

	/**
	 * Generate a double variable within range [0, 1).
	 *
	 * \return
	 */
	inline double uniform_rand()
	{
		uniform_real_distribution<double> distb(0., 1.);
		return distb(_gen);
	}

	/**
	 * Generate a double variable within range [val_min, val_max).
	 *
	 * \param val_min
	 * \param val_max
	 * \return
	 */
	inline double uniform_rand(const double val_min, const double val_max)
	{
		uniform_real_distribution<double> distb(val_min, val_max);
		return distb(_gen);
	}

	/**
	 * Generate a int variable within range [val_min, val_max].
	 *
	 * \param val_min
	 * \param val_max
	 * \return
	 */
	inline int uniform_rand(const int val_min, const int val_max)
	{
		uniform_int_distribution<int> distb(val_min, val_max);
		return distb(_gen);
	}

	template <typename T, typename A>
	inline void vec_uniform_rand_sample(const vector<T, A>& vec_in,
		const size_t num_elems, vector<T, A>& vec_out)
	{
		std::sample(vec_in.begin(), vec_in.end(),
			std::back_inserter(vec_out), num_elems, _gen);
	}

	/**
	 * Draw num_elems samples from a vector (without replacement)
	 * based on uniform distribution.
	 *
	 * \param vec_in
	 * \param num_elems
	 * \param vec_out
	 */
	template <typename T, typename A>
	inline void vec_uniform_rand_sample_exclude_1st_last(const vector<T, A>& vec_in,
		const size_t num_elems, vector<T, A>& vec_out)
	{
		std::sample(vec_in.begin() + 1, vec_in.end() - 1,
			std::back_inserter(vec_out), num_elems, _gen);
	}

	template <typename T, typename A>
	inline void vec_biased_sample_one_idx(const vector<T, A>& vec_in,
		const vector<double>& weights, size_t& idx)
	{
		discrete_distribution<size_t> distb{ weights.begin(), weights.end() };
		idx = distb(_gen);
	}

	template <typename T, typename A>
	inline void vec_biased_sample_one_elem(const vector<T, A>& vec_in,
		const vector<double>& weights, T& elem)
	{
		discrete_distribution<size_t> distb{ weights.begin(), weights.end() };
		elem = vec_in[distb(_gen)];
	}

	/**
	 * Shuffle the vector, but skip first n elements.
	 *
	 * \param num_skip
	 * \param vec
	 */
	template<typename T, typename A>
	inline void vec_shuffle_skip_first_n(const size_t num_skip, vector<T, A>& vec)
	{
		std::shuffle(vec.begin() + num_skip, vec.end(), _gen);
	}

private:
	mt19937 _gen;
};

#endif // !UTILS_H

