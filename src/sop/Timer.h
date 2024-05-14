/*****************************************************************//**
 * \file   Timer.h
 * \brief  Timer for counting algorithm execution time
 *
 * \author Qiuchen Qian
 * \date   August 2023
 *********************************************************************/
#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <assert.h>
using namespace std::chrono_literals;

class Timer
{
public:
	Timer() = default;
	~Timer() {}

	void tick()
	{
		_end = timep_t{};
		_start = std::chrono::steady_clock::now();
	}

	void tock() { _end = std::chrono::steady_clock::now(); }

	auto duration() const {
		// Use gsl_Expects if your project supports it.
		assert(_end != timep_t{} && "Timer must toc before reading the time");
		return std::chrono::duration_cast<std::chrono::milliseconds>(_end - _start);
	}

private:
	using timep_t = decltype(std::chrono::steady_clock::now());

	timep_t _start = std::chrono::steady_clock::now();
	timep_t _end = {};
};

#endif // !TIMER_H
