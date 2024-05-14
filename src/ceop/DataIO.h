#ifndef DATAIO_H
#define DATAIO_H

#include <string>
#include <iostream>
#include <filesystem>
#include <fstream>
#include "StnZone.h"

using std::string, std::to_string;
using std::cout;
using std::runtime_error;
namespace fs = std::filesystem;

namespace DataIO
{
	constexpr auto DOTDELIMITER = ".";
	typedef std::pair<std::size_t, std::pair<double, Point>> IdxValPt;

	void load_ceop_data(const fs::path data_dir, const string f_name,
		vector<Circle>& circs);

	void write_stn_zone_vtx(const fs::path data_fdir, const string f_name,
		const size_t t, const vector<StnZone>& stn_zones);

	void write_ceop_sol(const fs::path data_fdir, const string f_name,
		const time_t time_stamp, const string suffix, const double alg_time,
		const double cost, const double prize, const vector<Point>& pos);

	void write_ceop_sol(const fs::path data_fdir, const string f_name,
		const string suffix, const double cost, const double prize, 
		const vector<Point>& waypts);

	void write_sz_p_pt(const fs::path data_fdir, const string f_name,
		const vector<IdxValPt>& sz_p_pt);
}


#endif // !DATAIO_H

