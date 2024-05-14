/*****************************************************************//**
 * \file   DataIO.h
 * \brief  Data input and output management.
 * 
 * \author Qiuchen Qian
 * \date   October 2023
 *********************************************************************/
#ifndef DATAIO_H
#define DATAIO_H

#include <string>
#include <iostream>
#include <filesystem>
#include <fstream>
#include "StnZone.h"
#include "ACS.h"

using std::string;
namespace fs = std::filesystem;

namespace DataIO 
{
	constexpr auto DOTDELIMITER = ".";

	void load_tddp_data(const fs::path data_dir, const string f_name,
		vector<Circle>& circs);

	void write_stn_zone(const fs::path data_fdir, const string f_name,
		const vector<StnZone>& stn_zones);

	void write_stn_zone_vtx(const fs::path data_fdir, const string f_name,
		const vector<StnZone>& stn_zones);

	void write_points(const fs::path data_fdir, const size_t cnt, 
		const string suffix, const double cost, 
		const double prize, const vector<Point>& pts);

	void write_nbhd_c(const fs::path data_fdir, const size_t cnt,
		const string suffix, const vector<NbhdWayPt>& nbhd_waypts);

	void write_final_sol(const fs::path data_fdir, const string f_name,
		const time_t time_stamp, const string suffix, const double alg_time,
		const double cost, const double prize, const vector<Point>& pos);

	void write_ant_fit(const fs::path data_fdir, const string f_name,
		const size_t ptcl_idx, const vector<double>& gb_fits);
}


#endif // !DATAIO_H

