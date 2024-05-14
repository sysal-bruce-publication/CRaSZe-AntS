/*****************************************************************//**
 * \file   main.cpp
 * \brief  
 * 
 * \author Qiuchen Qian
 * \date   October 2023
 *********************************************************************/
#include "DataIO.h"
#include "json.hpp"
#include "Timer.h"
#include "StnZoneSrch.h"
#include "ACS.h"

#define ROOT_PATH fs::current_path()
#define PROBLEM_NAME "ceop"

using json = nlohmann::json;

string data_fname;
fs::path data_fdir;
unique_ptr<Timer> timer;
unique_ptr<StnZoneSrch> szs;
unique_ptr<ACS> acs;
double alg_time = 0;
bool print_result = false, set_timer = false;

void initialize(const string config_fpath, const bool print_config)
{
	// Load JSON configuration file.
	std::ifstream json_f(config_fpath);
	json json_data = json::parse(json_f);
	// Print configuration before close the file
	if (print_config) { cout << json_data.dump(4) << "\n"; } cout.flush();
	print_result = json_data["print_result"];
	set_timer = json_data["set_timer"];
	// file name and directory
	string data_dir_name = json_data["data_dir"];
	string problem_name  = json_data["problem"]["name"];
	double budget        = json_data["problem"]["budget"];
	data_fname          += "." + problem_name;
	data_fdir            = ROOT_PATH / data_dir_name / problem_name;
	// SZS related parameters
	size_t num_iters_szs = json_data["SZS"]["num_iters"];
	size_t max_circs     = json_data["SZS"]["max_circs"];
	size_t num_vtxs      = json_data["SZS"]["num_disc_pts"];
	// ACS related parameters
	size_t num_ants        = json_data["ACS"]["num_ants"];
	size_t num_max_no_impr = json_data["ACS"]["max_no_impr"];
	size_t num_iters_acs   = json_data["ACS"]["num_iters"];
	size_t num_iters_2opt  = json_data["2opt"]["num_iters"];
	size_t num_iters_arc   = json_data["arc"]["num_iters"];
	double q0              = json_data["ACS"]["q0"];
	double rho             = json_data["ACS"]["rho"];
	double alpha           = json_data["ACS"]["alpha"];
	double beta            = json_data["ACS"]["beta"];
	double tol_acs         = json_data["ACS"]["tol"];
	// Close file
	json_f.close();

	// Load data
	vector<Circle> all_circs;
	DataIO::load_ceop_data(data_fdir, data_fname, all_circs);
	// Initialize SZS
	szs = make_unique<StnZoneSrch>(num_iters_szs, max_circs, num_vtxs, all_circs);
	// Initialize ACS
	acs = make_unique<ACS>(num_ants, num_max_no_impr, num_iters_acs, num_iters_2opt, 
                        num_iters_arc, budget, q0, alpha, beta, rho, tol_acs);
	// Initialize timer
	if (set_timer) { timer = make_unique<Timer>(); }
}

int main(int argc, char* argv[])
{
	// string config_fname = "general.ceop.json";
	// data_fname = "bubbles1";
	// bool print_config = true;
	// bool write_vertex = true;
	string config_fname = argv[1];
	data_fname = argv[2];
	bool print_config = std::stoi(argv[3]);
	bool write_vertex = std::stoi(argv[4]);

	// ROOT_PATH should be located in `workspace/`
	fs::path config_fpath = ROOT_PATH / "configs" / config_fname;
	initialize(config_fpath.string(), print_config);
	timer->tick();
	vector<StnZone> stn_zones = szs->randomized_search();
	acs->init_colony(szs->start_pt(), szs->end_pt(), stn_zones);
	acs->evolve_until_stop_criterion();

	double final_prize = 0, final_cost = 0;
	// Arc search algorithm phase
	vector<Point> final_pos;
 	acs->get_gb_info(final_prize, final_cost, final_pos);
	acs->arc_search(final_prize, final_cost, final_pos);
	timer->tock();
	alg_time += timer->duration().count();

	time_t end_time = std::chrono::system_clock::to_time_t(
		std::chrono::system_clock::now());
	if (print_result) {
		cout << end_time << " ";
		if (set_timer) { cout << "Time " << alg_time << " ms, "; }
		cout << "Cost " << final_cost << ", Prize " << final_prize << std::endl;
	}
	if (write_vertex) {
		DataIO::write_sz_p_pt(data_fdir, data_fname, acs->sz_p_pt());
	}
	DataIO::write_ceop_sol(data_fdir, data_fname, end_time, ".ceop_sol",
		alg_time, final_cost, final_prize, final_pos);

	return 0;
}
