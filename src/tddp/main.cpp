/*****************************************************************//**
 * \file   main.cpp
 * \brief  Source code for solving TDDP with PSO and ACS.
 * 
 * \author Qiuchen Qian
 * \date   August 2023
 *********************************************************************/
#include "DataIO.h"
#include "json.hpp"
#include "Timer.h"
#include "StnZoneSrch.h"
#include "PSO.h"

using json = nlohmann::json;
using std::cout;
using std::runtime_error;

fs::path root_dir = fs::current_path();
string data_fname;
fs::path data_fdir;
unique_ptr<Timer> timer;
unique_ptr<StnZoneSrch> szs;
unique_ptr<PSO> pso;
double alg_time = 0;
bool print_result = false, set_timer = false;

/**
 * Load configuration file and initialize the algorithm.
 *
 * \param config_fpath
 * \param fname
 * \return
 */
void initialize(const string config_fpath, const bool print_config)
{
	// Load JSON configuration file.
	std::ifstream json_f(config_fpath);
	if (!json_f.is_open()) { throw runtime_error(config_fpath + " cannot open."); }
	json json_data = json::parse(json_f);
	// Print configuration before close the file
	if (print_config) { cout << json_data.dump(4) << "\n"; } cout.flush();
	print_result = json_data["print_result"];
	set_timer = json_data["set_timer"];
	// file name and directory
	string data_dir_name   = json_data["data_dir"];
	string problem_name    = json_data["problem"]["name"];
	data_fname            += "." + problem_name;
	data_fdir              = root_dir / data_dir_name / problem_name;
	// RSZS related parameters
	size_t num_iters_szs   = json_data["SZS"]["num_iters"];
	size_t max_circs       = json_data["SZS"]["max_circs"];
	double v_uav           = json_data["SZS"]["v_uav"];
	// PSO related parameters
	size_t num_ptcls       = json_data["PSO"]["num_ptcls"];
	size_t num_iters_pso   = json_data["PSO"]["num_iters"];
	size_t max_no_impr_pso = json_data["PSO"]["max_no_impr"];
	double c1              = json_data["PSO"]["c1"];
	double c2              = json_data["PSO"]["c2"];
	double w_max           = json_data["PSO"]["w_max"];
	double w_min           = json_data["PSO"]["w_min"];
	double tol_pso         = json_data["PSO"]["tol"];
	// ACS related parameters
	size_t num_ants        = json_data["ACS"]["num_ants"];
	size_t num_iters_acs   = json_data["ACS"]["num_iters"];
	size_t num_iters_2opt  = json_data["2opt"]["num_iters"];
	size_t max_no_impr_acs = json_data["ACS"]["max_no_impr"];
	double q0              = json_data["ACS"]["q0"];
	double rho             = json_data["ACS"]["rho"];
	double alpha           = json_data["ACS"]["alpha"];
	double beta            = json_data["ACS"]["beta"];
	double tol_acs         = json_data["ACS"]["tol"];
	double v_truck         = json_data["ACS"]["v_truck"];
	double budget          = json_data["ACS"]["budget"];
	// Close file
	json_f.close();

	// Load data
	vector<Circle> all_circs;
	DataIO::load_tddp_data(data_fdir, data_fname, all_circs);
	// Initialize SZS
	szs = make_unique<StnZoneSrch>(num_iters_szs, max_circs, v_uav, all_circs);
	pso = make_unique<PSO>(num_ptcls, max_no_impr_pso, 
		num_iters_pso, c1, c2, w_min, w_max, tol_pso);
	pso->create_swarm_with_acs(num_ants, num_iters_acs, num_iters_2opt, 
		max_no_impr_acs, v_truck, q0, alpha, beta, rho, tol_acs, budget);
	// Initialize timer
	if (set_timer) { timer = make_unique<Timer>(); }
}
 
int main(int argc, char* argv[])
{
	// string config_fname = "general.tddp.json";
	// data_fname = "bubbles1";
	// bool print_config = true;
	// bool write_vertex = true;
	string config_fname = argv[1];
	data_fname = argv[2];
	bool print_config = std::stoi(argv[3]);
	bool write_vertex = std::stoi(argv[4]);

	fs::path config_fpath = root_dir / "configs" / config_fname;
	initialize(config_fpath.string(), print_config);
	if (set_timer) { timer->tick(); }
	vector<StnZone> stn_zones = szs->randomized_search();
	// DataIO::write_stn_zone(data_fdir, data_fname, stn_zones);
	Point start_pt = szs->start_pt(), end_pt = szs->end_pt();
	pso->init_swarm_with_stn_zones(start_pt, end_pt, stn_zones);
	pso->evolve_until_stop_criterion();
	if (set_timer) { timer->tock(); alg_time += timer->duration().count(); }
	Fitness final_fit;
	vector<Point> final_pos;
	pso->get_output(final_fit, final_pos);
	time_t end_time = std::chrono::system_clock::to_time_t(
		std::chrono::system_clock::now());
	if (print_result) {
		cout << end_time << " ";
		if (set_timer) { cout << "Time " << alg_time << " ms, "; }
		cout << "Cost " << final_fit.cost << ", Prize " << final_fit.prize << "\n";
	}
	if (write_vertex) {
		DataIO::write_stn_zone_vtx(data_fdir, data_fname, pso->get_stn_zones());
	}
	DataIO::write_final_sol(data_fdir, data_fname, end_time, ".tddp_sol", alg_time,
		final_fit.cost, final_fit.prize, final_pos);
	return 1;
}
