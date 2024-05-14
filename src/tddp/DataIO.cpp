#include "DataIO.h"

using std::ifstream, std::ofstream;
using std::getline;
using std::runtime_error;
using std::unique_ptr, std::make_unique;
using std::to_string;
using std::cout;

void DataIO::load_tddp_data(const fs::path data_dir, const string f_name,
	vector<Circle>& circs)
{
	// Initialize 
	if (!circs.empty()) { circs.clear(); }
	// Complete full path data file
	string f_full_path = (data_dir / "instances" / f_name).string();
	ifstream file(f_full_path);
	if (!file.is_open()) { throw runtime_error(f_full_path + " can't open."); }

	size_t cnt = 0;
	string x, y, unused, r, p, w;
	while (getline(file, x, ' ')) {
		getline(file, y, ' ');
		getline(file, unused, ' ');
		getline(file, r, ' ');
		getline(file, p, ' ');
		getline(file, w);
		if (cnt == 0) { circs.push_back(Circle(cnt, Point(stod(x), stod(y)))); }
		else if (cnt == 1) { circs.push_back(Circle(cnt, Point(stod(x), stod(y)))); }
		else {
			circs.push_back(
				Circle(cnt, Point(stod(x), stod(y)), stod(r), stod(p), stod(w)));
		}
		cnt++;
	}
	file.close();
}

void DataIO::write_stn_zone(const fs::path data_fdir, const string f_name,
	const vector<StnZone>& stn_zones)
{
	string token = f_name.substr(0, f_name.find(DOTDELIMITER));
	string out_fdir = ((data_fdir / "sz-vtxs") / (
		token + ".sz")).string();
	ofstream out_file(out_fdir, std::ios::out);
	if (!out_file.is_open()) throw runtime_error(out_fdir + " cannot open.");

	for (const StnZone& sz : stn_zones) {
		out_file << sz.p() << " " << sz.ctr().x << " " << sz.ctr().y << " " 
			<< sz.num_vetxs() << " ";
		for (const Point& pt : sz.vs()) {
			out_file << pt.x << " " << pt.y << " ";
		}
		vector<Circle> cs_ = sz.cs();
		for (size_t i = 0; i < cs_.size(); i++) {
			out_file << cs_[i].w() << " " << cs_[i].ctr().x << " " << cs_[i].ctr().y;
			if (i < cs_.size() - 1) { out_file << " "; }
			else { out_file << "\n"; }
		}
	}
	out_file.close();
}

void DataIO::write_stn_zone_vtx(const fs::path data_fdir, const string f_name,
	const vector<StnZone>& stn_zones)
{
	string token = f_name.substr(0, f_name.find(DOTDELIMITER));
	string out_fdir = ((data_fdir / "sz-vtxs") / (
		token + ".szppt")).string();
	ofstream out_file(out_fdir, std::ios::out);
	try { if (!out_file.is_open()) throw runtime_error(out_fdir + " cannot open."); }
	catch (const std::exception& e) { cout << e.what() << '\n'; }

	for (const StnZone& sz : stn_zones) {
		vector<Point> vs = sz.vs();
		if (sz.is_singleton()) { out_file << sz.id() << " " << sz.p()
			<< " " << sz.ctr().x << " " << sz.ctr().y << "\n"; }
		else {
			for (size_t i = 0; i < vs.size(); i++) {
				out_file << sz.id() << " " << sz.p() 
					<< " " << vs[i].x << " " << vs[i].y << "\n";
			}
		}
	}
	out_file.close();
}

void DataIO::write_nbhd_c(const fs::path data_fdir, const size_t cnt,
	const string suffix, const vector<NbhdWayPt>& nbhd_waypts)
{
	string of_dir = ((data_fdir / "sols") / ("ptcl" + to_string(cnt) + suffix)).string();
	ofstream out_file(of_dir, std::ios::app);
	if (!out_file.is_open()) { throw std::runtime_error(of_dir + " cannot open."); }

	for (size_t i = 0; i < nbhd_waypts.size(); i++) {
		out_file << nbhd_waypts[i].c;
		if (i < nbhd_waypts.size() - 1) { out_file << " "; }
		else { out_file << "\n"; }
	}
}

void DataIO::write_points(const fs::path data_fdir, const size_t cnt, 
	const string suffix, const double cost, 
	const double prize, const vector<Point>& pts)
{
	string of_dir = ((data_fdir / "sols") / ("ptcl" + to_string(cnt) + suffix)).string();
	ofstream out_file(of_dir, std::ios::app);
	if (!out_file.is_open()) { throw std::runtime_error(of_dir + " cannot open."); }

	out_file << cost << " " << prize << " ";
	for (size_t i = 0; i < pts.size(); i++) {
		out_file << pts[i].x << " " << pts[i].y;
		if (i < pts.size() - 1) { out_file << " "; }
		else { out_file << "\n"; }
	}
	out_file.close();
}

void DataIO::write_final_sol(const fs::path data_fdir, const string f_name,
	const time_t time_stamp, const string suffix, const double alg_time, 
	const double cost, const double prize, const vector<Point>& pos)
{
	string token = f_name.substr(0, f_name.find(DOTDELIMITER));
	string of_dir = ((data_fdir / "sols") / (token + suffix)).string();
	std::ofstream file(of_dir, std::ios::app);
	if (!file.is_open()) { throw std::runtime_error(of_dir + " cannot open."); }

	file << time_stamp << " " << alg_time << " " << cost << " " << prize;
	for (size_t i = 0; i < pos.size(); i++) {
		file << " " << pos[i].x << " " << pos[i].y;
	}
	file << "\n";
	file.close();
}

void DataIO::write_ant_fit(const fs::path data_fdir, const string f_name,
	const size_t ptcl_idx, const vector<double>& gb_fits)
{
	string token = f_name.substr(0, f_name.find(DOTDELIMITER));
	string of_dir = ((data_fdir / "sols") / (token + ".ptcl" + to_string(ptcl_idx))).string();
	std::ofstream file(of_dir, std::ios::app);
	if (!file.is_open()) { throw std::runtime_error(of_dir + " cannot open."); }
	for (size_t i = 0; i < gb_fits.size(); i++) {
		file << i << " " << gb_fits[i] << "\n";
	}
	file.close();
}
