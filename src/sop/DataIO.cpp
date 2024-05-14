#include "DataIO.h"

using std::ifstream, std::ofstream;
using std::getline;
using std::to_string;

void DataIO::load_sop_data(const fs::path data_dir, const string f_name,
	vector<Circle>& circs)
{
	// Initialize 
	if (!circs.empty()) { circs.clear(); }
	// Complete full path data file
	string f_full_path = (data_dir / "instances" / f_name).string();
	ifstream file(f_full_path);
	if (!file.is_open()) { throw runtime_error(f_full_path + " can't open."); }

	size_t cnt = 0;
	string x, y, unused, r, p;
	while (getline(file, x, ' ')) {
		getline(file, y, ' ');
		getline(file, unused, ' ');
		getline(file, r, ' ');
		getline(file, p);
		if (cnt == 0) { circs.push_back(Circle(cnt, Point(stod(x), stod(y)))); }
		else if (cnt == 1) { circs.push_back(Circle(cnt, Point(stod(x), stod(y)))); }
		else {
			circs.push_back(Circle(cnt, Point(stod(x), stod(y)), stod(r), stod(p)));
		}
		cnt++;
	}
	file.close();
}

void DataIO::write_stn_zone_vtx(const fs::path data_fdir, const string f_name,
	const size_t t, const vector<StnZone>& stn_zones)
{
	string token = f_name.substr(0, f_name.find(DOTDELIMITER));
	string out_fdir = ((data_fdir / "sz-vtxs") / (
		token + ".iter" + to_string(t) + ".vert")).string();
	ofstream out_file(out_fdir, std::ios::out);
	try { if (!out_file.is_open()) throw runtime_error(out_fdir + " cannot open."); }
	catch (const std::exception& e) { cout << e.what() << '\n'; }

	for (const StnZone& sz : stn_zones) {
		vector<Point> vs = sz.vs();
		if (vs.empty()) { throw std::runtime_error("vertices empty."); }
		else {
			for (size_t i = 0; i < vs.size(); i++) {
				out_file << vs[i].x << " " << vs[i].y;
				if (i < vs.size() - 1) { out_file << " "; }
				else { out_file << "\n"; }
			}
		}
	}
	out_file.close();
}

void DataIO::write_sop_sol(const fs::path data_fdir, const string f_name,
	const time_t time_stamp, const string suffix, const double alg_time,
	const double cost, const double prize, const vector<Point>& pos)
{
	string token = f_name.substr(0, f_name.find(DOTDELIMITER));
	string of_dir = ((data_fdir / "sols") / (token + ".acs")).string();
	std::ofstream file(of_dir, std::ios::app);
	try { if (!file.is_open()) throw std::runtime_error(of_dir + " cannot open."); }
	catch (const std::exception& e) { std::cout << e.what() << '\n'; }

	file << time_stamp << " " << alg_time << " " << cost << " " << prize;
	for (size_t i = 0; i < pos.size(); i++) {
		file << " " << pos[i].x << " " << pos[i].y;
	}
	file << "\n";
	file.close();
}

void DataIO::write_sop_sol(const fs::path data_fdir, const string f_name,
	const string suffix, const double cost, const double prize,
	const vector<Point>& waypts)
{
	string token = f_name.substr(0, f_name.find(DOTDELIMITER));
	string of_dir = ((data_fdir / "data" / "sop" / "sols") / (token + suffix)).string();
	std::ofstream file(of_dir, std::ios::out);
	if (!file.is_open()) { throw runtime_error(of_dir + " cannot open."); }

	file << cost << " " << prize;
	for (size_t i = 0; i < waypts.size(); i++) {
		file << " " << waypts[i].x << " " << waypts[i].y;
	}
	file << "\n";
	file.close();
}

void DataIO::write_sz_p_pt(const fs::path data_fdir, const string f_name, 
	const vector<IdxValPt>& sz_p_pt)
{
	string token = f_name.substr(0, f_name.find(DOTDELIMITER));
	string of_dir = ((data_fdir / "sz-vtxs") / (token + ".szppt")).string();
	std::ofstream file(of_dir, std::ios::out);
	if (!file.is_open()) { throw runtime_error(of_dir + " cannot open."); }

	for (const IdxValPt& x : sz_p_pt) {
		file << x.first << " " << x.second.first << " " << x.second.second.x
			<< " " << x.second.second.y << "\n";
	}
	file.close();
}
