#include <string>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <cstdio>
#include <map>
#include <vector>
#include <limits>

#include "tomography.h"
#include "detector.h"
#include "event.h"
#include "Tray.h"
#include "datareader.h"

using std::string;
using std::cout;
using std::endl;
using std::flush;
using std::remove;
using std::map;
using std::vector;
using std::numeric_limits;

using boost::property_tree::ptree;

int main(int argc, char ** argv){
	if(argc<2){
		cout << "You must indicate a config file which contains the data run and ped run carac, you can specify a CosmicBench cosmic file (if you don't, default is \"config_default.json\")" << endl;
		return 1;
	}

	string config_file_wrapper = argv[1];
	ptree config_tree_wrapper;
	read_json(config_file_wrapper, config_tree_wrapper);
	ptree config_tree_bench;
	string config_file_CB = "config_default.json";
	if(argc>2){
		config_file_CB = argv[2];
	}
	read_json(config_file_CB, config_tree_bench);
	Tomography::Init(config_tree_bench);
	config_tree_bench.put<string>("Ped",config_tree_wrapper.get<string>("Ped"));
	config_tree_bench.put<string>("RMSPed",config_tree_wrapper.get<string>("RMSPed"));
	remove((config_tree_bench.get<string>("signal_file")).c_str());
	config_tree_bench.put<string>("data_file_basename",config_tree_wrapper.get<string>("pedrun_name"));
	config_tree_bench.put<int>("data_file_first",config_tree_wrapper.get<int>("pedrun_min"));
	config_tree_bench.put<int>("data_file_last",config_tree_wrapper.get<int>("pedrun_max"));
	DataReader blah_ped(config_tree_bench,true);
	blah_ped.process();
	blah_ped.compute_ped();
	blah_ped.read_ped();
	blah_ped.do_ped_sub();
	blah_ped.do_common_noise_sub();
	blah_ped.compute_RMSPed();
	remove((config_tree_bench.get<string>("signal_file")).c_str());
	config_tree_bench.put<string>("data_file_basename",config_tree_wrapper.get<string>("datarun_name"));
	config_tree_bench.put<int>("data_file_first",config_tree_wrapper.get<int>("datarun_min"));
	config_tree_bench.put<int>("data_file_last",config_tree_wrapper.get<int>("datarun_max"));

	Tray * rayFile = new Tray(config_tree_wrapper.get<string>("outTree"));
	DataReader blah(config_tree_bench,false);
	blah.read_ped();
	CosmicBench * bench = new CosmicBench(config_tree_bench);
	int total_det = bench->get_det_N_tot();
	double Z_Up = numeric_limits<double>::min();
	double Z_Down = numeric_limits<double>::max();
	for(int i=0;i<total_det;i++){
		double current_z = (bench->get_detector(i))->get_z();
		if(current_z>Z_Up) Z_Up = current_z;
		if(current_z<Z_Down) Z_Down = current_z;
	}
	long event_nb = 0;
	int Nevent = 0;
	double evttime = 0;
	while(!(blah.is_end())){
		if((event_nb%100) == 0) cout << "\r" << "event processed : " << event_nb << flush;
		blah.process_event();
		Nevent = blah.get_event_n();
		evttime = blah.get_evttime();
		blah.do_ped_CMN_sub_event();
		map<Tomography::det_type,vector<vector<vector<float> > > > current_data = blah.get_data();
		map<Tomography::det_type,vector<vector<vector<double> > > > current_data_d;
		for(map<Tomography::det_type,vector<vector<vector<float> > > >::iterator type_it = current_data.begin();type_it!=current_data.end();++type_it){
			current_data_d[type_it->first] = vector<vector<vector<double> > >((type_it->second).size());
			for(unsigned int i=0;i<(type_it->second).size();i++){
				current_data_d[type_it->first][i] = vector<vector<double> >((type_it->second)[i].size());
				for(unsigned int j=0;j<(type_it->second)[i].size();j++){
					current_data_d[type_it->first][i][j] = vector<double>((type_it->second)[i][j].size());
					for(unsigned int k=0;k<(type_it->second)[i][j].size();k++){
						current_data_d[type_it->first][i][j][k] = (type_it->second)[i][j][k];
					}
				}
			}
		}
		vector<Event*> events;
		for(int i=0;i<total_det;i++){
			Detector * det = bench->get_detector(i);
			events.push_back(det->build_event(current_data_d[det->get_type()][det->get_n_in_tree()],Nevent,evttime));
			(events.back())->MultiCluster();
		}
		CosmicBenchEvent CBEvent(bench,events);
		CBEvent.do_cuts();
		rayFile->fillTree(Nevent,evttime,CBEvent.get_absorption_rays(),Z_Up,Z_Down);
		for(vector<Event*>:: iterator event_it=events.begin();event_it!=events.end();++event_it){
			delete (*event_it);
		}
		event_nb++;
	}
	cout << "\r" << "event processed : " << event_nb << endl;
	rayFile->Write();
	rayFile->CloseFile();
	delete rayFile;
	delete bench;
	Tomography::Quit();
}