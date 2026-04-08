#define write_rays_task_cpp

#include "task/write_rays_task.h"

#include "Tray.h"
#include "event.h"
#include "cluster.h"

Write_Rays_Task::Write_Rays_Task(Tray * writer_, double z_up_, double z_down_, string ampl_file_name_): Output_Task<ray_data>(){
	writer = writer_;
	z_up = z_up_;
	z_down = z_down_;
	ampl_file_name = ampl_file_name_;
	last_timestamp = time(NULL);
}

Write_Rays_Task::Write_Rays_Task(Tray * writer_, double z_up_, double z_down_, string ampl_file_name_, Typed_Task<ray_data> * next_task_): Output_Task<ray_data>(next_task_){
	writer = writer_;
	z_up = z_up_;
	z_down = z_down_;
	ampl_file_name = ampl_file_name_;
	last_timestamp = time(NULL);
}

Write_Rays_Task::~Write_Rays_Task(){

}
bool Write_Rays_Task::do_task(){
	ray_data * current_data = get_next_data();
	if(current_data->CBevent == NULL){
		delete current_data;
		return false;
	}
	for(vector<Ray>::const_iterator ray_it=(current_data->rays).begin();ray_it!=(current_data->rays).end();++ray_it){
		const vector<Cluster*> current_clus = ray_it->get_clus();
		map<int,double> current_ampl;
		for(vector<Cluster*>::const_iterator clus_it=current_clus.begin();clus_it!=current_clus.end();++clus_it){
			current_ampl[(*clus_it)->get_n_in_tree()] = (*clus_it)->get_maxStripAmpl();
			delete *clus_it;
		}
		pthread_mutex_lock(&IO_mutex);
		ampl_history.push(current_ampl);
		pthread_mutex_unlock(&IO_mutex);
	}
	time_t current_time = time(NULL);
	pthread_mutex_lock(&IO_mutex);
	// MODIFIY THIS VALUE TO CHANGE FEEDBACK TIME
	if((current_time - last_timestamp) > 300){
		map<int,double> ampl_mean;
		map<int,int> ampl_number;
		while(!(ampl_history.empty())){
			map<int,double> current_ampl = ampl_history.front();
			ampl_history.pop();
			for(map<int,double>::const_iterator ampl_it=current_ampl.begin();ampl_it!=current_ampl.end();++ampl_it){
				if(ampl_mean.count(ampl_it->first)>0){
					ampl_mean[ampl_it->first] += ampl_it->second;
					ampl_number[ampl_it->first]++;
				}
				else{
					ampl_mean[ampl_it->first] = ampl_it->second;
					ampl_number[ampl_it->first] = 1;
				}
			}
		}
		ofstream ampl_file(ampl_file_name.c_str());
		ampl_file << current_time << "\n";
		for(map<int,double>::iterator ampl_it=ampl_mean.begin();ampl_it!=ampl_mean.end();++ampl_it){
			ampl_file << ampl_it->first << " " << (ampl_it->second)/(ampl_number[ampl_it->first]) << "\n";
		}

		last_timestamp = current_time;
	}
	writer->fillTree((current_data->CBevent)->get_evn(), (current_data->CBevent)->get_evttime(), current_data->rays, z_up, z_down);
	if((data_treated%1000) == 0) writer->Write();
	pthread_mutex_unlock(&IO_mutex);
	if(next_task==NULL) delete current_data;
	else next_task->push_next_data(current_data);
	return true;
}
bool Write_Rays_Task::can_exec() const{
	return (writer!=NULL && (!is_queue_empty()));
}
