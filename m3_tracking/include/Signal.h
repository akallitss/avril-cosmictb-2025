#ifndef Signal_h
#define Signal_h
#include "Tsignal_R.h"
#include "detector.h"
#include "tomography.h"

#include <string>
#include <map>

#include <boost/property_tree/ptree.hpp>

using std::string;
using std::map;

using boost::property_tree::ptree;

class TProfile;
class TFile;

class Signal: public Tsignal_R, public CosmicBench{
	public:
		Signal(string configFilePath);
		Signal(ptree config_tree);
		~Signal();
		//process all pedestal and common noise subtracted events in signal tree to make the clustering and output the analyse tree
		void MultiCluster();
		//process all raw events in signal tree to make the clustering and output the analyse tree
		void MultiCluster_raw();
		//void ElecToAnalyse();
		//void ElecToRays(string outFileName);
		//display the signal shape corresponding to the given correction (raw, pedestal subtracted or pedestal and common noise subtracted) of multiple events from evn_min to evn_max
		void EventDisplay(int evn_min = 0, int evn_max = 20, Tomography::signal_type signal_correction = Tomography::corr);
		//compute and display the houg cluster and tracks for the given event
		void HoughTracking(long event_nb);
		//test function to output caracteristics variables of tracks computed using the cluster made with the convolution method
		void ConvClusterTest();
		map<Tomography::det_type,map<int,TProfile*> > SignalOverNoise();
		void SignalOverNoiseDisplay();
		void SignalDispersion();
		void DebugHoles(long event_nb);
		void NoiseLevels();
	protected:
		string analyseTree;
		bool use_srf;
		Tomography::elec_type electronic_type;
		string data_file_basename;
		string signalName;
		string PedName;
		string RMSName;
		long max_event;
		int data_file_first;
		int data_file_last;
		TFile * fIn;
		bool exists;
};
#endif
