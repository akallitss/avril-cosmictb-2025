#ifndef Tanalyse_W_h
#define Tanalyse_W_h

#include <TTree.h>
#include <TFile.h>

#include <string>
#include <vector>
#include <map>

#include "tomography.h"
#include "detector.h"

using std::string;
using std::vector;

class Event;

class Tanalyse_W{
   public:
      TTree          *T;
      TFile          *saveFile;

      // Declaration of leaf types
      map<Tomography::det_type,unsigned short> det_N;

      Int_t           evn;
      Double_t        evttime;
      /*
      */
      map<Tomography::det_type,Int_t*> NClus;
      map<Tomography::det_type,Int_t*> Spark;
      map<Tomography::det_type,Double_t*> ClusAmpl;
      map<Tomography::det_type,Double_t*> ClusSize;
      map<Tomography::det_type,Double_t*> ClusPos;
      map<Tomography::det_type,Double_t*> ClusMaxStripAmpl;
      map<Tomography::det_type,Int_t*> ClusMaxStrip;
      map<Tomography::det_type,Double_t*> ClusMaxSample;
      map<Tomography::det_type,Double_t*> ClusTOT;
      map<Tomography::det_type,Double_t*> ClusT;
      map<Tomography::det_type,Double_t*> StripMaxAmpl;

      //Tanalyse(string saveFileName);
      Tanalyse_W(string saveFileName, map<Tomography::det_type,unsigned short> det_N_);
      ~Tanalyse_W();
      void Init();
      TTree * getTree() const;
      void Write();
      void CloseFile();
      void fillTree(int evn_, double evttime_, map<Tomography::det_type,vector<Event*> > events);
};

#endif

