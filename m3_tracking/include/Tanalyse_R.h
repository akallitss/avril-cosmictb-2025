//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 18 13:34:43 2014 by ROOT version 5.18/00b
// from TTree T/event
// found on file: run6_Pb_analyse.root
//////////////////////////////////////////////////////////

#ifndef Tanalyse_R_h
#define Tanalyse_R_h

#include <TChain.h>
#include <TFile.h>

#include "tomography.h"
#include "detector.h"
#include <map>

using std::map;

class CM_Detector;
class MG_Detector;
class MGv2_Detector;
class MGv3_Detector;
class CosmicBench;

class Tanalyse_R{
public:
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   long            current_entry;

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

   // List of branches
   TBranch        *b_evn;   //!
   TBranch        *b_evttime;
   /*
   */
   map<Tomography::det_type,TBranch*> b_NClus; //!
   map<Tomography::det_type,TBranch*> b_Spark;
   map<Tomography::det_type,TBranch*> b_ClusAmpl; //!
   map<Tomography::det_type,TBranch*> b_ClusSize; //!
   map<Tomography::det_type,TBranch*> b_ClusPos; //!
   map<Tomography::det_type,TBranch*> b_ClusMaxStripAmpl; //!
   map<Tomography::det_type,TBranch*> b_ClusMaxStrip; //!
   map<Tomography::det_type,TBranch*> b_ClusMaxSample; //!
   map<Tomography::det_type,TBranch*> b_ClusTOT; //!
   map<Tomography::det_type,TBranch*> b_ClusT; //!
   map<Tomography::det_type,TBranch*> b_StripMaxAmpl; //!

   Tanalyse_R();
   Tanalyse_R(TTree *tree, map<Tomography::det_type,unsigned short> det_N_);
   virtual ~Tanalyse_R();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, map<Tomography::det_type,unsigned short> det_N_);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual bool GetNext();
   virtual long GetCurrentEntry() const;
};

#endif

