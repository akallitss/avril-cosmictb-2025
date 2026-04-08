//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 15 13:53:10 2014 by ROOT version 5.34/10
// from TTree T/event
// found on file: ../MG2Dalone_200fC_shaping63_gasstarted_140818_11H23_025_signal.root
//////////////////////////////////////////////////////////

#ifndef Tsignal_R_h
#define Tsignal_R_h

#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <map>

#include "tomography.h"
#include "detector.h"

using std::vector;
using std::map;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Tsignal_R {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   long            current_entry;

   // Declaration of leaf types
   map<Tomography::det_type,unsigned short> det_N;
   Int_t           Nevent;
   Double_t        evttime;
   /*
   */
   map<Tomography::det_type,Float_t*> StripAmpl;
   map<Tomography::det_type,Float_t*> StripAmpl_ped;
   map<Tomography::det_type,Float_t*> StripAmpl_corr;

   // List of branches
   TBranch        *b_Nevent;   //!
   TBranch        *b_evttime;   //!
   /*
   */
   map<Tomography::det_type,TBranch*> b_StripAmpl; //!
   map<Tomography::det_type,TBranch*> b_StripAmpl_ped; //!
   map<Tomography::det_type,TBranch*> b_StripAmpl_corr; //!

   Tsignal_R(TTree *tree, map<Tomography::det_type,unsigned short> det_N_);
   Tsignal_R();
   virtual ~Tsignal_R();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, map<Tomography::det_type,unsigned short> det_N_);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual bool     GetNext();
   virtual long     GetCurrentEntry() const;
   template<typename T> vector<vector<T> > get_ampl(Tomography::det_type type_, unsigned short det_n_);
   template<typename T> vector<vector<T> > get_ampl_ped(Tomography::det_type type_, unsigned short det_n_);
   template<typename T> vector<vector<T> > get_ampl_raw(Tomography::det_type type_, unsigned short det_n_);

};

#endif
