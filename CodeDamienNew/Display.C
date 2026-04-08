#include "Riostream.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TLine.h"

//mode = 0 -> use StripAmpl 
//mode = 1 -> use StripAmpl_corrped

int Display(int drmin, int drmax, int mode, int nevent=10, int skipevent=0)
{
  gROOT->Reset();

  TTree *T = (TTree*)gDirectory->Get("T");

  TCanvas *cp = new TCanvas("cp","",800,600);
  //cp->Divide(4,2);
  //cp->Divide(2,1);
  cp->cd(1);
  int Nsample=32;
  double Ymin=-400;
  double Ymax=800;
  if(mode==0) {Ymin = 0, Ymax = 2500;}

  int dreammin=drmin;
  int dreammax=drmax;
  
  TH2F *h1_X = new TH2F("h1_X","EventDisplay;Sample;Amp",Nsample,0,Nsample,500,Ymin,Ymax);
  h1_X->SetStats(0);
  char cut[100];
  char char0[200];
  TLatex llx;
  llx.SetTextSize(0.08);
  llx.SetNDC();
  double isdata;
  for(int event=skipevent;event<nevent+skipevent;event++)
    {
      isdata=0;
      sprintf(cut,"Nevent==%d",event);
      
      int pmin=0*64;
      int pmax=pmin+64;
      
      cp->cd(1);
      for(int idream=dreammin;idream<dreammax;idream++)
	for(int strip=pmin;strip<pmax;strip++)
	  {
	    if(mode==1)sprintf(char0,"StripAmpl_corrped[%d][%d][]:TsampleNum[]>>h1_X",idream,strip);
	    else if(mode==0)sprintf(char0,"StripAmpl[%d][%d][]:TsampleNum[]>>h1_X",idream,strip);
	    if(strip==pmin && idream==dreammin) T->Draw(char0,cut,"l");
	    else T->Draw(char0,cut,"lsame");
	  }
      
      llx.DrawLatex(0.3,0.8, cut);
      cout << "Event number = " << event << endl;
      cp->Update();
      //gSystem->Sleep(2000);
    }
  return 1;
}
