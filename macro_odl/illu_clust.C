#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1.h"
#include "TProfile.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "TPaveText.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TCrown.h"
#include "TBox.h"
#include "TPolyLine.h"
#include "TText.h"
#include "TLine.h"
#include "TArc.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TStyle.h"
#include "TMarker.h"

#include <iomanip>
#include <sstream>


using namespace std;

void ReadMapping(string folderMapping, vector<vector<double> > &Map_P2){

   int Connector; 
   string PinName;
   int Channel; 
   int PadIndex;
   string PadName; 
   float X,Y,Radius,Phi,DeltaPhi;
   ifstream inrms; string line;
  
   for(int b=0; b<10; b++){

      char connf[20]; sprintf(connf, "connector_%d_v2.txt",b);
      string fileMapping = connf;
      string fileconn = folderMapping + fileMapping;
      cout << "open mapping file : " << fileconn << endl; 
      inrms.open(fileconn.c_str());
      
      if(inrms.is_open()) while(std::getline(inrms, line)){
         //std::cout << line << std::endl;
         int nv = sscanf(line.c_str(),"%d %*s %d %d %*s %f %f %f %f %f", &Connector, &Channel, &PadIndex, &X, &Y, &Radius, &Phi, &DeltaPhi);
         //cout << "["<< nv<<"]"<< Connector << " " << X << endl;
         if(nv>=7) {
            
 
            Map_P2[Channel + 128*Connector][0] = X;
            Map_P2[Channel + 128*Connector][1] = Y;
            Map_P2[Channel + 128*Connector][2] = Phi;
            Map_P2[Channel + 128*Connector][3] = Connector;
            Map_P2[Channel + 128*Connector][4] = Channel;
            
            //cout<< PadIndex <<" "<< Channel + 128*Connector <<" "<< X <<" "<< Y <<endl;
         }

   } else cout << "failed to open mapping file : " << fileconn << endl; 
      inrms.close();
   }
}

void Transformation(vector<double> &X, vector<double> &Y, double dX, double dY, double Theta){

   for(int k=0; k<X.size(); k++){

      double X0 = X[k];
      double Y0 = Y[k];

      X[k] = dX + X0*TMath::Cos(Theta) - Y0*TMath::Sin(Theta);
      Y[k] = dY + X0*TMath::Sin(Theta) + Y0*TMath::Cos(Theta);
   }  
}

void illu_clust() {

   
 const int NCHANNELS_P2 = 1280;
    vector<vector<double> > Map_P2(NCHANNELS_P2, vector<double>(5, 0));
 
    ReadMapping("mappingP2_v2/", Map_P2);











TCanvas *c4 = new TCanvas("c4","c4", 1000, 1000);
    c4->SetRightMargin(0.15);
   
    gStyle->SetPalette(55);
    TH2D *hP2_gain_pad = new TH2D("hP2_gain_pad","hP2_gain_pad",540,60,600,540,0,540);
   
    hP2_gain_pad->SetStats(kFALSE);

    hP2_gain_pad->GetZaxis()->SetTitle("Gain");


    TPolyLine **P2_Pad2 = new TPolyLine *[NCHANNELS_P2];

    hP2_gain_pad->Draw("COLZ");
    hP2_gain_pad->SetContour(100); 

    TLatex *tlChannel;

    for(int i=0; i<NCHANNELS_P2; i++){
        double b = 1.0;
      
        vector<double> Pad_coord_X; Pad_coord_X.push_back(-(5+b));Pad_coord_X.push_back( (5+b));Pad_coord_X.push_back( (5+b));Pad_coord_X.push_back(-(5+b));Pad_coord_X.push_back(-(5+b));

        vector<double> Pad_coord_Y;Pad_coord_Y.push_back(-(5+b));Pad_coord_Y.push_back(-(5+b));Pad_coord_Y.push_back( (5+b));Pad_coord_Y.push_back( (5+b));Pad_coord_Y.push_back(-(5+b));
      
        double dX = Map_P2[i][0]; 
        double dY = Map_P2[i][1]; 
        double Phi = Map_P2[i][2];
        double Connector = Map_P2[i][3];
        double Channel = Map_P2[i][4];
      
        Transformation(Pad_coord_X, Pad_coord_Y, dX, dY, Phi);
       
      
        P2_Pad2[i] = new TPolyLine(5, &Pad_coord_X[0], &Pad_coord_Y[0]);
      
        P2_Pad2[i]->SetLineColor(kBlack);
        P2_Pad2[i]->SetLineWidth(1);
     
        if(i == 516 ){
          double val = 0.8;
        int idx = int(val * (gStyle->GetNumberOfColors() - 1));
        int color = TColor::GetColorPalette(idx);
        P2_Pad2[i]->SetFillColor(color);
        }
          if(i == 591 ){
          double val = 0.86;
        int idx = int(val * (gStyle->GetNumberOfColors() - 1));
        int color = TColor::GetColorPalette(idx);
        P2_Pad2[i]->SetFillColor(color);
        }
          if(i == 608 ){
          double val = 0.86;
        int idx = int(val * (gStyle->GetNumberOfColors() - 1));
        int color = TColor::GetColorPalette(idx);
        P2_Pad2[i]->SetFillColor(color);
        }
          if(i == 625 ){
          double val = 0.95;
        int idx = int(val * (gStyle->GetNumberOfColors() - 1));
        int color = TColor::GetColorPalette(idx);
        P2_Pad2[i]->SetFillColor(color);
        }
     
       
        
        
        P2_Pad2[i]->Draw("fsame"); 
        P2_Pad2[i]->Draw("same");

      

        // string tlChannel_str = "";
        // std::ostringstream oss;
        // oss << int(i);
        // tlChannel_str += oss.str();

        // tlChannel = new TLatex(dX-4,dY,tlChannel_str.c_str());
        // tlChannel->SetTextSize(0.01);
        // tlChannel->Draw();
    }

}