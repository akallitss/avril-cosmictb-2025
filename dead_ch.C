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
         // if (Channel < 64){ Channel = 63 - Channel; }
         // else{ Channel = 191 - Channel;}
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

void dead_ch() {

   // --- Ouvrir le fichier ---
   TFile *f = new TFile("/data/cosmic_data/2026/W2/workdir_p2_m400v_d600v_cosmic_longrun_260205/outtree_aligned.root");
   if (!f || f->IsZombie()) {
      cout << "Impossible d'ouvrir outtree.root" << endl;
      return;
   }

   // --- Charger le TTree ---
   TTree *t = (TTree*)f->Get("Tout");
   if (!t) {
      cout << "Impossible de trouver le TTree 'Tout'" << endl;
      return;
   }
   const int NCHANNELS_P2 = 1280;

   // --- Branches ---
   double Chi2X, Chi2Y; 
   double RayX_corr, RayY_corr; 
   vector<int> *DetCh = 0;         
   vector<double> *DetX = 0;   
   vector<double> *DetY = 0;
   vector<double> *DetAmp = 0;
   vector<double> *rms = 0;

   t->SetBranchAddress("RayX_corr", &RayX_corr);
   t->SetBranchAddress("RayY_corr", &RayY_corr);
   t->SetBranchAddress("DetX", &DetX);
   t->SetBranchAddress("DetY", &DetY);
   t->SetBranchAddress("DetCh", &DetCh);
   t->SetBranchAddress("DetAmp", &DetAmp);  
   t->SetBranchAddress("Chi2X", &Chi2X);
   t->SetBranchAddress("Chi2Y", &Chi2Y);
   t->SetBranchAddress("DetRMS", &rms);

    vector<vector<double> > Map_P2(NCHANNELS_P2, vector<double>(5, 0));
 
    ReadMapping("mappingP2_v2/", Map_P2);

    map<int,int> passed;
    map<int,int> seen;

    Long64_t N = t->GetEntries();
    int rayon = 20;

    for (Long64_t iev = 0; iev < N; iev++) {

      t->GetEntry(iev);
            
      if (!(Chi2X  < 0.1))  continue;
      if (!(Chi2Y  < 0.1))  continue;

      double minDist = 999;
      int ch= -1;
      for(int i=0; i<NCHANNELS_P2; i++){
         double padx = Map_P2[i][0]; 
         double pady = Map_P2[i][1]; 
         double dist = sqrt((RayX_corr-padx)*(RayX_corr-padx) + (RayY_corr-pady)*(RayY_corr-pady));
         if (dist < minDist) {minDist = dist;ch=i;}
      }	
      if (minDist > 8) continue;
      //cout << minDist << endl;
      passed[ch]++;
      

      int nhits = DetX->size();

      for (int h = 0; h < nhits; h++) {
         if (DetAmp->at(h) <= 50) continue;
         
         double detx = DetX->at(h);
         double dety = DetY->at(h);
		   int detCh = DetCh->at(h);
    	   double distance = sqrt((RayX_corr-detx)*(RayX_corr-detx) + (RayY_corr-dety)*(RayY_corr-dety));
		   //if(distance<rayon) {seen[ch]++;break;}
         if(detCh == ch) {seen[ch]++;break;}
      }
   }
   map<int,double> amp;
   map<int,int> nb;
  
   for (Long64_t iev = 0; iev < N; iev++) {
    t->GetEntry(iev);  
    int nhits = DetAmp->size();
    for (int h = 0; h < nhits; h++) {
      int detCh = DetCh->at(h);
      amp[detCh]+= DetAmp->at(h);
      nb[detCh]++;
    }
   }

vector<double> amps;
TH1F *hamp= new TH1F("hamp", "amplitude", 100, 0, 5000);

   for(int i=0; i<NCHANNELS_P2; i++){
   double ampMoy = amp[i]/nb[i];
   amps.push_back(ampMoy);
   if(ampMoy < 25)cout << "dead ch by amp : " << i << endl;
   hamp->Fill(ampMoy);
}
    
    TCanvas *camp = new TCanvas("camp","camp",800,600);
    hamp->Draw();




int BorderPadsArray[] = {78,75,73,70,68,66,65,64,87,85,82,80,77,74,71,69,67,97,94,91,88,84,83,79,76,72,112,110,107,101,96,92,89,86,81,11,9,3,123,119,114,108,
   207,202,197,193,192,209,204,199,195,194,212,206,201,196,222,216,210,203,198,227,219,213,208,200,232,225,217,211,205,243,230,223,214,135,251,244,237,221,154,140,254,238,
   336,328,322,320,333,325,321,338,329,323,343,334,326,348,339,330,324,345,335,327,351,342,331,360,347,337,372,358,344,332,369,353,340,259,367,350,270,260,362,285,271,354,
   465,454,448,461,451,469,457,449,464,453,473,460,450,468,456,478,463,452,474,459,485,467,455,481,462,493,475,458,487,470,504,484,466,501,476,387,496,471,384,490,418,479,
   576,581,596,583,577,585,578,586,579,588,580,590,606,591,608,593,582,595,613,597,584,601,620,603,587,605,589,607,633,609,592,615,518,618,599,621,529,630,536,635,550,
   704,708,719,706,715,705,711,724,709,718,707,714,730,712,725,710,720,735,716,731,713,726,744,721,739,717,733,755,727,750,722,743,641,737,766,729,757,657,751,654,740,683,
   832,834,843,852,833,842,851,835,841,850,836,840,849,837,839,848,866,838,847,867,844,854,868,845,855,869,846,856,870,889,857,871,891,858,874,769,859,877,772,860,882,807,
   960,961,967,973,980,962,965,971,978,986,964,969,977,984,963,968,975,983,992,966,972,981,991,998,970,979,990,1004,1015,976,987,1005,1016,974,985,1006,1017,908,997,1007,1022,932,
1088,1089,1091,1095,1098,1102,1106,1110,1090,1092,1094,1097,1101,1105,1109,1113,1121,1093,1096,1100,1104,1108,1116,1122,1129,1134,1099,1103,1107,1117,1123,1131,1139,1147,1033,1112,1118,1126,1136,1027,1034,1060};
int NBorderPads = sizeof(BorderPadsArray) / sizeof(int);

vector<int> BorderPads(BorderPadsArray, BorderPadsArray + NBorderPads);

TCanvas *c4 = new TCanvas("c4","c4", 1000, 1000);
   c4->SetRightMargin(0.15);
   
   gStyle->SetPalette(55);
   TH2D *hP2_eff_pad = new TH2D("hP2_eff_pad","hP2_eff_pad",540,60,600,540,0,540);
   
   hP2_eff_pad->SetStats(kFALSE);
   hP2_eff_pad->SetMinimum(0);
   hP2_eff_pad->SetMaximum(1);
   hP2_eff_pad->GetZaxis()->SetTitle("Efficiency");


   TPolyLine **P2_Pad2 = new TPolyLine *[NCHANNELS_P2];

   hP2_eff_pad->SetBinContent(1, 1, 0.0001);
   hP2_eff_pad->Draw("COLZ");
   hP2_eff_pad->SetContour(100); 

   double valmoy = 0;
   int nbch =0;
   vector<double> effs;
   vector<int> chn;

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

         // if (Channel < 64){ Channel = 63 - Channel; }
         // else{ Channel = 191 - Channel;}
         
        Transformation(Pad_coord_X, Pad_coord_Y, dX, dY, Phi);

         bool isBorder = (find(BorderPads.begin(), BorderPads.end(), i) != BorderPads.end());
         bool BorderSide[4] = {false, false, true, false};
      
        P2_Pad2[i] = new TPolyLine(5, &Pad_coord_X[0], &Pad_coord_Y[0]);
      
        P2_Pad2[i]->SetLineColor(kBlack);
        P2_Pad2[i]->SetLineWidth(1);
     
	
	    double val = 0.0;

        map<int,int>::iterator itP = passed.find(i);
        map<int,int>::iterator itS = seen.find(i);

        if (itP != passed.end() && itP->second > 0) {
            int p = itP->second;
            int s = (itS != seen.end()) ? itS->second : 0;
            val = double(s) / double(p);
            
            valmoy += val;
            nbch++;
            if (Connector == 2){
                effs.push_back(val);
                chn.push_back(i);
                if(val < 0.3) cout <<"dead channel : " << Channel + 128*Connector << " ou " << Channel << " eff = " << val << endl;
                }
           
        }

   

        int idx = int(val * (gStyle->GetNumberOfColors() - 1));
        int color = TColor::GetColorPalette(idx);


        P2_Pad2[i]->SetFillColor(color);
        
        P2_Pad2[i]->Draw("fsame"); 
        P2_Pad2[i]->Draw("same");


         if(isBorder){
		      for (int s = 0; s < 4; s++) {
   		   if (!BorderSide[s]) continue;

   		   int p1 = s;
   		   int p2 = s + 1;

   		   TLine *edge = new TLine(
      		Pad_coord_X[p1], Pad_coord_Y[p1],
      		Pad_coord_X[p2], Pad_coord_Y[p2]
   		   );
   		   edge->SetLineColor(kBlack);
   		   edge->SetLineWidth(5);
   		   edge->Draw("same");
	         }
         }

      //   string tlChannel_str = "";
      //   std::ostringstream oss;
      //   oss << int(Channel + 128*Connector);
      //   tlChannel_str += oss.str();

      //   tlChannel = new TLatex(dX-4,dY,tlChannel_str.c_str());
      //   tlChannel->SetTextSize(0.01);
      //   tlChannel->Draw();
    }
    cout << endl;
    cout << "eff moy par pad = " << valmoy/nbch << endl;

    TH1F *heff = new TH1F("heff", "efficiency", 100, 0, 1);

    for(int v = 0; v < effs.size(); v++){
        heff->Fill(effs[v]);
    }
    
    TCanvas *ceff = new TCanvas("ceff","ceff",800,600);
    heff->Draw();


    TH1F *h = new TH1F("h", "RMS", 100, 0, 10);

    
    ifstream file("/data/cosmic_data/2025/W43/workdir_CosTb_P2_mapping_M400V_v2/RMSPed.dat");

    int col1, col2;
    float col3;

    while (file >> col1 >> col2 >> col3) {
      if(col1 == 4 || col1 == 5) {
        h->Fill(col3);
        
       cout << "channel : " << 64*col1 + col2 << " amp moy = " << amps[64*col1 + col2] << " rms = " << col3 << endl;
      }
    }

    file.close();

    
    TCanvas *c = new TCanvas("c","c",800,600);
    h->Draw();

   //  TGraph *g = new TGraph();
   //  g->SetTitle("efficiency vs channel;Channel;efficiency");
   //  int i = 0;

   //  for(int v = 0; v < effs.size(); v++){
   //      g->SetPoint(i, chn[v], effs[v]);
   //      i++;
   //  } 
 

   //  TCanvas *c1 = new TCanvas("c","eff vs Channel",900,600);

   //  g->SetMarkerStyle(20);
   //  g->SetMarkerSize(0.8);
   //  g->Draw("AP");
}
