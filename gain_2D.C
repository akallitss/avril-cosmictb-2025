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

void gain_2D() {

    // --- Ouvrir le fichier ---
    TFile *f = new TFile("outtree_fe55.root");
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
    vector<int> *DetCh = 0;         
    vector<double> *DetAmp = 0;
    vector<double> *ClustAmp = 0;
    vector<double> *ClustMaxStrip = 0;
  
    t->SetBranchAddress("DetCh", &DetCh);
    t->SetBranchAddress("DetAmp", &DetAmp);
    t->SetBranchAddress("ClustAmp", &ClustAmp);
    t->SetBranchAddress("ClustMaxStrip", &ClustMaxStrip);

    vector<vector<double> > Map_P2(NCHANNELS_P2, vector<double>(5, 0));
 
    ReadMapping("mappingP2_v2/", Map_P2);




//----------------histo for each pad---------
TH1::AddDirectory(kFALSE);

std::vector<TH1F*> hAmp(NCHANNELS_P2);

for (int ch = 0; ch < NCHANNELS_P2; ch++) {
    hAmp[ch] = new TH1F(
        Form("hAmp_%d", ch),
        Form("Channel %d;Amplitude;Counts", ch),
        150, 250, 5000
    );
}
Long64_t nentries = t->GetEntries();

for (Long64_t iev = 0; iev < nentries; iev++) {
    t->GetEntry(iev);

    int nhits = ClustAmp->size();
    for (int h = 0; h < nhits; h++) {
        int ch  = ClustMaxStrip->at(h);
        double a = ClustAmp->at(h);

        if (ch < 0 || ch >= NCHANNELS_P2) continue;
        if (a <= 250) continue;

        hAmp[ch]->Fill(a);
    }
}

// for (Long64_t iev = 0; iev < nentries; iev++) {
//     t->GetEntry(iev);

//     int nhits = DetCh->size();
//     for (int h = 0; h < nhits; h++) {
//         int ch  = DetCh->at(h);
//         double a = DetAmp->at(h);

//         if (ch < 0 || ch >= NCHANNELS_P2) continue;
//         if (a <= 100) continue;

//         hAmp[ch]->Fill(a);
//     }
// }




//----------------fit ---------
std::vector<double> gains(NCHANNELS_P2, 0.0);
vector<double> chi2s(NCHANNELS_P2, 0.0);
TH1F *hchi2 = new TH1F(
        "chi2","chi2;count",
        200, 0, 1000
    );

for (int ch = 0; ch < NCHANNELS_P2; ch++) {

    if (hAmp[ch]->GetEntries() < 10) continue;

   
    double mean = -1;
    double sigmaIni  = 0;
    double maxLoc = 0;
    int maxLocBin = -1;

    int nbins = hAmp[ch]->GetNbinsX();
    vector<int> localMaxBin;
    vector<double> localMax;
    vector<double> Means;

    for (int i = 6; i < nbins-3; i++) {
        double s_prev =  (hAmp[ch]->GetBinContent(i-2) + hAmp[ch]->GetBinContent(i-1) + hAmp[ch]->GetBinContent(i)) / 3.0 ;
        double s_curr =  (hAmp[ch]->GetBinContent(i-1) + hAmp[ch]->GetBinContent(i) + hAmp[ch]->GetBinContent(i+1)) / 3.0 ;
        double s_next =  (hAmp[ch]->GetBinContent(i) + hAmp[ch]->GetBinContent(i+1) + hAmp[ch]->GetBinContent(i+2)) / 3.0 ;

        if (s_curr > s_prev && s_curr > s_next) {
            localMaxBin.push_back(i);
            localMax.push_back(hAmp[ch]->GetBinContent(i));
            Means.push_back(hAmp[ch]->GetBinCenter(i));
        }
    }
    double newSigma = 3000;
    double m1 =0;
    double chi2 = -1;
    int maxTries = 10;
    int tries = 0;
    bool fitOk = false;

    while(fitOk == false && localMax.size() > 0 && tries < maxTries ){
        fitOk = true;
        tries++;
        double mean = -1;
        double sigmaIni  = 0;
        double maxLoc = -1;
        int maxLocBin = -1;
        int indice = -1;

        for(int maxInd = 0; maxInd < localMax.size(); maxInd++){
            if(localMax[maxInd]>maxLoc){
                maxLoc = localMax[maxInd];
                mean = Means[maxInd];
                maxLocBin = localMaxBin[maxInd];
                indice = maxInd;
            }
        }
       
        if(maxLoc < 10){
            fitOk = false;
            if (indice >= 0 && indice < localMax.size()) {
                localMax.erase(localMax.begin() + indice);
            }
            if (indice >= 0 && indice < localMaxBin.size()) {
                localMaxBin.erase(localMaxBin.begin() + indice);
            }
            if (indice >= 0 && indice < Means.size()) {
                Means.erase(Means.begin() + indice);
            }
            continue;
            }
        int j = maxLocBin; 
        while(hAmp[ch]->GetBinContent(j) > maxLoc/2){
            j++;
        }
        if (j >= nbins) { 
            if (indice >= 0 && indice < localMax.size()) {
                localMax.erase(localMax.begin() + indice);
            }
            if (indice >= 0 && indice < localMaxBin.size()) {
                localMaxBin.erase(localMaxBin.begin() + indice);
            }
            if (indice >= 0 && indice < Means.size()) {
                Means.erase(Means.begin() + indice);
            }
            fitOk = false;
            continue;
        }
        sigmaIni = (hAmp[ch]->GetBinCenter(j) - mean)/1.1774;

        int nmax = localMaxBin.size();

        TGraph *gMax = new TGraph(nmax);

        for (int k = 0; k < nmax; k++) {
            int bin = localMaxBin[k];
            double x = hAmp[ch]->GetBinCenter(bin);
            double y = hAmp[ch]->GetBinContent(bin);

            gMax->SetPoint(k, x, y);
        }
        gMax->SetMarkerStyle(20);     
        gMax->SetMarkerSize(2);
        gMax->SetMarkerColor(kRed);

        TMarker *mBest = new TMarker(mean, maxLoc, 29);
        mBest->SetMarkerColor(kBlue);
        mBest->SetMarkerSize(2);

        double xmin = hAmp[ch]->GetXaxis()->GetXmin();
        double xmax = hAmp[ch]->GetXaxis()->GetXmax();

        double fitMin = std::max(mean - 4*sigmaIni, xmin);
        double fitMax = std::min(mean + 4*sigmaIni, xmax);

        if (fitMax <= fitMin) { 
            if (indice >= 0 && indice < localMax.size()) {
                localMax.erase(localMax.begin() + indice);
            }
            if (indice >= 0 && indice < localMaxBin.size()) {
                localMaxBin.erase(localMaxBin.begin() + indice);
            }
            if (indice >= 0 && indice < Means.size()) {
                Means.erase(Means.begin() + indice);
            }
            fitOk = false;
            continue;
        }

        TF1 *gfit = new TF1(Form("gfit_%d", ch),"gaus",fitMin, fitMax);


        gfit->SetParameters(maxLoc, mean, sigmaIni);

        hAmp[ch]->Fit(gfit, "QR");

    //     if (ch % 10 == 0 || ch == 49) {
    //         cout << "ch=" << ch
    //     << " nLocalMax=" << localMaxBin.size()
    //     << " maxLoc=" << maxLoc
    //     << " mean ini =" << mean 
    //     << " mean after fit = " << gfit->GetParameter(1)
    //     << " sigma= " << gfit->GetParameter(2)  <<  endl;
    //     TCanvas *c = new TCanvas(Form("c_%d", ch),
    //                          Form("Channel %d", ch),
    //                          800, 600);
    //     hAmp[ch]->Draw();
    //     gfit->Draw("same");
    //     gMax->Draw("P same");
    //     mBest->Draw("same");
    //     c->Update();
    //     gPad->WaitPrimitive();

    // }
    
    m1 = gfit->GetParameter(1);
    newSigma = gfit->GetParameter(2);
    chi2 = gfit->GetChisquare();

    if (ch % 10 == 0 || ch == 49)cout << "chi2 = " << chi2 << endl;

    if(newSigma > 700 ||newSigma < 20 || m1 < 0 || indice < 0){fitOk = false;}

    delete gfit;
    delete gMax;
    delete mBest;

    if (indice >= 0 && indice < localMax.size()) {
        localMax.erase(localMax.begin() + indice);
        }
    if (indice >= 0 && indice < localMaxBin.size()) {
        localMaxBin.erase(localMaxBin.begin() + indice);
        }
    if (indice >= 0 && indice < Means.size()) {
        Means.erase(Means.begin() + indice);
        }
    }

    hchi2->Fill(chi2);
    gains[ch] = 1.033 * m1;
    if(!fitOk){ cout << "fit pas ok pour ch " << ch << endl;
    gains[ch]= -1;}
    

}
 TCanvas *c1 = new TCanvas("c1",
                             "chi2",
                             800, 600);
hchi2->Draw();
double max = gains[0];
int max_i = -1;
int min_i = -1;
double min = gains[0];
for (int i = 0; i < gains.size(); i++){
    if (gains[i]< min){min = gains[i]; min_i = i;}
    if (gains[i]>max){max = gains[i]; max_i = i;}
}

cout << "max : " << max << " for channel : " << max_i << " min : " << min << " for channel : "<< min_i << endl;


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
        bool isBorder = (find(BorderPads.begin(), BorderPads.end(), i) != BorderPads.end());
        bool BorderSide[4] = {false, false, true, false};
      
        P2_Pad2[i] = new TPolyLine(5, &Pad_coord_X[0], &Pad_coord_Y[0]);
      
        P2_Pad2[i]->SetLineColor(kBlack);
        P2_Pad2[i]->SetLineWidth(1);
     
        double val = (gains[i])/(max);
        
        if(val >0){
        int idx = int(val * (gStyle->GetNumberOfColors() - 1));
        int color = TColor::GetColorPalette(idx);
        P2_Pad2[i]->SetFillColor(color);
        }
        else P2_Pad2[i]->SetFillColor(kBlack);
       
        
        
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

        // string tlChannel_str = "";
        // std::ostringstream oss;
        // oss << int(i);
        // tlChannel_str += oss.str();

        // tlChannel = new TLatex(dX-4,dY,tlChannel_str.c_str());
        // tlChannel->SetTextSize(0.01);
        // tlChannel->Draw();
    }

}