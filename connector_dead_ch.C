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

void ReadMapping(string folderMapping, map<int,int> &Map_P2){

   int Pin; 
   int Channel; 
   
   ifstream inrms; string line;
  
    char connf[256]; sprintf(connf, "dream2pin.txt");
      string fileMapping = connf;
      string fileconn = folderMapping + fileMapping;
      cout << "open mapping file : " << fileconn << endl; 
      inrms.open(fileconn.c_str());
      int i = 3;
      
      if(inrms.is_open()) while(std::getline(inrms, line)){
         //std::cout << line << std::endl;
         int nv = sscanf(line.c_str(),"%*d %*d %*d %*d %d %*d %*d %d ", &Pin, &Channel);
         //cout << "["<< nv<<"]"<< Connector << " " << X << endl;
         if(nv>=2) {
            // if (Pin>68)Map_P2[Pin-4] = Channel-1;
            // else if(Pin>100)Map_P2[Pin-6] = Channel-1;
            // else Map_P2[Pin] = Channel-1;
            Map_P2[i] = Channel - 1;
            i ++;
            
        //cout<< Pin <<" "<< Channel <<endl;
         }

   } else cout << "failed to open mapping file : " << fileconn << endl; 
      inrms.close();

      cout << "dernier pin " << i << endl;
  
}

void ReadPed(string folderMapping, vector<vector<int> > &Map_Ch){
   int ct = 0;
   ifstream inrms; string line;
  
    char connf[256]; sprintf(connf, "selfTPOTFe_ped_fe55_feu3_50fC_3_pedthr_260218_09H40_000_03_thr.aux");
      string fileMapping = connf;
      string fileconn = folderMapping + fileMapping;
      cout << "open mapping file : " << fileconn << endl; 
      inrms.open(fileconn.c_str());
     
    if(inrms.is_open()) while(std::getline(inrms, line)){
      
         if (line.empty() || line[0] == '#') {
        cout << "ligne ignoree" << endl;
        continue;
    }

    if (line.find('D') == string::npos || line.find('C') == string::npos) {
        cout << "fin des lignes de mapping" << endl;
        break;
    }
        int index = atoi(line.c_str());
      if (ct == 0) cout << "index lu : " << index << endl;
      ct++;

        size_t posD = line.find('D');
        if (posD == std::string::npos) continue;
        int d = atoi(line.c_str() + posD + 1);
        if (ct == 0)cout << "dream lu" << endl;

        size_t posC = line.find('C', posD);
        if (posC == std::string::npos) continue;
        int c = atoi(line.c_str() + posC + 1);
        if (ct == 0)cout << "cannal lu" << endl;

        size_t posMax = line.find("Max");
        if (posMax == std::string::npos) continue;

        size_t posEq = line.find('=', posMax);
        if (posEq == std::string::npos) continue;

        int maxVal = atoi(line.c_str() + posEq + 1);
        if (ct == 0)cout << "max lu " <<endl;

        // std::cout << "index=" << index
        //           << " D=" << d
        //           << " C=" << c
        //           << " Max=" << maxVal << std::endl;

            if (index >= (int)Map_Ch.size()) {
            Map_Ch.resize(index + 1, vector<int>(3, 0));
        }

            Map_Ch[index][0] = d;
            Map_Ch[index][1] = c;
            Map_Ch[index][2] = maxVal;
        
   } else cout << "failed to open mapping file : " << fileconn << endl; 
      inrms.close();
  cout << "nb de ligne comptabilsees : " << ct << endl;
}


void connector_dead_ch(){

    const int NCHANNELS = 512;
    const int NPINS = 128;
    int connector = 2; 

    map<int,int> dream2pin;
    vector< vector<int> > Map_P2(NCHANNELS, vector<int>(3, 0));

    ReadMapping("",dream2pin);
    cout << " map oin size : " << dream2pin.size() << endl;
     for (std::map<int,int>::iterator it = dream2pin.begin(); it != dream2pin.end(); ++it) {
        int key = it->first;
        int value = it->second;
        //std::cout << "pin" << key << " -> channel " << value << std::endl;
    }
    ReadPed("/data/cosmic_data/2025/W43/test_P2/",Map_P2);
    cout << "vector size " << Map_P2.size() << endl;

    if (gROOT->FindObject("c3")) {
    delete gROOT->FindObject("c3");
}
    TCanvas *c3 = new TCanvas("c3","c3", 1000, 1000);
   
    TH2D *hp2_conn = new TH2D("hP2_conn","hP2_conn",140,0,140,190,0,190);

    hp2_conn->SetStats(kFALSE);
   
    hp2_conn->Draw();

    TPolyLine **P2_Conn = new TPolyLine *[NPINS];

    TPolyLine **pline = new TPolyLine *[4];

    TLatex *tlChannel;

    for (int conn = 0; conn < 4; conn++){

        double cX = 136;
        double cY = 174-50*conn;

        double X[5] = {cX-1,cX+1,cX+1,cX-1,cX-1};
        double Y[5] = {cY-2, cY-2, cY+2, cY+2,cY-2};

        pline[conn] = new TPolyLine(5,X,Y);

        pline[conn]->Draw("fsame");
        pline[conn]->Draw("same");

        string tlInd = "";
        std::ostringstream oss;
        oss << int(3);
        tlInd += oss.str();

        tlChannel = new TLatex(cX-0.5,cY,tlInd.c_str());
        tlChannel->SetTextSize(0.01);
        tlChannel->Draw();
        

        string tlConn_str = "Connector ";
        std::ostringstream oss2;
        oss2 << int(conn);
        tlConn_str += oss2.str();

        tlChannel = new TLatex(60,180-50*conn,tlConn_str.c_str());
        tlChannel->SetTextSize(0.03);
        tlChannel->Draw();
       
        int num = 1;
        for(int i = 0; i < NPINS; i++){
            
            int key = i+3; 
            if(dream2pin.find(key) == dream2pin.end()) continue;
            int ch = dream2pin[key];

            int idx = 128*conn + ch;
            if(idx >= (int)Map_P2.size()) continue;

            double coordX = (key%2==0) ? 135-i : 134-i;
            double coordY = (key%2==0) ? 160-50*conn : 170-50*conn;

            double x[5] = {coordX-1,coordX+1,coordX+1,coordX-1,coordX-1};
            double y[5] = {coordY-2, coordY-2, coordY+2, coordY+2,coordY-2};

            P2_Conn[i] = new TPolyLine(5,x,y);
            P2_Conn[i]->SetFillColor(38);

            if(Map_P2[idx][2]<600) P2_Conn[i]->SetFillColor(kRed);
            else if(Map_P2[idx][2]<1000) P2_Conn[i]->SetFillColor(kOrange);
            //else if(Map_P2[idx][2]<2000) P2_Conn[i]->SetFillColor(kBlue);

            P2_Conn[i]->Draw("fsame");
            P2_Conn[i]->Draw("same");

            string tlChannel_str = "";
            std::ostringstream oss;
            oss << int(num);
            if (key%2==0) num++;
            tlChannel_str += oss.str();

            tlChannel = new TLatex(coordX-0.5,coordY,tlChannel_str.c_str());
            tlChannel->SetTextSize(0.01);
            tlChannel->Draw();
            
                
        }
       
    }

    // for(int i=0;i<NPINS;i++){
    // delete P2_Conn[i];
    //     }
    // delete[] P2_Conn;
    // delete hp2_conn;
    c3->Update();
}