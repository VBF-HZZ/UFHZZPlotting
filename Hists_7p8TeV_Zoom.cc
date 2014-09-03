#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <map>
#include <utility>
#include <algorithm>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"
#include "Collection.h"
#include "MergedCollection.h"
#include "PlotHelper.h"

using namespace std;

typedef std::vector<double> dv;
typedef std::vector<TH1F*> hv;
typedef std::vector<THStack*> stackv;
typedef std::vector<ULong64_t> uv;
typedef std::vector<TString> tsv;
typedef std::vector<std::string> stringv;

double nevents_data(0);
double nevents_mc(0);

int Run = 8;

int main(int argc, char* argv[]){

  if (argc < 6){
    cout << "Usage: "<< argv[0] << " <fileDir> <plotXlow_zoom> <plotXhigh_zoom> <binSize_zoom> <Run>" << endl;
    return 1;
  }

  PlotHelper *helper = new PlotHelper();

  string outDir = (string)argv[1];
  TString ToutDir = (TString)outDir;
  TString epsPlotDir = ToutDir+"/eps";
  TString pngPlotDir = ToutDir+"/png";
  string txtDir     = outDir+"/txt";
 
  TString fileDir_7 = "rootfiles_7TeV_RMDUP";
  TString fileDir_8 = "rootfiles_8TeV_RMDUP";
  //TString fileDir_7 = "/scratch/osghpc/snowball/UF_HZZ4L_Analysis/FullAnalysis2012/Histogramming/rootFiles_Legacy_dataMC";
  //TString fileDir_8 = "/scratch/osghpc/snowball/UF_HZZ4L_Analysis/FullAnalysis2012/Histogramming_8TeV/rootFiles_Legacy_dataMC";
  //TString fileDir_8 = "/scratch/osghpc/snowball/UF_HZZ4L_Analysis/FullAnalysis2012/MEKD_Z4L/Histogramming/rootFiles";
  //TString fileDir_7 = "../Histogramming/rootFiles_44X";
  //TString fileDir_8 = "../Histogramming_8TeV/rootFiles_Moriond";

  TString treeName = "passedEvents_dataMC";
  TString treeName_7 = "passedEvents_dataMC";
  TString treeName_data = "AnaAfterHlt/passedEvents";
  TString treeName_data_7 = "passedEventsNoDuplicate";

  // ------------- Histograms -------------- //

  //lumi
  double lumi_7 = 5.051;  
  double lumi_8 = 19.712;  

  // run: 7, 8 or 78
  Run = atoi(argv[5]);
  std::cout << "Processing Run " << Run << std::endl;

  //low mass zoom
  Double_t binSize_zoom = atof(argv[4]);
  Double_t xLow4l_zoom = (float)atof(argv[2]);// 70 // range of tree
  Double_t xHigh4l_zoom = (float)atof(argv[3]);// 190 //range of tree
  Double_t plotXlow_zoom = atof(argv[2]);//range of plot 
  Double_t plotXhigh_zoom = atof(argv[3]);//range of plot
  Int_t nBins4l_zoom = (plotXhigh_zoom-plotXlow_zoom)/binSize_zoom;
  Double_t plotYmax_zoom = 36;
  TString yAxisLabel_zoom = (TString)"Events / "+Form("%.0f",binSize_zoom)+" GeV";
  TString xAxisLabel_m4l = "m_{4l} [GeV]";
  TString xAxisLabel_m4e = "m_{4e} [GeV]";
  TString xAxisLabel_m4mu = "m_{4#mu} [GeV]";
  TString xAxisLabel_m2e2mu = "m_{2e2#mu} [GeV]";

  //reducible background
  Double_t  nBgEvents_4mu = 3.67; //11.5;//11.54;//3.67;//0.58+3.09;
  Double_t  nBgEvents_4e = 7.46;//3.6;//3.67;//7.46;//1.37+6.09;
  Double_t  nBgEvents_2e2mu = 11.54; //7.4;//7.46;//11.54; //2.29+9.25;

  // reproducing wrong combination of the norms
  bool doWrongCombine=false;
  if (doWrongCombine)
  {
    nBgEvents_4mu = 11.54; // 2e2mu 
    nBgEvents_4e = 3.67; // 4mu
    nBgEvents_2e2mu = 7.46; // 4e
  }
  if (Run==7)
  {
    nBgEvents_4mu = 0.58;
    nBgEvents_4e = 1.37;
    nBgEvents_2e2mu = 2.29;
    if (doWrongCombine)
    {
      nBgEvents_4mu = 2.29; // 2e2mu
      nBgEvents_4e = 0.58; // 4mu
      nBgEvents_2e2mu = 1.37; // 4e
    }
  }
  else if (Run==8)
  {
    nBgEvents_4mu = 3.09;
    nBgEvents_4e = 6.09;
    nBgEvents_2e2mu = 9.25;
    if (doWrongCombine)
    {
      nBgEvents_4mu = 9.25; // 2e2mu
      nBgEvents_4e = 3.09; // 4mu
      nBgEvents_2e2mu = 6.09; // 4e
    }
  }

  Double_t  nBgEvents = nBgEvents_4mu+nBgEvents_4e+nBgEvents_2e2mu;
  Double_t  l_4e_p0 = 1.0; //norm
  Double_t  l_4e_p1 = 0.117024;  // normalisation of landau1  (i.e. 2P2F)
  Double_t  l_4e_p2 = 195.407;   // MPV of Landau1
  Double_t  l_4e_p3 = 38.9472;   // sigma of Landau
  Double_t  l_4e_p4 = 3.68476;    // a  in pol1 = a + bx
  Double_t  l_4e_p5 = -0.00580439;    // b  in pol1 = a + bx
  Double_t  l_4e_p6 = 2.57278;     // normalisation of Landau2  (i.e. 3P1F)
  Double_t  l_4e_p7 = 110.862;     // MPV of Landau2
  Double_t  l_4e_p8 = 9.59455;     // sigma of Laudau2
  Double_t  l_4mu_p0 = 1.0; // norm
  Double_t  l_4mu_p1 = 129.0;  // MPV
  Double_t  l_4mu_p2 = 15.0;  // sigma
  Double_t  l_2e2mu_p0 = 1.0;  // norm
  Double_t  l_2e2mu_p1 = 0.00439895;  // normalisation of landau1  (i.e. 2P2F 2mu2e)
  Double_t  l_2e2mu_p2 = 195.407;     // MPV of Landau1
  Double_t  l_2e2mu_p3 = 38.9472;     // sigma of Landau
  Double_t  l_2e2mu_p4 = 3.68476;     // a  in pol1 = a + bx
  Double_t  l_2e2mu_p5 = -0.00580439;     // b  in pol1 = a + bx
  Double_t  l_2e2mu_p6 = 1.92769;     // normalisation of Landau2  (i.e. 3P1F 2mu2e)
  Double_t  l_2e2mu_p7 = 110.862;     // MPV of Landau2
  Double_t  l_2e2mu_p8 = 9.59455;     // sigma of Laudau2
  Double_t  l_2e2mu_p9 = 1.;               // normalisation of Landau3  (i.e. 2e2mu), set to one by convention
  Double_t  l_2e2mu_p10 = 129;            // MPV of Landau3
  Double_t  l_2e2mu_p11 = 15.0;            // sigma of Landau3   
  
  //cuts
  double massZ1Cut = 40;
  double massZ2Cut = 12;
  double m4lCut = 0;
  double ptMuCut = 5;
  double ptElCut = 7;

  // ------------ Style --------------- //
  Color_t redBgColor = kGreen-5;
  Color_t ZZBgColor = kAzure-9;
  Color_t h126Color = kOrange+10;
  Color_t h350Color = kRed+1;
  gStyle->SetOptStat(0000);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);


  double scale126_4e = 1;
  double scale126_4mu = 1;
  double scale126_2e2mu = 1.;
  double scale126_4l = 1.;

  //scale126_4e = 0.6+0.780174*lumi_8/5.261;
  //scale126_4mu = 1.148+1.594*lumi_8/5.261;
  //scale126_2e2mu = 1.53+2.22227*lumi_8/5.261;
  //scale126_4l = scale126_4e+scale126_4mu+scale126_2e2mu;

  //if (Run==7) 
  //{
  //  scale126_4e *= lumi_7/(lumi_8+lumi_7);
  //  scale126_4mu *= lumi_7/(lumi_8+lumi_7);
  //  scale126_2e2mu *= lumi_7/(lumi_8+lumi_7);
  //  scale126_4l *= lumi_7/(lumi_8+lumi_7);
  //}
  //else if (Run==8)
  //{
  //  scale126_4e *= lumi_8/(lumi_8+lumi_7);
  //  scale126_4mu *= lumi_8/(lumi_8+lumi_7);
  //  scale126_2e2mu *= lumi_8/(lumi_8+lumi_7);
  //  scale126_4l *= lumi_8/(lumi_8+lumi_7);
  //}


  // 7 TeV
  //Collection *CggH126_7    = new Collection("mH126",fileDir_8+ "/mH_126_ggH.root"     ,treeName,lumi_7,false,massZ1Cut,massZ2Cut,m4lCut);
  //Collection *CVBF126_7    = new Collection("mH126",fileDir_8+ "/mH_126_VBF.root"     ,treeName,lumi_7,false,massZ1Cut,massZ2Cut,m4lCut);
  Collection *CggH350_7    = new Collection("mH350ggH",fileDir_7+ "/mH_350_ggH.root"  ,treeName_7,lumi_7,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *CVBF350_7    = new Collection("mH350VBF",fileDir_7+ "/mH_350_VBF.root"  ,treeName_7,lumi_7,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz4mu_7     = new Collection("ZZ4mu",fileDir_7+ "/ZZ_4mu.root"         ,treeName_7,lumi_7,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz4e_7      = new Collection("ZZ4e",fileDir_7+ "/ZZ_4e.root"           ,treeName_7,lumi_7,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz2e2mu_7   = new Collection("ZZ2e2mu",fileDir_7+ "/ZZ_2e2mu.root"     ,treeName_7,lumi_7,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz4tau_7    = new Collection("ZZ4tau",fileDir_7+ "/ZZto4tau.root"      ,treeName_7,lumi_7,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz2e2tau_7  = new Collection("ZZ2e2tau",fileDir_7+ "/ZZto2e2tau.root"  ,treeName_7,lumi_7,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz2mu2tau_7 = new Collection("ZZ2mu2tau",fileDir_7+ "/ZZto2mu2tau.root" ,treeName_7,lumi_7,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Cggzz4l_7    = new Collection("ggZZ4l",fileDir_7+ "/ggZZ_4l.root"       ,treeName_7,lumi_7,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Cggzz2l2l_7  = new Collection("ggZZ2l2l",fileDir_7+ "/ggZZ_2e2mu.root"  ,treeName_7,lumi_7,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Cdata_7      = new Collection("Data",fileDir_7+ "/Data_7TeV.root"            ,treeName_data_7,1,true,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);

  // 8 TeV
  Collection *CSMH126_8    = new Collection("mH126",fileDir_8+ "/mH_126_SMH.root"     ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *CggH126_8    = new Collection("mH126",fileDir_8+ "/mH_126_ggH_powheg15.root"     ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *CVBF126_8    = new Collection("mH126",fileDir_8+ "/mH_126_VBF.root"     ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *CTTH126_8    = new Collection("mH126",fileDir_8+ "/mH_126_TTH.root"     ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *CWH126_8    = new Collection("mH126",fileDir_8+ "/mH_126_WH.root"     ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *CZH126_8    = new Collection("mH126",fileDir_8+ "/mH_126_ZH.root"     ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *CggH350_8    = new Collection("mH350ggH",fileDir_8+ "/mH_350_ggH.root"  ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *CVBF350_8    = new Collection("mH350VBF",fileDir_8+ "/mH_350_VBF.root"  ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz4mu_8     = new Collection("ZZ4mu",fileDir_8+ "/ZZ_4mu.root"         ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz4e_8      = new Collection("ZZ4e",fileDir_8+ "/ZZ_4e.root"           ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz2e2mu_8   = new Collection("ZZ2e2mu",fileDir_8+ "/ZZ_2e2mu.root"     ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz4tau_8    = new Collection("ZZ4tau",fileDir_8+ "/ZZto4tau.root"      ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz2e2tau_8  = new Collection("ZZ2e2tau",fileDir_8+ "/ZZto2e2tau.root"  ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz2mu2tau_8 = new Collection("ZZ2mu2tau",fileDir_8+ "/ZZto2mu2tau.root" ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  //Collection *Cggzz4l_8    = new Collection("ggZZ4l",fileDir_7+ "/ggZZ_4l.root"       ,treeName_7,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  //Collection *Cggzz2l2l_8  = new Collection("ggZZ2l2l",fileDir_7+ "/ggZZ_2e2mu.root"  ,treeName_7,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  //Collection *Cggzz4l_8    = new Collection("ggZZ4l",fileDir_7+ "/ggZZ_4l.root"       ,treeName_7,lumi_8*1.18,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  //Collection *Cggzz2l2l_8  = new Collection("ggZZ2l2l",fileDir_7+ "/ggZZ_2e2mu.root"  ,treeName_7,lumi_8*1.18,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Cggzz4l_8    = new Collection("ggZZ4l",fileDir_8+ "/ggZZ_4l.root"       ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut,0.1);
  Collection *Cggzz2l2l_8  = new Collection("ggZZ2l2l",fileDir_8+ "/ggZZ_2e2mu.root"  ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut,0.1);

  //Collection *CZZJetsTo4L_8     = new Collection("ZZJetsTo4L",fileDir_8+ "/ZZJetsTo4L.root"         ,treeName,lumi_8,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Cdata_8      = new Collection("Data",fileDir_8+ "/Data_8TeV.root"            ,treeName_data,1,true,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);


  std::vector<Collection*> cvZZ, cvH350,cvggzz,cvH126,cvData;
  if (Run==7||Run==78) cvZZ.push_back(Czz4mu_7);
  if (Run==7||Run==78) cvZZ.push_back(Czz4e_7);
  if (Run==7||Run==78) cvZZ.push_back(Czz2e2mu_7);
  if (Run==8||Run==78) cvZZ.push_back(Czz4mu_8);
  if (Run==8||Run==78) cvZZ.push_back(Czz4e_8);
  if (Run==8||Run==78) cvZZ.push_back(Czz2e2mu_8);
  //MergedCollection *MCZZ_noTaus = new MergedCollection("ZZnoTaus",cvZZ);
  if (Run==7||Run==78) cvZZ.push_back(Czz4tau_7);
  if (Run==7||Run==78) cvZZ.push_back(Czz2e2tau_7);
  if (Run==7||Run==78) cvZZ.push_back(Czz2mu2tau_7);
  if (Run==8||Run==78) cvZZ.push_back(Czz4tau_8);
  if (Run==8||Run==78) cvZZ.push_back(Czz2e2tau_8);
  if (Run==8||Run==78) cvZZ.push_back(Czz2mu2tau_8);
  //MergedCollection *MCZZ_taus = new MergedCollection("ZZwithTaus",cvZZ);
  if (Run==7||Run==78) cvZZ.push_back(Cggzz4l_7);
  if (Run==7||Run==78) cvZZ.push_back(Cggzz2l2l_7);
  if (Run==8||Run==78) cvZZ.push_back(Cggzz4l_8);
  if (Run==8||Run==78) cvZZ.push_back(Cggzz2l2l_8);
  //if (Run==8||Run==78) cvZZ.push_back(CZZJetsTo4L_8);

  if (Run==7||Run==78) cvggzz.push_back(Cggzz4l_7);
  if (Run==7||Run==78) cvggzz.push_back(Cggzz2l2l_7);
  if (Run==8||Run==78) cvggzz.push_back(Cggzz4l_8);
  if (Run==8||Run==78) cvggzz.push_back(Cggzz2l2l_8);
  if (Run==7||Run==78) cvH350.push_back(CggH350_7);
  if (Run==7||Run==78) cvH350.push_back(CVBF350_7);
  if (Run==8||Run==78) cvH350.push_back(CggH350_8);
  if (Run==8||Run==78) cvH350.push_back(CVBF350_8);
  //if (Run==8||Run==78) cvH126.push_back(CggH126_8);
  if (Run==8||Run==78) cvH126.push_back(CSMH126_8);
  if (Run==8||Run==78) cvH126.push_back(CVBF126_8);
  if (Run==8||Run==78) cvH126.push_back(CTTH126_8);
  if (Run==8||Run==78) cvH126.push_back(CWH126_8);
  if (Run==8||Run==78) cvH126.push_back(CZH126_8);

  if (Run==8||Run==78) cvData.push_back(Cdata_8);
  if (Run==7||Run==78) cvData.push_back(Cdata_7);

  MergedCollection *MCZZ = new MergedCollection("ZZ",cvZZ);
  //MergedCollection *MCggZZ = new MergedCollection("ggZZ",cvggzz);
  MergedCollection *MCH350 = new MergedCollection("mH350",cvH350);
  MergedCollection *MCH126 = new MergedCollection("mH126",cvH126);
  MergedCollection *MCData = new MergedCollection("Data",cvData);

  // ------------ Legend ------------ //

  TLegend *leg = new TLegend(0.6,0.58,0.85,0.88);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TLatex *CP = new TLatex();
  CP->SetNDC(kTRUE);
  CP->SetTextSize(0.03);
  CP->SetTextAlign(31);
  CP->SetTextFont(42);
  CP->SetTextAlign(11);

  // ================================================================================================= //

  TF1* redBgFunc_4e = new TF1("redBgFunc_4e", "[0]*(landau(1) * (1 + exp( pol1(4))) + landau(6))", 10, 1000);
  redBgFunc_4e->FixParameter(0,l_4e_p0);
  redBgFunc_4e->FixParameter(1,l_4e_p1);
  redBgFunc_4e->FixParameter(2,l_4e_p2);
  redBgFunc_4e->FixParameter(3,l_4e_p3);
  redBgFunc_4e->FixParameter(4,l_4e_p4);
  redBgFunc_4e->FixParameter(5,l_4e_p5);
  redBgFunc_4e->FixParameter(6,l_4e_p6);
  redBgFunc_4e->FixParameter(7,l_4e_p7);
  redBgFunc_4e->FixParameter(8,l_4e_p8);

  TF1* redBgFunc_4mu = new TF1("redBgFunc_4mu", "landau", 10, 1000);
  redBgFunc_4mu->FixParameter(0,l_4mu_p0);
  redBgFunc_4mu->FixParameter(1,l_4mu_p1);
  redBgFunc_4mu->FixParameter(2,l_4mu_p2);

  TF1* redBgFunc_2e2mu = new TF1("redBgFunc_2e2mu", "[0]*(landau(1) * (1 + exp( pol1(4))) + landau(6) + landau(9))", 10, 1000);
  redBgFunc_2e2mu->FixParameter(0,l_2e2mu_p0);
  redBgFunc_2e2mu->FixParameter(1,l_2e2mu_p1);
  redBgFunc_2e2mu->FixParameter(2,l_2e2mu_p2);
  redBgFunc_2e2mu->FixParameter(3,l_2e2mu_p3);
  redBgFunc_2e2mu->FixParameter(4,l_2e2mu_p4);
  redBgFunc_2e2mu->FixParameter(5,l_2e2mu_p5);
  redBgFunc_2e2mu->FixParameter(6,l_2e2mu_p6);
  redBgFunc_2e2mu->FixParameter(7,l_2e2mu_p7);
  redBgFunc_2e2mu->FixParameter(8,l_2e2mu_p8);
  redBgFunc_2e2mu->FixParameter(9,l_2e2mu_p9);
  redBgFunc_2e2mu->FixParameter(10,l_2e2mu_p10);
  redBgFunc_2e2mu->FixParameter(11,l_2e2mu_p11);

  //norm
  l_4e_p0 = nBgEvents_4e/redBgFunc_4e->Integral(10, 2000);
  l_4mu_p0 = nBgEvents_4mu/redBgFunc_4mu->Integral(10, 2000);
  l_2e2mu_p0 = nBgEvents_2e2mu/redBgFunc_2e2mu->Integral(10, 2000);

  redBgFunc_4e->FixParameter(0,l_4e_p0);
  redBgFunc_4mu->FixParameter(0,l_4mu_p0);
  redBgFunc_2e2mu->FixParameter(0,l_2e2mu_p0);

  // 4l
  TF1* redBgFunc_4l = new TF1("redBgFunc_4l", "[0]*(landau(1) * (1 + exp( pol1(4))) + landau(6))+landau(9)+[12]*(landau(13) * (1 + exp( pol1(16))) + landau(18) + landau(21))", 10, 1000);
  redBgFunc_4l->FixParameter(0,l_4e_p0);
  redBgFunc_4l->FixParameter(1,l_4e_p1);
  redBgFunc_4l->FixParameter(2,l_4e_p2);
  redBgFunc_4l->FixParameter(3,l_4e_p3);
  redBgFunc_4l->FixParameter(4,l_4e_p4);
  redBgFunc_4l->FixParameter(5,l_4e_p5);
  redBgFunc_4l->FixParameter(6,l_4e_p6);
  redBgFunc_4l->FixParameter(7,l_4e_p7);
  redBgFunc_4l->FixParameter(8,l_4e_p8);
  redBgFunc_4l->FixParameter(9,l_4mu_p0);
  redBgFunc_4l->FixParameter(10,l_4mu_p1);
  redBgFunc_4l->FixParameter(11,l_4mu_p2);
  redBgFunc_4l->FixParameter(12,l_2e2mu_p0);
  redBgFunc_4l->FixParameter(13,l_2e2mu_p1);
  redBgFunc_4l->FixParameter(14,l_2e2mu_p2);
  redBgFunc_4l->FixParameter(15,l_2e2mu_p3);
  redBgFunc_4l->FixParameter(16,l_2e2mu_p4);
  redBgFunc_4l->FixParameter(17,l_2e2mu_p5);
  redBgFunc_4l->FixParameter(18,l_2e2mu_p6);
  redBgFunc_4l->FixParameter(19,l_2e2mu_p7);
  redBgFunc_4l->FixParameter(20,l_2e2mu_p8);
  redBgFunc_4l->FixParameter(21,l_2e2mu_p9);
  redBgFunc_4l->FixParameter(22,l_2e2mu_p10);
  redBgFunc_4l->FixParameter(23,l_2e2mu_p11);


  // ------------ m4l --------------- //
  TH1F* histm4l_h126_zoom = new TH1F("m4l_h126_zoom","Mass of 4 leptons; m_{4l} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4l_h350_zoom = new TH1F("m4l_h350_zoom","Mass of 4 leptons; m_{4l} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4l_ZZ_zoom = new TH1F("m4l_ZZ_zoom","Mass of 4 leptons; m_{4l} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4l_ZX_zoom = new TH1F("m4l_ZX_zoom","Mass of 4 leptons; m_{4l} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom); 
  TH1F* histm4l_data_zoom = new TH1F("m4l_data_zoom","Mass of 4 leptons; m_{4l} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);

  double histm4l_data_zoom_integral, histm4l_data_zoom_integralerr, histm4l_data_zoom_entries;
  double histm4l_ZZ_zoom_integral, histm4l_ZZ_zoom_integralerr, histm4l_ZZ_zoom_entries;
  double histm4l_ZX_zoom_integral, histm4l_ZX_zoom_integralerr;
  double histm4l_h126_zoom_integral, histm4l_h126_zoom_integralerr, histm4l_h126_zoom_entries;

  helper->getHistFromTF1(redBgFunc_4l,histm4l_ZX_zoom,-1);
  histm4l_data_zoom_entries = helper->fillHistFromVariable(MCData,histm4l_data_zoom,"mass4l");
  TGraphAsymmErrors* histm4l_dataE_zoom;
  helper->getAsymErr(histm4l_data_zoom,histm4l_dataE_zoom);
    
  histm4l_ZZ_zoom_entries = helper->fillHistFromVariable(MCZZ,histm4l_ZZ_zoom,"mass4l");
  helper->fillHistFromVariable(MCH350,histm4l_h350_zoom,"mass4l");
  histm4l_h126_zoom_entries = helper->fillHistFromVariable(MCH126,histm4l_h126_zoom,"mass4l");
  //histm4l_h126_zoom->Scale(scale126_4l/histm4l_h126_zoom->Integral());

  histm4l_data_zoom_integral = histm4l_data_zoom->IntegralAndError(1, histm4l_data_zoom->GetNbinsX(), histm4l_data_zoom_integralerr);
  histm4l_ZZ_zoom_integral = histm4l_ZZ_zoom->IntegralAndError(1, histm4l_ZZ_zoom->GetNbinsX(), histm4l_ZZ_zoom_integralerr);
  histm4l_ZX_zoom_integral = histm4l_ZX_zoom->IntegralAndError(1, histm4l_ZX_zoom->GetNbinsX(), histm4l_ZX_zoom_integralerr);
  histm4l_h126_zoom_integral = histm4l_h126_zoom->IntegralAndError(1, histm4l_h126_zoom->GetNbinsX(), histm4l_h126_zoom_integralerr);

  //histm4l_data_zoom_entries = histm4l_data_zoom->GetEffectiveEntries();
  //histm4l_ZZ_zoom_entries = histm4l_ZZ_zoom->GetEffectiveEntries();
  //histm4l_h126_zoom_entries = histm4l_h126_zoom->GetEffectiveEntries();

  // test give 2 events for ZZ background and 1 event more to Signal:
  //histm4l_ZZ_zoom->Scale((7.5+histm4l_ZZ_zoom->Integral())/histm4l_ZZ_zoom->Integral());
  //histm4l_h126_zoom->Scale((2.0+histm4l_h126_zoom->Integral())/histm4l_h126_zoom->Integral());

  nevents_data = histm4l_data_zoom->Integral();
  nevents_mc = histm4l_ZZ_zoom->Integral();
  nevents_mc += histm4l_ZX_zoom->Integral();
  nevents_mc += histm4l_h126_zoom->Integral();
  std::cout << "Integral Data: " <<  nevents_data << "; Integral MC: " << nevents_mc << std::endl;

  helper->setHistProperties(histm4l_ZZ_zoom,ZZBgColor,1001,helper->filled);
  helper->setHistProperties(histm4l_ZX_zoom,redBgColor,1001,helper->filled);
  helper->setHistProperties(histm4l_h126_zoom,h126Color,1001,helper->line);
  helper->setHistProperties(histm4l_h350_zoom,h350Color,1001,helper->line);
  helper->setHistProperties(histm4l_data_zoom,kBlack,1001,helper->markers,20);

  helper->addLegendEntry(leg,histm4l_data_zoom,"Data"    ,"P");
  helper->addLegendEntry(leg,histm4l_ZZ_zoom  ,"ZZ/Z#gamma*"       ,"F");
  helper->addLegendEntry(leg,histm4l_ZX_zoom ,"Z+X"       ,"F");
  helper->addLegendEntry(leg,histm4l_h126_zoom,"m_{H}=126","F");

  dv highestBin_m4lz;
  highestBin_m4lz.push_back(helper->findHighestBin(histm4l_ZX_zoom));
  highestBin_m4lz.push_back(helper->findHighestBin(histm4l_ZZ_zoom));
  highestBin_m4lz.push_back(helper->findHighestBin(histm4l_h126_zoom));
  highestBin_m4lz.push_back(helper->findHighestBin(histm4l_h350_zoom));
  highestBin_m4lz.push_back(helper->findHighestBin(histm4l_data_zoom));

  THStack *stackM4l_zoom = new THStack("stackM4l_zoom", "stackM4l_zoom");
  stackM4l_zoom->Add(histm4l_ZX_zoom);
  stackM4l_zoom->Add(histm4l_ZZ_zoom);
  stackM4l_zoom->Add(histm4l_h126_zoom);
  //stackM4l_zoom->Add(histm4l_h350_zoom);

  plotYmax_zoom = 35.5;
  helper->draw7p8Plot(stackM4l_zoom,histm4l_dataE_zoom,leg,CP,pngPlotDir+"/m4l_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".png",
	   xAxisLabel_m4l,yAxisLabel_zoom,lumi_7,lumi_8,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);
  helper->draw7p8Plot(stackM4l_zoom,histm4l_dataE_zoom,leg,CP,epsPlotDir+"/m4l_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".eps",
		   xAxisLabel_m4l,yAxisLabel_zoom,lumi_7,lumi_8,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);



  // ------------- m4e ------------ //
  TH1F* histm4e_h126_zoom = new TH1F("m4e_h126_zoom","Mass of 4 leptons; m_{4e} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4e_h350_zoom = new TH1F("m4e_h350_zoom","Mass of 4 leptons; m_{4e} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4e_ZZ_zoom = new TH1F("m4e_ZZ_zoom","Mass of 4 leptons; m_{4e} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4e_ZX_zoom = new TH1F("m4e_ZX_zoom","Mass of 4 leptons; m_{4e} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom); 
  TH1F* histm4e_data_zoom = new TH1F("m4e_data_zoom","Mass of 4 leptons; m_{4e} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);

  double histm4e_data_zoom_integral, histm4e_data_zoom_integralerr, histm4e_data_zoom_entries;
  double histm4e_ZZ_zoom_integral, histm4e_ZZ_zoom_integralerr, histm4e_ZZ_zoom_entries;
  double histm4e_ZX_zoom_integral, histm4e_ZX_zoom_integralerr;
  double histm4e_h126_zoom_integral, histm4e_h126_zoom_integralerr, histm4e_h126_zoom_entries;

  //helper->getHistFromTF1(redBgFunc_4e,histm4e_ZX_zoom, nBgEvents_4e);
  helper->getHistFromTF1(redBgFunc_4e,histm4e_ZX_zoom, -1);
  histm4e_data_zoom_entries = helper->fillHistFromVariable(MCData,histm4e_data_zoom,"mass4e");
  TGraphAsymmErrors* histm4e_dataE_zoom;
  helper->getAsymErr(histm4e_data_zoom,histm4e_dataE_zoom);
    
  histm4e_ZZ_zoom_entries = helper->fillHistFromVariable(MCZZ,histm4e_ZZ_zoom,"mass4e");
  helper->fillHistFromVariable(MCH350,histm4e_h350_zoom,"mass4e");
  histm4e_h126_zoom_entries = helper->fillHistFromVariable(MCH126,histm4e_h126_zoom,"mass4e");
  //histm4e_h126_zoom->Scale(scale126_4e/histm4e_h126_zoom->Integral());

  histm4e_data_zoom_integral = histm4e_data_zoom->IntegralAndError(1, histm4e_data_zoom->GetNbinsX(), histm4e_data_zoom_integralerr);
  histm4e_ZZ_zoom_integral = histm4e_ZZ_zoom->IntegralAndError(1, histm4e_ZZ_zoom->GetNbinsX(), histm4e_ZZ_zoom_integralerr);
  histm4e_ZX_zoom_integral = histm4e_ZX_zoom->IntegralAndError(1, histm4e_ZX_zoom->GetNbinsX(), histm4e_ZX_zoom_integralerr);
  histm4e_h126_zoom_integral = histm4e_h126_zoom->IntegralAndError(1, histm4e_h126_zoom->GetNbinsX(), histm4e_h126_zoom_integralerr);

  //histm4e_data_zoom_entries = histm4e_data_zoom->GetEffectiveEntries();
  //histm4e_ZZ_zoom_entries = histm4e_ZZ_zoom->GetEffectiveEntries();
  //histm4e_h126_zoom_entries = histm4e_h126_zoom->GetEffectiveEntries();

  helper->setHistProperties(histm4e_ZZ_zoom,ZZBgColor,1001,helper->filled);
  helper->setHistProperties(histm4e_ZX_zoom,redBgColor,1001,helper->filled);
  helper->setHistProperties(histm4e_h126_zoom,h126Color,1001,helper->line);
  helper->setHistProperties(histm4e_h350_zoom,h350Color,1001,helper->line);
  helper->setHistProperties(histm4e_data_zoom,kBlack,1001,helper->markers,20);

  dv highestBin_m4ez;
  highestBin_m4ez.push_back(helper->findHighestBin(histm4e_ZX_zoom));
  highestBin_m4ez.push_back(helper->findHighestBin(histm4e_ZZ_zoom));
  highestBin_m4ez.push_back(helper->findHighestBin(histm4e_h126_zoom));
  highestBin_m4ez.push_back(helper->findHighestBin(histm4e_h350_zoom));
  highestBin_m4ez.push_back(helper->findHighestBin(histm4e_data_zoom));

  THStack *stackM4e_zoom = new THStack();
  stackM4e_zoom->Add(histm4e_ZX_zoom);
  stackM4e_zoom->Add(histm4e_ZZ_zoom);
  stackM4e_zoom->Add(histm4e_h126_zoom);
  //stackM4e_zoom->Add(histm4e_h350_zoom);

  plotYmax_zoom = 29.5;
  helper->draw7p8Plot(stackM4e_zoom,histm4e_dataE_zoom,leg,CP,pngPlotDir+"/m4e_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".png",
		      xAxisLabel_m4e,yAxisLabel_zoom,lumi_7,lumi_8,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);
  helper->draw7p8Plot(stackM4e_zoom,histm4e_dataE_zoom,leg,CP,epsPlotDir+"/m4e_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".eps",
		      xAxisLabel_m4e,yAxisLabel_zoom,lumi_7,lumi_8,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);



  // ------------ m4mu ------------- //
  TH1F* histm4mu_h126_zoom = new TH1F("m4mu_h126_zoom","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4mu_h350_zoom = new TH1F("m4mu_h350_zoom","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4mu_ZZ_zoom = new TH1F("m4mu_ZZ_zoom","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4mu_ZX_zoom = new TH1F("m4mu_ZX_zoom","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom); 
  TH1F* histm4mu_data_zoom = new TH1F("m4mu_data_zoom","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);

  double histm4mu_data_zoom_integral, histm4mu_data_zoom_integralerr, histm4mu_data_zoom_entries;
  double histm4mu_ZZ_zoom_integral, histm4mu_ZZ_zoom_integralerr, histm4mu_ZZ_zoom_entries;
  double histm4mu_ZX_zoom_integral, histm4mu_ZX_zoom_integralerr;
  double histm4mu_h126_zoom_integral, histm4mu_h126_zoom_integralerr, histm4mu_h126_zoom_entries;

  //helper->getHistFromTF1(redBgFunc_4mu,histm4mu_ZX_zoom,nBgEvents_4mu);
  helper->getHistFromTF1(redBgFunc_4mu,histm4mu_ZX_zoom,-1);
  histm4mu_data_zoom_entries = helper->fillHistFromVariable(MCData,histm4mu_data_zoom,"mass4mu");
  TGraphAsymmErrors* histm4mu_dataE_zoom;
  helper->getAsymErr(histm4mu_data_zoom,histm4mu_dataE_zoom);
    
  histm4mu_ZZ_zoom_entries = helper->fillHistFromVariable(MCZZ,histm4mu_ZZ_zoom,"mass4mu");
  helper->fillHistFromVariable(MCH350,histm4mu_h350_zoom,"mass4mu");
  histm4mu_h126_zoom_entries = helper->fillHistFromVariable(MCH126,histm4mu_h126_zoom,"mass4mu");
  //histm4mu_h126_zoom->Scale(scale126_4mu/histm4mu_h126_zoom->Integral());


  histm4mu_data_zoom_integral = histm4mu_data_zoom->IntegralAndError(1, histm4mu_data_zoom->GetNbinsX(), histm4mu_data_zoom_integralerr);
  histm4mu_ZZ_zoom_integral = histm4mu_ZZ_zoom->IntegralAndError(1, histm4mu_ZZ_zoom->GetNbinsX(), histm4mu_ZZ_zoom_integralerr);
  histm4mu_ZX_zoom_integral = histm4mu_ZX_zoom->IntegralAndError(1, histm4mu_ZX_zoom->GetNbinsX(), histm4mu_ZX_zoom_integralerr);
  histm4mu_h126_zoom_integral = histm4mu_h126_zoom->IntegralAndError(1, histm4mu_h126_zoom->GetNbinsX(), histm4mu_h126_zoom_integralerr);

  //histm4mu_data_zoom_entries = histm4mu_data_zoom->GetEffectiveEntries();
  //histm4mu_ZZ_zoom_entries = histm4mu_ZZ_zoom->GetEffectiveEntries();
  //histm4mu_h126_zoom_entries = histm4mu_h126_zoom->GetEffectiveEntries();

  helper->setHistProperties(histm4mu_ZZ_zoom,ZZBgColor,1001,helper->filled);
  helper->setHistProperties(histm4mu_ZX_zoom,redBgColor,1001,helper->filled);
  helper->setHistProperties(histm4mu_h126_zoom,h126Color,1001,helper->line);
  helper->setHistProperties(histm4mu_h350_zoom,h350Color,1001,helper->line);
  helper->setHistProperties(histm4mu_data_zoom,kBlack,1001,helper->markers,20);

  dv highestBin_m4muz;
  highestBin_m4muz.push_back(helper->findHighestBin(histm4mu_ZX_zoom));
  highestBin_m4muz.push_back(helper->findHighestBin(histm4mu_ZZ_zoom));
  highestBin_m4muz.push_back(helper->findHighestBin(histm4mu_h126_zoom));
  highestBin_m4muz.push_back(helper->findHighestBin(histm4mu_h350_zoom));
  highestBin_m4muz.push_back(helper->findHighestBin(histm4mu_data_zoom));

  THStack *stackM4mu_zoom = new THStack();
  stackM4mu_zoom->Add(histm4mu_ZX_zoom);
  stackM4mu_zoom->Add(histm4mu_ZZ_zoom);
  stackM4mu_zoom->Add(histm4mu_h126_zoom);
  //stackM4mu_zoom->Add(histm4mu_h350_zoom);

  plotYmax_zoom = 29.5;
  helper->draw7p8Plot(stackM4mu_zoom,histm4mu_dataE_zoom,leg,CP,pngPlotDir+"/m4mu_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".png",
		      xAxisLabel_m4mu,yAxisLabel_zoom,lumi_7,lumi_8,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);
  helper->draw7p8Plot(stackM4mu_zoom,histm4mu_dataE_zoom,leg,CP,epsPlotDir+"/m4mu_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".eps",
		      xAxisLabel_m4mu,yAxisLabel_zoom,lumi_7,lumi_8,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);




  // ------------ m2e2mu ------------- //
  TH1F* histm2e2mu_h126_zoom = new TH1F("m2e2mu_h126_zoom","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm2e2mu_h350_zoom = new TH1F("m2e2mu_h350_zoom","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm2e2mu_ZZ_zoom = new TH1F("m2e2mu_ZZ_zoom","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm2e2mu_ZX_zoom = new TH1F("m2e2mu_ZX_zoom","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom); 
  TH1F* histm2e2mu_data_zoom = new TH1F("m2e2mu_data_zoom","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);

  double histm2e2mu_data_zoom_integral, histm2e2mu_data_zoom_integralerr, histm2e2mu_data_zoom_entries;
  double histm2e2mu_ZZ_zoom_integral, histm2e2mu_ZZ_zoom_integralerr, histm2e2mu_ZZ_zoom_entries;
  double histm2e2mu_ZX_zoom_integral, histm2e2mu_ZX_zoom_integralerr;
  double histm2e2mu_h126_zoom_integral, histm2e2mu_h126_zoom_integralerr, histm2e2mu_h126_zoom_entries;

  //helper->getHistFromTF1(redBgFunc_2e2mu,histm2e2mu_ZX_zoom,nBgEvents_2e2mu);
  helper->getHistFromTF1(redBgFunc_2e2mu,histm2e2mu_ZX_zoom,-1);
  histm2e2mu_data_zoom_entries = helper->fillHistFromVariable(MCData,histm2e2mu_data_zoom,"mass2e2mu");
  TGraphAsymmErrors* histm2e2mu_dataE_zoom;
  helper->getAsymErr(histm2e2mu_data_zoom,histm2e2mu_dataE_zoom);
    
  histm2e2mu_ZZ_zoom_entries = helper->fillHistFromVariable(MCZZ,histm2e2mu_ZZ_zoom,"mass2e2mu");
  helper->fillHistFromVariable(MCH350,histm2e2mu_h350_zoom,"mass2e2mu");
  histm2e2mu_h126_zoom_entries = helper->fillHistFromVariable(MCH126,histm2e2mu_h126_zoom,"mass2e2mu");
  //histm2e2mu_h126_zoom->Scale(scale126_2e2mu/histm2e2mu_h126_zoom->Integral());


  histm2e2mu_data_zoom_integral = histm2e2mu_data_zoom->IntegralAndError(1, histm2e2mu_data_zoom->GetNbinsX(), histm2e2mu_data_zoom_integralerr);
  histm2e2mu_ZZ_zoom_integral = histm2e2mu_ZZ_zoom->IntegralAndError(1, histm2e2mu_ZZ_zoom->GetNbinsX(), histm2e2mu_ZZ_zoom_integralerr);
  histm2e2mu_ZX_zoom_integral = histm2e2mu_ZX_zoom->IntegralAndError(1, histm2e2mu_ZX_zoom->GetNbinsX(), histm2e2mu_ZX_zoom_integralerr);
  histm2e2mu_h126_zoom_integral = histm2e2mu_h126_zoom->IntegralAndError(1, histm2e2mu_h126_zoom->GetNbinsX(), histm2e2mu_h126_zoom_integralerr);

  //histm2e2mu_data_zoom_entries = histm2e2mu_data_zoom->GetEffectiveEntries();
  //histm2e2mu_ZZ_zoom_entries = histm2e2mu_ZZ_zoom->GetEffectiveEntries();
  //histm2e2mu_h126_zoom_entries = histm2e2mu_h126_zoom->GetEffectiveEntries();

  helper->setHistProperties(histm2e2mu_ZZ_zoom,ZZBgColor,1001,helper->filled);
  helper->setHistProperties(histm2e2mu_ZX_zoom,redBgColor,1001,helper->filled);
  helper->setHistProperties(histm2e2mu_h126_zoom,h126Color,1001,helper->line);
  helper->setHistProperties(histm2e2mu_h350_zoom,h350Color,1001,helper->line);
  helper->setHistProperties(histm2e2mu_data_zoom,kBlack,1001,helper->markers,20);

  dv highestBin_m2e2muz;
  highestBin_m2e2muz.push_back(helper->findHighestBin(histm2e2mu_ZX_zoom));
  highestBin_m2e2muz.push_back(helper->findHighestBin(histm2e2mu_ZZ_zoom));
  highestBin_m2e2muz.push_back(helper->findHighestBin(histm2e2mu_h126_zoom));
  highestBin_m2e2muz.push_back(helper->findHighestBin(histm2e2mu_h350_zoom));
  highestBin_m2e2muz.push_back(helper->findHighestBin(histm2e2mu_data_zoom));

  THStack *stackM2e2mu_zoom = new THStack();
  stackM2e2mu_zoom->Add(histm2e2mu_ZX_zoom);
  stackM2e2mu_zoom->Add(histm2e2mu_ZZ_zoom);
  stackM2e2mu_zoom->Add(histm2e2mu_h126_zoom);
  //stackM2e2mu_zoom->Add(histm2e2mu_h350_zoom);

  plotYmax_zoom = 29.5; 
  helper->draw7p8Plot(stackM2e2mu_zoom,histm2e2mu_dataE_zoom,leg,CP,pngPlotDir+"/m2e2mu_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".png",
		      xAxisLabel_m2e2mu,yAxisLabel_zoom,lumi_7,lumi_8,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);
  helper->draw7p8Plot(stackM2e2mu_zoom,histm2e2mu_dataE_zoom,leg,CP,epsPlotDir+"/m2e2mu_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".eps",
		      xAxisLabel_m2e2mu,yAxisLabel_zoom,lumi_7,lumi_8,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);
  

  // ================================================================================================= //

  // print yields
  std::cout << " ======================================== " << std::endl;
  std::cout << "   Run " << Run << " " << plotXlow_zoom << " < m4l < " << plotXhigh_zoom << std::endl;
  std::cout << " ======================================== " << std::endl;
  std::cout << "4l :" << std::endl;
  std::cout << " Data:  " << histm4l_data_zoom_integral << " +/- " << histm4l_data_zoom_integralerr << "; Entries: " << histm4l_data_zoom_entries << std::endl;
  std::cout << " ZZ:    " << histm4l_ZZ_zoom_integral << " +/- " << histm4l_ZZ_zoom_integralerr << "; Entries: " << histm4l_ZZ_zoom_entries << std::endl;
  std::cout << " ZX:    " << histm4l_ZX_zoom_integral << " +/- " << histm4l_ZX_zoom_integralerr << std::endl;
  std::cout << " mH126: " << histm4l_h126_zoom_integral << " +/- " << histm4l_h126_zoom_integralerr << "; Entries: " << histm4l_h126_zoom_entries << std::endl;
  std::cout << " ======================================== " << std::endl;

  std::cout << " ======================================== " << std::endl;
  std::cout << "4e :" << std::endl;
  std::cout << " Data:  " << histm4e_data_zoom_integral << " +/- " << histm4e_data_zoom_integralerr << "; Entries: " << histm4e_data_zoom_entries << std::endl;
  std::cout << " ZZ:    " << histm4e_ZZ_zoom_integral << " +/- " << histm4e_ZZ_zoom_integralerr << "; Entries: " << histm4e_ZZ_zoom_entries << std::endl;
  std::cout << " ZX:    " << histm4e_ZX_zoom_integral << " +/- " << histm4e_ZX_zoom_integralerr << std::endl;
  std::cout << " mH126: " << histm4e_h126_zoom_integral << " +/- " << histm4e_h126_zoom_integralerr << "; Entries: " << histm4e_h126_zoom_entries << std::endl;
  std::cout << " ======================================== " << std::endl;

  std::cout << " ======================================== " << std::endl;
  std::cout << "4mu :" << std::endl;
  std::cout << " Data:  " << histm4mu_data_zoom_integral << " +/- " << histm4mu_data_zoom_integralerr << "; Entries: " << histm4mu_data_zoom_entries << std::endl;
  std::cout << " ZZ:    " << histm4mu_ZZ_zoom_integral << " +/- " << histm4mu_ZZ_zoom_integralerr << "; Entries: " << histm4mu_ZZ_zoom_entries << std::endl;
  std::cout << " ZX:    " << histm4mu_ZX_zoom_integral << " +/- " << histm4mu_ZX_zoom_integralerr << std::endl;
  std::cout << " mH126: " << histm4mu_h126_zoom_integral << " +/- " << histm4mu_h126_zoom_integralerr << "; Entries: " << histm4mu_h126_zoom_entries << std::endl;
  std::cout << " ======================================== " << std::endl;

  std::cout << " ======================================== " << std::endl;
  std::cout << "2e2mu :" << std::endl;
  std::cout << " Data:  " << histm2e2mu_data_zoom_integral << " +/- " << histm2e2mu_data_zoom_integralerr << "; Entries: " << histm2e2mu_data_zoom_entries << std::endl;
  std::cout << " ZZ:    " << histm2e2mu_ZZ_zoom_integral << " +/- " << histm2e2mu_ZZ_zoom_integralerr << "; Entries: " << histm2e2mu_ZZ_zoom_entries << std::endl;
  std::cout << " ZX:    " << histm2e2mu_ZX_zoom_integral << " +/- " << histm2e2mu_ZX_zoom_integralerr << std::endl;
  std::cout << " mH126: " << histm2e2mu_h126_zoom_integral << " +/- " << histm2e2mu_h126_zoom_integralerr << "; Entries: " << histm2e2mu_h126_zoom_entries << std::endl;
  std::cout << " ======================================== " << std::endl;

  return 0;

}

