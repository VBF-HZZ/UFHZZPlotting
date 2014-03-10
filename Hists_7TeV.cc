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

int main(int argc, char* argv[]){

  using namespace std;

  if (argc != 9){
    cout << "Usage: ./makeHists7TeV.exe <fileDir> <plotXlow> <plotXhigh> <binSize> <plotXlow_zoom> <plotXhigh_zoom> <binSize_zoom> <print data>" << endl;
    return 1;
  }

  PlotHelper *helper = new PlotHelper();

  string outDir = (string)argv[1];
  TString ToutDir = (TString)outDir;
  TString epsPlotDir = ToutDir+"/eps";
  TString pngPlotDir = ToutDir+"/png";
  string txtDir     = outDir+"/txt";

  TString fileDir = "../Histogramming/rootFiles_Legacy_dataMC";
  TString treeName = "passedEvents_dataMC";

  // ------------- Histograms -------------- //

  //lumi
  double lumi = 5.051;
  double sqrts = 7;

  //whole range
  Double_t binSize = atof(argv[4]);
  Double_t xLow4l = -5;//range of tree 
  Double_t xHigh4l = 1005;//range of tree
  Double_t plotXlow = atof(argv[2]);//range of plot
  Double_t plotXhigh = atof(argv[3]);//range of plot
  Int_t nBins4l = (xHigh4l-xLow4l)/binSize;
  //Double_t plotYmax = 25;
  TString yAxisLabel = (TString)"Events / "+Form("%.0f",binSize)+" GeV";
  TString xAxisLabel_m4l = "m_{4l} [GeV]";
  TString xAxisLabel_m4e = "m_{4e} [GeV]";
  TString xAxisLabel_m4mu = "m_{4#mu} [GeV]";
  TString xAxisLabel_m2e2mu = "m_{2e2#mu} [GeV]";
  
  //low mass zoom
  Double_t binSize_zoom = atof(argv[7]);
  Double_t xLow4l_zoom = 70;//range of tree
  Double_t xHigh4l_zoom = 190;//range of tree
  Double_t plotXlow_zoom = atof(argv[5]);//range of plot 
  Double_t plotXhigh_zoom = atof(argv[6]);//range of plot
  Int_t nBins4l_zoom = (xHigh4l_zoom-xLow4l_zoom)/binSize_zoom;
  //Double_t plotYmax_zoom = 20;
  TString yAxisLabel_zoom = (TString)"Events / "+Form("%.0f",binSize_zoom)+" GeV";

  //reducible background
  Double_t  nBgEvents_4mu = 1.0;
  Double_t  nBgEvents_4e = 1.6;
  Double_t  nBgEvents_2e2mu = 2.6;
  Double_t  nBgEvents = nBgEvents_4mu+nBgEvents_4e+nBgEvents_2e2mu;
  Double_t  l_mu = 145.2;
  Double_t  l_sigma = 17.8;
  
  //cuts
  double massZ1Cut = 40;
  double massZ2Cut = 4;
  double m4lCut = 0;
  double ptMuCut = 5;
  double ptElCut = 7;

  bool printData = (bool)atoi(argv[8]);
  double scale126_4e = 0.6;
  double scale126_4mu = 1.148;
  double scale126_2e2mu = 1.53;
  double scale126_4l = 3.278;
  double percentMaxBin = 1.5;


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

  Collection *CggH126    = new Collection("mH126ggH","../Histogramming_8TeV/rootFiles_Legacy_dataMC/mH_126_ggH.root" ,"AnaAfterHlt/passedEvents",lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut, m4lCut);
  Collection *CVBF126    = new Collection("mH126VBF","../Histogramming_8TeV/rootFiles_Legacy_dataMC/mH_126_VBF.root" ,"AnaAfterHlt/passedEvents",lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut, m4lCut);
  //Collection *CggH126    = new Collection("mH126ggH",fileDir+"/mH_126_ggH.root"     ,treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  //Collection *CVBF126    = new Collection("mH126VBF",fileDir+"/mH_126_VBF.root"     ,treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *CggH350    = new Collection("mH350ggH",fileDir+"/mH_350_ggH.root"  ,treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *CVBF350    = new Collection("mH350VBF",fileDir+"/mH_350_VBF.root"  ,treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz4mu     = new Collection("ZZ4mu",fileDir+"/ZZ_4mu.root"         ,treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz4e      = new Collection("ZZ4e",fileDir+"/ZZ_4e.root"           ,treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz2e2mu   = new Collection("ZZ2e2mu",fileDir+"/ZZ_2e2mu.root"     ,treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz4tau    = new Collection("ZZ4tau",fileDir+"/ZZto4tau.root"      ,treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz2e2tau  = new Collection("ZZ2e2tau",fileDir+"/ZZto2e2tau.root"  ,treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Czz2mu2tau = new Collection("ZZ2mu2tau",fileDir+"/ZZto2mu2tau.root",treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Cggzz4l    = new Collection("ggZZ4l",fileDir+"/ggZZ_4l.root",treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Cggzz2l2l  = new Collection("ggZZ2l2l",fileDir+"/ggZZ_2e2mu.root",treeName,lumi,false,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);
  Collection *Cdata      = new Collection("Data",fileDir+"/Data.root"            ,"Ana/passedEvents",1,true,massZ1Cut,massZ2Cut,ptMuCut, ptElCut,m4lCut);

  std::vector<Collection*> cvZZ, cvH350,cvggzz,cvH126;
  cvZZ.push_back(Czz4mu);
  cvZZ.push_back(Czz4e);
  cvZZ.push_back(Czz2e2mu);
  MergedCollection *MCZZ_noTaus = new MergedCollection("ZZnoTaus",cvZZ);
  cvZZ.push_back(Czz4tau);
  cvZZ.push_back(Czz2e2tau);
  cvZZ.push_back(Czz2mu2tau);
  MergedCollection *MCZZ_taus = new MergedCollection("ZZwithTaus",cvZZ);
  cvZZ.push_back(Cggzz4l);
  cvZZ.push_back(Cggzz2l2l);
  cvggzz.push_back(Cggzz4l);
  cvggzz.push_back(Cggzz2l2l);
  cvH350.push_back(CggH350);
  cvH350.push_back(CVBF350);
  cvH126.push_back(CggH126);
  cvH126.push_back(CVBF126);

  MergedCollection *MCZZ = new MergedCollection("ZZ",cvZZ);
  MergedCollection *MCggZZ = new MergedCollection("ggZZ",cvggzz);
  MergedCollection *MCH350 = new MergedCollection("mH350",cvH350);
  MergedCollection *MCH126 = new MergedCollection("mH126",cvH126);

  // ------------ Legend ------------ //
  TLegend *leg = new TLegend(0.6,0.62,0.85,0.93);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TLatex *CP = new TLatex();
  CP->SetNDC(kTRUE);
  CP->SetTextSize(0.035);
  CP->SetTextAlign(31);
  CP->SetTextFont(42);
  CP->SetTextAlign(11);



  // ================================================================================================= //
  double plotYmax, plotYmax_zoom;

  // ------------ m4l --------------- //
  TF1 *redBgFunc_4l = new TF1("redBgFunc_4l","TMath::Landau(x,[0],[1],0)", xLow4l,xHigh4l);
  redBgFunc_4l->FixParameter(0,l_mu);
  redBgFunc_4l->FixParameter(1,l_sigma);

  TH1F* histm4l_h126_zoom = new TH1F("m4l_h126_zoom","Mass of 4 leptons; m_{4l} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4l_h350_zoom = new TH1F("m4l_h350_zoom","Mass of 4 leptons; m_{4l} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4l_ZZ_zoom = new TH1F("m4l_ZZ_zoom","Mass of 4 leptons; m_{4l} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  
  TH1F* histm4l_ZX_zoom = new TH1F("m4l_ZX_zoom","Mass of 4 leptons; m_{4l} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom); 
  helper->getHistFromTF1(redBgFunc_4l,histm4l_ZX_zoom,nBgEvents);

  TH1F* histm4l_data_zoom = new TH1F("m4l_data_zoom","Mass of 4 leptons; m_{4l} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  helper->fillHistFromVariable(Cdata,histm4l_data_zoom,"mass4l");
  TGraphAsymmErrors* histm4l_dataE_zoom;
  helper->getAsymErr(histm4l_data_zoom,histm4l_dataE_zoom);
    
  helper->fillHistFromVariable(MCZZ,histm4l_ZZ_zoom,"mass4l");
  helper->fillHistFromVariable(MCH350,histm4l_h350_zoom,"mass4l");
  helper->fillHistFromVariable(MCH126,histm4l_h126_zoom,"mass4l");
  histm4l_h126_zoom->Scale(scale126_4l/histm4l_h126_zoom->Integral());

  helper->setHistProperties(histm4l_ZZ_zoom,ZZBgColor,1001,helper->filled);
  helper->setHistProperties(histm4l_ZX_zoom,redBgColor,1001,helper->filled);
  helper->setHistProperties(histm4l_h126_zoom,h126Color,1001,helper->line);
  helper->setHistProperties(histm4l_h350_zoom,h350Color,1001,helper->line);
  helper->setHistProperties(histm4l_data_zoom,kBlack,1001,helper->markers,20);


  helper->addLegendEntry(leg,histm4l_data_zoom,"Data"    ,"P");
  helper->addLegendEntry(leg,histm4l_ZZ_zoom  ,"ZZ"       ,"F");
  helper->addLegendEntry(leg,histm4l_ZX_zoom ,"Z+X"       ,"F");
  helper->addLegendEntry(leg,histm4l_h126_zoom,"m_{H}=126","F");


  dv highestBin_m4lz;
  highestBin_m4lz.push_back(helper->findHighestBin(histm4l_ZX_zoom));
  highestBin_m4lz.push_back(helper->findHighestBin(histm4l_ZZ_zoom));
  highestBin_m4lz.push_back(helper->findHighestBin(histm4l_h126_zoom));
  highestBin_m4lz.push_back(helper->findHighestBin(histm4l_h350_zoom));
  highestBin_m4lz.push_back(helper->findHighestBin(histm4l_data_zoom));

  for(unsigned int k = 0; k < highestBin_m4lz.size(); k++)
    {
      if(percentMaxBin*highestBin_m4lz[k] > plotYmax_zoom) plotYmax_zoom = percentMaxBin*highestBin_m4lz[k];
    }


  THStack *stackM4l_zoom = new THStack();
  stackM4l_zoom->Add(histm4l_ZX_zoom);
  stackM4l_zoom->Add(histm4l_ZZ_zoom);
  stackM4l_zoom->Add(histm4l_h126_zoom);
  //stackM4l_zoom->Add(histm4l_h350_zoom);

  helper->drawPlot(stackM4l_zoom,histm4l_dataE_zoom,leg,CP,pngPlotDir+"/m4l_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".png",
		   xAxisLabel_m4l,yAxisLabel_zoom,sqrts,lumi,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);
  helper->drawPlot(stackM4l_zoom,histm4l_dataE_zoom,leg,CP,epsPlotDir+"/m4l_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".eps",
		   xAxisLabel_m4l,yAxisLabel_zoom,sqrts,lumi,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);
  plotYmax_zoom = 0;
  // ------------- m4e ------------ //
  TH1F* histm4e_h126_zoom = new TH1F("m4e_h126_zoom","Mass of 4 leptons; m_{4e} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4e_h350_zoom = new TH1F("m4e_h350_zoom","Mass of 4 leptons; m_{4e} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4e_ZZ_zoom = new TH1F("m4e_ZZ_zoom","Mass of 4 leptons; m_{4e} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  
  TH1F* histm4e_ZX_zoom = new TH1F("m4e_ZX_zoom","Mass of 4 leptons; m_{4e} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom); 
  helper->getHistFromTF1(redBgFunc_4l,histm4e_ZX_zoom,nBgEvents_4e);

  TH1F* histm4e_data_zoom = new TH1F("m4e_data_zoom","Mass of 4 leptons; m_{4e} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  helper->fillHistFromVariable(Cdata,histm4e_data_zoom,"mass4e");
  TGraphAsymmErrors* histm4e_dataE_zoom;
  helper->getAsymErr(histm4e_data_zoom,histm4e_dataE_zoom);
    
  helper->fillHistFromVariable(MCZZ,histm4e_ZZ_zoom,"mass4e");
  helper->fillHistFromVariable(MCH350,histm4e_h350_zoom,"mass4e");
  helper->fillHistFromVariable(MCH126,histm4e_h126_zoom,"mass4e");
  histm4e_h126_zoom->Scale(scale126_4e/histm4e_h126_zoom->Integral());

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

  for(unsigned int k = 0; k < highestBin_m4ez.size(); k++)
    {
      if(percentMaxBin*highestBin_m4ez[k] > plotYmax_zoom) plotYmax_zoom = percentMaxBin*highestBin_m4ez[k];
    }


  THStack *stackM4e_zoom = new THStack();
  stackM4e_zoom->Add(histm4e_ZX_zoom);
  stackM4e_zoom->Add(histm4e_ZZ_zoom);
  stackM4e_zoom->Add(histm4e_h126_zoom);
  //stackM4e_zoom->Add(histm4e_h350_zoom);

  helper->drawPlot(stackM4e_zoom,histm4e_dataE_zoom,leg,CP,pngPlotDir+"/m4e_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".png",
		   xAxisLabel_m4e,yAxisLabel_zoom,sqrts,lumi,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);
  helper->drawPlot(stackM4e_zoom,histm4e_dataE_zoom,leg,CP,epsPlotDir+"/m4e_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".eps",
		   xAxisLabel_m4e,yAxisLabel_zoom,sqrts,lumi,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);

  plotYmax_zoom = 0;
  // ------------ m4mu ------------- //
  TH1F* histm4mu_h126_zoom = new TH1F("m4mu_h126_zoom","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4mu_h350_zoom = new TH1F("m4mu_h350_zoom","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm4mu_ZZ_zoom = new TH1F("m4mu_ZZ_zoom","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  
  TH1F* histm4mu_ZX_zoom = new TH1F("m4mu_ZX_zoom","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom); 
  helper->getHistFromTF1(redBgFunc_4l,histm4mu_ZX_zoom,nBgEvents_4mu);

  TH1F* histm4mu_data_zoom = new TH1F("m4mu_data_zoom","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  helper->fillHistFromVariable(Cdata,histm4mu_data_zoom,"mass4mu");
  TGraphAsymmErrors* histm4mu_dataE_zoom;
  helper->getAsymErr(histm4mu_data_zoom,histm4mu_dataE_zoom);
    
  helper->fillHistFromVariable(MCZZ,histm4mu_ZZ_zoom,"mass4mu");
  helper->fillHistFromVariable(MCH350,histm4mu_h350_zoom,"mass4mu");
  helper->fillHistFromVariable(MCH126,histm4mu_h126_zoom,"mass4mu");
  histm4mu_h126_zoom->Scale(scale126_4mu/histm4mu_h126_zoom->Integral());
 

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

  for(unsigned int k = 0; k < highestBin_m4muz.size(); k++)
    {
      if(percentMaxBin*highestBin_m4muz[k] > plotYmax_zoom) plotYmax_zoom = percentMaxBin*highestBin_m4muz[k];
    }



  THStack *stackM4mu_zoom = new THStack();
  stackM4mu_zoom->Add(histm4mu_ZX_zoom);
  stackM4mu_zoom->Add(histm4mu_ZZ_zoom);
  stackM4mu_zoom->Add(histm4mu_h126_zoom);
  //stackM4mu_zoom->Add(histm4mu_h350_zoom);

  helper->drawPlot(stackM4mu_zoom,histm4mu_dataE_zoom,leg,CP,pngPlotDir+"/m4mu_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".png",
		   xAxisLabel_m4mu,yAxisLabel_zoom,sqrts,lumi,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);
  helper->drawPlot(stackM4mu_zoom,histm4mu_dataE_zoom,leg,CP,epsPlotDir+"/m4mu_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".eps",
		   xAxisLabel_m4mu,yAxisLabel_zoom,sqrts,lumi,plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);

  plotYmax_zoom = 0;
  // ------------ m2e2mu ------------- //
  TH1F* histm2e2mu_h126_zoom = new TH1F("m2e2mu_h126_zoom","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm2e2mu_h350_zoom = new TH1F("m2e2mu_h350_zoom","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  TH1F* histm2e2mu_ZZ_zoom = new TH1F("m2e2mu_ZZ_zoom","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  
  TH1F* histm2e2mu_ZX_zoom = new TH1F("m2e2mu_ZX_zoom","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom); 
  helper->getHistFromTF1(redBgFunc_4l,histm2e2mu_ZX_zoom,nBgEvents_2e2mu);

  TH1F* histm2e2mu_data_zoom = new TH1F("m2e2mu_data_zoom","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 3 GeV",nBins4l_zoom,xLow4l_zoom,xHigh4l_zoom);
  helper->fillHistFromVariable(Cdata,histm2e2mu_data_zoom,"mass2e2mu");
  TGraphAsymmErrors* histm2e2mu_dataE_zoom;
  helper->getAsymErr(histm2e2mu_data_zoom,histm2e2mu_dataE_zoom);
    
  helper->fillHistFromVariable(MCZZ,histm2e2mu_ZZ_zoom,"mass2e2mu");
  helper->fillHistFromVariable(MCH350,histm2e2mu_h350_zoom,"mass2e2mu");
  helper->fillHistFromVariable(MCH126,histm2e2mu_h126_zoom,"mass2e2mu");
  histm2e2mu_h126_zoom->Scale(scale126_2e2mu/histm2e2mu_h126_zoom->Integral());
 
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

  for(unsigned int k = 0; k < highestBin_m2e2muz.size(); k++)
    {
      if(percentMaxBin*highestBin_m2e2muz[k] > plotYmax_zoom) plotYmax_zoom = percentMaxBin*highestBin_m2e2muz[k];
    }




  THStack *stackM2e2mu_zoom = new THStack();
  stackM2e2mu_zoom->Add(histm2e2mu_ZX_zoom);
  stackM2e2mu_zoom->Add(histm2e2mu_ZZ_zoom);
  stackM2e2mu_zoom->Add(histm2e2mu_h126_zoom);
  //stackM2e2mu_zoom->Add(histm2e2mu_h350_zoom);

  helper->drawPlot(stackM2e2mu_zoom,histm2e2mu_dataE_zoom,leg,CP,pngPlotDir+"/m2e2mu_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".png",
		   xAxisLabel_m2e2mu,yAxisLabel_zoom, sqrts, lumi, plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);
  helper->drawPlot(stackM2e2mu_zoom,histm2e2mu_dataE_zoom,leg,CP,epsPlotDir+"/m2e2mu_"+Form("%.0f",plotXlow_zoom)+"_"+Form("%.0f",plotXhigh_zoom)+".eps",
		   xAxisLabel_m2e2mu,yAxisLabel_zoom, sqrts, lumi, plotXlow_zoom,plotXhigh_zoom,plotYmax_zoom);

  plotYmax_zoom = 0;

  // --------------- m4l ---------------- //
  TH1F* histm4l_h126 = new TH1F("m4l_h126","Mass of 4 leptons; m_{4l} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  TH1F* histm4l_h350 = new TH1F("m4l_h350","Mass of 4 leptons; m_{4l} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  TH1F* histm4l_ZZ = new TH1F("m4l_ZZ","Mass of 4 leptons; m_{4l} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  
  TH1F* histm4l_ZX = new TH1F("m4l_ZX","Mass of 4 leptons; m_{4l} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l); 
  helper->getHistFromTF1(redBgFunc_4l,histm4l_ZX,nBgEvents);

  TH1F* histm4l_data = new TH1F("m4l_data","Mass of 4 leptons; m_{4l} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  helper->fillHistFromVariable(Cdata,histm4l_data,"mass4l");
  TGraphAsymmErrors* histm4l_dataE;
  helper->getAsymErr(histm4l_data,histm4l_dataE);
    
  helper->fillHistFromVariable(MCZZ,histm4l_ZZ,"mass4l");
  helper->fillHistFromVariable(MCH350,histm4l_h350,"mass4l");
  helper->fillHistFromVariable(MCH126,histm4l_h126,"mass4l");
  histm4l_h126->Scale(scale126_4l/histm4l_h126->Integral());

  helper->setHistProperties(histm4l_ZZ,ZZBgColor,1001,helper->filled);
  helper->setHistProperties(histm4l_ZX,redBgColor,1001,helper->filled);
  helper->setHistProperties(histm4l_h126,h126Color,1001,helper->line);
  helper->setHistProperties(histm4l_h350,h350Color,1001,helper->line);
  helper->setHistProperties(histm4l_data,kBlack,1001,helper->markers,20);

  //helper->addLegendEntry(leg,histm4l_h350,"m_{H}=350","F");

  dv highestBin_m4l;
  highestBin_m4l.push_back(helper->findHighestBin(histm4l_ZX));
  highestBin_m4l.push_back(helper->findHighestBin(histm4l_ZZ));
  highestBin_m4l.push_back(helper->findHighestBin(histm4l_h126));
  highestBin_m4l.push_back(helper->findHighestBin(histm4l_h350));
  highestBin_m4l.push_back(helper->findHighestBin(histm4l_data));

  for(unsigned int k = 0; k < highestBin_m4l.size(); k++)
    {
      if(percentMaxBin*highestBin_m4l[k] > plotYmax) plotYmax = percentMaxBin*highestBin_m4l[k];
    }




  THStack *stackM4l = new THStack();
  stackM4l->Add(histm4l_ZX);
  stackM4l->Add(histm4l_ZZ);
  stackM4l->Add(histm4l_h126);
  //stackM4l->Add(histm4l_h350);


  helper->drawPlot(stackM4l,histm4l_dataE,leg,CP,pngPlotDir+"/m4l_"+Form("%.0f",plotXlow)+"_"+Form("%.0f",plotXhigh)+".png",
		   xAxisLabel_m4l,yAxisLabel,sqrts,lumi,plotXlow,plotXhigh,plotYmax);
  helper->drawPlot(stackM4l,histm4l_dataE,leg,CP,epsPlotDir+"/m4l_"+Form("%.0f",plotXlow)+"_"+Form("%.0f",plotXhigh)+".eps",
		   xAxisLabel_m4l,yAxisLabel,sqrts,lumi,plotXlow,plotXhigh,plotYmax);

  plotYmax = 0;
  // ------------- m4e ------------ //
  TH1F* histm4e_h126 = new TH1F("m4e_h126","Mass of 4 leptons; m_{4e} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  TH1F* histm4e_h350 = new TH1F("m4e_h350","Mass of 4 leptons; m_{4e} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  TH1F* histm4e_ZZ = new TH1F("m4e_ZZ","Mass of 4 leptons; m_{4e} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  
  TH1F* histm4e_ZX = new TH1F("m4e_ZX","Mass of 4 leptons; m_{4e} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l); 
  helper->getHistFromTF1(redBgFunc_4l,histm4e_ZX,nBgEvents_4e);

  TH1F* histm4e_data = new TH1F("m4e_data","Mass of 4 leptons; m_{4e} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  helper->fillHistFromVariable(Cdata,histm4e_data,"mass4e");
  TGraphAsymmErrors* histm4e_dataE;
  helper->getAsymErr(histm4e_data,histm4e_dataE);
    
  helper->fillHistFromVariable(MCZZ,histm4e_ZZ,"mass4e");
  helper->fillHistFromVariable(MCH350,histm4e_h350,"mass4e");
  helper->fillHistFromVariable(MCH126,histm4e_h126,"mass4e");
  histm4e_h126->Scale(scale126_4e/histm4e_h126->Integral());

  helper->setHistProperties(histm4e_ZZ,ZZBgColor,1001,helper->filled);
  helper->setHistProperties(histm4e_ZX,redBgColor,1001,helper->filled);
  helper->setHistProperties(histm4e_h126,h126Color,1001,helper->line);
  helper->setHistProperties(histm4e_h350,h350Color,1001,helper->line);
  helper->setHistProperties(histm4e_data,kBlack,1001,helper->markers,20);


  dv highestBin_m4e;
  highestBin_m4e.push_back(helper->findHighestBin(histm4e_ZX));
  highestBin_m4e.push_back(helper->findHighestBin(histm4e_ZZ));
  highestBin_m4e.push_back(helper->findHighestBin(histm4e_h126));
  highestBin_m4e.push_back(helper->findHighestBin(histm4e_h350));
  highestBin_m4e.push_back(helper->findHighestBin(histm4e_data));

  for(unsigned int k = 0; k < highestBin_m4e.size(); k++)
    {
      if(percentMaxBin*highestBin_m4e[k] > plotYmax) plotYmax = percentMaxBin*highestBin_m4e[k];
    }


  THStack *stackM4e = new THStack();
  stackM4e->Add(histm4e_ZX);
  stackM4e->Add(histm4e_ZZ);
  stackM4e->Add(histm4e_h126);
  //stackM4e->Add(histm4e_h350);

  plotYmax = 15;


  helper->drawPlot(stackM4e,histm4e_dataE,leg,CP,pngPlotDir+"/m4e_"+Form("%.0f",plotXlow)+"_"+Form("%.0f",plotXhigh)+".png",
		   xAxisLabel_m4e,yAxisLabel,sqrts,lumi,plotXlow,plotXhigh,plotYmax);
  helper->drawPlot(stackM4e,histm4e_dataE,leg,CP,epsPlotDir+"/m4e_"+Form("%.0f",plotXlow)+"_"+Form("%.0f",plotXhigh)+".eps",
		   xAxisLabel_m4e,yAxisLabel,sqrts,lumi,plotXlow,plotXhigh,plotYmax);

  plotYmax = 0;

  // ------------ m4mu ------------- //
  TH1F* histm4mu_h126 = new TH1F("m4mu_h126","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  TH1F* histm4mu_h350 = new TH1F("m4mu_h350","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  TH1F* histm4mu_ZZ = new TH1F("m4mu_ZZ","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  
  TH1F* histm4mu_ZX = new TH1F("m4mu_ZX","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l); 
  helper->getHistFromTF1(redBgFunc_4l,histm4mu_ZX,nBgEvents_4mu);

  TH1F* histm4mu_data = new TH1F("m4mu_data","Mass of 4 leptons; m_{4#mu} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  helper->fillHistFromVariable(Cdata,histm4mu_data,"mass4mu");
  TGraphAsymmErrors* histm4mu_dataE;
  helper->getAsymErr(histm4mu_data,histm4mu_dataE);
    
  helper->fillHistFromVariable(MCZZ,histm4mu_ZZ,"mass4mu");
  helper->fillHistFromVariable(MCH350,histm4mu_h350,"mass4mu");
  helper->fillHistFromVariable(MCH126,histm4mu_h126,"mass4mu");
  histm4mu_h126->Scale(scale126_4mu/histm4mu_h126->Integral());

  helper->setHistProperties(histm4mu_ZZ,ZZBgColor,1001,helper->filled);
  helper->setHistProperties(histm4mu_ZX,redBgColor,1001,helper->filled);
  helper->setHistProperties(histm4mu_h126,h126Color,1001,helper->line);
  helper->setHistProperties(histm4mu_h350,h350Color,1001,helper->line);
  helper->setHistProperties(histm4mu_data,kBlack,1001,helper->markers,20);



  dv highestBin_m4mu;
  highestBin_m4mu.push_back(helper->findHighestBin(histm4mu_ZX));
  highestBin_m4mu.push_back(helper->findHighestBin(histm4mu_ZZ));
  highestBin_m4mu.push_back(helper->findHighestBin(histm4mu_h126));
  highestBin_m4mu.push_back(helper->findHighestBin(histm4mu_h350));
  highestBin_m4mu.push_back(helper->findHighestBin(histm4mu_data));

  for(unsigned int k = 0; k < highestBin_m4mu.size(); k++)
    {
      if(percentMaxBin*highestBin_m4mu[k] > plotYmax) plotYmax = percentMaxBin*highestBin_m4mu[k];
    }




  THStack *stackM4mu = new THStack();
  stackM4mu->Add(histm4mu_ZX);
  stackM4mu->Add(histm4mu_ZZ);
  stackM4mu->Add(histm4mu_h126);
  //stackM4mu->Add(histm4mu_h350);

  plotYmax = 15;


  helper->drawPlot(stackM4mu,histm4mu_dataE,leg,CP,pngPlotDir+"/m4mu_"+Form("%.0f",plotXlow)+"_"+Form("%.0f",plotXhigh)+".png",
		   xAxisLabel_m4mu,yAxisLabel,sqrts,lumi,plotXlow,plotXhigh,plotYmax);
  helper->drawPlot(stackM4mu,histm4mu_dataE,leg,CP,epsPlotDir+"/m4mu_"+Form("%.0f",plotXlow)+"_"+Form("%.0f",plotXhigh)+".eps",
		   xAxisLabel_m4mu,yAxisLabel,sqrts,lumi,plotXlow,plotXhigh,plotYmax);

  plotYmax = 0;

  // ------------ m2e2mu ------------- //
  TH1F* histm2e2mu_h126 = new TH1F("m2e2mu_h126","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  TH1F* histm2e2mu_h350 = new TH1F("m2e2mu_h350","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  TH1F* histm2e2mu_ZZ = new TH1F("m2e2mu_ZZ","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  
  TH1F* histm2e2mu_ZX = new TH1F("m2e2mu_ZX","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l); 
  helper->getHistFromTF1(redBgFunc_4l,histm2e2mu_ZX,nBgEvents_2e2mu);

  TH1F* histm2e2mu_data = new TH1F("m2e2mu_data","Mass of 4 leptons; m_{2e2#mu} [GeV]; Events / 10 GeV",nBins4l,xLow4l,xHigh4l);
  helper->fillHistFromVariable(Cdata,histm2e2mu_data,"mass2e2mu");
  TGraphAsymmErrors* histm2e2mu_dataE;
  helper->getAsymErr(histm2e2mu_data,histm2e2mu_dataE);
    
  helper->fillHistFromVariable(MCZZ,histm2e2mu_ZZ,"mass2e2mu");
  helper->fillHistFromVariable(MCH350,histm2e2mu_h350,"mass2e2mu");
  helper->fillHistFromVariable(MCH126,histm2e2mu_h126,"mass2e2mu");
  histm2e2mu_h126->Scale(scale126_2e2mu/histm2e2mu_h126->Integral());

  helper->setHistProperties(histm2e2mu_ZZ,ZZBgColor,1001,helper->filled);
  helper->setHistProperties(histm2e2mu_ZX,redBgColor,1001,helper->filled);
  helper->setHistProperties(histm2e2mu_h126,h126Color,1001,helper->line);
  helper->setHistProperties(histm2e2mu_h350,h350Color,1001,helper->line);
  helper->setHistProperties(histm2e2mu_data,kBlack,1001,helper->markers,20);


  dv highestBin_m2e2mu;
  highestBin_m2e2mu.push_back(helper->findHighestBin(histm2e2mu_ZX));
  highestBin_m2e2mu.push_back(helper->findHighestBin(histm2e2mu_ZZ));
  highestBin_m2e2mu.push_back(helper->findHighestBin(histm2e2mu_h126));
  highestBin_m2e2mu.push_back(helper->findHighestBin(histm2e2mu_h350));
  highestBin_m2e2mu.push_back(helper->findHighestBin(histm2e2mu_data));

  for(unsigned int k = 0; k < highestBin_m2e2mu.size(); k++)
    {
      if(percentMaxBin*highestBin_m2e2mu[k] > plotYmax) plotYmax = percentMaxBin*highestBin_m2e2mu[k];
    }


  THStack *stackM2e2mu = new THStack();
  stackM2e2mu->Add(histm2e2mu_ZX);
  stackM2e2mu->Add(histm2e2mu_ZZ);
  stackM2e2mu->Add(histm2e2mu_h126);
  //stackM2e2mu->Add(histm2e2mu_h350);

  plotYmax = 15;


  helper->drawPlot(stackM2e2mu,histm2e2mu_dataE,leg,CP,pngPlotDir+"/m2e2mu_"+Form("%.0f",plotXlow)+"_"+Form("%.0f",plotXhigh)+".png",
		   xAxisLabel_m2e2mu,yAxisLabel,sqrts,lumi,plotXlow,plotXhigh,plotYmax);
  helper->drawPlot(stackM2e2mu,histm2e2mu_dataE,leg,CP,epsPlotDir+"/m2e2mu_"+Form("%.0f",plotXlow)+"_"+Form("%.0f",plotXhigh)+".eps",
		   xAxisLabel_m2e2mu,yAxisLabel,sqrts,lumi,plotXlow,plotXhigh,plotYmax);



  // ----------------- Print Data ------------- //
  if(printData)
    {
      treeName = "Ana/passedEvents";
      helper->printData(fileDir+"/Data.root",treeName,       txtDir+"/SChannelEvents_7TeV.txt",      40.0, 4.0, ptMuCut,ptElCut, 80.0,  100.0);
      helper->printData(fileDir+"/Data.root",treeName,      txtDir+"/PassedEvents_7TeV.txt",         40.0, 12.0,ptMuCut,ptElCut, 100.0, 1000.0);
      helper->printDataShort(fileDir+"/Data.root",treeName, txtDir+"/PassedEvents_short_7TeV.txt",   40.0, 12.0,ptMuCut,ptElCut, 100.0, 1000.0);
      helper->printDataTex(fileDir+"/Data.root",treeName,       txtDir+"/PassedEvents_7TeV.tex",     40.0, 12.0,ptMuCut,ptElCut, 100.0, 1000.0);
      
      //helper->printDataErr(fileDir+"/Data_new.root",treeName,       txtDir+"/EventErrors.txt",       40.0, 12.0, 70.0, 170.0);
      //helper->printFourMomenta(fileDir+"/Data_new.root",treeName,       txtDir+"/EventFourMomenta_8TeV.txt",  40.0, 12.0, 115.0, 135.0);
      //helper->printData(fileDir+"/Data_new.root","AnaAfterHlt/passedEvents",      txtDir+"/PassedEvents_forTiz.txt",       40.0, 12.0, 121.5, 130.5);
    }

  // ---------------- Print Yields ------------- //
  ofstream out;
  string yieldFile_100_800 = txtDir+"/HistoYields_100_1000.txt";
  out.open(yieldFile_100_800.c_str());
  helper->printCollectionYield(MCZZ_noTaus, out, 40., 12., 100., 800.);
  helper->printCollectionYield(MCZZ_taus, out, 40., 12., 100., 800.);
  helper->printCollectionYield(MCggZZ, out, 40., 12., 100., 800.);
  helper->printCollectionYield(MCH126, out, 40., 12., 100., 800.);
  helper->printCollectionYield(MCH350, out, 40., 12., 100., 800.);
  out << "----------------" << endl;
  out << "Z4L Fit"          << endl;
  helper->printCollectionYield(MCZZ_taus, out, 40., 4., 100., 110.);
  helper->printCollectionYield(MCZZ_taus, out, 40., 4., 70., 100.);
  out << "----------------" << endl;
  out << "Z4L Fit MZ2 > 12"          << endl;
  helper->printCollectionYield(MCZZ_taus, out, 40., 12., 100., 110.);
  helper->printCollectionYield(MCZZ_taus, out, 40., 12., 70., 100.);

  out.close();


  return 0;

}

