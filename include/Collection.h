#ifndef COLLECTION_H
#define COLLECTION_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <map>

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

//#include "Hists.h"

using namespace std;

typedef vector<double> dv;
typedef vector<TH1F*> hv;
typedef vector<THStack*> stackv;
typedef vector<ULong64_t> uv;
typedef vector<TString> tsv;
typedef vector<string> stringv;


class Collection
{

 public:
  Collection(std::string collectionName, TString fileName, TString treeName, double lumi,bool isData, double mZ1LowCut, double mZ2LowCut, double ptMuCut, double ptElCut, double m4lCut);
  ~Collection();

  dv weight;
  double nEvents;
  ULong64_t Run, Event, LumiSect;
  double mass4l, pT4l, massZ1, massZ2;
  double massError, mass4mu, mass2e2mu, mass4e;
  double FSRPhot1_Pt,FSRPhot2_Pt,FSRPhot1_eta,FSRPhot2_eta,FSRPhot1_phi,FSRPhot2_phi;
  double scaleWeight, eventWeight, dataMCweight;
  int idL1, idL2, idL3, idL4;
  double ptL1, ptL2, ptL3, ptL4;
  bool passedFullSelection;

  void fillVariables(TString fileName,TString treeName);
  std::string getName();

  pair<TString,dv> mass4lPair;
  pair<TString,dv> mass4muPair;
  pair<TString,dv> mass4ePair;
  pair<TString,dv> mass2e2muPair;
  pair<TString,dv> massZ1Pair;
  pair<TString,dv> massZ2Pair;
  pair<TString,dv> massErrPair;

 private:

  double globalMZ1Low;
  double globalMZ2Low;
  double globalM4lCut;
  double globalPtMuCut;
  double globalPtElCut;

  TString globalFileName;
  TString globalTreeName;
  double globalLumi;
  bool globalIsData;
  std::string globalName;

};


#endif


#ifndef COLLECTION_CC
#define COLLECTION_CC

Collection::Collection(std::string collectionName, TString fileName, TString treeName, double lumi,bool isData=false, double mZ1LowCut=40, double mZ2LowCut=4, double ptMuCut=5, double ptElCut=7, double m4lCut=0)
{

  globalName = collectionName;
  globalFileName = fileName;
  globalTreeName = treeName;
  globalLumi = lumi;
  globalIsData = isData;
  globalMZ1Low = mZ1LowCut;
  globalMZ2Low = mZ2LowCut;
  globalM4lCut = m4lCut;
  globalPtMuCut = ptMuCut;
  globalPtElCut = ptElCut;
  fillVariables(globalFileName,globalTreeName);

  mass4lPair.first = "mass4l";
  mass4ePair.first = "mass4e";
  mass4muPair.first = "mass4mu";
  mass2e2muPair.first = "mass2e2mu";
  massZ1Pair.first = "massZ1";
  massZ2Pair.first = "massZ2";
  massErrPair.first = "massErr";

  cout << globalName << "   " << globalIsData << endl;

}



Collection::~Collection()
{
  
  // --- do nothing here
  
}



void Collection::fillVariables(TString fileName,TString treeName)
{
  using namespace std;

  TFile *file = new TFile(fileName,"READ");
  if(!file){
    cout << "Cannot find file " << fileName << endl;
    exit(1);
  }
  TTree *tree = (TTree*)file->Get(treeName);
  if(!tree){
    cout << "Cannot find tree " << treeName << " in file " << fileName << endl;
    exit(1);
  }

  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("pTL1",&ptL1);
  tree->SetBranchAddress("pTL2",&ptL2);
  tree->SetBranchAddress("pTL3",&ptL3);
  tree->SetBranchAddress("pTL4",&ptL4);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("mass4e",&mass4e);
  tree->SetBranchAddress("mass4mu",&mass4mu);
  tree->SetBranchAddress("mass2e2mu",&mass2e2mu);
  tree->SetBranchAddress("scaleWeight",&scaleWeight);
  tree->SetBranchAddress("eventWeight",&eventWeight);
  tree->SetBranchAddress("dataMC_weight",&dataMCweight);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);
  tree->SetBranchAddress("massErrorUFCorr",&massError);
  tree->SetBranchAddress("FSRPhot1_Pt",&FSRPhot1_Pt);
  tree->SetBranchAddress("FSRPhot2_Pt",&FSRPhot2_Pt);
  tree->SetBranchAddress("FSRPhot1_eta",&FSRPhot1_eta);
  tree->SetBranchAddress("FSRPhot2_eta",&FSRPhot2_eta);
  tree->SetBranchAddress("FSRPhot1_phi",&FSRPhot1_phi);
  tree->SetBranchAddress("FSRPhot2_phi",&FSRPhot2_phi);
  tree->SetBranchAddress("passedFullSelection",&passedFullSelection);


  TH1F *hist = (TH1F*)file->Get("AnaAfterHlt/nEvents");
  TH1F *hist2;
  if(!hist){
    hist2 = (TH1F*)file->Get("Ana/nEvents");
    if(!hist2){
      cout << "Cannot find AnaAfterHlt/nEvents!" << endl;
      exit(1);
    }
  }

  double normEvents;
  if(!hist) normEvents = hist2->GetBinContent(1);
  else normEvents = hist->GetBinContent(1);

  nEvents = normEvents;
  
  cout<<"normEvents: "<<normEvents<<endl;
  
  for (int i = 0; i < tree->GetEntries(); i++)
    {
      tree->GetEntry(i);

      /*
      if( massZ1 < globalMZ1Low ) continue;
      if( massZ2 < globalMZ2Low ) continue;
      if( mass4l < globalM4lCut ) continue;

      if( abs(idL1) == 13 && ptL1 < globalPtMuCut ) continue;
      if( abs(idL2) == 13 && ptL2 < globalPtMuCut ) continue;
      if( abs(idL3) == 13 && ptL3 < globalPtMuCut ) continue;
      if( abs(idL4) == 13 && ptL4 < globalPtMuCut ) continue;

      if( abs(idL1) == 11 && ptL1 < globalPtElCut ) continue;
      if( abs(idL2) == 11 && ptL2 < globalPtElCut ) continue;
      if( abs(idL3) == 11 && ptL3 < globalPtElCut ) continue;
      if( abs(idL4) == 11 && ptL4 < globalPtElCut ) continue;
      */

      if (!passedFullSelection) continue;

      mass4lPair.second.push_back(mass4l);
      mass4muPair.second.push_back(mass4mu);
      mass4ePair.second.push_back(mass4e);
      mass2e2muPair.second.push_back(mass2e2mu);
      massZ1Pair.second.push_back(massZ1);
      massZ2Pair.second.push_back(massZ2);
      massErrPair.second.push_back(massError);
      if(globalIsData) weight.push_back(1);
      else weight.push_back(scaleWeight*eventWeight*dataMCweight*(globalLumi*1000)/normEvents);
      //else weight.push_back(scaleWeight*(globalLumi*1000)/normEvents);
      
    }
 
  file->Close();
    
}


std::string Collection::getName()
{
  return globalName;
}

#endif
