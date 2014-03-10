#ifndef MERGEDCOLLECTION_H
#define MERGEDCOLLECTION_H

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

#include "Collection.h"


using namespace std;

typedef vector<double> dv;
typedef vector<TH1F*> hv;
typedef vector<THStack*> stackv;
typedef vector<ULong64_t> uv;
typedef vector<TString> tsv;
typedef vector<string> stringv;


class MergedCollection
{

 public:
  MergedCollection(std::string collectionName, Collection *myCol1, Collection *myCol2);
  MergedCollection(std::string collectionName, std::vector<Collection*> myColV);
  ~MergedCollection();

  std::string getName();

  dv weight;

  pair<TString,dv> mass4lPair;
  pair<TString,dv> mass4muPair;
  pair<TString,dv> mass4ePair;
  pair<TString,dv> mass2e2muPair;
  pair<TString,dv> massZ1Pair;
  pair<TString,dv> massZ2Pair;
  pair<TString,dv> massErrPair;

 private:
  
  void MergeVariables(Collection *myCol1, Collection *myCol2);
  void MergeVariablesVec(std::vector<Collection*> myColV);
  std::string globalName;


};


#endif


#ifndef MERGEDCOLLECTION_CC
#define MERGEDCOLLECTION_CC

MergedCollection::MergedCollection(std::string collectionName,Collection *myCol1, Collection *myCol2)
{

  mass4lPair.first = "mass4l";
  mass4ePair.first = "mass4e";
  mass4muPair.first = "mass4mu";
  mass2e2muPair.first = "mass2e2mu";
  massZ1Pair.first = "massZ1";
  massZ2Pair.first = "massZ2";
  massErrPair.first = "massErr";
  MergeVariables(myCol1,myCol2);
  globalName = collectionName;

}

MergedCollection::MergedCollection(std::string collectionName, std::vector<Collection*> myColV)
{

  mass4lPair.first = "mass4l";
  mass4ePair.first = "mass4e";
  mass4muPair.first = "mass4mu";
  mass2e2muPair.first = "mass2e2mu";
  massZ1Pair.first = "massZ1";
  massZ2Pair.first = "massZ2";
  massErrPair.first = "massErr";
  MergeVariablesVec(myColV);
  globalName = collectionName;

}


MergedCollection::~MergedCollection()
{
  
  // --- do nothing here
  
}


void MergedCollection::MergeVariables(Collection *myCol1, Collection *myCol2)
{
  
  for (unsigned int i = 0; i < myCol1->mass4lPair.second.size(); i++)
    {
      mass4lPair.second.push_back(myCol1->mass4lPair.second[i]);
      mass4muPair.second.push_back(myCol1->mass4muPair.second[i]);
      mass4ePair.second.push_back(myCol1->mass4ePair.second[i]);
      mass2e2muPair.second.push_back(myCol1->mass2e2muPair.second[i]);
      massZ1Pair.second.push_back(myCol1->massZ1Pair.second[i]);
      massZ2Pair.second.push_back(myCol1->massZ2Pair.second[i]);
      massErrPair.second.push_back(myCol1->massErrPair.second[i]);
      weight.push_back(myCol1->weight[i]);
    }

  for (unsigned int i = 0; i < myCol2->mass4lPair.second.size(); i++)
    {
      mass4lPair.second.push_back(myCol2->mass4lPair.second[i]);
      mass4muPair.second.push_back(myCol2->mass4muPair.second[i]);
      mass4ePair.second.push_back(myCol2->mass4ePair.second[i]);
      mass2e2muPair.second.push_back(myCol2->mass2e2muPair.second[i]);
      massZ1Pair.second.push_back(myCol2->massZ1Pair.second[i]);
      massZ2Pair.second.push_back(myCol2->massZ2Pair.second[i]);
      massErrPair.second.push_back(myCol2->massErrPair.second[i]);
      weight.push_back(myCol2->weight[i]);
    }

 
}

void MergedCollection::MergeVariablesVec(std::vector<Collection*> myColV)
{
  for(unsigned int i = 0; i < myColV.size(); i++)
    {
      for (unsigned int j = 0; j < myColV[i]->mass4lPair.second.size(); j++)
	{
	  mass4lPair.second.push_back(myColV[i]->mass4lPair.second[j]);
	  mass4muPair.second.push_back(myColV[i]->mass4muPair.second[j]);
	  mass4ePair.second.push_back(myColV[i]->mass4ePair.second[j]);
	  mass2e2muPair.second.push_back(myColV[i]->mass2e2muPair.second[j]);
	  massZ1Pair.second.push_back(myColV[i]->massZ1Pair.second[j]);
	  massZ2Pair.second.push_back(myColV[i]->massZ2Pair.second[j]);
	  massErrPair.second.push_back(myColV[i]->massErrPair.second[j]);
	  weight.push_back(myColV[i]->weight[j]);
	}
    }
 
}

std::string MergedCollection::getName()
{
  return globalName;

}


#endif
