#ifndef PLOTHELPER_H
#define PLOTHELPER_H

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
#include "TRandom.h"

#include "Collection.h"
#include "MergedCollection.h"

typedef std::vector<TH1F*> hv;
typedef std::vector<THStack*> sv;
typedef std::vector<ULong64_t> uv;

class PlotHelper
{


 public:
  PlotHelper();
  ~PlotHelper();
  enum DrawStyle{line,filled,markers};


  int fillHistFromVariable(Collection *myCollection,TH1F* &hist,TString myVariable, double newWeight);
  int fillHistFromVariable(MergedCollection *myCollection,TH1F* &hist,TString myVariable, double newWeight);
  void setHistProperties(TH1F* &hist,Color_t color,Style_t fillStyle,DrawStyle draw,Style_t markerStyle=20);
  void setStackStyle(THStack* hist,TString xTitle, TString yTitle);
  void clearHist(TH1F* &hist);
  std::vector<TGraphAsymmErrors*> getAsymErr(std::vector<TH1F*> h);
  void getHistFromTF1(TF1* &funcFit, TH1F* &histFit, double norm);
  void addLegendEntry(TLegend* &leg, TH1F* &hist, TString label, TString style);
  void getAsymErr(TH1F* h,TGraphAsymmErrors* &gr);

  void printDataList(TString file, TString treeName, std::string outFile, double massZ1Cut, double massZ2Cut, double m4lLowCut, double m4lHighCut);
  void printFourMomenta(TString file, TString treeName, std::string outFile, double massZ1Cut, double massZ2Cut, double m4lLowCut, double m4lHighCut);
  void printData(TString file, TString treeName, std::string outFile, double massZ1Cut, double massZ2Cut, double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut);
  void printDataSync(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut);
  void printDataTex(TString file, TString treeName, std::string outFile, double massZ1Cut, double massZ2Cut, double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut);
  void printVsRun(TString file, TString treeName, std::string outFile, double massZ1Cut, double massZ2Cut, double m4lLowCut, double m4lHighCut);
  void printDataShort(TString file, TString treeName, std::string outFile, double massZ1Cut, double massZ2Cut, double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut);
  void printDataErr(TString file, TString treeName, std::string outFile, double massZ1Cut, double massZ2Cut, double m4lLowCut, double m4lHighCut);
  
  void printStepPlot(std::vector<TH1F*> step, std::vector<std::string> sampleName, std::string fileName);
  void printHistoYield(std::string name,std::vector<TH1F*> hist_4l,std::vector<TH1F*> hist_4mu,std::vector<TH1F*> hist_4e,std::vector<TH1F*> hist_2e2mu, std::vector<TString> histName, double xlow, double xhigh);

  void printCollectionYield(Collection *myCol, ofstream &out, double mZ1Low, double mZ2Low, double m4lLow, double m4lHigh);
  void printCollectionYield(MergedCollection *myCol, ofstream &out, double mZ1Low, double mZ2Low, double m4lLow, double m4lHigh);

  std::vector<double> getCollectionYield(Collection *myCol,  double mZ1Low, double mZ2Low, double m4lLow, double m4lHigh);
  std::vector<double> getCollectionYield(MergedCollection *myCol,  double mZ1Low, double mZ2Low, double m4lLow, double m4lHigh);
  
  void rescale(TH1F *hist);
  double correctionFactor(double mH, int channel);
  void rescale(TH1F *hist, int channel, double n_tot, double n_4mu, double n_4e, double n_2e2mu);

  void drawPlot(THStack *stack, TGraphAsymmErrors *gr, TLegend *leg, TLatex *latex,TString saveName,TString xTitle, TString yTitle, double sqrts, 
		double lumi,double xLow=999, double xHigh=999, double yMax=999);
  void draw7p8Plot(THStack *stack, TGraphAsymmErrors *gr, TLegend *leg, TLatex *latex,TString saveName,TString xTitle, TString yTitle, 
		   double lumi_7, double lumi_8, double xLow, double xHigh, double yMax);

  void draw7p8Plot(THStack *stack,TH1F* gr, TLegend *leg, TLatex *latex,TString saveName,TString xTitle, TString yTitle,
		   double lumi_7, double lumi_8, double xLow, double xHigh, double yMax);

  void makeArrayFromVector(std::vector<double> myVec, double myArray[]);

  void PaintOverflow(TH1F* h,TH1F* &htmp,Color_t color,Style_t fillStyle,DrawStyle draw, Style_t markerStyle);
  TH1F* setOverflowBin(TH1F* histo, float xmax);
  double findHighestBin(TH1F* histo);


 private:
  
  //runvecs
  uv runVec, lumiVec, eventVec;
  uv runVec2, lumiVec2, eventVec2;
  uv runVec3, lumiVec3, eventVec3;
  uv runVec4, lumiVec4, eventVec4;
  uv runVec5, lumiVec5, eventVec5;
  uv runVec6, lumiVec6, eventVec6;
  uv runVec7, lumiVec7, eventVec7;
  uv runVec8, lumiVec8, eventVec8;

  double canvasX, canvasY;

};


#endif

#ifndef HISTMAKER_CC
#define HISTMAKER_CC


PlotHelper::PlotHelper()
{
  
  canvasX = 600;
  canvasY = 600;

}


PlotHelper::~PlotHelper(){

  // do nothing 

}


int
PlotHelper::fillHistFromVariable(Collection *myCollection,TH1F* &hist,TString myVariable, double newWeight=1.0)
{

  hist->Sumw2();
  using namespace std;
  vector<double> myVec;
  vector<double> weight = myCollection->weight;
  for (int i=0; i<(int)weight.size(); i++) weight.at(i) *= newWeight;  

  double total=0.0;

  if(myVariable == myCollection->mass4lPair.first)
    {
      myVec = myCollection->mass4lPair.second;
    }
  if(myVariable == myCollection->mass4ePair.first)
    {
      myVec = myCollection->mass4ePair.second;
    }
  if(myVariable == myCollection->mass4muPair.first)
    {
      myVec = myCollection->mass4muPair.second;
    }
  if(myVariable == myCollection->mass2e2muPair.first) 
    {
      myVec = myCollection->mass2e2muPair.second;
    }

  int Entries = 0;
  double min = hist->GetXaxis()->GetXmin();
  double max = hist->GetXaxis()->GetXmax();  

  for( unsigned int i = 0; i < myVec.size(); i++ )
    {
      hist->Fill(myVec[i],(double)weight[i]);
      total += (double)weight[i];
      if (myVec[i]>min&&myVec[i]<max) Entries++;
    }

  std::cout<<hist->GetName()<<" total: "<<total<<std::endl;

  return Entries;

}

int
PlotHelper::fillHistFromVariable(MergedCollection *myCollection,TH1F* &hist,TString myVariable, double newWeight=1.0)
{
  hist->Sumw2();

  using namespace std;
  vector<double> myVec;
  vector<double> weight = myCollection->weight;
  for (int i=0; i<(int)weight.size(); i++) weight.at(i) *= newWeight;

  if(myVariable == myCollection->mass4lPair.first)
    {
      myVec = myCollection->mass4lPair.second;
    }
  if(myVariable == myCollection->mass4ePair.first)
    {
      myVec = myCollection->mass4ePair.second;
    }
  if(myVariable == myCollection->mass4muPair.first)
    {
      myVec = myCollection->mass4muPair.second;
    }
  if(myVariable == myCollection->mass2e2muPair.first) 
    {
      myVec = myCollection->mass2e2muPair.second;
    }

  
  int Entries=0;
  double min = hist->GetXaxis()->GetXmin();
  double max = hist->GetXaxis()->GetXmax();
  for( unsigned int i = 0; i < myVec.size(); i++ )
    {
      hist->Fill(myVec[i],(double)weight[i]);
      if (myVec[i]>min&&myVec[i]<max) Entries++;
    }

  return Entries;

}



void 
PlotHelper::clearHist(TH1F* &hist)
{
  for(int i = 0; i <= hist->GetNbinsX(); i++)
    {
      hist->SetBinContent(i,0);
    }
}


std::vector<TGraphAsymmErrors*> PlotHelper::getAsymErr(std::vector<TH1F*> h) 
{
  
  using namespace std;
  vector<TGraphAsymmErrors*> gr1;
  
  for( unsigned int j = 0; j < h.size(); j++)
    {
      const Int_t n = h[j]->GetNbinsX();
      Double_t x[n],y[n],exl[n],eyl[n],exh[n],eyh[n];

      double q = (1-0.6827)/2.;

      for (Int_t i=0;i<n;i++)
        {
          x[i] = h[j]->GetXaxis()->GetBinCenter(i+1);
          double N = h[j]->GetBinContent(i+1);

          exl[i] = 0; exh[i] = 0;

          if( N == 0 ){eyl[i] = 0; eyh[i] = 0; y[i] = -1;continue;}

          y[i] = N;
          eyl[i] = (N==0)?0:(N-ROOT::Math::chisquared_quantile_c(1-q,2*N)/2.);
          eyh[i] = ROOT::Math::chisquared_quantile_c(q,2*(N+1))/2.-N;

        }

      TGraphAsymmErrors *gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
      gr->SetMarkerStyle(20);
      gr1.push_back(gr);

    }

  return gr1;
}



void 
PlotHelper::getAsymErr(TH1F* h,TGraphAsymmErrors* &gr)
{
  
  using namespace std;
  
  const Int_t n = h->GetNbinsX();
  Double_t x[n],y[n],exl[n],eyl[n],exh[n],eyh[n];
  
  double q = (1-0.6827)/2.;
  
  for (Int_t i=0;i<n;i++)
    {
      x[i] = h->GetXaxis()->GetBinCenter(i+1);
      double N = h->GetBinContent(i+1);
      
      exl[i] = 0; exh[i] = 0;
      
      if( N == 0 ){eyl[i] = 0; eyh[i] = 0; y[i] = -1;continue;}
      
      y[i] = N;
      eyl[i] = (N==0)?0:(N-ROOT::Math::chisquared_quantile_c(1-q,2*N)/2.);
      eyh[i] = ROOT::Math::chisquared_quantile_c(q,2*(N+1))/2.-N;
      
    }
  
  gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlack);
  gr->SetMarkerSize(1.1);

}




void 
PlotHelper::getHistFromTF1(TF1* &funcFit, TH1F* &hist, double norm) 
{
  
  double xwidth = hist->GetXaxis()->GetBinWidth(2);
  double rangeX1 = hist->GetXaxis()->GetXmin();
  for (Int_t i = 1; i <= hist->GetNbinsX(); i++) 
    {
      if (norm>0.0) hist->SetBinContent(i,funcFit->Integral(rangeX1+(i-1)*xwidth,rangeX1+(i)*xwidth)/xwidth);
      else hist->SetBinContent(i,funcFit->Integral(rangeX1+(i-1)*xwidth,rangeX1+(i)*xwidth));
    
    }

  if (norm>0.0) hist->Scale(norm/hist->Integral());

}



void 
PlotHelper::setHistProperties(TH1F* &hist,Color_t color,Style_t fillStyle,DrawStyle draw,Style_t markerStyle)
{
  double lineWidth = 2;

  if(draw == markers)
    {
      hist->SetMarkerStyle(20);
      hist->SetMarkerColor(color);
      hist->SetMarkerSize(1.1);
    }
  else if(draw == line)
    {
      hist->SetLineColor(color);
      hist->SetLineWidth(lineWidth);
    }
  else if(draw == filled){
    hist->SetFillStyle(fillStyle);
    hist->SetFillColor(color);
  }

}

void 
PlotHelper::setStackStyle(THStack* hist,TString xTitle, TString yTitle)
{
  double ytitleOffset = 1.36;
  double xtitleOffset = 1.18;
  double labelSize = 0.05;
  double titleSize = 0.05;

  hist->GetXaxis()->SetTitleOffset(xtitleOffset);
  hist->GetYaxis()->SetTitleOffset(ytitleOffset);
  hist->GetXaxis()->SetLabelSize(labelSize);
  hist->GetYaxis()->SetLabelSize(labelSize);
  hist->GetXaxis()->SetTitleSize(titleSize);
  hist->GetYaxis()->SetTitleSize(titleSize);
  hist->GetXaxis()->SetTitle(xTitle);
  hist->GetYaxis()->SetTitle(yTitle);


}




void 
PlotHelper::addLegendEntry(TLegend* &leg, TH1F* &hist, TString label, TString style)
{

  leg->AddEntry(hist,label,style);

}


void 
PlotHelper::printDataList(TString file,TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;

  ULong64_t Run, Event, LumiSect;
  double mass4l, pT4l, massZ1, massZ2;
  double eta4l, phi4l;
  int idL1, idL2, idL3, idL4;
  double EL1, EL2, EL3, EL4;
  double EZ1, EZ2;
  double pTL1, pTL2, pTL3, pTL4;
  double pXL1, pXL2, pXL3, pXL4;
  double pYL1, pYL2, pYL3, pYL4;
  double pZL1, pZL2, pZL3, pZL4;
  double pTZ1, pTZ2;
  double pXZ1, pXZ2;
  double pYZ1, pYZ2;
  double pZZ1, pZZ2;
  double etaL1, etaL2, etaL3, etaL4;
  double phiL1, phiL2, phiL3, phiL4;
  double isoNHL1, isoNHL2, isoNHL3, isoNHL4;
  double isoCHL1, isoCHL2, isoCHL3, isoCHL4;
  double isoPhotL1, isoPhotL2, isoPhotL3, isoPhotL4;
  double RelIsoL1, RelIsoL2, RelIsoL3, RelIsoL4;
  double muRhoCor, elRhoCor;
  double worstIso, worstSIP;
  double cosTheta1, cosTheta2, Phi, cosThetaStar, phiStar1, phiStar2, Phi1, Phi2;
  int nJets, nVtx, nPhotons;
  double metVal, rawRho;
  double massError;

  TFile *f = new TFile(file,"READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("pT4l",&pT4l);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);
  tree->SetBranchAddress("eta4l",&eta4l);
  tree->SetBranchAddress("phi4l",&phi4l);
  tree->SetBranchAddress("EL1",&EL1);
  tree->SetBranchAddress("EL2",&EL2);
  tree->SetBranchAddress("EL3",&EL3);
  tree->SetBranchAddress("EL4",&EL4);
  tree->SetBranchAddress("EZ1",&EZ1);
  tree->SetBranchAddress("EZ2",&EZ2);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("pTL1",&pTL1);
  tree->SetBranchAddress("pTL2",&pTL2);
  tree->SetBranchAddress("pTL3",&pTL3);
  tree->SetBranchAddress("pTL4",&pTL4);
  tree->SetBranchAddress("pXL1",&pXL1);
  tree->SetBranchAddress("pXL2",&pXL2);
  tree->SetBranchAddress("pXL3",&pXL3);
  tree->SetBranchAddress("pXL4",&pXL4);
  tree->SetBranchAddress("pYL1",&pYL1);
  tree->SetBranchAddress("pYL2",&pYL2);
  tree->SetBranchAddress("pYL3",&pYL3);
  tree->SetBranchAddress("pYL4",&pYL4);
  tree->SetBranchAddress("pZL1",&pZL1);
  tree->SetBranchAddress("pZL2",&pZL2);
  tree->SetBranchAddress("pZL3",&pZL3);
  tree->SetBranchAddress("pZL4",&pZL4);
  tree->SetBranchAddress("pTZ1",&pTZ1);
  tree->SetBranchAddress("pTZ2",&pTZ2);
  tree->SetBranchAddress("pXZ1",&pXZ1);
  tree->SetBranchAddress("pXZ2",&pXZ2);
  tree->SetBranchAddress("pYZ1",&pYZ1);
  tree->SetBranchAddress("pYZ2",&pYZ2);
  tree->SetBranchAddress("pZZ1",&pZZ1);
  tree->SetBranchAddress("pZZ2",&pZZ2);
  tree->SetBranchAddress("etaL1",&etaL1);
  tree->SetBranchAddress("etaL2",&etaL2);
  tree->SetBranchAddress("etaL3",&etaL3);
  tree->SetBranchAddress("etaL4",&etaL4);
  tree->SetBranchAddress("phiL1",&phiL1);
  tree->SetBranchAddress("phiL2",&phiL2);
  tree->SetBranchAddress("phiL3",&phiL3);
  tree->SetBranchAddress("phiL4",&phiL4);
  tree->SetBranchAddress("isoNHL1",&isoNHL1);
  tree->SetBranchAddress("isoNHL2",&isoNHL2);
  tree->SetBranchAddress("isoNHL3",&isoNHL3);
  tree->SetBranchAddress("isoNHL4",&isoNHL4);
  tree->SetBranchAddress("isoCHL1",&isoCHL1);
  tree->SetBranchAddress("isoCHL2",&isoCHL2);
  tree->SetBranchAddress("isoCHL3",&isoCHL3);
  tree->SetBranchAddress("isoCHL4",&isoCHL4);
  tree->SetBranchAddress("isoPhotL1",&isoPhotL1);
  tree->SetBranchAddress("isoPhotL2",&isoPhotL2);
  tree->SetBranchAddress("isoPhotL3",&isoPhotL3);
  tree->SetBranchAddress("isoPhotL4",&isoPhotL4);
  tree->SetBranchAddress("RelIsoL1",&RelIsoL1);
  tree->SetBranchAddress("RelIsoL2",&RelIsoL2);
  tree->SetBranchAddress("RelIsoL3",&RelIsoL3);
  tree->SetBranchAddress("RelIsoL4",&RelIsoL4);
  tree->SetBranchAddress("muRhoCor",&muRhoCor);
  tree->SetBranchAddress("elRhoCor",&elRhoCor);
  tree->SetBranchAddress("worstIso",&worstIso);
  tree->SetBranchAddress("worstSIP",&worstSIP);
  tree->SetBranchAddress("cosTheta1",&cosTheta1);
  tree->SetBranchAddress("cosTheta2",&cosTheta2);
  tree->SetBranchAddress("Phi",&Phi);
  tree->SetBranchAddress("cosThetaStar",&cosThetaStar);
  tree->SetBranchAddress("phiStar1",&phiStar1);
  tree->SetBranchAddress("phiStar2",&phiStar2);
  tree->SetBranchAddress("Phi1",&Phi1);
  tree->SetBranchAddress("Phi2",&Phi2);
  tree->SetBranchAddress("nJets",&nJets);
  tree->SetBranchAddress("nVtx",&nVtx);
  tree->SetBranchAddress("nPhotons",&nPhotons);
  tree->SetBranchAddress("metVal",&metVal);
  tree->SetBranchAddress("rawRho",&rawRho);
  tree->SetBranchAddress("massErrorUFCorr",&massError);
  

  ofstream out;
  out.open(outFile.c_str());

  TString PD;
  
  TString twenty10 = "2010";
  TString Aug = "Aug05ReReco";
  TString Oct = "Oct03ReReco";
  TString Jul = "Jul05ReReco";
  TString Pv1 = "Promptv1_2011B";

  TString ev;

  for (int i = 0; i < tree->GetEntries(); i++) 
    {

      tree->GetEntry(i);
      bool notDuplicateEvent = true;

      ULong64_t runId = Run;
      ULong64_t eventId = Event;
      ULong64_t lumiId = LumiSect;

    
      for (unsigned int n = 0; n < runVec3.size(); n++)
	{
	  if(runId == runVec3[n] && lumiId == lumiVec3[n] && eventId == eventVec3[n]){notDuplicateEvent = false;}
	}
    

      if (notDuplicateEvent)
	{
	  if(massZ1 < massZ1Cut) continue;
	  if(massZ2 < massZ2Cut) continue;
	  if(mass4l < m4lLowCut) continue;
	  if(mass4l > m4lHighCut) continue;

	  if( Run < 16404 ){ PD = twenty10; }
	  else if( Run > 160000 && Run <= 167913 ){ PD = Jul; }
	  else if( Run >= 170249 && Run <=172619 ){ PD = Aug; }
	  else if( Run >= 172620 && Run <= 175770){ PD = Oct; }
	  else if( Run >= 175832 ){ PD = Pv1; }
	  else{ PD = "error"; }
	  
	  runVec3.push_back(runId);
	  lumiVec3.push_back(lumiId);
	  eventVec3.push_back(eventId);
	  
          //counter++;
          if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 13 && abs(idL4) == 13 )
            {
              if( PD == Jul ){ PD = "/DoubleMu/Run2011A-05Jul2011ReReco-ECAL-v1/AOD"; }
              else if (PD == Oct){ PD = "/DoubleMu/Run2011A-03Oct2011-v1/AOD"; }
              else if (PD == Aug){ PD = "/DoubleMu/Run2011A-05Aug2011-v1/AOD"; }
              else if (PD = Pv1){ PD = "/DoubleMu/Run2011B-PromptReco-v1/AOD"; }

            }
          if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 11 && abs(idL4) == 11 )
            {
              if( PD == Jul ){ PD = "/DoubleElectron/Run2011A-05Jul2011ReReco-ECAL-v1/AOD"; }
              else if (PD == Oct){ PD ="/DoubleElectron/Run2011A-03Oct2011-v1/AOD"; }
              else if (PD == Aug){ PD ="/DoubleElectron/Run2011A-05Aug2011-v1/AOD"; }
              else if (PD = Pv1){ PD = "/DoubleElectron/Run2011B-PromptReco-v1/AOD"; }
            }
          if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 11 && abs(idL4) == 11 )
            {
              if( PD == Jul ){ PD = "/DoubleMu/Run2011A-05Jul2011ReReco-ECAL-v1/AOD"; }
              else if (PD == Oct){ PD ="/DoubleMu/Run2011A-03Oct2011-v1/AOD"; }
              else if (PD == Aug){ PD ="/DoubleMu/Run2011A-05Aug2011-v1/AOD"; }
              else if (PD = Pv1){ PD = "/DoubleMu/Run2011B-PromptReco-v1/AOD"; }
            }
          if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 13 && abs(idL4) == 13 )
            {
              if( PD == Jul ){ PD = "/DoubleElectron/Run2011A-05Jul2011ReReco-ECAL-v1/AOD"; }
              else if (PD == Oct){ PD ="/DoubleElectron/Run2011A-03Oct2011-v1/AOD"; }
              else if (PD == Aug){ PD ="/DoubleElectron/Run2011A-05Aug2011-v1/AOD"; }
              else if (PD = Pv1){ PD = "/DoubleElectron/Run2011B-PromptReco-v1/AOD"; }
            }


          out << PD << "_" << Run << ":" << LumiSect << ":" << Event <<  endl;

	}
    
    }

  runVec3.clear();
  lumiVec3.clear();
  eventVec3.clear();

  
  out.close();
}


void 
PlotHelper::printFourMomenta(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;
  double mass4l;
  double EL1, EL2, EL3, EL4;
  double pXL1, pXL2, pXL3, pXL4;
  double pYL1, pYL2, pYL3, pYL4;
  double pZL1, pZL2, pZL3, pZL4;
  ULong64_t Run, Event, LumiSect;
  double massZ1, massZ2;
  int idL1,idL2,idL3,idL4;

  TFile *f = new TFile(file,"READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("EL1",&EL1);
  tree->SetBranchAddress("EL2",&EL2);
  tree->SetBranchAddress("EL3",&EL3);
  tree->SetBranchAddress("EL4",&EL4);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("pXL1",&pXL1);
  tree->SetBranchAddress("pXL2",&pXL2);
  tree->SetBranchAddress("pXL3",&pXL3);
  tree->SetBranchAddress("pXL4",&pXL4);
  tree->SetBranchAddress("pYL1",&pYL1);
  tree->SetBranchAddress("pYL2",&pYL2);
  tree->SetBranchAddress("pYL3",&pYL3);
  tree->SetBranchAddress("pYL4",&pYL4);
  tree->SetBranchAddress("pZL1",&pZL1);
  tree->SetBranchAddress("pZL2",&pZL2);
  tree->SetBranchAddress("pZL3",&pZL3);
  tree->SetBranchAddress("pZL4",&pZL4);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);


  ofstream out;
  out.open(outFile.c_str());

  for (int i = 0; i < tree->GetEntries(); i++) 
    {

      tree->GetEntry(i);
      bool notDuplicateEvent = true;

      ULong64_t runId = Run;
      ULong64_t eventId = Event;
      ULong64_t lumiId = LumiSect;

    
      for (unsigned int n = 0; n < runVec4.size(); n++)
	{
	  if(runId == runVec4[n] && lumiId == lumiVec4[n] && eventId == eventVec4[n]){notDuplicateEvent = false;}
	}
    

      if (notDuplicateEvent)
	{
          runVec4.push_back(runId);
          lumiVec4.push_back(lumiId);
          eventVec4.push_back(eventId);

          if(massZ1 < massZ1Cut) continue;
          if(massZ2 < massZ2Cut) continue;
          if(mass4l < m4lLowCut) continue;
          if(mass4l > m4lHighCut) continue;


	  out << idL1 << "  " << idL2 << "  " << idL3 << "  " << idL4 << "  " 
	      << EL1 << "  " << pXL1 << "  " << pYL1 << "  " << pZL1 << "  " 
	      << EL2 << "  " << pXL2 << "  " << pYL2 << "  " << pZL2 << "  " 
	      << EL3 << "  " << pXL3 << "  " << pYL3 << "  " << pZL3 << "  "
	      << EL4 << "  " << pXL4 << "  " << pYL4 << "  " << pZL4 << endl;
	}
      
    }

  runVec4.clear();
  lumiVec4.clear();
  eventVec4.clear();  
  out.close();
}



void 
PlotHelper::printData(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;

  ULong64_t Run, Event, LumiSect;
  double mass4l, mass4lNoFSR, pT4l, massZ1, massZ2;
  double eta4l, phi4l;
  int idL1, idL2, idL3, idL4;
  double IPL1, IPL2, IPL3, IPL4;
  double dIPL1, dIPL2, dIPL3, dIPL4;
  double EL1, EL2, EL3, EL4;
  double EZ1, EZ2;
  double pTVBFJet1, pTVBFJet2, VBFDiJetMass, VBFDeltaEta, FisherDiscrim; 
  double pTL1, pTL2, pTL3, pTL4;
  double pXL1, pXL2, pXL3, pXL4;
  double pYL1, pYL2, pYL3, pYL4;
  double pZL1, pZL2, pZL3, pZL4;
  double pTZ1, pTZ2;
  double pXZ1, pXZ2;
  double pYZ1, pYZ2;
  double pZZ1, pZZ2;
  double etaL1, etaL2, etaL3, etaL4;
  double phiL1, phiL2, phiL3, phiL4;
  double isoCHL1, isoCHL2, isoCHL3, isoCHL4;
  double isoNHL1, isoNHL2, isoNHL3, isoNHL4;
  double isoPhotL1, isoPhotL2, isoPhotL3, isoPhotL4;
  double RelIsoL1, RelIsoL2, RelIsoL3, RelIsoL4;
  double RelIsoUCL1, RelIsoUCL2, RelIsoUCL3, RelIsoUCL4;
  double muRhoCor, elRhoCor;
  double worstIso, worstSIP;
  double cosTheta1, cosTheta2, Phi, cosThetaStar, phiStar1, phiStar2, Phi1, Phi2;
  int nJets, nVtx, nPhotons;
  double metVal, rawRho;
  int extraLep_id[10];
  double extraLep_pT[10], extraLep_iso[10], extraLep_e[10];
  double extraLep_pX[10], extraLep_pY[10], extraLep_pZ[10];
  double extraLep_eta[10], extraLep_phi[10], extraLep_sip[10];
  double extraLep_chIso[10], extraLep_nhIso[10], extraLep_phIso[10];
  double massError;
  double minM3l, minDeltR, maxP;
  double melaLD;
  double FSRPhot1_Pt,FSRPhot2_Pt,FSRPhot1_eta,FSRPhot2_eta,FSRPhot1_phi,FSRPhot2_phi;
  bool FSR_Z1,FSR_Z2;

  TBranch *b_extraLep_id;
  TBranch *b_extraLep_pT, *b_extraLep_iso, *b_extraLep_e;
  TBranch *b_extraLep_pX, *b_extraLep_pY, *b_extraLep_pZ;
  TBranch *b_extraLep_eta, *b_extraLep_phi, *b_extraLep_sip;
  TBranch *b_extraLep_nhIso, *b_extraLep_chIso, *b_extraLep_phIso;

  double thetaPhoton_deg, theta12_deg, theta13_deg, theta14_deg, minMass2Lep, maxMass2Lep;
  double thetaPhotonZ_deg;

  TFile *f = new TFile(file,"READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("mass4lNoFSR",&mass4lNoFSR);
  tree->SetBranchAddress("minMass3l",&minM3l);
  tree->SetBranchAddress("minDeltaR",&minDeltR);
  tree->SetBranchAddress("Z4l_maxP",&maxP);
  tree->SetBranchAddress("pT4l",&pT4l);
  tree->SetBranchAddress("Z4l_thetaPhoton_deg",&thetaPhoton_deg);
  tree->SetBranchAddress("Z4l_thetaPhotonZ_deg",&thetaPhotonZ_deg);
  tree->SetBranchAddress("Z4l_theta12_deg",&theta12_deg);
  tree->SetBranchAddress("Z4l_theta13_deg",&theta13_deg);
  tree->SetBranchAddress("Z4l_theta14_deg",&theta14_deg);
  tree->SetBranchAddress("Z4l_minMass2l",&minMass2Lep);
  tree->SetBranchAddress("Z4l_maxMass2l",&maxMass2Lep);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);
  tree->SetBranchAddress("eta4l",&eta4l);
  tree->SetBranchAddress("phi4l",&phi4l);
  tree->SetBranchAddress("EL1",&EL1);
  tree->SetBranchAddress("EL2",&EL2);
  tree->SetBranchAddress("EL3",&EL3);
  tree->SetBranchAddress("EL4",&EL4);
  tree->SetBranchAddress("EZ1",&EZ1);
  tree->SetBranchAddress("EZ2",&EZ2);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("pTL1",&pTL1);
  tree->SetBranchAddress("pTL2",&pTL2);
  tree->SetBranchAddress("pTL3",&pTL3);
  tree->SetBranchAddress("pTL4",&pTL4);
  tree->SetBranchAddress("pXL1",&pXL1);
  tree->SetBranchAddress("pXL2",&pXL2);
  tree->SetBranchAddress("pXL3",&pXL3);
  tree->SetBranchAddress("pXL4",&pXL4);
  tree->SetBranchAddress("pYL1",&pYL1);
  tree->SetBranchAddress("pYL2",&pYL2);
  tree->SetBranchAddress("pYL3",&pYL3);
  tree->SetBranchAddress("pYL4",&pYL4);
  tree->SetBranchAddress("pZL1",&pZL1);
  tree->SetBranchAddress("pZL2",&pZL2);
  tree->SetBranchAddress("pZL3",&pZL3);
  tree->SetBranchAddress("pZL4",&pZL4);
  tree->SetBranchAddress("pTZ1",&pTZ1);
  tree->SetBranchAddress("pTZ2",&pTZ2);
  tree->SetBranchAddress("pXZ1",&pXZ1);
  tree->SetBranchAddress("pXZ2",&pXZ2);
  tree->SetBranchAddress("pYZ1",&pYZ1);
  tree->SetBranchAddress("pYZ2",&pYZ2);
  tree->SetBranchAddress("pZZ1",&pZZ1);
  tree->SetBranchAddress("pZZ2",&pZZ2);
  tree->SetBranchAddress("etaL1",&etaL1);
  tree->SetBranchAddress("etaL2",&etaL2);
  tree->SetBranchAddress("etaL3",&etaL3);
  tree->SetBranchAddress("etaL4",&etaL4);
  tree->SetBranchAddress("phiL1",&phiL1);
  tree->SetBranchAddress("phiL2",&phiL2);
  tree->SetBranchAddress("phiL3",&phiL3);
  tree->SetBranchAddress("phiL4",&phiL4);
  tree->SetBranchAddress("IPL1",&IPL1);
  tree->SetBranchAddress("IPL2",&IPL2);
  tree->SetBranchAddress("IPL3",&IPL3);
  tree->SetBranchAddress("IPL4",&IPL4);
  tree->SetBranchAddress("dIPL1",&dIPL1);
  tree->SetBranchAddress("dIPL2",&dIPL2);
  tree->SetBranchAddress("dIPL3",&dIPL3);
  tree->SetBranchAddress("dIPL4",&dIPL4);
  tree->SetBranchAddress("isoNHL1",&isoNHL1);
  tree->SetBranchAddress("isoNHL2",&isoNHL2);
  tree->SetBranchAddress("isoNHL3",&isoNHL3);
  tree->SetBranchAddress("isoNHL4",&isoNHL4);
  tree->SetBranchAddress("isoCHL1",&isoCHL1);
  tree->SetBranchAddress("isoCHL2",&isoCHL2);
  tree->SetBranchAddress("isoCHL3",&isoCHL3);
  tree->SetBranchAddress("isoCHL4",&isoCHL4);
  tree->SetBranchAddress("isoPhotL1",&isoPhotL1);
  tree->SetBranchAddress("isoPhotL2",&isoPhotL2);
  tree->SetBranchAddress("isoPhotL3",&isoPhotL3);
  tree->SetBranchAddress("isoPhotL4",&isoPhotL4);
  tree->SetBranchAddress("RelIsoL1",&RelIsoL1);
  tree->SetBranchAddress("RelIsoL2",&RelIsoL2);
  tree->SetBranchAddress("RelIsoL3",&RelIsoL3);
  tree->SetBranchAddress("RelIsoL4",&RelIsoL4);
  tree->SetBranchAddress("RelIsoUCL1",&RelIsoUCL1);
  tree->SetBranchAddress("RelIsoUCL2",&RelIsoUCL2);
  tree->SetBranchAddress("RelIsoUCL3",&RelIsoUCL3);
  tree->SetBranchAddress("RelIsoUCL4",&RelIsoUCL4);
  tree->SetBranchAddress("muRhoCor",&muRhoCor);
  tree->SetBranchAddress("elRhoCor",&elRhoCor);
  tree->SetBranchAddress("worstIso",&worstIso);
  tree->SetBranchAddress("worstSIP",&worstSIP);
  tree->SetBranchAddress("cosTheta1",&cosTheta1);
  tree->SetBranchAddress("cosTheta2",&cosTheta2);
  tree->SetBranchAddress("Phi",&Phi);
  tree->SetBranchAddress("cosThetaStar",&cosThetaStar);
  tree->SetBranchAddress("phiStar1",&phiStar1);
  tree->SetBranchAddress("phiStar2",&phiStar2);
  tree->SetBranchAddress("Phi1",&Phi1);
  tree->SetBranchAddress("Phi2",&Phi2);
  tree->SetBranchAddress("nJets",&nJets);
  tree->SetBranchAddress("nVtx",&nVtx);
  tree->SetBranchAddress("nPhotons",&nPhotons);
  tree->SetBranchAddress("metVal",&metVal);
  tree->SetBranchAddress("rawRho",&rawRho);
  tree->SetBranchAddress("extraLep_id", extraLep_id, &b_extraLep_id);
  tree->SetBranchAddress("extraLep_pT", extraLep_pT, &b_extraLep_pT);
  tree->SetBranchAddress("extraLep_pX", extraLep_pX, &b_extraLep_pX);
  tree->SetBranchAddress("extraLep_pY", extraLep_pY, &b_extraLep_pY);
  tree->SetBranchAddress("extraLep_pZ", extraLep_pZ, &b_extraLep_pZ);
  tree->SetBranchAddress("extraLep_e",extraLep_e, &b_extraLep_e);
  tree->SetBranchAddress("extraLep_iso", extraLep_iso, &b_extraLep_iso);
  tree->SetBranchAddress("extraLep_chIso",extraLep_chIso, &b_extraLep_chIso);
  tree->SetBranchAddress("extraLep_nhIso", extraLep_nhIso, &b_extraLep_nhIso);
  tree->SetBranchAddress("extraLep_phIso",extraLep_phIso, &b_extraLep_phIso);
  tree->SetBranchAddress("extraLep_eta",extraLep_eta, &b_extraLep_eta);
  tree->SetBranchAddress("extraLep_phi",extraLep_phi, &b_extraLep_phi);
  tree->SetBranchAddress("extraLep_sip",extraLep_sip, &b_extraLep_sip);
  tree->SetBranchAddress("massErrorUFCorr",&massError);
  tree->SetBranchAddress("melaLD",&melaLD);
  tree->SetBranchAddress("FSRPhot1_Pt",&FSRPhot1_Pt);
  tree->SetBranchAddress("FSRPhot2_Pt",&FSRPhot2_Pt);
  tree->SetBranchAddress("FSRPhot1_eta",&FSRPhot1_eta);
  tree->SetBranchAddress("FSRPhot2_eta",&FSRPhot2_eta);
  tree->SetBranchAddress("FSRPhot1_phi",&FSRPhot1_phi);
  tree->SetBranchAddress("FSRPhot2_phi",&FSRPhot2_phi);
  tree->SetBranchAddress("FSR_Z1",&FSR_Z1);
  tree->SetBranchAddress("FSR_Z2",&FSR_Z2);

  tree->SetBranchAddress("pTVBFJet1",&pTVBFJet1);
  tree->SetBranchAddress("pTVBFJet2",&pTVBFJet2);
  tree->SetBranchAddress("VBFDiJetMass",&VBFDiJetMass);
  tree->SetBranchAddress("VBFDeltaEta",&VBFDeltaEta);
  tree->SetBranchAddress("FisherDiscrim",&FisherDiscrim);




  ofstream out;
  out.open(outFile.c_str());

  int counter = 0, counter_4mu = 0, counter_4e = 0, counter_2e2mu = 0, VBFCounter = 0;

  tree->BuildIndex("Run");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();

  for( int i = 0; i < index->GetN(); i++ )
    {

      Long64_t local = tree->LoadTree( index->GetIndex()[i] );
      tree->GetEntry(local);

      bool notDuplicateEvent = true;

      ULong64_t runId = Run;
      ULong64_t eventId = Event;
      ULong64_t lumiId = LumiSect;

    
      for (unsigned int n = 0; n < runVec.size(); n++)
	{
	  if(runId == runVec[n] && lumiId == lumiVec[n] && eventId == eventVec[n]){notDuplicateEvent = false;}
	}
    

      if (notDuplicateEvent)
	{

	  if(massZ1 < massZ1Cut) continue;
	  if(massZ2 < massZ2Cut) continue;
	  if(mass4l < m4lLowCut) continue;
	  if(mass4l > m4lHighCut) continue;

	  if( abs(idL1) == 13 && pTL1 < ptMuCut ) continue;
	  if( abs(idL2) == 13 && pTL2 < ptMuCut ) continue;
	  if( abs(idL3) == 13 && pTL3 < ptMuCut ) continue;
	  if( abs(idL4) == 13 && pTL4 < ptMuCut ) continue;

	  if( abs(idL1) == 11 && pTL1 < ptElCut ) continue;
	  if( abs(idL2) == 11 && pTL2 < ptElCut ) continue;
	  if( abs(idL3) == 11 && pTL3 < ptElCut ) continue;
	  if( abs(idL4) == 11 && pTL4 < ptElCut ) continue;
	  

	  runVec.push_back(runId);
	  lumiVec.push_back(lumiId);
	  eventVec.push_back(eventId);
	  
	  counter++;
	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 13 && abs(idL4) == 13 ){counter_4mu++;}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 11 && abs(idL4) == 11 ){counter_4e++;}
	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 11 && abs(idL4) == 11 ){counter_2e2mu++;}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 13 && abs(idL4) == 13 ){counter_2e2mu++;}


		  

	  out << "########" << endl
	      << counter << endl
	      << "########" << endl;

	      out.precision(7);
	      out << left
		  << "Run:  " << Run          << "    Event: " << Event    << "   LumiSec: " << LumiSect << endl
		  << "m4l:  " << mass4l   << "   m4lError: " << massError   << "      mZ1: " << massZ1   << "       mZ2: " << massZ2  << endl
		  << "FSR_Z1  " << FSR_Z1 << " pt: " << FSRPhot1_Pt << "    FSR_Z2: " << FSR_Z2 << " pt: " << FSRPhot2_Pt << endl  
		  << "pTVBFJet1: " << pTVBFJet1 << "  pTVBFJet2: " << pTVBFJet2 << "  jjmass/deta: " << VBFDiJetMass <<"/"<<VBFDeltaEta
		  << "  Fisher: " << FisherDiscrim << endl 
		  << "cos(theta1): " << cosTheta1 << "  cos(theta2): " << cosTheta2 << "  cos(theta*): " << cosThetaStar << endl
		  << "phi: " << Phi << "  Phi1: " << Phi1 << "  Phi2: " << Phi2 << "  Phi*1: " << phiStar1 << "  Phi*2: " << phiStar2 << endl
		  << " MuonRhoCorr: " << muRhoCor << "      ElecRhoCor: " << elRhoCor << " Rho: " << rawRho << endl
		  << "MELALD: " << melaLD << "   m4lNoFSR: " << mass4lNoFSR << endl
		  << "minM3l:  " << minM3l << "   maxP:  " << maxP << "  minDeltaR:  " << minDeltR << endl
		  << "thetaPhotonZ:  " << thetaPhotonZ_deg << "  thetaPhoton:  " << thetaPhoton_deg << endl 
		  << "maxMass2L:  " << maxMass2Lep << "  minMass2L:  " << minMass2Lep << endl
		  << "pT4l: " << pT4l         << "    eta4l: " << eta4l    << "     phi4l: " << phi4l   << endl
		  << "worstIso: " << worstIso << "    worstSIP: " << worstSIP << endl
		  << "nJets: " << nJets << "    nPhotons: " << nPhotons << "   nVtx: " << nVtx << endl
		  << "MET: " << metVal << endl
		  << "L1: " << idL1         << "       pT: " << pTL1 << "  eta: " << etaL1 << "  phi: " << phiL1
		  << "  iso: CH " << isoCHL1 << " NH " << isoNHL1 << " Ph " << isoPhotL1 << "  RelIso: " << RelIsoL1 << "," << RelIsoUCL1 << endl
		  << "L2: " << idL2         << "       pT: " << pTL2 << "  eta: " << etaL2 << "  phi: " << phiL2
		  << "  iso: CH " << isoCHL2 << " NH " << isoNHL2 << " Ph " << isoPhotL2 << "  RelIso: " << RelIsoL2 << "," << RelIsoUCL2 << endl
		  << "L3: " << idL3         << "       pT: " << pTL3 << "  eta: " << etaL3 << "  phi: " << phiL3
		  << "  iso: CH " << isoCHL3 << " NH " << isoNHL3 << " Ph " << isoPhotL3 << "  RelIso: " << RelIsoL3 << "," << RelIsoUCL3 << endl
		  << "L4: " << idL4         << "       pT: " << pTL4 << "  eta: " << etaL4 << "  phi: " << phiL4
		  << "  iso: CH " << isoCHL4 << " NH " << isoNHL4 << " Ph " << isoPhotL4 << "  RelIso: " << RelIsoL4 << "," << RelIsoUCL4 << endl
		  << "  ID             px             py              pz             e        IP       deltaIP  "  << endl
		  << "-------------------------------------------------------------------------------------------" << endl
		  << right
		  << "  Z1    " << setprecision(4) << setw(12) << pXZ1 << "   " << setprecision(4) << setw(12) << pYZ1 << "   " 
		  << setprecision(4) << setw(12) << pZZ1 << "  " << setprecision(4) << setw(12) <<  EZ1 << endl
		  << "  Z2    " << setprecision(4) << setw(12) << pXZ2 << "   " << setprecision(4) << setw(12) <<  pYZ2 << "   " 
		  << setprecision(4) << setw(12) << pZZ2 << "  " << setprecision(4) << setw(12) << EZ2 << endl
		  << "  L1    " << setprecision(4) << setw(12) << pXL1 << "   " << setprecision(4) << setw(12) << pYL1 << "   " 
		  << setprecision(4) << setw(12) << pZL1 << "  " << setprecision(4) << setw(12) << EL1 
		  << "   " << setprecision(4) << setw(12) << IPL1 << "  " << setprecision(4) << setw(12) << dIPL1 << endl
		  << "  L2    " << setprecision(4) << setw(12) << pXL2 << "   " << setprecision(4) << setw(12) << pYL2 << "   " 
		  << setprecision(4) << setw(12) << pZL2 << "  " << setprecision(4) << setw(12) << EL2 
                  << "   " << setprecision(4) << setw(12) << IPL2 << "  " << setprecision(4) << setw(12) << dIPL2 << endl
		  << "  L3    " << setprecision(4) << setw(12) << pXL3 << "   " << setprecision(4) << setw(12) << pYL3 << "   " 
		  << setprecision(4) << setw(12) << pZL3 << "  " << setprecision(4) << setw(12) << EL3 
                  << "   " << setprecision(4) << setw(12) << IPL3 << "  " << setprecision(4) << setw(12) << dIPL3 << endl
		  << "  L4    " << setprecision(4) << setw(12) << pXL4 << "   " << setprecision(4) << setw(12) << pYL4 << "   " 
		  << setprecision(4) << setw(12) << pZL4 << "  " << setprecision(4) << setw(12) << EL4 
                  << "   " << setprecision(4) << setw(12) << IPL4 << "  " << setprecision(4) << setw(12) << dIPL4 << endl;
	      

	      out << "  Other Leptons" << endl
		  << "---------------" << endl;

	      for( int k = 0; k < 10; k++ )
		{
		  if( abs(extraLep_id[k]) == 11 || abs(extraLep_id[k]) == 13 )
		    {
		      if( extraLep_pT[k] != pTL1 && extraLep_pT[k] != pTL2 && extraLep_pT[k] != pTL3 && extraLep_pT[k] != pTL4 )
			{
		      
		      out << " " << extraLep_id[k] << setprecision(4) << setw(12) << extraLep_pX[k] << "   " << setprecision(4) << setw(12) << extraLep_pY[k] << "   "
			  << setprecision(4) << setw(12) << extraLep_pZ[k] << "  " << setprecision(4) << setw(12) <<  extraLep_e[k] << "  "  << endl
			  << "      pT: " << extraLep_pT[k] << "   eta: " << extraLep_eta[k] << "   phi: " << extraLep_phi[k]  
			  << "  SIP: " << extraLep_sip[k] << "  iso: " << extraLep_iso[k]   
			  << "  CH " << extraLep_chIso[k] << " NH " << extraLep_nhIso[k] << " h " << extraLep_phIso[k] << endl;
			}
		    }
		  
		}

	      
	      out << endl << endl;
	      
	      if(pTVBFJet1 > 0 && pTVBFJet2 > 0 && VBFDiJetMass > 300 && VBFDeltaEta > 3) VBFCounter++;
	      
	}
    
    }

  out << counter << " event candidates found" << endl;
  out << counter_4mu << " 4mu events " << endl;
  out << counter_4e << " 4e events " << endl;
  out << counter_2e2mu << " 2e2mu events" << endl;
  out << VBFCounter << " VBF tagged events" << endl;

  runVec.clear();
  lumiVec.clear();
  eventVec.clear();  
  
  out.close();
}



void 
PlotHelper::printDataShort(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;

  ULong64_t Run, Event, LumiSect;
  double mass4l,  massZ1, massZ2;
  int idL1, idL2, idL3, idL4;
  double pTL1, pTL2, pTL3, pTL4;
  double melaLD;
  string type;

  TFile *f = new TFile(file,"READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("pTL1",&pTL1);
  tree->SetBranchAddress("pTL2",&pTL2);
  tree->SetBranchAddress("pTL3",&pTL3);
  tree->SetBranchAddress("pTL4",&pTL4);
  tree->SetBranchAddress("melaLD",&melaLD);


  ofstream out;
  out.open(outFile.c_str());

  int counter = 0, counter_4mu = 0, counter_4e = 0, counter_2e2mu = 0;

  tree->BuildIndex("Run");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();

  for( int i = 0; i < index->GetN(); i++ )
    {

      Long64_t local = tree->LoadTree( index->GetIndex()[i] );
      tree->GetEntry(local);

      bool notDuplicateEvent = true;

      ULong64_t runId = Run;
      ULong64_t eventId = Event;
      ULong64_t lumiId = LumiSect;
  
      for (unsigned int n = 0; n < runVec6.size(); n++)
	{
	  if(runId == runVec6[n] && lumiId == lumiVec6[n] && eventId == eventVec6[n]){notDuplicateEvent = false;}
	}
    

      if (notDuplicateEvent)
	{

	  if(massZ1 < massZ1Cut) continue;
	  if(massZ2 < massZ2Cut) continue;
	  if(mass4l < m4lLowCut) continue;
	  if(mass4l > m4lHighCut) continue;

	  if( abs(idL1) == 13 && pTL1 < ptMuCut ) continue;
	  if( abs(idL2) == 13 && pTL2 < ptMuCut ) continue;
	  if( abs(idL3) == 13 && pTL3 < ptMuCut ) continue;
	  if( abs(idL4) == 13 && pTL4 < ptMuCut ) continue;

	  if( abs(idL1) == 11 && pTL1 < ptElCut ) continue;
	  if( abs(idL2) == 11 && pTL2 < ptElCut ) continue;
	  if( abs(idL3) == 11 && pTL3 < ptElCut ) continue;
	  if( abs(idL4) == 11 && pTL4 < ptElCut ) continue;

	  runVec6.push_back(runId);
	  lumiVec6.push_back(lumiId);
	  eventVec6.push_back(eventId);
	  
	  counter++;
	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 13 && abs(idL4) == 13 ){counter_4mu++;type="4mu";}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 11 && abs(idL4) == 11 ){counter_4e++;type="4e";}
	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 11 && abs(idL4) == 11 ){counter_2e2mu++;type="2mu2e";}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 13 && abs(idL4) == 13 ){counter_2e2mu++;type="2e2mu";}
		  
	      out.precision(7);
	      out << left
		  << setw(10) << Run  <<  "   " << setw(10) <<  Event  << "  " << setw(10) << LumiSect << "  " << setw(5) << type  
		  << "  " << setw(8) <<  mass4l   << "  " << setw(8) <<  massZ1  << "  " << setw(8) <<  massZ2  << "  " << melaLD << endl;
		
	}
    
    }

  out << endl << endl;
  out << counter << " event candidates found" << endl;
  out << counter_4mu << " 4mu events " << endl;
  out << counter_4e << " 4e events " << endl;
  out << counter_2e2mu << " 2e2mu events" << endl;

  runVec6.clear();
  lumiVec6.clear();
  eventVec6.clear();  
  
  out.close();
}


void 
PlotHelper::printDataSync(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;

  ULong64_t Run, Event, LumiSect;
  double mass4l,  massZ1, massZ2,massErrRaw,massErrCorr;
  int idL1, idL2, idL3, idL4;
  double pTL1, pTL2, pTL3, pTL4;
  double melaLD,pT4l;
  string type;
  bool VBFJet1, VBFJet2;
  double pTVBFJet1, pTVBFJet2;
  double JHUKD_H_qqZZ_noPDF,JHUKD_H_h0M_noPDF,JHUKD_H_h0P_noPDF;
  double JHUKD_H_h1P_noPDF,JHUKD_H_h1M_noPDF,JHUKD_H_ggh2P_noPDF;
  double JHUKD_H_qqh2P_noPDF;

  double VBFDiJetMass, VBFDeltaEta, FisherDiscrim;

  int nJets = 0;

  TFile *f = new TFile(file,"READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("pTL1",&pTL1);
  tree->SetBranchAddress("pTL2",&pTL2);
  tree->SetBranchAddress("pTL3",&pTL3);
  tree->SetBranchAddress("pTL4",&pTL4);
  tree->SetBranchAddress("pT4l",&pT4l);
  tree->SetBranchAddress("melaLD",&melaLD);
  tree->SetBranchAddress("massErrorUF",&massErrRaw);
  tree->SetBranchAddress("massErrorUFCorr",&massErrCorr);

  tree->SetBranchAddress("VBFJet1",&VBFJet1);
  tree->SetBranchAddress("pTVBFJet1",&pTVBFJet1);
  tree->SetBranchAddress("VBFJet2",&VBFJet2);
  tree->SetBranchAddress("pTVBFJet2",&pTVBFJet2);
  tree->SetBranchAddress("VBFDiJetMass",&VBFDiJetMass);
  tree->SetBranchAddress("VBFDeltaEta",&VBFDeltaEta);
  tree->SetBranchAddress("FisherDiscrim",&FisherDiscrim);

  tree->SetBranchAddress("JHUKD_H_qqZZ_noPDF",&JHUKD_H_qqZZ_noPDF);
  tree->SetBranchAddress("JHUKD_H_h0M_noPDF",&JHUKD_H_h0M_noPDF);
  tree->SetBranchAddress("JHUKD_H_h0P_noPDF",&JHUKD_H_h0P_noPDF);
  tree->SetBranchAddress("JHUKD_H_h1P_noPDF",&JHUKD_H_h1P_noPDF);
  tree->SetBranchAddress("JHUKD_H_h1M_noPDF",&JHUKD_H_h1M_noPDF);
  tree->SetBranchAddress("JHUKD_H_ggh2P_noPDF",&JHUKD_H_ggh2P_noPDF);
  tree->SetBranchAddress("JHUKD_H_qqh2P_noPDF",&JHUKD_H_qqh2P_noPDF);
	  

  ofstream out;
  out.open(outFile.c_str());

  tree->BuildIndex("Run");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();

  for( int i = 0; i < index->GetN(); i++ )
    {

      Long64_t local = tree->LoadTree( index->GetIndex()[i] );
      tree->GetEntry(local);

      bool notDuplicateEvent = true;

      ULong64_t runId = Run;
      ULong64_t eventId = Event;
      ULong64_t lumiId = LumiSect;
  
      for (unsigned int n = 0; n < runVec8.size(); n++)
	{
	  if(runId == runVec8[n] && lumiId == lumiVec8[n] && eventId == eventVec8[n]){notDuplicateEvent = false;}
	}
    

      if (notDuplicateEvent)
	{

	  if(massZ1 < massZ1Cut) continue;
	  if(massZ2 < massZ2Cut) continue;
	  if(mass4l < m4lLowCut) continue;
	  if(mass4l > m4lHighCut) continue;

	  //if(mass4l > 110 && mass4l < 140) continue;
	  //if(mass4l > 300) continue;
	  
	  if( abs(idL1) == 13 && pTL1 < ptMuCut ) continue;
	  if( abs(idL2) == 13 && pTL2 < ptMuCut ) continue;
	  if( abs(idL3) == 13 && pTL3 < ptMuCut ) continue;
	  if( abs(idL4) == 13 && pTL4 < ptMuCut ) continue;

	  if( abs(idL1) == 11 && pTL1 < ptElCut ) continue;
	  if( abs(idL2) == 11 && pTL2 < ptElCut ) continue;
	  if( abs(idL3) == 11 && pTL3 < ptElCut ) continue;
	  if( abs(idL4) == 11 && pTL4 < ptElCut ) continue;

	  runVec8.push_back(runId);
	  lumiVec8.push_back(lumiId);
	  eventVec8.push_back(eventId);

	  //out <<Run<<":"<<LumiSect<<":"<<Event<<":";
	  //out.setf(ios::fixed,ios::floatfield);
	  //out.precision(2);
	  //out << mass4l <<":"<< massZ1 <<":"<< massZ2 <<":"
	  //   << massErrRaw <<":"<< massErrCorr <<":";
	  //out.setf(ios::fixed,ios::floatfield);
	  //out.precision(3);
	  //out << melaLD <<endl;
	  nJets = 0;

	  if((VBFJet1 && !VBFJet2) || (VBFJet2 && !VBFJet1)) nJets = 1;
	  if(!VBFJet2 && !VBFJet1) nJets = 0;
	  if(VBFJet1 && VBFJet2) nJets = 2;

	  if(!VBFJet1){pTVBFJet1 = 0; VBFDiJetMass = 0; VBFDeltaEta = 0;}
	  if(!VBFJet2){pTVBFJet2 = 0; VBFDiJetMass = 0; VBFDeltaEta = 0;}

	  out <<Run<<":"<<LumiSect<<":"<<Event<<":";
	  out.setf(ios::fixed,ios::floatfield);
	  out.precision(2);
	  out << mass4l <<":"<< massZ1 <<":"<< massZ2 <<":"
	      << massErrRaw <<":"<< massErrCorr <<":";
	  out.setf(ios::fixed,ios::floatfield);
	  out.precision(3);
	  out << JHUKD_H_qqZZ_noPDF<<":";
	  out.precision(2);
	  out << pT4l << ":";
	  out.precision(0);
	  out << nJets << ":";
	  out.precision(2);
	  out << pTVBFJet1 << ":" << pTVBFJet2 << ":" << VBFDiJetMass << ":";
	  out.precision(3);
	  out << VBFDeltaEta << ":" << FisherDiscrim << ":" << JHUKD_H_h0M_noPDF << ":" << JHUKD_H_h0P_noPDF << ":" <<JHUKD_H_h1P_noPDF 
	      << ":" << JHUKD_H_h1M_noPDF << ":" << JHUKD_H_ggh2P_noPDF << ":" << JHUKD_H_qqh2P_noPDF <<endl;


	}
    
    }


  runVec8.clear();
  lumiVec8.clear();
  eventVec8.clear();  
  
  out.close();
}






void 
PlotHelper::printDataErr(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;

  TRandom *pRnd = new TRandom();



  ULong64_t Run, Event, LumiSect;
  double mass4l,m4lerr;
  double massZ1,massZ2;
  int idL1, idL2, idL3, idL4;
  string type;

  TFile *f = new TFile(file,"READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("massErrorUFCorr",&m4lerr);

  ofstream out;
  out.open(outFile.c_str());

  int counter = 0, counter_4mu = 0, counter_4e = 0, counter_2e2mu = 0;

  tree->BuildIndex("Run");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();

  for( int i = 0; i < index->GetN(); i++ )
    {

      Long64_t local = tree->LoadTree( index->GetIndex()[i] );
      tree->GetEntry(local);

      bool notDuplicateEvent = true;

      ULong64_t runId = Run;
      ULong64_t eventId = Event;
      ULong64_t lumiId = LumiSect;
  
      for (unsigned int n = 0; n < runVec7.size(); n++)
	{
	  if(runId == runVec7[n] && lumiId == lumiVec7[n] && eventId == eventVec7[n]){notDuplicateEvent = false;}
	}
    

      if (notDuplicateEvent)
	{

          if(massZ1 < massZ1Cut) continue;
          if(massZ2 < massZ2Cut) continue;
          if(mass4l < m4lLowCut) continue;
          if(mass4l > m4lHighCut) continue;

	  runVec7.push_back(runId);
	  lumiVec7.push_back(lumiId);
	  eventVec7.push_back(eventId);

	  counter++;
	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 13 && abs(idL4) == 13 ){counter_4mu++;type="4mu";}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 11 && abs(idL4) == 11 ){counter_4e++;type="4e";}
	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 11 && abs(idL4) == 11 ){counter_2e2mu++;type="2mu2e";}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 13 && abs(idL4) == 13 ){counter_2e2mu++;type="2e2mu";}
		  
	      out.precision(7);
	      //out << left << setw(5) << type << "  " << setw(8) << Run << "  " << setw(10) << Event << "   " << setw(8) << 
	      out << left << setw(5) << type << "  " << setw(8) << mass4l   << "  " << setw(8) << m4lerr << "   " << pRnd->Uniform(0,5) << endl;
		
	}
    
    }

  out << endl << endl;
  out << counter << " event candidates found" << endl;
  out << counter_4mu << " 4mu events " << endl;
  out << counter_4e << " 4e events " << endl;
  out << counter_2e2mu << " 2e2mu events" << endl;

  runVec7.clear();
  lumiVec7.clear();
  eventVec7.clear();  
  
  out.close();
}


void 
PlotHelper::printDataTex(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut,double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;


  ULong64_t Run, Event, LumiSect;
  double mass4l, pT4l, massZ1, massZ2;
  double eta4l, phi4l;
  int idL1, idL2, idL3, idL4;
  double EL1, EL2, EL3, EL4;
  double EZ1, EZ2;
  double pTL1, pTL2, pTL3, pTL4;
  double pXL1, pXL2, pXL3, pXL4;
  double pYL1, pYL2, pYL3, pYL4;
  double pZL1, pZL2, pZL3, pZL4;
  double pTZ1, pTZ2;
  double pXZ1, pXZ2;
  double pYZ1, pYZ2;
  double pZZ1, pZZ2;
  double etaL1, etaL2, etaL3, etaL4;
  double phiL1, phiL2, phiL3, phiL4;
  double isoCHL1, isoCHL2, isoCHL3, isoCHL4;
  double isoNHL1, isoNHL2, isoNHL3, isoNHL4;
  double isoPhotL1, isoPhotL2, isoPhotL3, isoPhotL4;
  double RelIsoL1, RelIsoL2, RelIsoL3, RelIsoL4;
  double muRhoCor, elRhoCor;
  double worstIso, worstSIP;
  double cosTheta1, cosTheta2, Phi, cosThetaStar, phiStar1, phiStar2, Phi1, Phi2;
  int nJets, nVtx, nPhotons;
  double metVal, rawRho;
  double massError;
  double melaLD;

  TFile *f = new TFile(file,"READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("pT4l",&pT4l);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);
  tree->SetBranchAddress("eta4l",&eta4l);
  tree->SetBranchAddress("phi4l",&phi4l);
  tree->SetBranchAddress("EL1",&EL1);
  tree->SetBranchAddress("EL2",&EL2);
  tree->SetBranchAddress("EL3",&EL3);
  tree->SetBranchAddress("EL4",&EL4);
  tree->SetBranchAddress("EZ1",&EZ1);
  tree->SetBranchAddress("EZ2",&EZ2);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("pTL1",&pTL1);
  tree->SetBranchAddress("pTL2",&pTL2);
  tree->SetBranchAddress("pTL3",&pTL3);
  tree->SetBranchAddress("pTL4",&pTL4);
  tree->SetBranchAddress("pXL1",&pXL1);
  tree->SetBranchAddress("pXL2",&pXL2);
  tree->SetBranchAddress("pXL3",&pXL3);
  tree->SetBranchAddress("pXL4",&pXL4);
  tree->SetBranchAddress("pYL1",&pYL1);
  tree->SetBranchAddress("pYL2",&pYL2);
  tree->SetBranchAddress("pYL3",&pYL3);
  tree->SetBranchAddress("pYL4",&pYL4);
  tree->SetBranchAddress("pZL1",&pZL1);
  tree->SetBranchAddress("pZL2",&pZL2);
  tree->SetBranchAddress("pZL3",&pZL3);
  tree->SetBranchAddress("pZL4",&pZL4);
  tree->SetBranchAddress("pTZ1",&pTZ1);
  tree->SetBranchAddress("pTZ2",&pTZ2);
  tree->SetBranchAddress("pXZ1",&pXZ1);
  tree->SetBranchAddress("pXZ2",&pXZ2);
  tree->SetBranchAddress("pYZ1",&pYZ1);
  tree->SetBranchAddress("pYZ2",&pYZ2);
  tree->SetBranchAddress("pZZ1",&pZZ1);
  tree->SetBranchAddress("pZZ2",&pZZ2);
  tree->SetBranchAddress("etaL1",&etaL1);
  tree->SetBranchAddress("etaL2",&etaL2);
  tree->SetBranchAddress("etaL3",&etaL3);
  tree->SetBranchAddress("etaL4",&etaL4);
  tree->SetBranchAddress("phiL1",&phiL1);
  tree->SetBranchAddress("phiL2",&phiL2);
  tree->SetBranchAddress("phiL3",&phiL3);
  tree->SetBranchAddress("phiL4",&phiL4);
  tree->SetBranchAddress("isoNHL1",&isoNHL1);
  tree->SetBranchAddress("isoNHL2",&isoNHL2);
  tree->SetBranchAddress("isoNHL3",&isoNHL3);
  tree->SetBranchAddress("isoNHL4",&isoNHL4);
  tree->SetBranchAddress("isoCHL1",&isoCHL1);
  tree->SetBranchAddress("isoCHL2",&isoCHL2);
  tree->SetBranchAddress("isoCHL3",&isoCHL3);
  tree->SetBranchAddress("isoCHL4",&isoCHL4);
  tree->SetBranchAddress("isoPhotL1",&isoPhotL1);
  tree->SetBranchAddress("isoPhotL2",&isoPhotL2);
  tree->SetBranchAddress("isoPhotL3",&isoPhotL3);
  tree->SetBranchAddress("isoPhotL4",&isoPhotL4);
  tree->SetBranchAddress("RelIsoL1",&RelIsoL1);
  tree->SetBranchAddress("RelIsoL2",&RelIsoL2);
  tree->SetBranchAddress("RelIsoL3",&RelIsoL3);
  tree->SetBranchAddress("RelIsoL4",&RelIsoL4);
  tree->SetBranchAddress("muRhoCor",&muRhoCor);
  tree->SetBranchAddress("elRhoCor",&elRhoCor);
  tree->SetBranchAddress("worstIso",&worstIso);
  tree->SetBranchAddress("worstSIP",&worstSIP);
  tree->SetBranchAddress("cosTheta1",&cosTheta1);
  tree->SetBranchAddress("cosTheta2",&cosTheta2);
  tree->SetBranchAddress("Phi",&Phi);
  tree->SetBranchAddress("cosThetaStar",&cosThetaStar);
  tree->SetBranchAddress("phiStar1",&phiStar1);
  tree->SetBranchAddress("phiStar2",&phiStar2);
  tree->SetBranchAddress("Phi1",&Phi1);
  tree->SetBranchAddress("Phi2",&Phi2);
  tree->SetBranchAddress("nJets",&nJets);
  tree->SetBranchAddress("nVtx",&nVtx);
  tree->SetBranchAddress("nPhotons",&nPhotons);
  tree->SetBranchAddress("metVal",&metVal);
  tree->SetBranchAddress("rawRho",&rawRho);
  tree->SetBranchAddress("massErrorUFCorr",&massError);
  tree->SetBranchAddress("melaLD",&melaLD);


  ofstream out;
  out.open(outFile.c_str());

  string ev = "", z1 = "", z2 = "";


  tree->BuildIndex("mass4l");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();

  for( int i = 0; i < index->GetN(); i++ )
    {

      Long64_t local = tree->LoadTree( index->GetIndex()[i] );
      tree->GetEntry(local);

      bool notDuplicateEvent = true;
      
      ULong64_t runId = Run;
      ULong64_t eventId = Event;
      ULong64_t lumiId = LumiSect;
      
      
      for (unsigned int n = 0; n < runVec2.size(); n++)
	{
	  if(runId == runVec2[n] && lumiId == lumiVec2[n] && eventId == eventVec2[n]){notDuplicateEvent = false;}
	}
      
    
      
      if (notDuplicateEvent)
	{
	  
	  runVec2.push_back(runId);
	  lumiVec2.push_back(lumiId);
	  eventVec2.push_back(eventId);
	  
          if(massZ1 < massZ1Cut) continue;
          if(massZ2 < massZ2Cut) continue;
          if(mass4l < m4lLowCut) continue;
          if(mass4l > m4lHighCut) continue;


	  if( abs(idL1) == 13 && pTL1 < ptMuCut ) continue;
	  if( abs(idL2) == 13 && pTL2 < ptMuCut ) continue;
	  if( abs(idL3) == 13 && pTL3 < ptMuCut ) continue;
	  if( abs(idL4) == 13 && pTL4 < ptMuCut ) continue;

	  if( abs(idL1) == 11 && pTL1 < ptElCut ) continue;
	  if( abs(idL2) == 11 && pTL2 < ptElCut ) continue;
	  if( abs(idL3) == 11 && pTL3 < ptElCut ) continue;
	  if( abs(idL4) == 11 && pTL4 < ptElCut ) continue;


	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 13 && abs(idL4) == 13 ){ev = "$4\\mu$"; z1 = "($\\mu\\mu$)"; z2 = "($\\mu\\mu$)";}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 11 && abs(idL4) == 11 ){ev = "$4e$"; z1 = "($ee$)"; z2 = "($ee$)";}
	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 11 && abs(idL4) == 11 ){ev = "$2e2\\mu$"; z1 = "($\\mu\\mu$)"; z2 = "($ee$)";}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 13 && abs(idL4) == 13 ){ev = "$2e2\\mu$"; z1 = "($ee$)"; z2 = "($\\mu\\mu$)";}
	  
	    
	      out.precision(4);
	      out << left << ev 
		  << "  &  " << setprecision(5) <<  mass4l  << setprecision(3) << "  &  " << massError 
		  << "  &  " << setprecision(4) <<  pT4l << "  &  " <<  setprecision(4) << eta4l 
		  << "  &  " << setprecision(4) << melaLD 
		  << "  &  " << setprecision(4) << massZ1 << " " << z1 << "  &  " << setprecision(4) << massZ2 << " " << z2 
		  << "  &  " << pTL1 << "  &  " << etaL1 << "  &  " << phiL1  
		  << "  &  " << pTL2 << "  &  " << etaL2 << "  &  " << phiL2
                  << "  &  " << pTL3 << "  &  " << etaL3 << "  &  " << phiL3
                  << "  &  " << pTL4 << "  &  " << etaL4 << "  &  " << phiL4
		  << "  &  " << Run <<  "  &  "  << Event << " \\\\ " <<  endl;
	      
	      ev = "";
	      z1 = "";
	      z2 = "";

	}
      
    }


  runVec2.clear();
  lumiVec2.clear();
  eventVec2.clear();
  
  out.close();
}


void 
PlotHelper::printVsRun(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;

  ULong64_t Run, Event, LumiSect;
  double mass4l;

  TFile *f = new TFile(file,"READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("mass4l",&mass4l);


  ofstream out;
  out.open(outFile.c_str());

  tree->BuildIndex("Run");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();

  for( int i = 0; i < index->GetN(); i++ )
    {

      Long64_t local = tree->LoadTree( index->GetIndex()[i] );
      tree->GetEntry(local);

      bool notDuplicateEvent = true;

      ULong64_t runId = Run;
      ULong64_t eventId = Event;
      ULong64_t lumiId = LumiSect;

    
      for (unsigned int n = 0; n < runVec5.size(); n++)
	{
	  if(runId == runVec5[n] && lumiId == lumiVec5[n] && eventId == eventVec5[n]){notDuplicateEvent = false;}
	}
    

      if (notDuplicateEvent)
	{

	  runVec5.push_back(runId);
	  lumiVec5.push_back(lumiId);
	  eventVec5.push_back(eventId);

	  out << eventId << "  " << lumiId << "  " << runId << "  " << mass4l << endl;

	}
    }


  runVec5.clear();
  lumiVec5.clear();
  eventVec5.clear();

  out.close();

}

void 
PlotHelper::printHistoYield(std::string name,std::vector<TH1F*> hist_4l,std::vector<TH1F*> hist_4mu,std::vector<TH1F*> hist_4e,std::vector<TH1F*> hist_2e2mu, std::vector<TString> histName,double xlow, double xhigh)
{

  using namespace std;

  ofstream out;
  out.open(name.c_str());

  
  for( unsigned int i = 0; i < hist_4l.size(); i++)
    {
      double nEvents_4l = hist_4l[i]->Integral(hist_4l[i]->FindBin(xlow),hist_4l[i]->FindBin(xhigh));
      double nEvents_4mu = hist_4mu[i]->Integral(hist_4mu[i]->FindBin(xlow),hist_4mu[i]->FindBin(xhigh));
      double nEvents_4e = hist_4e[i]->Integral(hist_4e[i]->FindBin(xlow),hist_4e[i]->FindBin(xhigh));
      double nEvents_2e2mu = hist_2e2mu[i]->Integral(hist_2e2mu[i]->FindBin(xlow),hist_2e2mu[i]->FindBin(xhigh));
      
      out << histName[i] << " [" << xlow << "," << xhigh << "] " << endl;
      out << "4l:     " << nEvents_4l << endl;
      out << "4mu:    " << nEvents_4mu << endl;
      out << "4e:     " << nEvents_4e << endl;
      out << "2e2mu:  " << nEvents_2e2mu << endl;
      out << nEvents_4l << "/ " << nEvents_4e << "/ " << nEvents_4mu << "/ " << nEvents_2e2mu << endl;
      out << "------------------ " << endl;
    }
  out.close();
  
}


void 
PlotHelper::printCollectionYield(Collection *myCol, ofstream &out, double mZ1Low, double mZ2Low, double m4lLow, double m4lHigh)
{


  double nEvents_4l=0,nEvents_4e=0,nEvents_4mu=0,nEvents_2e2mu=0;
  double Entries_4l=0,Entries_4e=0,Entries_4mu=0,Entries_2e2mu=0;

  for(unsigned int i = 0; i < myCol->mass4lPair.second.size(); i++)
    {
      if(myCol->massZ1Pair.second[i] < mZ1Low) continue;
      if(myCol->massZ2Pair.second[i] < mZ2Low) continue;
      if(myCol->mass4lPair.second[i] < m4lLow) continue;
      if(myCol->mass4lPair.second[i] > m4lHigh) continue;
      
      nEvents_4l += myCol->weight[i]; Entries_4l+=1.0;
      if(myCol->mass4ePair.second[i] > 0    ) { nEvents_4e    += myCol->weight[i]; Entries_4e+=1.; }
      if(myCol->mass4muPair.second[i] > 0   ) { nEvents_4mu   += myCol->weight[i]; Entries_4mu+=1.; }
      if(myCol->mass2e2muPair.second[i] > 0 ) { nEvents_2e2mu += myCol->weight[i]; Entries_2e2mu+=1.; }

    }

  out << myCol->getName() << " [" << m4lLow << "," << m4lHigh << "] " << endl;
  out << "4l:     " << nEvents_4l << "; Entries " << Entries_4l << endl;
  out << "4mu:    " << nEvents_4mu << "; Entries " << Entries_4mu << endl;
  out << "4e:     " << nEvents_4e << "; Entries " << Entries_4e << endl;
  out << "2e2mu:  " << nEvents_2e2mu << "; Entries " << Entries_2e2mu << endl;
  out << nEvents_4l << "/ " << nEvents_4e << "/ " << nEvents_4mu << "/ " << nEvents_2e2mu << endl;
  out << "------------------ " << endl;
  
  out.close();

}

void 
PlotHelper::printCollectionYield(MergedCollection *myCol, ofstream &out, double mZ1Low, double mZ2Low, double m4lLow, double m4lHigh)
{

  double nEvents_4l=0,nEvents_4e=0,nEvents_4mu=0,nEvents_2e2mu=0;
  double Entries_4l=0,Entries_4e=0,Entries_4mu=0,Entries_2e2mu=0;

  for(unsigned int i = 0; i < myCol->mass4lPair.second.size(); i++)
    {
      if(myCol->massZ1Pair.second[i] < mZ1Low) continue;
      if(myCol->massZ2Pair.second[i] < mZ2Low) continue;
      if(myCol->mass4lPair.second[i] < m4lLow) continue;
      if(myCol->mass4lPair.second[i] > m4lHigh) continue;
      
      nEvents_4l += myCol->weight[i]; Entries_4l+=1.0;
      if(myCol->mass4ePair.second[i] > 0    ) { nEvents_4e    += myCol->weight[i]; Entries_4e+=1.; }
      if(myCol->mass4muPair.second[i] > 0   ) { nEvents_4mu   += myCol->weight[i]; Entries_4mu+=1.; }
      if(myCol->mass2e2muPair.second[i] > 0 ) { nEvents_2e2mu += myCol->weight[i]; Entries_2e2mu+=1.; }

    }

  out << myCol->getName() << " [" << m4lLow << "," << m4lHigh << "] " << endl;
  out << "4l:     " << nEvents_4l << "; Entries " << Entries_4l << endl;
  out << "4mu:    " << nEvents_4mu << "; Entries " << Entries_4mu << endl;
  out << "4e:     " << nEvents_4e << "; Entries " << Entries_4e << endl;
  out << "2e2mu:  " << nEvents_2e2mu << "; Entries " << Entries_2e2mu << endl;
  out << nEvents_4l << "/ " << nEvents_4e << "/ " << nEvents_4mu << "/ " << nEvents_2e2mu << endl;
  out << "------------------ " << endl;

}



std::vector<double> 
PlotHelper::getCollectionYield(Collection *myCol, double mZ1Low, double mZ2Low, double m4lLow, double m4lHigh)
{
  double nEvents_4l=0,nEvents_4e=0,nEvents_4mu=0,nEvents_2e2mu=0;
  double Entries_4l=0,Entries_4e=0,Entries_4mu=0,Entries_2e2mu=0;

  for(unsigned int i = 0; i < myCol->mass4lPair.second.size(); i++)
    {
      if(myCol->massZ1Pair.second[i] < mZ1Low) continue;
      if(myCol->massZ2Pair.second[i] < mZ2Low) continue;
      if(myCol->mass4lPair.second[i] < m4lLow) continue;
      if(myCol->mass4lPair.second[i] > m4lHigh) continue;
      
      nEvents_4l += myCol->weight[i]; Entries_4l+=1.0;
      if(myCol->mass4ePair.second[i] > 0    ) { nEvents_4e    += myCol->weight[i]; Entries_4e+=1.; }
      if(myCol->mass4muPair.second[i] > 0   ) { nEvents_4mu   += myCol->weight[i]; Entries_4mu+=1.; }
      if(myCol->mass2e2muPair.second[i] > 0 ) { nEvents_2e2mu += myCol->weight[i]; Entries_2e2mu+=1.; }

    }


  std::vector<double> mv;
  mv.push_back(nEvents_4l);
  mv.push_back(nEvents_4mu);
  mv.push_back(nEvents_4e);
  mv.push_back(nEvents_2e2mu);
  mv.push_back(Entries_4l);
  mv.push_back(Entries_4mu);
  mv.push_back(Entries_4e);
  mv.push_back(Entries_2e2mu);

  return mv;

  /*
  out << myCol->getName() << " [" << m4lLow << "," << m4lHigh << "] " << endl;
  out << "4l:     " << nEvents_4l << endl;
  out << "4mu:    " << nEvents_4mu << endl;
  out << "4e:     " << nEvents_4e << endl;
  out << "2e2mu:  " << nEvents_2e2mu << endl;
  out << nEvents_4l << "/ " << nEvents_4e << "/ " << nEvents_4mu << "/ " << nEvents_2e2mu << endl;
  out << "------------------ " << endl;
  
  out.close();
  */
}

std::vector<double>
PlotHelper::getCollectionYield(MergedCollection *myCol, double mZ1Low, double mZ2Low, double m4lLow, double m4lHigh)
{

  double nEvents_4l=0,nEvents_4e=0,nEvents_4mu=0,nEvents_2e2mu=0;
  double Entries_4l=0,Entries_4e=0,Entries_4mu=0,Entries_2e2mu=0;

  for(unsigned int i = 0; i < myCol->mass4lPair.second.size(); i++)
    {
      if(myCol->massZ1Pair.second[i] < mZ1Low) continue;
      if(myCol->massZ2Pair.second[i] < mZ2Low) continue;
      if(myCol->mass4lPair.second[i] < m4lLow) continue;
      if(myCol->mass4lPair.second[i] > m4lHigh) continue;
      
      nEvents_4l += myCol->weight[i]; Entries_4l+=1.0;
      if(myCol->mass4ePair.second[i] > 0    ) { nEvents_4e    += myCol->weight[i]; Entries_4e+=1.; }
      if(myCol->mass4muPair.second[i] > 0   ) { nEvents_4mu   += myCol->weight[i]; Entries_4mu+=1.; }
      if(myCol->mass2e2muPair.second[i] > 0 ) { nEvents_2e2mu += myCol->weight[i]; Entries_2e2mu+=1.; }

    }


  std::vector<double> mv;
  mv.push_back(nEvents_4l);
  mv.push_back(nEvents_4mu);
  mv.push_back(nEvents_4e);
  mv.push_back(nEvents_2e2mu);
  mv.push_back(Entries_4l);
  mv.push_back(Entries_4mu);
  mv.push_back(Entries_4e);
  mv.push_back(Entries_2e2mu);

  /*
  out << myCol->getName() << " [" << m4lLow << "," << m4lHigh << "] " << endl;
  out << "4l:     " << nEvents_4l << endl;
  out << "4mu:    " << nEvents_4mu << endl;
  out << "4e:     " << nEvents_4e << endl;
  out << "2e2mu:  " << nEvents_2e2mu << endl;
  out << nEvents_4l << "/ " << nEvents_4e << "/ " << nEvents_4mu << "/ " << nEvents_2e2mu << endl;
  out << "------------------ " << endl;
  */

  return mv;

}





void 
PlotHelper::rescale(TH1F *hist)
{


   for( int i = 1; i < hist->GetNbinsX(); i++ )
    {

      double y = hist->GetBinContent(i);
      double x = hist->GetBinCenter(i);

      //std::cout << x << "  " << y << std::endl;

      double scaling  = 1;

      //scaling = (0.8 + 3.95*exp( -pow(x-91.2,2)/pow(2.6,2)))*(1+( 60/x + 0.8 * exp(-(x-92)*(x-92)/196) )*( 1 - TMath::Erf((x-178)/10) )) / (1 + (60/x + 0.8 + exp(-pow(x-92,2)/14/14) )*(1 - TMath::Erf( (x-178)/10 ) ));

      scaling = 1+( 60/x + 0.8 * exp(-(x-92)*(x-92)/196) )*( 1 - TMath::Erf((x-178)/10) ); 


      std::cout << scaling << std::endl;

      hist->SetBinContent(i,y*scaling);
      
    }



}





double 
PlotHelper::correctionFactor(double mH, int channel)
{

  using namespace std;

  double A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0;


  if( channel == 1 )
    {
      A = 0.623194;
      B = 0.00762586;
      C = -6.6463e-05;
      D = 3.11407e-07;
      E = -8.47919e-10;
      F = 1.34617e-12;
      G = -1.15716e-15;
      H = 4.16432e-19;
    }
  else if( channel == 2 )
    {
      A = 0.232153;
      B = 0.0143873;
      C = -0.000108814;
      D = 4.18026e-07;
      E = -8.67393e-10;
      F = 9.2644e-13;
      G = -3.99454e-16;
      H = 0;
    }
  else if( channel == 3 )
    {
      A = 0.627437;
      B = 0.00762152;
      C = -6.6484e-05;
      D = 3.11375e-07;
      E = -8.47949e-10;
      F = 1.34617e-12;
      G = -1.15707e-15;
      H = 4.16702e-19;
    }
  else{

    cout << "-----------------------------------------" << endl
         << "In correctionFactor(double mH, int channel): channel " << channel << " not found!" << endl
         << "-----------------------------------------" << endl;
  }

  double eff = A + B*mH + C*mH*mH + D*mH*mH*mH + E*pow(mH,4) + F*pow(mH,5) + G*pow(mH,6) + H*pow(mH,7);


  return eff;

}



void 
PlotHelper::rescale(TH1F *hist, int channel, double n_tot, double n_4mu, double n_4e, double n_2e2mu)
{


  for( int i = 1; i < hist->GetNbinsX(); i++ )
    {

      double y = hist->GetBinContent(i);
      double x = hist->GetBinCenter(i);

      //std::cout << x << "  " << y << std::endl;

      double scaling  = 1;

      //4l
      if( channel == 0 )
	{
	  double frac_4mu = (double)n_4mu/n_tot;//(n_4mu+n_4e+n_2e2mu);
	  double frac_4e = (double)n_4e/n_tot;//(n_4mu+n_4e+n_2e2mu);
          double frac_2e2mu = (double)n_2e2mu/n_tot;//(n_4mu+n_4e+n_2e2mu);

	  scaling = frac_4mu*correctionFactor(x,1) + frac_4e*correctionFactor(x,2) + frac_2e2mu*correctionFactor(x,3);
	  //if(isZZ)scaling *= 1+( 60/x + 0.8 * exp(-(x-92)*(x-92)/196) )*( 1 - TMath::Erf((x-178)/10) ); 


	}
      else{ 
	scaling = correctionFactor(x,channel); 
	
	    
      }
      
      //cout << i << "  " << y << "  " << scaling << "  " << y*scaling << endl;
      hist->SetBinContent(i,y*scaling);
      
    }
  
}


void 
PlotHelper::drawPlot(THStack *stack, TGraphAsymmErrors *gr, TLegend *leg, TLatex *latex,TString saveName,TString xTitle, TString yTitle, double sqrts, 
		     double lumi, double xLow, double xHigh, double yMax)
{

  TCanvas *c = new TCanvas("c","c",canvasX,canvasY);
  c->cd();
  stack->Draw();
  setStackStyle(stack,xTitle,yTitle);
  gr->Draw("pSAME");
  if(yMax != 999) stack->SetMaximum(yMax);
  if(xLow != 999 && xHigh != 999) stack->GetXaxis()->SetRangeUser(xLow,xHigh);
  leg->Draw("SAME");
  gPad->SetTicks(1,1);
  c->Update();
  latex->DrawLatex(0.145,0.96,Form("CMS Preliminary                          #sqrt{s} = %.0f TeV, L = %.2f fb^{-1}",sqrts,lumi));  
  c->SaveAs(saveName);

  delete c;

}


void 
PlotHelper::draw7p8Plot(THStack *stack, TGraphAsymmErrors *gr, TLegend *leg, TLatex *latex,TString saveName,TString xTitle, TString yTitle, 
			double lumi_7, double lumi_8, double xLow, double xHigh, double yMax)
{
  TCanvas *c = new TCanvas("c","c",canvasX,canvasY);
  c->cd();
  if(yMax != 999) stack->SetMaximum(yMax);
  stack->Draw();
  setStackStyle(stack,xTitle,yTitle);
  if(xLow != 999 && xHigh != 999) stack->GetXaxis()->SetRangeUser(xLow,xHigh);
  gr->Draw("pSAME");
  leg->Draw("SAME");
  gPad->SetTicks(1,1);
  c->Update();
  //latex->DrawLatex(0.145,0.95,Form("                                                   #sqrt{s} = 7 TeV, L = %.2f fb^{-1}",lumi_7));
  latex->DrawLatex(0.145,0.965,Form("CMS Preliminary   #sqrt{s} = 7 TeV, L = %.2f fb^{-1}  #sqrt{s} = 8 TeV, L = %.1f fb^{-1}",lumi_7,lumi_8));
  c->SaveAs(saveName);

  delete c;

}

void 
PlotHelper::draw7p8Plot(THStack *stack,TH1F* gr, TLegend *leg, TLatex *latex,TString saveName,TString xTitle, TString yTitle, 
			double lumi_7, double lumi_8, double xLow, double xHigh, double yMax)
{

  TCanvas *c = new TCanvas("c","c",canvasX,canvasY);
  c->cd();
  if(yMax != 999) stack->SetMaximum(yMax);
  stack->Draw();
  if(xLow != 999 && xHigh != 999) stack->GetXaxis()->SetRangeUser(xLow,xHigh);
  setStackStyle(stack,xTitle,yTitle);
  gr->Draw("e1pSAME");
  leg->Draw("SAME");
  gPad->SetTicks(1,1);
  c->Update();
  //latex->DrawLatex(0.145,0.95,Form("                                                   #sqrt{s} = 7 TeV, L = %.2f fb^{-1}",lumi_7));
  latex->DrawLatex(0.145,0.965,Form("CMS Preliminary   #sqrt{s} = 7 TeV, L = %.2f fb^{-1}  #sqrt{s} = 8 TeV, L = %.1f fb^{-1}",lumi_7,lumi_8));
  c->SaveAs(saveName);

  delete c;

}




void PlotHelper::makeArrayFromVector(std::vector<double> myVec, double myArray[])
{

  for( unsigned int i = 0; i < myVec.size(); i++)
    {
      myArray[i] = myVec[i];
    }

}


void PlotHelper::PaintOverflow(TH1F* h,TH1F* &htmp,Color_t color,Style_t fillStyle,DrawStyle draw, Style_t markerStyle)
{

  // This function paint the histogram h with an extra bin for overflows
  TString name  = (TString)h->GetName()+"new";
  TString title = (TString)h->GetTitle();
  Int_t nx    = h->GetNbinsX()+1;
  Double_t x1 = h->GetBinLowEdge(1);
  Double_t bw = h->GetBinWidth(nx);
  Double_t x2 = h->GetBinLowEdge(nx)+bw;  // Book a temporary histogram having ab extra bin for overflows
  htmp  = new TH1F(name, title, nx, x1, x2);    // Fill the new hitogram including the extra bin for overflows
  setHistProperties(htmp,color,fillStyle,draw,markerStyle);
 
 for (Int_t i=1; i<=nx; i++) {
    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
  } 
 // Fill the underflows
 htmp->Fill(x1-1, h->GetBinContent(0));    // Restore the number of entries
 htmp->SetEntries(h->GetEntries());    // Draw the temporary histogram
 


  //htmp->Draw();
  //TText *t = new TText(x2-bw/2,h->GetBinContent(nx),"Overflow");
  //t->SetTextAngle(90);
  //t->SetTextAlign(12);
  //t->SetTextSize(0.03);
  //t->Draw();
}


TH1F* PlotHelper::setOverflowBin(TH1F* histo, float xmax) 
{ 
  TH1F* h = (TH1F*)histo->Clone();
  int ovfBin = h->FindBin(xmax);
  
  float ovfTot=0; 
  float ovfTotErr2=0;
  for (int i=ovfBin; i<h->GetNbinsX()+1; ++i) {
    ovfTot+=h->GetBinContent(i); 
    ovfTotErr2 += pow(h->GetBinError(i),2);
    h->SetBinContent(i,0);
    //h->SetBinError(i,0);
  }
  h->SetBinContent(ovfBin,ovfTot);
  //h->SetBinError(ovfBin,sqrt(ovfTotErr2));
  return h;
}


double PlotHelper::findHighestBin(TH1F* histo)
{

  double highestN=0;

  for (int i=1; i<histo->GetNbinsX(); i++) 
    {
      double binCont = histo->GetBinContent(i);
      if(binCont>highestN) highestN=binCont;
      //cout << binCont << endl;
    }

  return highestN;

}


#endif
