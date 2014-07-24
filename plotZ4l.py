#!/usr/bin/python
import sys, os, pwd, commands
import optparse, shlex, re
import math
from ROOT import *
import ROOT
from array import array


  

def parseOptions():
  
  usage = ('usage: %prog [options] \n'
           + '%prog -h for help')
  parser = optparse.OptionParser(usage)
  
  parser.add_option('-b', action='store_true', dest='noX', default=True ,help='no X11 windows')
  parser.add_option('-d', '--dir', dest='dir', type='string', default='Zto4LPlots' ,help='output directory')
  parser.add_option('-s', '--sqrts', dest='sqrts', type='int', default=0 ,help='sqrts = 7,8,0')
  parser.add_option('-l', '--xlow', dest='xLow', type='float', default=70.5 ,help='Low mass cut off')
  parser.add_option('-t', '--xhigh', dest='xHigh', type='float', default=181.5 ,help='High mass cut off')
  parser.add_option('--dataMC', action='store_true', dest='dataMCWeight', default=True ,help='Include data/MC scaling')
  parser.add_option('-z', '--massz2', dest='MASSZ2', type='float', default=12 ,help='Mass of Z2')
  parser.add_option('--ptmu', dest='PTMU', type='float', default=5 ,help='pT of muons')
  parser.add_option('--ptel', dest='PTEL', type='float', default=7 ,help='pT of electrons')
  
  # store options and arguments as global variables
  global opt, args
  (opt, args) = parser.parse_args()
  


def getVectors(file,isMC,lumi,mass4mu,weight4mu,mass4e,weight4e,mass2e2mu,weight2e2mu,mass2mu2e,weight2mu2e,mass2e2L,weight2e2L,mass2mu2L,weight2mu2L,mass2L2mu,weight2L2mu,mass2L2e,weight2L2e):
  global opt, args
  
  tfile = ROOT.TFile(file,"READ")
  if isMC:
    tree = tfile.Get("passedEvents_dataMC")
  else:
    if lumi == 5.051:
      tree = tfile.Get("AnaAfterHlt/passedEvents")
    else:
      tree = tfile.Get("AnaAfterHlt/passedEvents")
            
      
    
  if isMC:
    normHist = tfile.Get("AnaAfterHlt/nEvents")
    norm = normHist.GetBinContent(1)
  else:
    norm = 1

  runVec = []
  lumiVec = []
  eventVec = []

  for i in range( tree.GetEntries() ):
    tree.GetEntry(i)

    notDuplicateEvent = True

    if not isMC:
      for n in xrange(len(runVec)):
        if tree.Run == runVec[n] and tree.LumiSect == lumiVec[n] and tree.Event == eventVec[n]:
          notDuplicateEvent = False

    if not notDuplicateEvent: continue
    runVec.append(tree.Run)
    lumiVec.append(tree.LumiSect)
    eventVec.append(tree.Event)
        
    #if not tree.passedZ4lSelection: continue
    if not tree.massZ2 > opt.MASSZ2: continue
    if not tree.massZ1 > 40: continue
    if abs(tree.idL1) == 13 and tree.pTL1 < opt.PTMU: continue
    if abs(tree.idL2) == 13 and tree.pTL2 < opt.PTMU: continue
    if abs(tree.idL3) == 13 and tree.pTL3 < opt.PTMU: continue
    if abs(tree.idL4) == 13 and tree.pTL4 < opt.PTMU: continue

    if abs(tree.idL1) == 11 and tree.pTL1 < opt.PTEL: continue
    if abs(tree.idL2) == 11 and tree.pTL2 < opt.PTEL: continue
    if abs(tree.idL3) == 11 and tree.pTL3 < opt.PTEL: continue
    if abs(tree.idL4) == 11 and tree.pTL4 < opt.PTEL: continue



    weight = 1
    if isMC:
      if opt.dataMCWeight:
        weight = tree.scaleWeight*tree.eventWeight*tree.dataMC_weight*lumi*1000/norm        
      else:
        weight = tree.scaleWeight*tree.eventWeight*lumi*1000/norm
        
        
    if abs(tree.idL1) == 11 and abs(tree.idL3) == 11:
      mass4e.append(tree.mass4l)
      weight4e.append(weight)
    if abs(tree.idL1) == 13 and abs(tree.idL3) == 13:
      mass4mu.append(tree.mass4l)
      weight4mu.append(weight)
    if abs(tree.idL1) == 11 and abs(tree.idL3) == 13:
      mass2e2mu.append(tree.mass4l)
      weight2e2mu.append(weight)
    if abs(tree.idL1) == 13 and abs(tree.idL3) == 11:
      mass2mu2e.append(tree.mass4l)
      weight2mu2e.append(weight)
    if abs(tree.idL1) == 13 and (abs(tree.idL3) == 11 or abs(tree.idL3) == 13):
      mass2mu2L.append(tree.mass4l)
      weight2mu2L.append(weight)
    if abs(tree.idL1) == 11 and (abs(tree.idL3) == 11 or abs(tree.idL3) == 13):
      mass2e2L.append(tree.mass4l)
      weight2e2L.append(weight)
    if abs(tree.idL3) == 13 and (abs(tree.idL1) == 11 or abs(tree.idL1) == 13):
      mass2L2mu.append(tree.mass4l)
      weight2L2mu.append(weight)
    if abs(tree.idL3) == 11 and (abs(tree.idL1) == 11 or abs(tree.idL1) == 13):
      mass2L2e.append(tree.mass4l)
      weight2L2e.append(weight)

      


def getCount(listMass, listWeight, xLow, xHigh, isMC):

  count = 0
  
  for i in xrange(len(listMass)):
    if listMass[i] < xLow: continue
    if listMass[i] > xHigh: continue

    if isMC:
      count+=listWeight[i]
    else:
      count+=1

  return count



def getAsymErr(h,gr):

  n = h.GetNbinsX()
  x = []
  y = []
  exl = []
  eyl = []
  exh = []
  eyh = []
  
  q = (1-0.6827)/2
  
  for i in xrange(n):

    x.append(h.GetXaxis().GetBinCenter(i+1))
    N = h.GetBinContent(i+1)
    exl.append(0)
    exh.append(0)

    if N == 0:
      eyl.append(0)
      eyh.append(0)
      y.append(-1)
      continue
    
    y.append(N)
    eyl.append(N-ROOT.Math.chisquared_quantile_c(1-q,2*N)/2)
    eyh.append(ROOT.Math.chisquared_quantile_c(q,2*(N+1))/2-N)


  X = array( 'f', x)
  Y = array( 'f', y)
  Exl = array( 'f', exl)
  Eyl = array( 'f', eyl)
  Exh = array( 'f', exh)
  Eyh = array( 'f', eyh)
  
  gr = ROOT.TGraphAsymmErrors(n,X,Y,Exl,Exh,Eyl,Eyh)
  gr.SetMarkerStyle(20)
  gr.SetMarkerColor(kBlack)
  gr.SetMarkerSize(1.1)






def makePlot(listMassMC, listWeightMC, listData,chan):
  global opt, args

  binwidth = 3
  nBins = (opt.xHigh-opt.xLow)/binwidth

  scaleBG = 1
  if opt.MASSZ2 > 8:
    scaleBG = 0.8
  if opt.MASSZ2 > 10:
    scaleBG = 0.65


  chanPlot = ''
  if chan == '4mu':
    chanPlot = '4#mu'
    nBgEvents_7 = 0.8
    nBgEvents_8 = 3.0
  if chan == '4e':
    chanPlot = '4e'
    nBgEvents_7 = 1.6
    nBgEvents_8 = 6.2
  if chan == '2e2mu':
    chanPlot = '2e2#mu'
    nBgEvents_7 = 1.0
    nBgEvents_8 = 3.0
  if chan == '2mu2e':
    chanPlot = '2#mu2e'
    nBgEvents_7 = 1.9
    nBgEvents_8 = 6.0
  if chan == '2e2muc':
    chanPlot = '2e2#mu'
    nBgEvents_7 = 2.9
    nBgEvents_8 = 9.0
  if chan == '2e2L':
    chanPlot = '2e2l'
    nBgEvents_7 = 2.6
    nBgEvents_8 = 9.2
  if chan == '2mu2L':
    chanPlot = '2#mu2l'
    nBgEvents_7 = 2.7
    nBgEvents_8 = 9.0
  if chan == '4L':
    chanPlot = '4l'
    nBgEvents_7 = 5.1
    nBgEvents_8 = 18.2
  if chan == '2L2mu':
    chanPlot = '2l2#mu'
    nBgEvents_7 = 1.8
    nBgEvents_8 = 6.0
  if chan == '2L2e':
    chanPlot = '2l2e'
    nBgEvents_7 = 2.5
    nBgEvents_8 = 10.0




  #### Z+X ####
  l_mu_7 = 143
  l_sigma_7 = 19
  l_mu_8 = 144
  l_sigma_8 = 19


  dummyLandau7 = ROOT.TH1F("dummy7","",10000,0,1000);
  dummyLandau = dummyLandau7.Clone()
  for i in xrange(10000):
    dummyLandau7.Fill(ROOT.gRandom.Landau(l_mu_7,l_sigma_7))

  dummyLandau8 = ROOT.TH1F("dummy8","",10000,0,1000);
  for i in xrange(10000):
    dummyLandau8.Fill(ROOT.gRandom.Landau(l_mu_8,l_sigma_8))


  dummyLandau7.Scale(5.051/dummyLandau7.Integral("width"))
  dummyLandau8.Scale(19.8/dummyLandau8.Integral("width"))

  if opt.sqrts == 7:
    dummyLandau.Add(dummyLandau7)
    dummyLandau.Scale(1/dummyLandau.Integral("width"))
  if opt.sqrts == 8:
    dummyLandau.Add(dummyLandau8)
    dummyLandau.Scale(1/dummyLandau.Integral("width"))
  if opt.sqrts == 0:
    dummyLandau.Add(dummyLandau7)
    dummyLandau.Add(dummyLandau8)
    dummyLandau.Scale(1/dummyLandau.Integral("width"))

  scaleFactorRB = dummyLandau.Integral(dummyLandau.GetXaxis().FindBin(opt.xLow),dummyLandau.GetXaxis().FindBin(opt.xHigh))
  print "ScaleFact ",scaleFactorRB
  
  hMC = ROOT.TH1F("histMC",";mass_{"+chanPlot+"} (GeV); N Events / 2 GeV",nBins,opt.xLow,opt.xHigh)
  hRB7 = hMC.Clone()
  hRB8 = hMC.Clone()
  hRB = hMC.Clone()
  hData = ROOT.TH1F("histData",";mass_{"+chanPlot+"} (GeV); N Events / 2 GeV",nBins,opt.xLow,opt.xHigh)
  hDataG = ROOT.TGraphAsymmErrors()
  hData.SetMarkerStyle(20)
  hData.SetMarkerSize(1.1)
  

  for i in xrange(10000):
    hRB7.Fill(ROOT.gRandom.Landau(l_mu_7,l_sigma_7))
    
  hRB7.Scale(nBgEvents_7*scaleFactorRB*scaleBG/hRB7.Integral())
  hRB7.SetFillColor(kGreen-5)
  hRB7.SetLineWidth(2)

  for i in xrange(10000):
    hRB8.Fill(ROOT.gRandom.Landau(l_mu_8,l_sigma_8))
    
  hRB8.Scale(nBgEvents_8*scaleFactorRB*scaleBG/hRB8.Integral())
  hRB8.SetFillColor(kGreen-5)
  hRB8.SetLineWidth(2)

  if opt.sqrts == 7:
    hRB.Add(hRB7)
  if opt.sqrts == 8:
    hRB.Add(hRB8)
  if opt.sqrts == 0:
    hRB.Add(hRB7)
    hRB.Add(hRB8)

  hRB.SetFillColor(kGreen-5)
  hRB.SetLineWidth(2)
  hRB.SetLineColor(kBlack)

  hMC.SetFillColor(kAzure-9)
  hMC.SetLineWidth(2)
  hMC.SetLineColor(kBlack)

  ROOT.gStyle.SetOptStat(0000)
  ROOT.gStyle.SetPadLeftMargin(0.14)
  ROOT.gStyle.SetOptTitle(0)
  ROOT.gStyle.SetPadTopMargin(0.05)
  ROOT.gStyle.SetPadBottomMargin(0.13)

  for i in range(len(listMassMC)):
    hMC.Fill(listMassMC[i],listWeightMC[i])
  for i in range(len(listData)):
    hData.Fill(listData[i],1)


  maxY = 0

  for i in range(hMC.GetNbinsX()):
    if hMC.GetBinContent(i+1) > maxY:
      maxY = hMC.GetBinContent(i+1)

  for i in range(hData.GetNbinsX()):
    if hData.GetBinContent(i+1) > maxY:
      maxY = hData.GetBinContent(i+1)

  maxY = maxY*1.25


  if len(listData) > 0:
    getAsymErr(hData,hDataG)

  stack = ROOT.THStack()
  stack.Add(hRB)
  stack.Add(hMC)
  #######################
  countSmall = getCount(listMassMC, listWeightMC, 86.2, 96.2, True)
  countSmallData = getCount(listData, listData, 86.2, 96.2, False)
  countLarge = getCount(listMassMC, listWeightMC, 81.2, 101.2, True)
  countLargeData = getCount(listData, listData, 81.2, 101.2, False)

  tmpRB = hRB.Integral(hRB.GetXaxis().FindBin(86.2),hRB.GetXaxis().FindBin(96.2))
  countSmall+=tmpRB
  tmpRB = hRB.Integral(hRB.GetXaxis().FindBin(81.2),hRB.GetXaxis().FindBin(101.2))
  countLarge+=tmpRB
  #######################

  leg = ROOT.TLegend(0.2,0.6,0.4,0.8)
  leg.SetTextSize(0.04)
  leg.SetTextFont(42)
  leg.SetFillColor(kWhite)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)

  leg.AddEntry(hMC,"ZZ/Z#gamma*","F")
  leg.AddEntry(hRB,"Z+X","F")
  leg.AddEntry(hData,"Data","P")

  CP = ROOT.TLatex()
  CP.SetNDC(kTRUE)
  CP.SetTextSize(0.035)
  if opt.sqrts == 0:
    CP.SetTextSize(0.03)
  #CP.SetTextAlign(31)
  CP.SetTextFont(42)
  CP.SetTextAlign(11)


  c = ROOT.TCanvas("c","c",750,750)
  c.cd()
  gPad.SetTicks(1,1)
  stack.Draw()
  ytitleOffset = 1.36
  xtitleOffset = 1.18
  labelSize = 0.05
  titleSize = 0.05
  stack.GetXaxis().SetTitleOffset(xtitleOffset)
  stack.GetYaxis().SetTitleOffset(ytitleOffset)
  stack.GetXaxis().SetLabelSize(labelSize)
  stack.GetYaxis().SetLabelSize(labelSize)
  stack.GetXaxis().SetTitleSize(titleSize)
  stack.GetYaxis().SetTitleSize(titleSize)
  stack.GetXaxis().SetTitle("m_{"+chanPlot+"} (GeV)")
  stack.GetYaxis().SetTitle("N Events / "+str(binwidth)+" GeV")
  hData.Draw("e1pSAME")
  leg.Draw()
  
  if chan == '4e':
    stack.SetMaximum(20)
  elif chan == '2mu2e':
    stack.SetMaximum(35)
  elif chan == '2L2mu':
    stack.SetMaximum(350)
  elif chan == '2L2e':
    stack.SetMaximum(50)
  elif chan == '2mu2L':
    stack.SetMaximum(260)
  elif chan == '2e2L':
    stack.SetMaximum(140)
  else:
    stack.SetMaximum(maxY)


  c.Update()
  if opt.sqrts == 7:
    CP.DrawLatex(0.145,0.965,"CMS Preliminary                          #sqrt{s} = 7 TeV, L = 5.1 fb^{-1}")
  if opt.sqrts == 8:
    CP.DrawLatex(0.145,0.965,"CMS Preliminary                          #sqrt{s} = 8 TeV, L = 19.8 fb^{-1}")
  if opt.sqrts == 0:
    CP.DrawLatex(0.145,0.965,"CMS Preliminary   #sqrt{s} = 7 TeV, L = 5.1 fb^{-1}  #sqrt{s} = 8 TeV, L = 19.8 fb^{-1}")

  CP.DrawLatex(0.2,0.55,"Exp[86.2,96.2]: {0:.2f}".format(countSmall))
  CP.DrawLatex(0.2,0.5,"Obs[86.2,96.2]: {0}".format(countSmallData))
  CP.DrawLatex(0.2,0.45,"Exp[81.2,101.2]: {0:.2f}".format(countLarge))
  CP.DrawLatex(0.2,0.4,"Obs[81.2,101.2]: {0}".format(countLargeData))

  ss = opt.sqrts
  if ss == 0:
    ss = '7p8'
  c.SaveAs(opt.dir+"/Z4l_Mass_"+str(opt.xLow)+"-"+str(opt.xHigh)+"_"+chan+"_"+str(ss)+"TeV.eps")
  c.SaveAs(opt.dir+"/Z4l_Mass_"+str(opt.xLow)+"-"+str(opt.xHigh)+"_"+chan+"_"+str(ss)+"TeV.png")
  c.SaveAs(opt.dir+"/Z4l_Mass_"+str(opt.xLow)+"-"+str(opt.xHigh)+"_"+chan+"_"+str(ss)+"TeV.pdf")
  c.SaveAs(opt.dir+"/Z4l_Mass_"+str(opt.xLow)+"-"+str(opt.xHigh)+"_"+chan+"_"+str(ss)+"TeV.C")
  c.SaveAs(opt.dir+"/Z4l_Mass_"+str(opt.xLow)+"-"+str(opt.xHigh)+"_"+chan+"_"+str(ss)+"TeV.root")

  
  
  
if __name__ == "__main__":
  parseOptions()
  ROOT.gROOT.ProcessLine(".L include/tdrstyle.cc")
  
  path7 = "/scratch/osghpc/snowball/UF_HZZ4L_Analysis/FullAnalysis2012/Histogramming/rootFiles_Legacy_dataMC/"
  path8 = "/scratch/osghpc/snowball/UF_HZZ4L_Analysis/FullAnalysis2012/Histogramming_8TeV/rootFiles_Legacy_dataMC/"
  #path7 = "../Histogramming/Moriond/rootFiles_Moriond_dataMC/"
  #path8 = "../Histogramming_8TeV/Moriond/rootFiles_Moriond_dataMC/"

  if opt.sqrts == 7:
    files = [path7+"ZZ_4mu.root",path7+"ZZ_4e.root",path7+"ZZ_2e2mu.root",path7+"ZZto2e2tau.root",path7+"ZZto2mu2tau.root",path7+"ZZto4tau.root"]
    lumi = [5.051, 5.051, 5.051, 5.051, 5.051, 5.051]
  if opt.sqrts == 8:
    files = [path8+"ZZ_4mu.root",path8+"ZZ_4e.root",path8+"ZZ_2e2mu.root",path8+"ZZto2e2tau.root",path8+"ZZto2mu2tau.root",path8+"ZZto4tau.root"]
    lumi = [19.8, 19.8, 19.8, 19.8, 19.8, 19.8]
  if opt.sqrts == 0:
    files = [path7+"ZZ_4mu.root",path7+"ZZ_4e.root",path7+"ZZ_2e2mu.root",path7+"ZZto2e2tau.root",path7+"ZZto2mu2tau.root",path7+"ZZto4tau.root",
             path8+"ZZ_4mu.root",path8+"ZZ_4e.root",path8+"ZZ_2e2mu.root",path8+"ZZto2e2tau.root",path8+"ZZto2mu2tau.root",path8+"ZZto4tau.root"]
    lumi = [5.051, 5.051, 5.051, 5.051, 5.051, 5.051, 19.8, 19.8, 19.8, 19.8, 19.8, 19.8]


  mass4mu = []
  mass4e = []
  mass2e2mu = []
  mass2mu2e = []
  mass2e2L = []
  mass2mu2L = []
  mass2L2e = []
  mass2L2mu = []
  weight4mu = []
  weight4e = []
  weight2e2mu = []
  weight2mu2e = []
  weight2e2L = []
  weight2mu2L = []
  weight2L2e = []
  weight2L2mu = []

  data4mu = []
  data4e = []
  data2e2mu = []
  data2mu2e = []
  data2e2L = []
  data2mu2L = []
  data2L2mu = []
  data2L2e = []
  dweight4mu = []
  dweight4e = []
  dweight2e2mu = []
  dweight2mu2e = []
  dweight2e2L = []
  dweight2mu2L = []
  dweight2L2e = []
  dweight2L2mu = []

  
  counter = 0
  for file in files:
    getVectors(file,True,lumi[counter],mass4mu,weight4mu,mass4e,weight4e,mass2e2mu,weight2e2mu,mass2mu2e,weight2mu2e,mass2e2L,weight2e2L,mass2mu2L,weight2mu2L,mass2L2mu,weight2L2mu,mass2L2e,weight2L2e)
    counter+=1
      

  if opt.sqrts == 7:
    getVectors(path7+"Data.root",False,5.051,data4mu,dweight4mu,data4e,dweight4e,data2e2mu,dweight2e2mu,data2mu2e,dweight2mu2e,data2e2L,dweight2e2L,data2mu2L,dweight2mu2L,data2L2mu,dweight2L2mu,data2L2e,dweight2L2e)
  if opt.sqrts == 8:
    getVectors(path8+"Data.root",False,1,data4mu,dweight4mu,data4e,dweight4e,data2e2mu,dweight2e2mu,data2mu2e,dweight2mu2e,data2e2L,dweight2e2L,data2mu2L,dweight2mu2L,data2L2mu,dweight2L2mu,data2L2e,dweight2L2e)
  if opt.sqrts == 0:
    getVectors(path7+"Data.root",False,5.051,data4mu,dweight4mu,data4e,dweight4e,data2e2mu,dweight2e2mu,data2mu2e,dweight2mu2e,data2e2L,dweight2e2L,data2mu2L,dweight2mu2L,data2L2mu,dweight2L2mu,data2L2e,dweight2L2e)
    getVectors(path8+"Data.root",False,1,data4mu,dweight4mu,data4e,dweight4e,data2e2mu,dweight2e2mu,data2mu2e,dweight2mu2e,data2e2L,dweight2e2L,data2mu2L,dweight2mu2L,data2L2mu,dweight2L2mu,data2L2e,dweight2L2e)
    

    
  
  makePlot(mass4mu, weight4mu, data4mu,'4mu')
  makePlot(mass4e, weight4e, data4e,'4e')
  makePlot(mass2e2mu, weight2e2mu, data2e2mu,'2e2mu')
  makePlot(mass2mu2e, weight2mu2e, data2mu2e,'2mu2e')
  makePlot(mass2mu2L, weight2mu2L, data2mu2L,'2mu2L')
  makePlot(mass2e2L, weight2e2L, data2e2L,'2e2L')
  makePlot(mass2L2mu, weight2L2mu, data2L2mu,'2L2mu')
  makePlot(mass2L2e, weight2L2e, data2L2e,'2L2e')

  mass4l = []
  weight4l = []
  data4l = []

  for i in xrange(len(mass4mu)):
    mass4l.append(mass4mu[i])
    weight4l.append(weight4mu[i])
  for i in xrange(len(mass4e)):
    mass4l.append(mass4e[i])
    weight4l.append(weight4e[i])
  for i in xrange(len(mass2e2mu)):
    mass4l.append(mass2e2mu[i])
    weight4l.append(weight2e2mu[i])
  for i in xrange(len(mass2mu2e)):
    mass4l.append(mass2mu2e[i])
    weight4l.append(weight2mu2e[i])

  for i in xrange(len(data4mu)):
    data4l.append(data4mu[i])
  for i in xrange(len(data4e)):
    data4l.append(data4e[i])
  for i in xrange(len(data2e2mu)):
    data4l.append(data2e2mu[i])
  for i in xrange(len(data2mu2e)):
    data4l.append(data2mu2e[i])
  
  makePlot(mass4l, weight4l, data4l,'4L')
  
