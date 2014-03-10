#!/usr/bin/python
import sys, os, pwd, commands
import optparse, shlex, re
import math
from ROOT import *
import ROOT
from array import array


def parseOptions():
  
  usage = ('usage: %prog [options] datasetList\n'
           + '%prog -h for help')
  parser = optparse.OptionParser(usage)
  
  parser.add_option('-b', action='store_true', dest='noX', default=True ,help='no X11 windows')
  parser.add_option('--xlow', dest='xLow', type='float', default=86 ,help='Low mass cut off')
  parser.add_option('--xhigh', dest='xHigh', type='float', default=96 ,help='High mass cut off')
  
  # store options and arguments as global variables
  global opt, args
  (opt, args) = parser.parse_args()



def getYield(fileName, treeName, histName, lumi, isData):
  global opt, args
    
  file = ROOT.TFile(fileName,"READ")
  tree = file.Get(treeName)
    
  normHist = file.Get(histName)
  normEvents = normHist.GetBinContent(1)

  yield_4l = 0
  yield_4mu = 0
  yield_4e = 0
  yield_2e2mu = 0
  
  for i in range( tree.GetEntries() ):
    tree.GetEntry(i)
    
    if not tree.passedZ4lSelection: continue
    if tree.massZ2 < 4: continue
    if tree.mass4l < opt.xLow or tree.mass4l > opt.xHigh: continue
#    if not tree.passedQCDcut: continue
    if tree.pTL1 < 5 and math.fabs(tree.idL1) == 13: continue
    if tree.pTL2 < 5 and math.fabs(tree.idL2) == 13: continue
    if tree.pTL3 < 5 and math.fabs(tree.idL3) == 13: continue
    if tree.pTL4 < 5 and math.fabs(tree.idL4) == 13: continue

    if tree.pTL1 < 7 and math.fabs(tree.idL1) == 11: continue
    if tree.pTL2 < 7 and math.fabs(tree.idL2) == 11: continue
    if tree.pTL3 < 7 and math.fabs(tree.idL3) == 11: continue
    if tree.pTL4 < 7 and math.fabs(tree.idL4) == 11: continue

    
    if isData:
      weight = 1
    else:
      weight = tree.scaleWeight*tree.eventWeight*tree.dataMC_weight*(lumi*1000)/normEvents
      
    yield_4l += weight
    if math.fabs(tree.idL1) == 13 and math.fabs(tree.idL3) == 13: yield_4mu += weight
    if math.fabs(tree.idL1) == 11 and math.fabs(tree.idL3) == 11: yield_4e += weight
    if (math.fabs(tree.idL1) == 11 and math.fabs(tree.idL3) == 13) or (math.fabs(tree.idL1) == 13 and math.fabs(tree.idL3) == 11) : yield_2e2mu += weight
    
    
  yields = []
  yields.append(yield_4l)
  yields.append(yield_4mu)
  yields.append(yield_4e)
  yields.append(yield_2e2mu)
  
  return yields


if __name__ == "__main__":
  global opt, args
  parseOptions()
  

  path = '/scratch/osghpc/snowball/UF_HZZ4L_Analysis/FullAnalysis2012/MEKD_Z4L/Histogramming/rootFiles'
#  files = ['ZZ4e_tchan.root','ZZ2e2mu_tchan.root','ZZ4mu_tchan.root']
  files = ['ZZTo2e2mu.root','ZZTo4mu.root','ZZTo4e.root','ZZTo2e2tau.root','ZZTo2mu2tau.root','ZZTo4tau.root']
#  files = ['ggH_91.2.root']
  
  myYield_4l = 0
  myYield_4mu = 0
  myYield_4e = 0
  myYield_2e2mu = 0

  for i in range(len(files)):
    print files[i]
    myYields = getYield(path+'/'+files[i], "passedEvents_dataMC","AnaAfterHlt/nEvents", 19.7, False)
    myYield_4l += myYields[0]
    myYield_4mu += myYields[1]
    myYield_4e += myYields[2]
    myYield_2e2mu += myYields[3]

  
#   path = '../Histogramming/rootFiles_Legacy_dataMC'
#   files = ['ZZ_2e2mu.root','ZZ_4e.root','ZZ_4mu.root','ZZto2e2tau.root','ZZto2mu2tau.root','ZZto4tau.root']

#   myYield_4l = 0
#   myYield_4mu = 0
#   myYield_4e = 0
#   myYield_2e2mu = 0

#   for i in range(len(files)):
#     print files[i]
#     myYields = getYield(path+'/'+files[i], "passedEvents_dataMC","AnaAfterHlt/nEvents", 5.051, False)
#     myYield_4l += myYields[0]
#     myYield_4mu += myYields[1]
#     myYield_4e += myYields[2]
#     myYield_2e2mu += myYields[3]

#   path = '../Histogramming_8TeV/rootFiles_Legacy_dataMC'
#   files = ['ZZ_2e2mu.root','ZZ_4e.root','ZZ_4mu.root','ZZto2e2tau.root','ZZto2mu2tau.root','ZZto4tau.root']

#   for i in range(len(files)):
#     print files[i]
#     myYields = getYield(path+'/'+files[i],"passedEvents_dataMC","AnaAfterHlt/nEvents", 19.63, False)
#     myYield_4l += myYields[0]
#     myYield_4mu += myYields[1]
#     myYield_4e += myYields[2]
#     myYield_2e2mu += myYields[3]
    

  myYieldData_4l = 0
  myYieldData_4mu = 0
  myYieldData_4e = 0
  myYieldData_2e2mu = 0
  
#   path = '../Histogramming/rootFiles_Legacy'
#   myYields = getYield(path+'/Data.root',"Ana/passedEvents","Ana/nEvents", 1, True)
#   print myYields
#   myYieldData_4l += myYields[0]
#   myYieldData_4mu += myYields[1]
#   myYieldData_4e += myYields[2]
#   myYieldData_2e2mu += myYields[3]

#   path = '../Histogramming_8TeV/rootFiles_Legacy'
#   myYields = getYield(path+'/Data.root',"AnaAfterHlt/passedEvents","AnaAfterHlt/nEvents", 1, True)
#   myYieldData_4l += myYields[0]
#   myYieldData_4mu += myYields[1]
#   myYieldData_4e += myYields[2]
#   myYieldData_2e2mu += myYields[3]


  print "["+str(opt.xLow)+","+str(opt.xHigh)+"]"
  print "               MC Yields              "
  print "--------------------------------------"
  print "4l: ", myYield_4l
  print "4mu: ", myYield_4mu
  print "4e: ", myYield_4e
  print "2e2mu: ", myYield_2e2mu
  print "              Data Yields             "
  print "--------------------------------------"
  print "4l: ", myYieldData_4l
  print "4mu: ", myYieldData_4mu
  print "4e: ", myYieldData_4e
  print "2e2mu: ", myYieldData_2e2mu

  
