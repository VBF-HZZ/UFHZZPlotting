#!/usr/bin/python
import sys, os, pwd, commands
import optparse, shlex
import math
from array import array

def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-s', '--sqrts',   dest='sqrts', type='string', default="8", help='Sqrts = 7, 8, or 7p8')
    parser.add_option('-d', '--dir',   dest='dir', type='string', default="plots_test", help='Output Directory')
    parser.add_option('-l', '--lists',  action='store_true', dest='lists',  default=False, help='Make Lists only')


    
    global opt, args
    (opt, args) = parser.parse_args()

    if (opt.sqrts != '7' and opt.sqrts != '8' and opt.sqrts != '7p8'):
        print "Sqrts ", opt.sqrts, " not an option, please choose 7, 8, or 7p8"
        sys.exit()

    makeDirectory(opt.dir)
    makeDirectory(opt.dir+'/eps')
    makeDirectory(opt.dir+'/png')
    makeDirectory(opt.dir+'/txt')    
        
    if not os.path.exists(opt.dir):
        print opt.dir, " not found!"
        sys.exit()


def processCmd(cmd):
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        sys.exit()
    return status


def makeDirectory(subDirName):
    if (not os.path.exists(subDirName)):
        cmd = 'mkdir -p '+subDirName
        status, output = commands.getstatusoutput(cmd)
        if status !=0:
            print 'Error in creating submission dir '+subDirName+'. Exiting...'
            sys.exit()
    else:
        print 'Directory '+subDirName+' already exists. Exiting...'
        sys.exit()


def makeHists():
    global opt, args
    parseOptions()
    
    myFileDir = opt.dir

    ####### Main Settings ######
    myPlotXlow =      [60,  60]
    myPlotXhigh =     [600, 800]
    myBinSize =       [10,  10]
    myPlotXlowZoom =  [70,  70.5]
    myPlotXhighZoom = [170, 181.5]
    myBinSizeZoom =   [3,   3]
    ############################

    mainName = ''
    compileName = ''
    
    if opt.sqrts == '7':
        if opt.lists:
            mainName = 'makeLists7TeV.exe'
            compileName = 'Lists_7TeV'
        else:
            mainName = 'makeHists7TeV.exe'
            compileName = 'Hists_7TeV'
    elif opt.sqrts == '8':
        if opt.lists:
            mainName = 'makeLists8TeV.exe'
            compileName = 'Lists_8TeV'
        else:
            mainName = 'makeHists8TeV.exe'
            compileName = 'Hists_8TeV'
    elif opt.sqrts == '7p8':
        mainName = 'makeHists7p8TeV.exe'
        compileName = 'Hists_7p8TeV'
    
    cmd = "make "+compileName
    status = processCmd(cmd)
    if status != 0:
        print compileName," could not be compiled, please check!"
        sys.exit()


    if opt.lists:
        cmd = "./"+mainName+" "+myFileDir
        status = processCmd(cmd)
        if status != 0:
            print cmd, " failed, fix your bug!"
            sys.exit()
    else:
        printData = 0
        for i in range(len(myPlotXlow)):
            if i == 0:
                printData = 1
            else: 
                printData = 0
            cmd = "./"+mainName+" "+myFileDir+" "+str(myPlotXlow[i])+" "+str(myPlotXhigh[i])+" "+str(myBinSize[i])+" "+str(myPlotXlowZoom[i])+" "+str(myPlotXhighZoom[i])+" "+str(myBinSizeZoom[i])+" "+str(printData)+" > log.txt"
            status = processCmd(cmd)
            if status != 0:
                print cmd, " failed, fix your bug!"
                sys.exit()
        

    sys.exit()


if __name__ == "__main__":
    makeHists()

