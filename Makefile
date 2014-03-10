CXX=g++
CC=gcc
CFLAGS=-O2 -Wall
LDFLAGS=-Llib
PRFFLAGS=-lProof
THRDFLAGS=-lThread
INS=-I$(ROOTSYS)/include/root
INS2=-I$(ROOFITSYS)/include
INSS=-I./include

LD1=-L$(ROOFITSYS)/lib

CFLAGS += `root-config --cflags`
LIBS += `root-config --glibs`

LDa=-lRooFitCore
LDb=-lRooFit

.PHONY: clean all main test

all: Hists_7TeV Hists_8TeV Hists_7p8TeV 

HistMakerFromTree: HistMakerFromTree.o
	$(CXX) -o makeHists.exe HistMakerFromTree.o $(LIBS)

Hists_7TeV: Hists_7TeV.o
	$(CXX) -o makeHists7TeV.exe Hists_7TeV.o $(LIBS)

Hists_8TeV: Hists_8TeV.o
	$(CXX) -o makeHists8TeV.exe Hists_8TeV.o $(LIBS)

Lists_8TeV: Lists_8TeV.o
	$(CXX) -o makeLists8TeV.exe Lists_8TeV.o $(LIBS)

Lists_7TeV: Lists_7TeV.o
	$(CXX) -o makeLists7TeV.exe Lists_7TeV.o $(LIBS)

Hists_7p8TeV: Hists_7p8TeV.o
	$(CXX) -o makeHists7p8TeV.exe Hists_7p8TeV.o $(LIBS)

Hists_m4lerr_7p8TeV: Hists_m4lerr_7p8TeV.o
	$(CXX) -o makeHistErrors_7p8TeV.exe Hists_m4lerr_7p8TeV.o $(LIBS)

Desert: Desert.o
	$(CXX) -o Desert.exe Desert.o $(LIBS)

ZpeakComposition_7TeV: ZpeakComposition_7TeV.o
	$(CXX) -o ZpeakComposition_7TeV.exe ZpeakComposition_7TeV.o $(LIBS)

ZpeakComposition_8TeV: ZpeakComposition_8TeV.o
	$(CXX) -o ZpeakComposition_8TeV.exe ZpeakComposition_8TeV.o $(LIBS)

test: test.o
	$(CXX) -o test.exe test.o $(LIBS)

clean:
	@rm *.o *.exe *~ 


##############RULES##############
.cc.o:
	$(CXX) $(CFLAGS) $(INS) $(INSS) -c $<
.cpp.o:
	$(CXX) $(CFLAGS) $(INS) $(INSS) -c $<


