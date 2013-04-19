ROOUNFOLDLIB=/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/libRooUnfold.so
ROOUNFOLDINC=/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/src/

CPPFLAGS=`root-config --libs --cflags` -I$(ROOUNFOLDINC) -L$(ROOUNFOLDLIB) -lRooUnfold

.PHONY: all
all: BaseAnalysis.a

.PHONY: BaseAnalysis
BaseAnalysis: BaseAnalysis.o ReadParameters.o

BaseAnalysis.o: BaseAnalysis.cc BaseAnalysis.h
	g++ $(CPPFLAGS) -c BaseAnalysis.cc

ReadParameters.o: BaseAnalysis.h ReadParameters.cc
	g++ $(CPPFLAGS) -c ReadParameters.cc

BaseAnalysis.a:BaseAnalysis.o ReadParameters.o
	ar rcs BaseAnalysis.a BaseAnalysis.o ReadParameters.o

.PHONY: clean
clean:
	[ -f BaseAnalysis.o ] && rm BaseAnalysis.o
	[ -f BaseAnalysis.a ] && rm BaseAnalysis.a
	[ -f ReadParameters.o ] && rm ReadParameters.o
