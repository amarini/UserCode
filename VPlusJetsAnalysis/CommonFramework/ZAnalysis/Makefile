ROOUNFOLDLIB=/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/
ROOUNFOLDINC=/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/src/

CPPFLAGS=`root-config --libs --cflags` -I$(ROOUNFOLDINC) -L$(ROOUNFOLDLIB) -lRooUnfold

.PHONY:all
all: ZAnalysis.a

.PHONY: ZAnalysis
ZAnalysis: ZAnalysis.o

ZAnalysis.o: ../BaseAnalysis/BaseAnalysis.h ZAnalysis.h ZAnalysis.cc
	g++ $(CPPFLAGS) -c ZAnalysis.cc

.PHONY: BaseAnalysis
BaseAnalysis:
	cd .. && make BaseAnalysis

ZAnalysis.a: ZAnalysis BaseAnalysis
	ar  rcs ZAnalysis.a ZAnalysis.o ../BaseAnalysis/BaseAnalysis.o ../BaseAnalysis/ReadParameters.o

zexe: zexe.cc ZAnalysis ZAnalysis.a
	g++ $(CPPFLAGS) ZAnalysis.a zexe.cc -o zexe

.PHONY:clean
clean:
	[ -f ZAnalysis.o ] && rm ZAnalysis.o
	[ -f ZAnalysis.a ] && rm ZAnalysis.a
	[ -f zexe ] && rm zexe
