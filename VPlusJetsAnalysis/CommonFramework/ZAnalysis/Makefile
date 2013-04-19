ROOUNFOLDLIB=/Users/andreamarini/Downloads/RooUnfold-1.1.1/
ROOUNFOLDINC=/Users/andreamarini/Downloads/RooUnfold-1.1.1/src/

CPPFLAGS=`root-config --libs --cflags` -I$(ROOUNFOLDINC) -L$(ROOUNFOLDLIB) -lRooUnfold

.PHONY: ZAnalysis
ZAnalysis: ../BaseAnalysis/BaseAnalysis.h ZAnalysis.h ZAnalysis.cc
	g++ $(CPPFLAGS) -c ZAnalysis.cc
ZAnalysis.o: ZAnalysis

ZAnalysis.a: ZAnalysis
	ar  rcs ZAnalysis.a ZAnalysis.o ../BaseAnalysis/BaseAnalysis.o ../BaseAnalysis/ReadParameters.o

zexe: zexe.cc ZAnalysis ZAnalysis.a
	g++ $(CPPFLAGS) ZAnalysis.a zexe.cc -o zexe

.PHONY:clean
clean:
	[ -f ZAnalysis.o ] && rm ZAnalysis.o
	[ -f ZAnalysis.a ] && rm ZAnalysis.a
	[ -f zexe ] && rm zexe
