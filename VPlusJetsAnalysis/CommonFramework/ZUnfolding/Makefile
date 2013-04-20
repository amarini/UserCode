ROOUNFOLDLIB=/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/
ROOUNFOLDINC=/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/src/

CPPFLAGS=`root-config --libs --cflags` -I$(ROOUNFOLDINC) -L$(ROOUNFOLDLIB) -lRooUnfold

ZUnfolding: ZUnfolding.cc ZAnalysis
	g++ $(CPPFLAGS) ZUnfolding.cc ../ZAnalysis/ZAnalysis.a  -o ZUnfolding

.PHONY:ZAnalysis
ZAnalysis:
	cd .. && make ZAnalysis
