ROOUNFOLDLIB=/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/
ROOUNFOLDINC=/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/src/

CPPFLAGS=`root-config --libs --cflags`

.PHONY:all
all: QGDumper.a qgexe

.PHONY: QGDumper
QGDumper: QGDumper.o

QGDumper.o: ../BaseAnalysis/BaseAnalysis.h QGDumper.h QGDumper.cc
	g++ $(CPPFLAGS) -c QGDumper.cc

.PHONY: BaseAnalysis
BaseAnalysis:
	cd .. && make BaseAnalysis

QGDumper.a: QGDumper BaseAnalysis
	ar  rcs QGDumper.a QGDumper.o ../BaseAnalysis/BaseAnalysis.o ../BaseAnalysis/ReadParameters.o

qgexe: qgexe.cc QGDumper QGDumper.a
	g++ $(CPPFLAGS) qgexe.cc QGDumper.a -o qgexe
#	echo "export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:/afs/cern.ch/work/a/amarini/RooUnfold-1.1.1/"

.PHONY:clean
clean:
	[ -f QGDumper.o ] && rm QGDumper.o
	[ -f QGDumper.a ] && rm QGDumper.a
	[ -f qgexe ] && rm qgexe
