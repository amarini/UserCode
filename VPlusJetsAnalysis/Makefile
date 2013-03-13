
EXTRAFLAGS=

ifneq (${ROOFITSYS}, )
EXTRAFLAGS+=-I${ROOFITSYS}/include -L${ROOFITSYS}/lib
endif


.PHONY: all
all: ValidationPlot ValidationHisto simlink

#different behaviour in mac and linux
.PHONY: simlink
simlink: 
	[ -h validation/src ] || ln -s ${PWD}/src validation/src
	[ -h susy/src ] || ln -s ${PWD}/src susy/src
	[ -h quarkflavor/src ] || ln -s ${PWD}/src quarkflavor/src
	[ -h ValidationBatch.sh ] || ln -s validation/ValidationBatch.sh
	[ -h SusyBatch.sh ] || ln -s susy/SusyBatch.sh

.PHONY: Validation
Validation:ValidationPlot ValidationHisto

ValidationPlot: src/ReadParameters.C validation/ValidationPlot.C simlink
	g++ -D STANDALONE `root-config --libs --cflags` validation/ValidationPlot.C -o ValidationPlot
ValidationHisto: src/ReadParameters.C validation/ValidationHisto.C simlink
	g++ -D STANDALONE `root-config --libs --cflags` validation/ValidationHisto.C -o ValidationHisto

.PHONY: susy 
susy: SusyPlot SusyHisto
SusyPlot: src/ReadParameters.C susy/SusyPlot.C simlink
	g++ -D STANDALONE `root-config --libs --cflags` susy/SusyPlot.C -o SusyPlot
SusyHisto: src/ReadParameters.C susy/SusyHisto.C simlink
	g++ -D STANDALONE `root-config --libs --cflags` susy/SusyHisto.C -o SusyHisto

.PHONY: QG 
QG: QGPlot QGHisto QGAnalysis
QGPlot: src/ReadParameters.C quarkflavor/QGPlot.C simlink
	g++ -D STANDALONE `root-config --libs --cflags` quarkflavor/QGPlot.C -o QGPlot
QGHisto: src/ReadParameters.C quarkflavor/QGHisto.C simlink
	g++ -D STANDALONE `root-config --libs --cflags` quarkflavor/QGHisto.C -o QGHisto
QGAnalysis: src/ReadParameters.C quarkflavor/QGAnalysis.C simlink
	g++ -D STANDALONE `root-config --libs --cflags` $(EXTRAFLAGS) -lRooFit -lRooFitCore -lMinuit quarkflavor/QGAnalysis.C -o QGAnalysis
