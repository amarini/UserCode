
.PHONY: all
all: ValidationPlot ValidationHisto simlink

#different behaviour in mac and linux
.PHONY: simlink
simlink: 
	[ -h validation/src ] || ln -s ${PWD}/src validation/src

ValidationPlot: src/ReadParameters.C validation/ValidationPlot.C
	g++ -D STANDALONE `root-config --libs --cflags` validation/ValidationPlot.C -o ValidationPlot
ValidationHisto: src/ReadParameters.C validation/ValidationHisto.C
	g++ -D STANDALONE `root-config --libs --cflags` validation/ValidationHisto.C -o ValidationHisto

.PHONY: susy 
susy: SusyPlot SusyHisto
SusyPlot: src/ReadParameters.C susy/SusyPlot.C
	g++ -D STANDALONE `root-config --libs --cflags` susy/SusyPlot.C -o SusyPlot
SusyHisto: src/ReadParameters.C susy/SusyHisto.C
	g++ -D STANDALONE `root-config --libs --cflags` susy/SusyHisto.C -o SusyHisto

.PHONY: QG 
QG: QGPlot QGHisto QGAnalysis
QGPlot: src/ReadParameters.C quarkflavor/QGPlot.C
	g++ -D STANDALONE `root-config --libs --cflags` quarkflavor/QGPlot.C -o QGPlot
QGHisto: src/ReadParameters.C quarkflavor/QGHisto.C
	g++ -D STANDALONE `root-config --libs --cflags` quarkflavor/QGHisto.C -o QGHisto
QGAnalysis: src/ReadParameters.C quarkflavor/QGAnalysis.C
	g++ -D STANDALONE `root-config --libs --cflags` -lRooFit -lRooFitCore -lMinuit quarkflavor/QGAnalysis.C -o QGAnalysis
