
EXTRAFLAGS=

ifneq (${ROOFITSYS}, )
EXTRAFLAGS+=-I${ROOFITSYS}/include -L${ROOFITSYS}/lib
endif


.PHONY: all
all: ValidationPlot ValidationHisto symlink

#different behaviour in mac and linux
.PHONY: symlink
symlink: 
	[ -h validation/src ] || ln -s "${PWD}/src" validation/src
	[ -h susy/src ] || ln -s "${PWD}/src" susy/src
	[ -h quarkflavor/src ] || ln -s "${PWD}/src" quarkflavor/src
	[ -h ValidationBatch.sh ] || ln -s validation/ValidationBatch.sh
	[ -h SusyBatch.sh ] || ln -s susy/SusyBatch.sh

.PHONY: Validation
Validation:ValidationPlot ValidationHisto

ValidationPlot: src/ReadParameters.C validation/ValidationPlot.C symlink
	g++ -D STANDALONE `root-config --libs --cflags` validation/ValidationPlot.C -o ValidationPlot
ValidationHisto: src/ReadParameters.C validation/ValidationHisto.C symlink
	g++ -D STANDALONE `root-config --libs --cflags` validation/ValidationHisto.C -o ValidationHisto

.PHONY: susy 
susy: SusyPlot SusyHisto
SusyPlot: src/ReadParameters.C susy/SusyPlot.C symlink
	g++ -D STANDALONE `root-config --libs --cflags` susy/SusyPlot.C -o SusyPlot
SusyHisto: src/ReadParameters.C susy/SusyHisto.C symlink
	g++ -D STANDALONE `root-config --libs --cflags` susy/SusyHisto.C -o SusyHisto

.PHONY: QG 
QG: QGPlot QGHisto QGAnalysis
QGPlot: src/ReadParameters.C quarkflavor/QGPlot.C symlink
	g++ -D STANDALONE `root-config --libs --cflags` quarkflavor/QGPlot.C -o QGPlot
QGHisto: src/ReadParameters.C quarkflavor/QGHisto.C symlink
	g++ -D STANDALONE `root-config --libs --cflags` quarkflavor/QGHisto.C -o QGHisto
QGAnalysis: src/ReadParameters.C quarkflavor/QGAnalysis.C symlink
	g++ -D STANDALONE `root-config --libs --cflags` $(EXTRAFLAGS) -lRooFit -lRooFitCore -lMinuit quarkflavor/QGAnalysis.C -o QGAnalysis

.PHONY: ValidationQG
ValidationQG:ValidationQGPlot ValidationQGHisto

symlink_vQG:
	[ -h validationQG/src ] || ln -s "${PWD}/src" validationQG/src

ValidationPlotQG: src/ReadParameters.C validationQG/ValidationPlotQG.C symlink_vQG
	g++ -D STANDALONE `root-config --libs --cflags` validationQG/ValidationPlotQG.C -o ValidationPlotQG
ValidationHistoQG: src/ReadParameters.C validationQG/ValidationHistoQG.C symlink_vQG
	g++ -D STANDALONE `root-config --libs --cflags` validationQG/ValidationHistoQG.C -o ValidationHistoQG
