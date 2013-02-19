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
