.PHONY: all
all: ValidationPlot ValidationHisto

ValidationPlot: src/ReadParameters.C validation/ValidationPlot.C
	g++ -D STANDALONE `root-config --libs --cflags` validation/ValidationPlot.C -o ValidationPlot
ValidationHisto: src/ReadParameters.C validation/ValidationHisto.C
	g++ -D STANDALONE `root-config --libs --cflags` validation/ValidationHisto.C -o ValidationHisto

