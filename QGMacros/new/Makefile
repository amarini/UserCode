CFLAGS=--ansi -fPIC
all: PtBins.o Plot.o SplitTree FitN.o PtMeans.o QQPlot.o CMSLatex.o MinChiSquare.o

clean: 
	rm PtBins.o Plot.o SplitTree FitN.o PtMeans.o

PtBins.o: PtBins.C PtBins.h Makefile
	g++ ${CFLAGS} `root-config --libs --cflags` -c PtBins.C

Plot.o: Plot.C Plot.h Makefile
	g++ ${CFLAGS} `root-config --libs --cflags` -c Plot.C

SplitTree: SplitTree.C Makefile
	g++ ${CFLAGS} `root-config --libs --cflags` -o SplitTree SplitTree.C

FitN.o: FitN.C FitN.h Makefile
	g++ ${CFLAGS} `root-config --libs --cflags` -c FitN.C

PtMeans.o: PtMeans.C PtMeans.h PtBins.o Makefile
	g++ ${CFLAGS} `root-config --libs --cflags` -c PtMeans.C

QQPlot.o: QQPlot.C QQPlot.h Makefile
	g++ ${CFLAGS} `root-config --libs --cflags` -c QQPlot.C

CMSLatex.o: CMSLatex.C CMSLatex.h Makefile
	g++ ${CFLAGS} `root-config --libs --cflags` -c CMSLatex.C

MinChiSquare.o: MinChiSquare.C MinChiSquare.h Makefile
	g++ ${CFLAGS} `root-config --libs --cflags` -c MinChiSquare.C
