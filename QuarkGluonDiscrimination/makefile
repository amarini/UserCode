CPP=g++
CFLAGS=`root-config --cflags` -fPIC
#LIBS=`root-config --libs`
LIBS= -L/afs/cern.ch/cms/slc5_ia32_gcc434/lcg/root/5.22.00d-cms18/lib -lCore -lCint -lRIO -lNet -lHist -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lz -pthread -lm -ldl -rdynamic

.PHONY:all
all: DensityFunction.o IndipendentDensityFunction.o QuarkGluonDiscrimination.a QuarkGluonDiscrimination.so
.PHONY:clean
clean:
	[ -f DensityFunction.o ] && rm DensityFunction.o || true
	[ -f IndipendentDensityFunction.o ] && rm IndipendentDensityFunction.o || true
	[ -f QuarkGluonDiscrimination.a ] && rm QuarkGluonDiscrimination.a || true
	[ -f libQuarkGluonDiscrimination.so.1.1 ] && rm libQuarkGluonDiscrimination.so.1.1 || true
	[ -h libQuarkGluonDiscrimination.so.1 ] && rm libQuarkGluonDiscrimination.so.1 || true
	[ -h libQuarkGluonDiscrimination.so ] && rm libQuarkGluonDiscrimination.so || true
DensityFunction.o: DensityFunction.cc QuarkGluonDiscrimination.h
	${CPP} ${LIBS} ${CFLAGS} -c  DensityFunction.cc

IndipendentDensityFunction.o: QuarkGluonDiscrimination.h IndipendentDensityFunction.cc
	${CPP} ${LIBS} ${CFLAGS} -c IndipendentDensityFunction.cc

QuarkGluonDiscrimination.a: IndipendentDensityFunction.o DensityFunction.o
	ar rcs QuarkGluonDiscrimination.a IndipendentDensityFunction.o DensityFunction.o
QuarkGluonDiscrimination.so: QuarkGluonDiscrimination.a
	${CPP} ${LIBS} ${CFLAGS} -shared -Wl,-soname,libQuarkGluonDiscrimination.so.1 QuarkGluonDiscrimination.a -o libQuarkGluonDiscrimination.so.1.1
	[ -h libQuarkGluonDiscrimination.so.1 ] && rm libQuarkGluonDiscrimination.so.1 || true
	[ -h libQuarkGluonDiscrimination.so ] && rm libQuarkGluonDiscrimination.so || true
	ln -s libQuarkGluonDiscrimination.so.1.1 libQuarkGluonDiscrimination.so.1
	ln -s libQuarkGluonDiscrimination.so.1 libQuarkGluonDiscrimination.so

