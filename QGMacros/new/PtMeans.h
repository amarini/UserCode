#include <stdio.h>
#include "PtBins.h"
#ifndef MEANS_H
#define MEANS_H
class Means : public Bins
{
public:
Means();
~Means();

double *PtMeans;
double *RhoMeans;
int nPtMeans;
int nRhoMeans;
int ReadTxt(const char *fileName);
int WriteTxt(const char *fileName,const char*outFileName,const char *treeName="omog");
};
#endif
