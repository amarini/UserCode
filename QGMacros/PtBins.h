
int getBins(double  *Bins,int nBins,double MinBin=15.0,double MaxBin=1000.,bool log=false);
int getBin(int nBins,double  *Bins,double value,double*x0=0,double*x1=0);


int getBins(double  *Bins,int nBins,double MinBin,double MaxBin,bool log)
{
double incr;
if(log)
	{
	incr=TMath::Power(MaxBin/MinBin,1.0/double(nBins));
	Bins[0]=MinBin;
	Bins[nBins]=MaxBin;
	for(int i=1;i<nBins;i++)
		Bins[i]=Bins[i-1]*incr;
	}
else
	{
	incr=(MaxBin-MinBin)/nBins;
	Bins[0]=MinBin;
	Bins[nBins+1]=MaxBin;
	for(int i=1; i<nBins+1;i++)
		Bins[i]=Bins[i-1]+incr;
	}
return 0;
}
int getBin(int nBins,double  Bins[],double value,double *x0,double *x1)
{
int R=0;
//int nBins=sizeof(Bins)/sizeof(double);//?
if(value <Bins[0])return -1;
if(value >Bins[nBins])return -1;
for(R=0;R<nBins;R++)
	{
	if(Bins[R]>value)break;	
	}
R--;
if(x0) *x0=Bins[R];
if(x1) *x1=Bins[R+1];
return R;	
}
