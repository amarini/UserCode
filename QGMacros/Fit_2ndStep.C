/* Author: andrea.carlo.marini@cern.ch
 * Date:   07/03/2012
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "TH1F.h"
#include "TH2F.h"
#include "PtBins.h"


char GetChar(FILE *fr);

int Fit_2ndStep(const char fileName[]="nCharged.txt",int parameter=0)
{

FILE *fr=fopen(fileName,"r");
char c; //testing character

//[quark]
//line
//ptmin ptmax rhomin rhomax 2+Npar xmin xmax par0 par1 ..
//...
//[gluon]
//line


//getting binnig
double RhoBins[25];int nRhoBins=20;
double PtBins[25];int nPtBins=18;
getBins_int(18,PtBins,20,1000,true);
PtBins[18]=3500;
getBins_int(21,RhoBins,0,20,false);
//---
TH2F *quark=new TH2F("quark","quark",nPtBins,PtBins,nRhoBins,RhoBins);
TH2F *gluon=new TH2F("gluon","gluon",nPtBins,PtBins,nRhoBins,RhoBins);

//fscanf(fr,"[quark]\n");//move into the file
while (( GetChar(fr) == '[') || (GetChar(fr)== '{'))
	{
	c='\0';while(c!='\n'){fscanf(fr,"%c",&c);fprintf(stderr,"%c",c);}
	} //[quark] {bla}

while (( GetChar(fr) != '[') && (GetChar(fr)!= EOF) )
{
	//read a formatted line
	float ptmin, ptmax, rhomin, rhomax,par[10];int nPar;
	fscanf(fr,"%f %f %f %f %d %*f %*f",&ptmin,&ptmax,&rhomin,&rhomax,&nPar);nPar-=2;
	fprintf(stderr,"%f %f %f %f %d %d\n",ptmin,ptmax,rhomin,rhomax,nPar,parameter);	

	for(int i=0;i<nPar;++i)fscanf(fr,"%f",&par[i]);
	c='\0';while(c!='\n'){fscanf(fr,"%c",&c);} //go to the end of line
	//check
	if(nPar<parameter){perror("I do not have that parameter\n");break;}
	//filling the histogram
	quark->SetBinContent(quark->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parameter] );
}//end of while: loop on the lines

//skip lines that begin with [ or { -> useless
while (( GetChar(fr) == '[') || (GetChar(fr)== '{'))
	{
	c='\0';while(c!='\n'){fscanf(fr,"%c",&c);fprintf(stderr,"%c",c);}
	} //[gluon] {bla}

while (( GetChar(fr) != '[') && (GetChar(fr)!= EOF))
{
	//read a formatted line
	float ptmin, ptmax, rhomin, rhomax,par[10];int nPar;
	fscanf(fr,"%f %f %f %f %d %*f %*f",&ptmin,&ptmax,&rhomin,&rhomax,&nPar);nPar-=2;
	fprintf(stderr,"%f %f %f %f %d %d\n",ptmin,ptmax,rhomin,rhomax,nPar,parameter);	

	for(int i=0;i<nPar;++i)fscanf(fr,"%f",&par[i]);
	c='\0';while(c!='\n'){fscanf(fr,"%c",&c);} //go to the end of line
	//check
	if(nPar<parameter){perror("I do not have that parameter\n");break;}
	//filling the histogram
	gluon->SetBinContent(gluon->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parameter] );
}//end of while: loop on the lines
fclose(fr);

//END OF READ FILE
float ChosenPt=90.;
TH1D *q=(TH1D*)quark->ProjectionY("_py", quark->GetXaxis()->FindBin( ChosenPt),quark->GetXaxis()->FindBin(ChosenPt ) );
TH1D *g=(TH1D*)gluon->ProjectionY("_py", gluon->GetXaxis()->FindBin( ChosenPt),gluon->GetXaxis()->FindBin(ChosenPt ) );
TCanvas *c1=new TCanvas("c1","c1",800,800);
g->SetMarkerColor(kRed);
g->SetLineColor(kRed);
q->SetMarkerStyle(20);
g->SetMarkerStyle(20);
q->Draw("AXIS");
q->Draw("AXIS X+ Y+ SAME");
q->Draw("P SAME");
g->Draw("P SAME");

}//enf of Fit_2ndStep function



char GetChar(FILE *fr)
{
fpos_t pos;//position in the file
char c;
//get position of the line to be analized;
fgetpos(fr,&pos);
c=fgetc(fr);
fsetpos(fr,&pos); //moving back to the beginning of the line
return c;
}
