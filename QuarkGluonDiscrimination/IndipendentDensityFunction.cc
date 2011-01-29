#include "QuarkGluonDiscrimination.h"
#include <TH1F.h>
#include <TFile.h>

Double_t Discrimination::IndipendentDensityFunction(TH1F *h1, 
						    TH1F*h2,
						    TH1F*h3,
						    TH1F*h4,
						    const Double_t value1, 
						    const Double_t value2, 
						    const Double_t value3,
 						    const Double_t value4
						) const
{
 h1->Scale(1./h1->Integral("width"));
 h2->Scale(1./h2->Integral("width"));
 h3->Scale(1./h3->Integral("width"));
 h4->Scale(1./h4->Integral("width"));
 return h1->GetBinContent(h1->FindBin(value1))*
        h2->GetBinContent(h2->FindBin(value2))*
        h3->GetBinContent(h3->FindBin(value3))*
        h4->GetBinContent(h4->FindBin(value4));
	
}

Double_t Discrimination::IndipendentDensityFunction(const char* filename,
						const Double_t value1, 
						const Double_t value2, 
						const Double_t value3, 
						const Double_t value4
						) const
{
TFile *f=new TFile(filename);//open file
if(f==NULL)return -1;//check file
TH1F *h1=(TH1F*)f->Get("ncharged")->Clone("ncharged1");
TH1F *h2=(TH1F*)f->Get("nneutral")->Clone("nneutral1");
TH1F *h3=(TH1F*)f->Get("PtD")->Clone("PtD1");
TH1F *h4=(TH1F*)f->Get("rRMS")->Clone("rRMS1");
Double_t R= IndipendentDensityFunction(h1,h2,h3,h4,value1,value2,value3,value4);
f->Close();//close file
delete f;
return R;
}

Double_t Discrimination::IndipendentLikelihood(const char* filename,
                                                const Double_t value1,
                                                const Double_t value2,
                                                const Double_t value3,
                                                const Double_t value4
                                                ) const
{
TFile *f=new TFile(filename);//open file
if(f==NULL)return -1;//check file
TH1F *h1=(TH1F*)f->Get("ncharged_gluon")->Clone("ncharged1");
TH1F *h2=(TH1F*)f->Get("nneutral_gluon")->Clone("nneutral1");
TH1F *h3=(TH1F*)f->Get("PtD_gluon")->Clone("PtD1");
TH1F *h4=(TH1F*)f->Get("rRMS_gluon")->Clone("rRMS1");
Double_t R= IndipendentDensityFunction(h1,h2,h3,h4,value1,value2,value3,value4);
h1=(TH1F*)f->Get("ncharged_quark")->Clone("ncharged2");
h2=(TH1F*)f->Get("nneutral_quark")->Clone("nneutral2");
h3=(TH1F*)f->Get("PtD_quark")->Clone("PtD2");
h4=(TH1F*)f->Get("rRMS_quark")->Clone("rRMS2");
//R= G/G+Q
R=R/(R+IndipendentDensityFunction(h1,h2,h3,h4,value1,value2,value3,value4));
f->Close();//close file
delete f;
return R;
}

Double_t Discrimination::IndipendentLikelihood(const Double_t jtpt,
						const Double_t value1,
                                                const Double_t value2,
                                                const Double_t value3,
                                                const Double_t value4,
						const char*filename
						)const
{
#include "PtBins.h"

char  filename2[255];
for(int i=0;i<NBins-1;i++)
	{
	if((PtBins[i]<jtpt) &&( jtpt<PtBins[i+1] )) 
		{
		sprintf(filename2,filename,PtBins[i],PtBins[i+1]);
		return IndipendentLikelihood(filename2,value1,value2,value3,value4); 
		}
	} 
return -1.0;
}

