#include "QuarkGluonDiscrimination.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

Double_t Discrimination::DensityFunction(TH1F *h,Double_t value) const
{
return 0.0;
}

int Discrimination::Populate(const char *filename,const char type,const char*treename)
{
TFile *f=new TFile(filename);//open file
TTree *t=(TTree*)f->Get(treename); //tree in file
//variable Declaration
Int_t ncharged,nneutral,pdgID;
Double_t PtD,rRMS;
//Setting Branch Address
t->SetBranchAddress("PtD",&PtD);
t->SetBranchAddress("r",&rRMS);
t->SetBranchAddress("ncharged",&ncharged);
t->SetBranchAddress("nneutral",&nneutral);
t->SetBranchAddress("pdgID",&pdgID);
for(long int i=1; i<t->GetEntries();i++)
	{
	t->GetEntry(i);
	switch(type)
		{
		case 'Q': if(pdgID==21)continue;
		     break;
		case 'G': if(pdgID!=21)continue;
		     break;
		default: ;//nothing
		}
	SpaceStructure *nodo = new SpaceStructure;
	nodo->data[0]=ncharged;
	nodo->data[1]=nneutral;
	nodo->data[2]=PtD;
	nodo->data[3]=rRMS;
	Space->push_back(nodo);
	}

f->Close();//close file
delete f;//delete pointer
return 0;
}
