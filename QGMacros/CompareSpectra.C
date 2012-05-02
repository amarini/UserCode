// Macro to Compare Pt Spectra between data ph+jet + dijet and MC

int CompareSpectra(int PtMin=80,int PtMax=100,int RhoMin=4,int RhoMax=6,
			const char *outFileName="",
			const char *fileName_Ph_MC="Omog_QGStudies_*_Summer11.root",
			const char *fileName_Di_MC="Omog_DiJet_QCD_HT_Summer11.root",
			const char *fileName_Ph_data="Omog_QGStudies_Photon_Run2011_FULL.root",
			const char *fileName_Di_data="Omog_DiJet_HT_Run2011_FULL.root"
		)
{
TChain *Ph_MC=new TChain("omog");
TChain *Di_MC=new TChain("omog");
TChain *Ph_data=new TChain("omog");
TChain *Di_data=new TChain("omog");

string Directory=string("/shome/amarini/CMSSW_4_2_8_patch7/src/UserCode/amarini/Omog/");
Ph_MC->Add( (Directory+fileName_Ph_MC).c_str());
Di_MC->Add( (Directory+fileName_Di_MC).c_str());
Ph_data->Add( (Directory+fileName_Ph_data).c_str());
Di_data->Add( (Directory+fileName_Di_data).c_str());

TString Selection=TString::Format("%d <ptJet0 && ptJet0<%d && %d <rhoPF && rhoPF<%d",PtMin,PtMax,RhoMin,RhoMax);

TH1F* h_Ph_MC  =new TH1F("Ph_MC"  ,"Ph_MC"  ,40,PtMin,PtMax);
TH1F* h_Di_MC  =new TH1F("Di_MC"  ,"Di_MC"  ,40,PtMin,PtMax);
TH1F* h_Ph_data=new TH1F("Ph_data","Ph_data",40,PtMin,PtMax);
TH1F* h_Di_data=new TH1F("Di_data","Di_data",40,PtMin,PtMax);

//Drawing
Ph_MC->Draw("ptJet0>>Ph_MC",("eventWeight*("+Selection+")").Data(),"NORM GOFF");
Di_MC->Draw("ptJet0>>Di_MC",("eventWeight*("+Selection+")").Data(),"NORM GOFF");
Ph_data->Draw("ptJet0>>Ph_data",Selection.Data(),"NORM GOFF");
Di_data->Draw("ptJet0>>Di_data",Selection.Data(),"NORM GOFF");

//Setting Colors and stuff
h_Ph_MC->SetLineColor(kBlue+2);
h_Ph_MC->SetFillColor(kBlue-4);
h_Ph_MC->SetFillStyle(3004);
h_Ph_data->SetMarkerColor(kBlue+2);
h_Ph_data->SetMarkerStyle(20);
h_Di_MC->SetLineColor(kRed+2);
h_Di_MC->SetFillColor(kRed-4);
h_Di_MC->SetFillStyle(3004);
h_Di_data->SetMarkerColor(kRed+2);
h_Di_data->SetMarkerStyle(20);
//Draw
TCanvas *c1=new TCanvas();
h_Ph_MC->Draw("HIST");
h_Di_MC->Draw("HIST SAME");
h_Ph_data->Draw("P SAME");
h_Di_data->Draw("P SAME");

if(outFileName[0]!='\0') c1->SaveAs(outFileName);

return 0;
}
