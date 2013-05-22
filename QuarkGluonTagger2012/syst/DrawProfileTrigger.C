{

//.L MyMath.C+
gSystem->Load("MyMath_C.so");
TChain *mc=new TChain("Hbb/events");
TChain *data=new TChain("Hbb/events");
mc->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_flatQCD_P6_Dijets.root");
data  ->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_data2012_JetMon_Dijets.root");

TChain *mch=new TChain("Hbb/events");
mch->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_flatQCD_Hpp_Dijets.root");

mc->SetLineColor(kRed); mch->SetLineColor(kGreen);
data->SetMarkerStyle(20);

//fast --DEBUG
//mc->SetMaxEntryLoop(10000);
//data->SetMaxEntryLoop(10000);
//mch->SetMaxEntryLoop(10000);

gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

string selection="((jetPt[2]*2.0)/(jetPt[0]+jetPt[1])<.3) && abs(MyMath::DeltaPhi(jetPhi[0],jetPhi[1]) )>2.5 && (jetPt[0]>20 && jetPt[0]<100) && abs(jetEta[0])<2.0 ";

//-----------------------------ptD-------------------------------------//

TCanvas *c=new TCanvas("c1","c1",800,800);
TPad *p1=new TPad("p1","p1",.4,.12,.64,.4);
c->cd();

TProfile *Pdata=new TProfile("Pdata","Pdata",100,20,100);
TProfile *PP6  =new TProfile("PP6","PP6",100,20,100);
TProfile *PHpp =new TProfile("PHpp","PHpp",100,20,100);
TProfile *PP6_noTrigger  =new TProfile("PP6_noTrigger","PP6_noTrigger",100,20,100);
TProfile *PHpp_noTrigger =new TProfile("PHpp_noTrigger","PHpp_noTrigger",100,20,100);

printf("plot data\n");
data->Draw("jetPtD_QC[0]:jetPt[0]>>Pdata", (selection+"&&triggerResult[1]").c_str(),"prof");
printf("plot pythia\n");
mc->Draw("jetPtD_QC[0]:jetPt[0]>>PP6",(selection+"&&triggerResult[1]").c_str(),"prof SAME");
printf("plot herwig\n");
mch->Draw("jetPtD_QC[0]:jetPt[0]>>PHpp",(selection+"&&triggerResult[1]").c_str(),"prof SAME");

mc->Draw("jetPtD_QC[0]:jetPt[0]>>PP6_noTrigger",(selection).c_str(),"prof SAME");
mch->Draw("jetPtD_QC[0]:jetPt[0]>>PHpp_noTrigger",(selection).c_str(),"prof SAME");

PP6->SetLineColor(kRed);   
PHpp->SetLineColor(kGreen);

PP6_noTrigger->SetLineColor(kRed+2);   PP6_noTrigger ->SetLineWidth(2); 
PHpp_noTrigger->SetLineColor(kGreen+3);PHpp_noTrigger->SetLineWidth(2); 

TLegend *L=new TLegend(0.65,.12,.89,.40,"P_{T}D QC");
L->AddEntry("Pdata","data");
L->AddEntry("PP6","P6");
L->AddEntry("PHpp","Hpp");
L->AddEntry("PP6_noTrigger","P6 (noTrigger)");
L->AddEntry("PHpp_noTrigger","Hpp (noTrigger)");
L->Draw();

Pdata->GetXaxis()->SetTitle("P_{T}^{Jet}");
Pdata->GetYaxis()->SetTitle("P_{T}D_{QC} mean");

p1->Draw(); p1->SetBottomMargin(0.05);p1->SetLeftMargin(0.05);p1->SetTopMargin(0);p1->SetRightMargin(0.01);
p1->cd();
TH1F *h1data=new TH1F("h1data","h1data",50,0,1.0001);
TH1F *h1P6=new TH1F("h1P6","h1P6",50,0,1.0001);
TH1F *h1Hpp=new TH1F("h1Hpp","h1Hpp",50,0,1.0001);
TH1F *h1P6_noTrigger=new TH1F("h1P6_noTrigger","h1P6_noTrigger",50,0,1.0001);
TH1F *h1Hpp_noTrigger=new TH1F("h1Hpp_noTrigger","h1Hpp_noTrigger",50,0,1.0001);
data->Draw("jetPtD_QC[0]>>h1data",(selection+string("&&triggerResult[1] && 40 < jetPt[0] && jetPt[0] < 45 ")).c_str(),"P NORM");
mc->Draw("jetPtD_QC[0]>>h1P6",(selection+string("&&triggerResult[1] && 40<jetPt[0] && jetPt[0]<45 ")).c_str(),"SAME NORM");
mch->Draw("jetPtD_QC[0]>>h1Hpp",(selection+string("&&triggerResult[1] && 40<jetPt[0] && jetPt[0]<45 ")).c_str(),"SAME NORM");
mc->Draw("jetPtD_QC[0]>>h1P6_noTrigger",(selection+string("&& 40<jetPt[0] && jetPt[0]<45 ")).c_str(),"SAME NORM");
mch->Draw("jetPtD_QC[0]>>h1Hpp_noTrigger",(selection+string(" && 40<jetPt[0] && jetPt[0]<45 ")).c_str(),"SAME NORM");

h1data->SetMarkerStyle(20);
h1P6->SetLineColor(kRed);   
h1Hpp->SetLineColor(kGreen);   
h1P6_noTrigger->SetLineColor(kRed+2);    h1P6_noTrigger ->SetLineWidth(2);
h1Hpp_noTrigger->SetLineColor(kGreen+3); h1Hpp_noTrigger->SetLineWidth(2); 
TLatex *l=new TLatex; l->SetNDC(); l->DrawLatex(0.5,.89,"40<P_{T}[GeV]<45");

//-----------------------------nChg-------------------------------------//
TCanvas *c2=new TCanvas("c2","c2",800,800);
TPad *p2=new TPad("p2","p2",.4,.12,.64,.4);
c2->cd();

TProfile *P2data=new TProfile("P2data","Pdata",100,20,100);
TProfile *P2P6  =new TProfile("P2P6","PP6",100,20,100);
TProfile *P2Hpp =new TProfile("P2Hpp","PHpp",100,20,100);
TProfile *P2P6_noTrigger  =new TProfile("P2P6_noTrigger","PP6_noTrigger",100,20,100);
TProfile *P2Hpp_noTrigger =new TProfile("P2Hpp_noTrigger","PHpp_noTrigger",100,20,100);

printf("plot data\n");
data->Draw("jetChgPart_QC[0]:jetPt[0]>>P2data",(selection+"&&triggerResult[1]").c_str(),"prof");
printf("plot pythia\n");
mc->Draw("jetChgPart_QC[0]:jetPt[0]>>P2P6",(selection+"&&triggerResult[1]").c_str(),"prof SAME");
printf("plot herwig\n");
mch->Draw("jetChgPart_QC[0]:jetPt[0]>>P2Hpp",(selection+"&&triggerResult[1]").c_str(),"prof SAME");

mc->Draw("jetChgPart_QC[0]:jetPt[0]>>P2P6_noTrigger",(selection).c_str(),"prof SAME");
mch->Draw("jetChgPart_QC[0]:jetPt[0]>>P2Hpp_noTrigger",(selection).c_str(),"prof SAME");

P2P6->SetLineColor(kRed);   
P2Hpp->SetLineColor(kGreen);   
P2P6_noTrigger->SetLineColor(kRed+2);    P2P6_noTrigger ->SetLineWidth(2);
P2Hpp_noTrigger->SetLineColor(kGreen+3); P2Hpp_noTrigger->SetLineWidth(2); 

TLegend *L=new TLegend(0.65,.12,.89,.4,"Chg Mult QC");
L->AddEntry("P2data","data");
L->AddEntry("P2P6","P6");
L->AddEntry("P2Hpp","Hpp");
L->AddEntry("P2P6_noTrigger","P6 (noTrigger)");
L->AddEntry("P2Hpp_noTrigger","Hpp (noTrigger)");
L->Draw();

P2data->GetXaxis()->SetTitle("P_{T}^{Jet}");
P2data->GetYaxis()->SetTitle("Chg Mult. mean");

p2->Draw();p2->SetBottomMargin(0.05);p2->SetLeftMargin(0.05);p2->SetTopMargin(0);p2->SetRightMargin(0.01);
p2->cd();
TH1F *h2data=new TH1F("h2data","h2data",20,0,20);
TH1F *h2P6=new TH1F("h2P6","h2P6",20,0,20);
TH1F *h2Hpp=new TH1F("h2Hpp","h2Hpp",20,0,20);
TH1F *h2P6_noTrigger=new TH1F("h2P6_noTrigger","h2P6_noTrigger",20,0,20);
TH1F *h2Hpp_noTrigger=new TH1F("h2Hpp_noTrigger","h2Hpp_noTrigger",20,0,20);
data->Draw("jetChgPart_QC[0]>>h2data",(selection+string("&&triggerResult[1] && 40<jetPt[0] && jetPt[0]<45")).c_str(),"P NORM");
mc->Draw("jetChgPart_QC[0]>>h2P6",(selection+string("&&triggerResult[1] && 40<jetPt[0] && jetPt[0]<45")).c_str(),"SAME NORM");
mch->Draw("jetChgPart_QC[0]>>h2Hpp",(selection+string("&&triggerResult[1] && 40<jetPt[0] && jetPt[0]<45")).c_str(),"SAME NORM");
mc->Draw("jetChgPart_QC[0]>>h2P6_noTrigger",(selection+string("&& 40<jetPt[0] && jetPt[0]<45")).c_str(),"SAME NORM");
mch->Draw("jetChgPart_QC[0]>>h2Hpp_noTrigger",(selection+string(" && 40<jetPt[0] && jetPt[0]<45")).c_str(),"SAME NORM");

h2data->SetMarkerStyle(20);
h2P6->SetLineColor(kRed);   
h2Hpp->SetLineColor(kGreen);   
h2P6_noTrigger->SetLineColor(kRed+2);    h2P6_noTrigger ->SetLineWidth(2);
h2Hpp_noTrigger->SetLineColor(kGreen+3); h2Hpp_noTrigger->SetLineWidth(2); 
TLatex *l=new TLatex; l->SetNDC(); l->DrawLatex(0.5,.89,"40<P_{T}[GeV]<45");


//-----------------------------nNeutral-------------------------------------//
TCanvas *c3=new TCanvas("c3","c3",800,800);
TPad *p3=new TPad("p3","p3",.4,.12,.64,.4);
c3->cd();

TProfile *P3data=new TProfile("P3data","Pdata",100,20,100);
TProfile *P3P6  =new TProfile("P3P6","PP6",100,20,100);
TProfile *P3Hpp =new TProfile("P3Hpp","PHpp",100,20,100);
TProfile *P3P6_noTrigger  =new TProfile("P3P6_noTrigger","PP6_noTrigger",100,20,100);
TProfile *P3Hpp_noTrigger =new TProfile("P3Hpp_noTrigger","PHpp_noTrigger",100,20,100);

printf("plot data\n");
data->Draw("jetNeutralPart_ptcut[0]:jetPt[0]>>P3data",(selection+"&&triggerResult[1]").c_str(),"prof");
printf("plot pythia\n");
mc->Draw("jetNeutralPart_ptcut[0]:jetPt[0]>>P3P6",(selection+"&&triggerResult[1]").c_str(),"prof SAME");
printf("plot herwig\n");
mch->Draw("jetNeutralPart_ptcut[0]:jetPt[0]>>P3Hpp",(selection+"&&triggerResult[1]").c_str(),"prof SAME");

mc->Draw("jetNeutralPart_ptcut[0]:jetPt[0]>>P3P6_noTrigger",(selection).c_str(),"prof SAME");
mch->Draw("jetNeutralPart_ptcut[0]:jetPt[0]>>P3Hpp_noTrigger",(selection).c_str(),"prof SAME");

P3P6->SetLineColor(kRed);   
P3Hpp->SetLineColor(kGreen);   
P3P6_noTrigger->SetLineColor(kRed+2);    P3P6_noTrigger ->SetLineWidth(2);
P3Hpp_noTrigger->SetLineColor(kGreen+3); P3Hpp_noTrigger->SetLineWidth(2); 

TLegend *L=new TLegend(0.65,.12,.89,.4,"Neutral Mult QC");
L->AddEntry("P3data","data");
L->AddEntry("P3P6","P6");
L->AddEntry("P3Hpp","Hpp");
L->AddEntry("P3P6_noTrigger","P6 (noTrigger)");
L->AddEntry("P3Hpp_noTrigger","Hpp (noTrigger)");
L->Draw();

P3data->GetXaxis()->SetTitle("P_{T}^{Jet}");
P3data->GetYaxis()->SetTitle("Neutral Mult. mean");

p3->Draw();p3->SetBottomMargin(0.05);p3->SetLeftMargin(0.05);p3->SetTopMargin(0);p3->SetRightMargin(0.01);
p3->cd();
TH1F *h3data=new TH1F("h3data","h3data",15,0,15);
TH1F *h3P6=new TH1F("h3P6","h3P6",15,0,15);
TH1F *h3Hpp=new TH1F("h3Hpp","h3Hpp",15,0,15);
TH1F *h3P6_noTrigger=new TH1F("h3P6_noTrigger","h3P6_noTrigger",15,0,15);
TH1F *h3Hpp_noTrigger=new TH1F("h3Hpp_noTrigger","h3Hpp_noTrigger",15,0,15);
data->Draw("jetNeutralPart_ptcut[0]>>h3data",(selection+string("&&triggerResult[1] && 40<jetPt[0] && jetPt[0]<45")).c_str(),"P NORM");
mc->Draw("jetNeutralPart_ptcut[0]>>h3P6",(selection+string("&&triggerResult[1] && 40<jetPt[0] && jetPt[0]<45")).c_str(),"SAME NORM");
mch->Draw("jetNeutralPart_ptcut[0]>>h3Hpp",(selection+string("&&triggerResult[1] && 40<jetPt[0] && jetPt[0]<45")).c_str(),"SAME NORM");
mc->Draw("jetNeutralPart_ptcut[0]>>h3P6_noTrigger",(selection+string("&& 40<jetPt[0] && jetPt[0]<45")).c_str(),"SAME NORM");
mch->Draw("jetNeutralPart_ptcut[0]>>h3Hpp_noTrigger",(selection+string(" && 40<jetPt[0] && jetPt[0]<45")).c_str(),"SAME NORM");

h3data->SetMarkerStyle(20);
h3P6->SetLineColor(kRed);   
h3Hpp->SetLineColor(kGreen);   
h3P6_noTrigger->SetLineColor(kRed+2);    h3P6_noTrigger ->SetLineWidth(2);
h3Hpp_noTrigger->SetLineColor(kGreen+3); h3Hpp_noTrigger->SetLineWidth(2); 
TLatex *l=new TLatex; l->SetNDC(); l->DrawLatex(0.5,.89,"40<P_{T}[GeV]<45");


//-----------------------------axis1-------------------------------------//
TCanvas *c4=new TCanvas("c4","c4",800,800);
TPad *p4=new TPad("p4","p4",.4,.6,.64,.89);
c4->cd();

TProfile *P4data=new TProfile("P4data","Pdata",100,20,100);
TProfile *P4P6  =new TProfile("P4P6","PP6",100,20,100);
TProfile *P4Hpp =new TProfile("P4Hpp","PHpp",100,20,100);
TProfile *P4P6_noTrigger  =new TProfile("P4P6_noTrigger","PP6_noTrigger",100,20,100);
TProfile *P4Hpp_noTrigger =new TProfile("P4Hpp_noTrigger","PHpp_noTrigger",100,20,100);

printf("plot data\n");
data->Draw("jetAxis_QC[1][0]:jetPt[0]>>P4data",(selection+string("&&triggerResult[1] && jetAxis_QC[1][0]>=0")).c_str(),"prof");
printf("plot pythia\n");
mc->Draw("jetAxis_QC[1][0]:jetPt[0]>>P4P6",(selection+string("&&triggerResult[1] && jetAxis_QC[1][0]>=0")).c_str(),"prof SAME");
printf("plot herwig\n");
mch->Draw("jetAxis_QC[1][0]:jetPt[0]>>P4Hpp",(selection+string("&&triggerResult[1] && jetAxis_QC[1][0]>=0")).c_str(),"prof SAME");

mc->Draw("jetAxis_QC[1][0]:jetPt[0]>>P4P6_noTrigger",(selection+string("&& jetAxis_QC[1][0]>=0")).c_str(),"prof SAME");
mch->Draw("jetAxis_QC[1][0]:jetPt[0]>>P4Hpp_noTrigger",(selection+string("&& jetAxis_QC[1][0]>=0")).c_str(),"prof SAME");

P4P6->SetLineColor(kRed);   
P4Hpp->SetLineColor(kGreen);   
P4P6_noTrigger->SetLineColor(kRed+2);    P4P6_noTrigger ->SetLineWidth(2);
P4Hpp_noTrigger->SetLineColor(kGreen+3); P4Hpp_noTrigger->SetLineWidth(2); 

TLegend *L=new TLegend(0.65,.6,.89,.89,"Axis2 QC");
L->AddEntry("P4data","data");
L->AddEntry("P4P6","P6");
L->AddEntry("P4Hpp","Hpp");
L->AddEntry("P4P6_noTrigger","P6 (noTrigger)");
L->AddEntry("P4Hpp_noTrigger","Hpp (noTrigger)");
L->Draw();

P4data->GetXaxis()->SetTitle("P_{T}^{Jet}");
P4data->GetYaxis()->SetTitle("Axis2 mean");

p4->Draw();p4->SetBottomMargin(0.05);p4->SetLeftMargin(0.05);p4->SetTopMargin(0);p4->SetRightMargin(0.01);
p4->cd();
TH1F *h4data=new TH1F("h4data","h4data",30,0,.2);
TH1F *h4P6=new TH1F("h4P6","h4P6",30,0,.2);
TH1F *h4Hpp=new TH1F("h4Hpp","h4Hpp",30,0,.2);
TH1F *h4P6_noTrigger=new TH1F("h4P6_noTrigger","h4P6_noTrigger",30,0,.2);
TH1F *h4Hpp_noTrigger=new TH1F("h4Hpp_noTrigger","h4Hpp_noTrigger",30,0,.2);
data->Draw("jetAxis_QC[1][0]>>h4data",(selection+string("&&triggerResult[1] && 40<jetPt[0] && jetPt[0]<45 && jetAxis_QC[1][0]>=0")).c_str(),"P NORM");
mc->Draw("jetAxis_QC[1][0]>>h4P6",(selection+string("&&triggerResult[1] && 40<jetPt[0] && jetPt[0]<45 && jetAxis_QC[1][0]>=0")).c_str(),"SAME NORM");
mch->Draw("jetAxis_QC[1][0]>>h4Hpp",(selection+string("&&triggerResult[1] && 40<jetPt[0] && jetPt[0]<45 && jetAxis_QC[1][0]>=0")).c_str(),"SAME NORM");
mc->Draw("jetAxis_QC[1][0]>>h4P6_noTrigger",(selection+string("&& 40<jetPt[0] && jetPt[0]<45 && jetAxis_QC[1][0]>=0")).c_str(),"SAME NORM");
mch->Draw("jetAxis_QC[1][0]>>h4Hpp_noTrigger",(selection+string(" && 40<jetPt[0] && jetPt[0]<45 && jetAxis_QC[1][0]>=0")).c_str(),"SAME NORM");

h4data->SetMarkerStyle(20);
h4P6->SetLineColor(kRed);   
h4Hpp->SetLineColor(kGreen);   
h4P6_noTrigger->SetLineColor(kRed+2);    h4P6_noTrigger ->SetLineWidth(2);
h4Hpp_noTrigger->SetLineColor(kGreen+3); h4Hpp_noTrigger->SetLineWidth(2); 
TLatex *l=new TLatex; l->SetNDC(); l->DrawLatex(0.5,.89,"40<P_{T}[GeV]<45");


}
