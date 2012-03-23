#include "TLatex.h"
#include "CMSLatex.h"
void CMSLatex::DrawSimulation(){
TLatex *lat=new TLatex();
lat->SetNDC();//normalized coordinate syst ref to pad
lat->SetTextSize(0.04);
lat->SetTextAlign(23);
lat->DrawLatex(0.25,.89,"CMS Simulation");
lat->DrawLatex(0.25,.84,"#sqrt{s} = 7 TeV");
}
void CMSLatex::DrawPreliminary(){
TLatex *lat=new TLatex();
lat->SetNDC();//normalized coordinate syst ref to pad
lat->SetTextSize(0.04);
lat->SetTextAlign(23);
lat->DrawLatex(0.25,.89,"CMS Preliminary");
lat->DrawLatex(0.25,.84,"#sqrt{s} = 7 TeV");
}
