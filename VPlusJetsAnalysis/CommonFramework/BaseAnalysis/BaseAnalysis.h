#include <map>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"

#ifndef BASE_ANALYSIS_H
#define BASE_ANALYSIS_H
using namespace std;

class BaseAnalysis
{
public:
	BaseAnalysis();
	virtual ~BaseAnalysis();
	int Loop(int type,int j=-1,int NJobs=-1);
	virtual int CreateHistos();
	int Write(string outputFile="");
	int Sumw2();
	//virtual int SetActive(); //read list of active Branches from a configuration file
	virtual int FillHistos(){};
	virtual int FillHistosGEN(){};
	virtual int OpenFiles();
	virtual bool BaseSelection(){return true;};
	int SetConfigFile(string a,string b="");
	int BookFile(string a);
	virtual int Smear();
	char* Read(const char*fileName,const char * parName,char*R=NULL);
	char* ReadMult(const char *fileName1,const char*fileName2,const char*parName);
	int SetSmearType(int smeartype,string label) {SmearType=smeartype;extraLabel=label;}
	int SetCuts();
	int GetChId(){return -CHID2;}
protected:	
	//tree - var
	double weight;
	int nVtx;
	int nLeptons;
        int nPhotons;
        int nJets;
        int nJetsVeto;
        float rho;
        float rho25;
        float llM;
        float llPt;
        float llPhi;
        float llY;
        float llDPhi;
        float llEta;
        vector<float> *lepPt       ;
        vector<float> *lepEta      ;
        vector<float> *lepPhi      ;
        vector<float> *lepE        ;
        vector<int>   *lepChId     ;
        vector<float> *jetPt       ;
        vector<float> *jetEta      ;
        vector<float> *jetPhi      ;
        vector<float> *jetE        ;
        vector<float> *jetUNC      ;
        vector<float> *jetQGL      ;
        vector<float> *jetBtag     ;
        vector<float> *jetBeta     ;
        vector<int>   *jetVeto     ;
        vector<float> *photonPt;
        vector<float> *photonEta;
        vector<float> *photonPhi;
        vector<float> *photonE;
        double PUWeight;
	
	///GEN
	float llMGEN;
	float llPtGEN;
	float llYGEN;
	vector<float> *lepPtGEN; 		
	vector<float> *lepEtaGEN; 		
	vector<float> *lepPhiGEN; 		
	vector<float> *lepEGEN; 		
	vector<int> *lepChIdGEN; 		
	vector<float> *jetPtGEN; 		
	vector<float> *jetEtaGEN; 		
	vector<float> *jetPhiGEN; 		
	vector<float> *jetEGEN; 		
//protected:
	//histograms
	map<string,TH1F*> histos;
	map<string,TH2F*> histos2;
	map<string,RooUnfoldResponse*> response;
	string extraLabel; //useful for syst smear
	int SmearType;
	//
	vector<int> jet;
	vector<int> jet_BS;
	
	string configFileBase;
	string configFileNew;
	vector<string> bookedFiles;
	TChain *t;

	
	int debug;
	//cut
	int CHID2;
	float llMCut;
	float JetPtCut;
	float JetEtaCut;
	float JetDRCut;
	float LepPtCut;
	float LepEtaCut;
	float ZYCut;
	
	int JetSel; // 3 = Z | 4 = gamma |
	
};
#endif
