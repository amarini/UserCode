#include "../BaseAnalysis/BaseAnalysis.h"


class QGDumper: public BaseAnalysis
{
public:
	QGDumper();
	~QGDumper();

	int CreateHistos(); //create TTree

	int FillHistos(); //FillTTree
	
	int SetActive();// SetBranches

	bool BaseSelection();
	bool BaseSelectionGEN();

//	int FillHistosGEN();

	int Smear();
	
	int Write(string outFile="");
	

	string treeName;
	string outFileName;

protected:
	bool selGen;
	TTree *t_;
	TFile *outFile;
	
	float out_QGLMLP;
	float out_QGLHisto;
	float out_ptJet0;
	float out_etaJet0;
	float out_phiJet0;
	float out_deltaPhi;
	float out_PUReWeight;
	int   out_pdgIdPartJet0;
};

