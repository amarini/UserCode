#include "../BaseAnalysis/BaseAnalysis.h"


class ZAnalysis: public BaseAnalysis
{
public:
	ZAnalysis();
	~ZAnalysis();

	int CreateHistos();
	//int SetActive(); //read list of active Branches from a configuration file
	int FillHistos();
//	int OpenFiles();
	bool BaseSelection();
	bool BaseSelectionGEN();
	int FillMiss();
	int FillHistosGEN();
	int Smear();
	
protected:
	bool selGen;
};

