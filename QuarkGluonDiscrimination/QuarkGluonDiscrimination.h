#ifndef QGDISCR_H
#define QGDISCR_H
#include "TH1F.h"
#include <vector>
class SpaceStructure
{
public:
	const static int n=4;
	Double_t data[n];
	

};
//Declaration of const static members
const int SpaceStructure::n;

class Discrimination
{
public:
	//Constructor
	Discrimination();
	//Destructor
	~Discrimination();
/********************************
 * OTHER METHOD TO EVAL DENSITY *
 ********************************/
	//Evaluate Density in one point - work in progress
	Double_t DensityFunction(TH1F*h,Double_t value) const;
	//Populate Density in data structure - work in progress
	int Populate(const char *filename,const char type='A',const char*treename="demo/t");
		//type can be 'A' for all, 'Q' for quarks, 'G' for gluons. Need to be set pdgid

/********************************
 * INDIPENDENT VARIABLES APPROX *
 ********************************/

	//Evaluate the product of four density -> works only if they are indipendent
	Double_t IndipendentDensityFunction(TH1F *h1, TH1F*h2,TH1F*h3,TH1F*h4,
						const Double_t value1,
						const Double_t value2, 
						const Double_t value3, 
						const Double_t value4
						) const;
	Double_t IndipendentDensityFunction(const char* filename,
						const Double_t value1, 
						const Double_t value2, 
						const Double_t value3, 
						const Double_t value4
						) const;
	//Evaluate Likelihood as for indipendent variables
	Double_t IndipendentLikelihood(const char* filename,
                                                const Double_t value1,
                                                const Double_t value2,
                                                const Double_t value3,
                                                const Double_t value4
                                                ) const;
	Double_t IndipendentLikelihood(const Double_t jtpt,
						const Double_t value1,
                                                const Double_t value2,
                                                const Double_t value3,
                                                const Double_t value4,
						const char*filename="VarDistribution/VarDistribution_%03.0lf_%03.0lf.root"
						)const;
protected:
private:
	std::vector<SpaceStructure*> *Space;
};

Discrimination::Discrimination()
	{
	Space=new std::vector<SpaceStructure*> ;
	}
Discrimination::~Discrimination()
	{
	while (!Space->empty())
		{
		delete Space->at(Space->size()-1);
		Space->pop_back();
		}
	delete Space;
	}

#endif
