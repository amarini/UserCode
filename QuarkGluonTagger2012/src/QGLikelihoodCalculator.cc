#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "../interface/QGLikelihoodCalculator.h"
#include <FWCore/Utilities/interface/Exception.h>

#include "TMath.h"

using namespace std;

//--------- Constructor ----------------------------------------------------------------------------------------
QGLikelihoodCalculator::QGLikelihoodCalculator( const std::string& dirName ) {

  //map<string,JetCorrectorParameters*> JCP;
  //map<string,SimpleJetCorrector*> SJC;
	JCP.clear();
	SJC.clear();
	names.clear();
	
	names.push_back("ptD_QC");
	names.push_back("nPFCand_QC_ptCut");
	names.push_back("axis1_QC");
	names.push_back("axis2_QC");
 
	
	for(vector<string>::iterator iStr=names.begin();iStr!=names.end();iStr++){
		string path=edm::FileInPath( (dirName+ (*iStr)+"Jet0.txt").c_str() ).fullPath();
		JCP[(*iStr)+".quark"] = new JetCorrectorParameters(path,"quark");
		JCP[(*iStr)+".gluon"] = new JetCorrectorParameters(path,"gluon");
		
		path=edm::FileInPath( (dirName+ (*iStr)+"Jet0_F.txt").c_str() ).fullPath();
		JCP[(*iStr)+".F.quark"] = new JetCorrectorParameters(path,"quark");
		JCP[(*iStr)+".F.gluon"] = new JetCorrectorParameters(path,"gluon");
	}
  //check that provided files are for correct variables:
	
	for(vector<string>::iterator iStr=names.begin();iStr!=names.end();iStr++){
	
		if(JCP[(*iStr)+".quark"]->definitions().level() != (string("QGL")+(*iStr)+"_quark").c_str() )
			throw cms::Exception("QuarkGluonTagger Config File Error")<<"quark section of file\'"<< *iStr <<"\' is not of the proper format. Check input files.";
		if(JCP[(*iStr)+".gluon"]->definitions().level() != (string("QGL")+(*iStr)+"_gluon").c_str() )
			throw cms::Exception("QuarkGluonTagger Config File Error")<<"gluon section of file\'"<< *iStr <<"\' is not of the proper format. Check input files.";
		if(JCP[(*iStr)+".F.quark"]->definitions().level() != (string("QGL")+(*iStr)+"_F_quark").c_str() )
			throw cms::Exception("QuarkGluonTagger Config File Error")<<"quark section of file\'"<< *iStr <<"_F\' is not of the proper format. Check input files.";
		if(JCP[(*iStr)+".F.gluon"]->definitions().level() != (string("QGL")+(*iStr)+"_F_gluon").c_str() )
			throw cms::Exception("QuarkGluonTagger Config File Error")<<"gluon section of file\'"<< *iStr <<"_F\' is not of the proper format. Check input files.";
	}

	for(vector<string>::iterator iStr=names.begin();iStr!=names.end();iStr++){
		SJC[(*iStr)+".quark"] = new SimpleJetCorrector(*JCP[(*iStr)+".quark"]);
		SJC[(*iStr)+".gluon"] = new SimpleJetCorrector(*JCP[(*iStr)+".gluon"]);
		
		SJC[(*iStr)+".F.quark"] = new SimpleJetCorrector((*iStr)+".F.quark");
		SJC[(*iStr)+".F.gluon"] = new SimpleJetCorrector((*iStr)+".F.gluon");
	}

}

//--------- Destructor ----------------------------------------------------------------------------------------
QGLikelihoodCalculator::~QGLikelihoodCalculator() {

  for(map<string,JetCorrectorParameters*>::iterator it= JCP.begin(); it!=JCP.end(); it++)
	delete it->second;
  for(map<string,SimpleJetCorrector*>::iterator it= SJC.begin();it!=SJC.end();it++)
	delete it->second;
	
 JCP.clear();
 SJC.clear();

}

float QGLikelihoodCalculator::computeQGLikelihood( float pt, float rhoPF, float eta, int nPFCand_QC_ptCut, float ptD_QC,float axis1_QC,float axis2_QC ) {


  std::vector<float> v_pt_rho;
  v_pt_rho.push_back( pt );
  v_pt_rho.push_back( rhoPF );

  std::vector<float> v_nPFCand_QC_ptCut;
  v_nPFCand_QC_ptCut.push_back( (float)nPFCand_QC_ptCut);

  if(axis1_QC<=0) return -2;
  std::vector<float> v_axis1_QC;
  v_axis1_QC.push_back( TMath::Log((float)axis1_QC) );

  if(axis2_QC<=0) return -2.1;
  std::vector<float> v_axis2_QC;
  v_axis2_QC.push_back( TMath::Log((float)axis2_QC) );

  std::vector<float> v_ptD_QC;
  v_ptD_QC.push_back( (float)ptD_QC );

	float qProb=1.0;
	float gProb=1.0;
	
	if(fabs(eta)<2.5){ //CENTRAL REGION
	qProb*=SJC[string("ptD_QC_quark")]->correction(v_pt_rho,v_ptD_QC);
	gProb*=SJC[string("ptD_QC_gluon")]->correction(v_pt_rho,v_ptD_QC);
	
	qProb*=SJC[string("axis1_QC_quark")]->correction(v_pt_rho,v_axis1_QC);
	gProb*=SJC[string("axis1_QC_gluon")]->correction(v_pt_rho,v_axis1_QC);
	
	qProb*=SJC[string("axis2_QC_quark")]->correction(v_pt_rho,v_axis2_QC);
	gProb*=SJC[string("axis2_QC_gluon")]->correction(v_pt_rho,v_axis2_QC);

	qProb*=SJC[string("nPFCand_QC_ptCut_quark")]->correction(v_pt_rho,v_nPFCand_QC_ptCut);
	gProb*=SJC[string("nPFCand_QC_ptCut_gluon")]->correction(v_pt_rho,v_nPFCand_QC_ptCut);
	}
	else if(fabs(eta)>=2.5) {
	qProb*=SJC[string("ptD_QC_F_quark")]->correction(v_pt_rho,v_ptD_QC);
	gProb*=SJC[string("ptD_QC_F_gluon")]->correction(v_pt_rho,v_ptD_QC);
	
	qProb*=SJC[string("axis1_QC_F_quark")]->correction(v_pt_rho,v_axis1_QC);
	gProb*=SJC[string("axis1_QC_F_gluon")]->correction(v_pt_rho,v_axis1_QC);
	
	qProb*=SJC[string("axis2_QC_F_quark")]->correction(v_pt_rho,v_axis2_QC);
	gProb*=SJC[string("axis2_QC_F_gluon")]->correction(v_pt_rho,v_axis2_QC);

	qProb*=SJC[string("nPFCand_QC_F_ptCut_quark")]->correction(v_pt_rho,v_nPFCand_QC_ptCut);
	gProb*=SJC[string("nPFCand_QC_F_ptCut_gluon")]->correction(v_pt_rho,v_nPFCand_QC_ptCut);
	}



  float QGLikelihood = ( (qProb+gProb) > 0) ? qProb/(qProb+gProb) : -1 ;

  return QGLikelihood;

}



