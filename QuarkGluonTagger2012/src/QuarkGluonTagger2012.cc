
/*
 * Original Author: Andrea Carlo Marini
 * Email: andrea.carlo.marini@cern.ch
 * Date: Jan 2013
 *
 *
 *
 */
#include <memory>
#include <iostream>
#include <vector>

#include "../interface/QuarkGluonTagger2012.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TMath.h"



QuarkGluonTagger2012::QuarkGluonTagger2012(const edm::ParameterSet& iConfig)
{
        src_        = iConfig.getParameter<edm::InputTag> ("jets");
        srcRho_     = iConfig.getParameter<edm::InputTag> ("rho");
        srcRho2_     = iConfig.getParameter<edm::InputTag> ("rho2");
        jecService_ = iConfig.getParameter<std::string>   ("jec");
	isPatJet_ = iConfig.existsAs<bool>("isPatJet") ? iConfig.getParameter<bool>("isPatJet") : false ; 
       
        produces<edm::ValueMap<float> >().setBranchAlias("qg");
        qglikeli_ = new QGLikelihoodCalculator(); 
	debug=0;
}

QuarkGluonTagger2012::~QuarkGluonTagger2012(){}

void QuarkGluonTagger2012::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
using namespace edm;
using namespace reco;
using namespace std;

//------ Jets -----------------------------------------------------------------------------
edm::Handle<reco::PFJetCollection> pfjets;
edm::Handle<vector<pat::Jet> > patjets;
if(!isPatJet_)
	iEvent.getByLabel(src_,pfjets);
else
	iEvent.getByLabel(src_,patjets);

vector<float> values;
if(!isPatJet_) 	values.reserve(pfjets->size());
else  		values.reserve(patjets->size());
//------ Rho -----------------------------------------------------------------------------
edm::Handle<double> rho;
iEvent.getByLabel(srcRho_,rho);

//------ JEC -----------------------------------------------------------------------------
JEC_ = JetCorrector::getJetCorrector(jecService_,iSetup);
//------ Vtxs -----------------------------------------------------------------------------
edm::Handle<reco::VertexCollection> recVtxs;
iEvent.getByLabel("goodOfflinePrimaryVertices",recVtxs);


//------ Variables -----------------------------------------------------------------------------
float QGL;

//RECO
if(!isPatJet_)
for(reco::PFJetCollection::const_iterator ijet=pfjets->begin();ijet!=pfjets->end();ijet++)
	{
	//JEC
		if(debug>1) std::cout<<"BEGIN RECO JET"<<endl;
	double cor = JEC_->correction(*ijet,iEvent,iSetup);
	double corPt = cor*ijet->pt();
	
//	nChg = ijet->getChargedHadronMultiplicity();
//	nNeutral = ijet->neutralHadronMultiplicity()+ ijet->photonMultiplicity();
//	ptD= ijet->constituentPtDistribution();
	
    //----- calculations based on the constituents -----
//----- calculations based on the constituents -----
    std::vector<PFCandidatePtr> pfConst(ijet->getPFConstituents());
    int n_pf = pfConst.size();                                                                                          
    float phiJet = ijet->phi();                                                                                         
    float etaJet = ijet->eta();                                                                                         
    float deta,dphi,weight,weight2,sumW(0.0),sumW2(0.0),sum_deta(0.0),sum_dphi(0.0),sum_deta2(0.0),sum_dphi2(0.0),sum_detadphi(0.0);
    //UNUSED---   float dR;
    //UNUSED---   float Teta(0),Tphi(0),Ttheta(-9),jetPtMax(0),axis1(-999),axis2(-999),tana(-999),ptD(-999);                      
    
    float sumW_QC(0.0),sumW2_QC(0.0),sum_deta_QC(0.0),sum_dphi_QC(0.0),sum_deta2_QC(0.0),sum_dphi2_QC(0.0),sum_detadphi_QC(0.0);
    float axis1_QC(-999),axis2_QC(-999),ptD_QC(-999);   

    //UNUSED---   float ave_deta(0.0),ave_dphi(0.0),ave_deta2(0.0),ave_dphi2(0.0);
    float ave_deta_QC(0.0),ave_dphi_QC(0.0),ave_deta2_QC(0.0),ave_dphi2_QC(0.0);
    //UNUSED---   float pull(0.0);
    //UNUSED---   float pull_QC(0.0);

    float pTMax(0.0),pTMaxChg(0.0),pTMaxNeutral(0.0),pTMaxChg_QC(0.0);
    int nChg_QC(0),nChg_ptCut(0),nChg_ptCut_QC(0),nNeutral_ptCut(0);
    std::vector<float> jetPart_pt,jetPart_deta,jetPart_dphi;
    std::vector<bool> jetPart_forMult,jetPart_forAxis;
     
    for(int j=0;j<n_pf;j++) {                                                                 
      reco::TrackRef itrk ;
      reco::PFCandidatePtr  part = ijet->getPFConstituent(j);
      if (part.isNonnull())
        itrk = (*part).trackRef();
      if (part->pt() > pTMax) 
        pTMax = part->pt();
      if (itrk.isNonnull() && part->pt() > pTMaxChg) 
        pTMaxChg = part->pt();
      if (!itrk.isNonnull() && part->pt() > pTMaxNeutral) 
        pTMaxNeutral = part->pt();
      if (!itrk.isNonnull() && part->pt() > 1.0) 
        nNeutral_ptCut++;
 
      bool trkForAxis = false;
      bool trkForMult = false;

      //-----matching with vertex tracks-------
      if (!itrk.isNonnull()) { 
        trkForMult = true;
        trkForAxis = true;
      }
      else {
        if (part->pt() > 1.0)
          nChg_ptCut++;
        float dZmin = 999;
        int index_min = 999;
        reco::VertexCollection::const_iterator vtxClose;
        for(unsigned ivtx = 0;ivtx < recVtxs->size();ivtx++) {
          float dZ_cut = fabs(itrk->dz((*recVtxs)[ivtx].position()));
          float sumpT = 0;
          for(reco::Vertex::trackRef_iterator itk = (*recVtxs)[ivtx].tracks_begin();itk!=(*recVtxs)[ivtx].tracks_end(); ++itk) {
            sumpT = sumpT + ((*itk)->pt())*((*itk)->pt());
          }
          if (dZ_cut < dZmin) {
            dZmin = dZ_cut;
            index_min = ivtx;
              //  std::cout<<"dz=="<<dZ_cut<<std::endl;
          }
        }//Loop over vertices 
        if (index_min == 0) {
          float dz = itrk->dz((*recVtxs)[0].position());
          float d0 = itrk->dxy((*recVtxs)[0].position());
          float vtx_xError = (*recVtxs)[0].xError();
          float vtx_yError = (*recVtxs)[0].yError();
          float vtx_zError = (*recVtxs)[0].zError();
          float d0_sigma=sqrt(pow(itrk->d0Error(),2) + pow(vtx_xError,2) + pow(vtx_yError,2));
          float dz_sigma=sqrt(pow(itrk->dzError(),2) + pow(vtx_zError,2));
          if (itrk->quality(reco::TrackBase::qualityByName("highPurity")) && fabs(dz/dz_sigma) < 5.) {
            trkForAxis = true;
            if (fabs(d0/d0_sigma) < 5.)
              trkForMult = true;
          }//
        }
        if (trkForMult)
          nChg_QC++;
        if (itrk.isNonnull() && trkForMult && part->pt() > 1.0)
          nChg_ptCut_QC++;
        if (part->pt() > pTMaxChg_QC && trkForAxis) 
          pTMaxChg_QC = part->pt();
      }// for charged particles only

      //-----------Store part info-----------------------
      jetPart_pt.push_back(pfConst[j]->pt());
      jetPart_forMult.push_back(trkForMult);
      jetPart_forAxis.push_back(trkForAxis);
  
      deta = pfConst[j]->eta() - etaJet;                                                                                               
      dphi = 2*atan(tan((pfConst[j]->phi()-phiJet)/2));                                                                              
      //UNUSED---   dR = sqrt(deta*deta + dphi*dphi);                                                                                              
      weight = pfConst[j]->pt(); // used for the thrust and ptD variables
      weight2 = weight * weight; // used for the jet axis variables                                                          
      sumW += weight;
      sumW2 += weight2;                                                                                                   
      //UNUSED---   Teta += weight * dR * deta;
      //UNUSED---   Tphi += weight * dR * dphi;
      sum_deta += deta*weight2;
      sum_dphi += dphi*weight2;
      sum_deta2 += deta*deta*weight2;
      sum_detadphi += deta*dphi*weight2;
      sum_dphi2 += dphi*dphi*weight2;
      //-----Axis using charged tracks with quality cuts--- 
      if (trkForAxis) {
        sumW2_QC += weight2;
        sumW_QC += weight;
        sum_deta_QC += deta*weight2;
        sum_dphi_QC += dphi*weight2;
        sum_deta2_QC += deta*deta*weight2;
        sum_detadphi_QC += deta*dphi*weight2;
        sum_dphi2_QC += dphi*dphi*weight2;
      }
      jetPart_deta.push_back(deta);
      jetPart_dphi.push_back(dphi);
      //-----------------------------------      
      //UNUSED---   if (fabs(pfConst[j]->charge()) > 0) {
        //UNUSED---   jetPtMax = TMath::Max(jetPtMax,float(pfConst[j]->pt()));
      //UNUSED---   }
    }// loop over the constituents
    //UNUSED---  if (sumW > 0) {
    //UNUSED---    Teta = Teta/sumW;
    //UNUSED---    Tphi = Tphi/sumW;
    //UNUSED---    if (Teta != 0 && Tphi !=0 ) {
    //UNUSED---      Ttheta = atan2(Tphi,Teta);
    //UNUSED---    }
    //UNUSED---    ptD = sqrt(sumW2)/sumW;
    //UNUSED---    ave_deta = sum_deta/sumW2;
    //UNUSED---    ave_dphi = sum_dphi/sumW2;
    //UNUSED---    ave_deta2 = sum_deta2/sumW2;
    //UNUSED---    ave_dphi2 = sum_dphi2/sumW2;
    //UNUSED---    float a = ave_deta2-ave_deta*ave_deta;
    //UNUSED---    float b = ave_dphi2-ave_dphi*ave_dphi;
    //UNUSED---    float c = -(sum_detadphi/sumW2-ave_deta*ave_dphi);
    //UNUSED---    float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
    //UNUSED---    if (a+b+delta > 0) {
    //UNUSED---      axis1 = sqrt(0.5*(a+b+delta));
    //UNUSED---    }
    //UNUSED---    if (a+b-delta > 0) {  
    //UNUSED---      axis2 = sqrt(0.5*(a+b-delta));
    //UNUSED---    }
    //UNUSED---    if (c != 0) {
    //UNUSED---      tana = 0.5*(b-a+delta)/c;
    //UNUSED---    }
    //UNUSED---  }
    if (sumW_QC > 0) {
      ave_deta_QC = sum_deta_QC/sumW2_QC;
      ave_dphi_QC = sum_dphi_QC/sumW2_QC;
      ave_deta2_QC = sum_deta2_QC/sumW2_QC;
      ave_dphi2_QC = sum_dphi2_QC/sumW2_QC;
      ptD_QC = sqrt(sumW2_QC)/sumW_QC;

      float a_QC = ave_deta2_QC-ave_deta_QC*ave_deta_QC;
      float b_QC = ave_dphi2_QC-ave_dphi_QC*ave_dphi_QC;
      float c_QC = -(sum_detadphi_QC/sumW2_QC-ave_deta_QC*ave_dphi_QC);
      float delta_QC = sqrt(fabs((a_QC-b_QC)*(a_QC-b_QC)+4*c_QC*c_QC));

      if (a_QC+b_QC+delta_QC > 0) {
        axis1_QC = sqrt(0.5*(a_QC+b_QC+delta_QC));
      }
      if (a_QC+b_QC-delta_QC > 0) {
        axis2_QC = sqrt(0.5*(a_QC+b_QC-delta_QC));
      }
    }
    //UNUSED---   int nChg  = ijet->chargedMultiplicity();
    //UNUSED---   int nNeutral  = ijet->neutralMultiplicity();
	
	//------------variables defined compute QGL
	QGL=-1.0;
	int nPFCand_QC_ptCut=nChg_QC+nNeutral_ptCut;
		if(debug>1) std::cout<<"--computeQGL"<<endl;
	QGL=qglikeli_->computeQGLikelihood(corPt,*rho,ijet->eta(),nPFCand_QC_ptCut, ptD_QC ,axis1_QC, axis2_QC);
	values.push_back(QGL);
	}
if(isPatJet_)
for(vector<pat::Jet>::const_iterator ijet=patjets->begin();ijet!=patjets->end();ijet++)
{
//----- calculations based on the constituents -----
		if(debug>1) std::cout<<"BEGIN PAT JET"<<endl;
    std::vector<PFCandidatePtr> pfConst(ijet->getPFConstituents());
    int n_pf = pfConst.size();                                                                                          
    float phiJet = ijet->phi();                                                                                         
    float etaJet = ijet->eta();                                                                                         
    float deta,dphi,weight,weight2,sumW(0.0),sumW2(0.0),sum_deta(0.0),sum_dphi(0.0),sum_deta2(0.0),sum_dphi2(0.0),sum_detadphi(0.0);              
    //UNUSED---   float dR;
    //UNUSED---   float Teta(0),Tphi(0),Ttheta(-9),jetPtMax(0),axis1(-999),axis2(-999),tana(-999),ptD(-999);                      
    
    float sumW_QC(0.0),sumW2_QC(0.0),sum_deta_QC(0.0),sum_dphi_QC(0.0),sum_deta2_QC(0.0),sum_dphi2_QC(0.0),sum_detadphi_QC(0.0);
    float axis1_QC(-999),axis2_QC(-999),ptD_QC(-999);   

    //UNUSED---   float ave_deta(0.0),ave_dphi(0.0),ave_deta2(0.0),ave_dphi2(0.0);
    float ave_deta_QC(0.0),ave_dphi_QC(0.0),ave_deta2_QC(0.0),ave_dphi2_QC(0.0);
    //UNUSED---   float pull(0.0);
    //UNUSED---   float pull_QC(0.0);

    float pTMax(0.0),pTMaxChg(0.0),pTMaxNeutral(0.0),pTMaxChg_QC(0.0);
    int nChg_QC(0),nChg_ptCut(0),nChg_ptCut_QC(0),nNeutral_ptCut(0);
    std::vector<float> jetPart_pt,jetPart_deta,jetPart_dphi;
    std::vector<bool> jetPart_forMult,jetPart_forAxis;
     
    for(int j=0;j<n_pf;j++) {                                                                 
      reco::TrackRef itrk ;
      reco::PFCandidatePtr  part = ijet->getPFConstituent(j);
      if (part.isNonnull())
        itrk = (*part).trackRef();
      if (part->pt() > pTMax) 
        pTMax = part->pt();
      if (itrk.isNonnull() && part->pt() > pTMaxChg) 
        pTMaxChg = part->pt();
      if (!itrk.isNonnull() && part->pt() > pTMaxNeutral) 
        pTMaxNeutral = part->pt();
      if (!itrk.isNonnull() && part->pt() > 1.0) 
        nNeutral_ptCut++;
 
      bool trkForAxis = false;
      bool trkForMult = false;

      //-----matching with vertex tracks-------
      if (!itrk.isNonnull()) { 
        trkForMult = true;
        trkForAxis = true;
      }
      else {
        if (part->pt() > 1.0)
          nChg_ptCut++;
        float dZmin = 999;
        int index_min = 999;
        reco::VertexCollection::const_iterator vtxClose;
        for(unsigned ivtx = 0;ivtx < recVtxs->size();ivtx++) {
          float dZ_cut = fabs(itrk->dz((*recVtxs)[ivtx].position()));
          float sumpT = 0;
          for(reco::Vertex::trackRef_iterator itk = (*recVtxs)[ivtx].tracks_begin();itk!=(*recVtxs)[ivtx].tracks_end(); ++itk) {
            sumpT = sumpT + ((*itk)->pt())*((*itk)->pt());
          }
          if (dZ_cut < dZmin) {
            dZmin = dZ_cut;
            index_min = ivtx;
              //  std::cout<<"dz=="<<dZ_cut<<std::endl;
          }
        }//Loop over vertices 
        if (index_min == 0) {
          float dz = itrk->dz((*recVtxs)[0].position());
          float d0 = itrk->dxy((*recVtxs)[0].position());
          float vtx_xError = (*recVtxs)[0].xError();
          float vtx_yError = (*recVtxs)[0].yError();
          float vtx_zError = (*recVtxs)[0].zError();
          float d0_sigma=sqrt(pow(itrk->d0Error(),2) + pow(vtx_xError,2) + pow(vtx_yError,2));
          float dz_sigma=sqrt(pow(itrk->dzError(),2) + pow(vtx_zError,2));
          if (itrk->quality(reco::TrackBase::qualityByName("highPurity")) && fabs(dz/dz_sigma) < 5.) {
            trkForAxis = true;
            if (fabs(d0/d0_sigma) < 5.)
              trkForMult = true;
          }//
        }
        if (trkForMult)
          nChg_QC++;
        if (itrk.isNonnull() && trkForMult && part->pt() > 1.0)
          nChg_ptCut_QC++;
        if (part->pt() > pTMaxChg_QC && trkForAxis) 
          pTMaxChg_QC = part->pt();
      }// for charged particles only

      //-----------Store part info-----------------------
      jetPart_pt.push_back(pfConst[j]->pt());
      jetPart_forMult.push_back(trkForMult);
      jetPart_forAxis.push_back(trkForAxis);
  
      deta = pfConst[j]->eta() - etaJet;                                                                                               
      dphi = 2*atan(tan((pfConst[j]->phi()-phiJet)/2));                                                                              
      //UNUSED---   dR = sqrt(deta*deta + dphi*dphi);                                                                                              
      weight = pfConst[j]->pt(); // used for the thrust and ptD variables
      weight2 = weight * weight; // used for the jet axis variables                                                          
      sumW += weight;
      sumW2 += weight2;                                                                                                   
      //UNUSED---   Teta += weight * dR * deta;
      //UNUSED---   Tphi += weight * dR * dphi;
      sum_deta += deta*weight2;
      sum_dphi += dphi*weight2;
      sum_deta2 += deta*deta*weight2;
      sum_detadphi += deta*dphi*weight2;
      sum_dphi2 += dphi*dphi*weight2;
      //-----Axis using charged tracks with quality cuts--- 
      if (trkForAxis) {
        sumW2_QC += weight2;
        sumW_QC += weight;
        sum_deta_QC += deta*weight2;
        sum_dphi_QC += dphi*weight2;
        sum_deta2_QC += deta*deta*weight2;
        sum_detadphi_QC += deta*dphi*weight2;
        sum_dphi2_QC += dphi*dphi*weight2;
      }
      jetPart_deta.push_back(deta);
      jetPart_dphi.push_back(dphi);
      //-----------------------------------      
      //UNUSED---   if (fabs(pfConst[j]->charge()) > 0) {
        //UNUSED---   jetPtMax = TMath::Max(jetPtMax,float(pfConst[j]->pt()));
      //UNUSED---   }
    }// loop over the constituents
    //UNUSED---   if (sumW > 0) {
      //UNUSED---   Teta = Teta/sumW;
      //UNUSED---   Tphi = Tphi/sumW;
      //UNUSED---   if (Teta != 0 && Tphi !=0 ) {
      //UNUSED---   Ttheta = atan2(Tphi,Teta);
      //UNUSED---   }
      //UNUSED---   ptD = sqrt(sumW2)/sumW;
      //UNUSED---   ave_deta = sum_deta/sumW2;
      //UNUSED---   ave_dphi = sum_dphi/sumW2;
      //UNUSED---   ave_deta2 = sum_deta2/sumW2;
      //UNUSED---   ave_dphi2 = sum_dphi2/sumW2;
      //UNUSED---   float a = ave_deta2-ave_deta*ave_deta;
      //UNUSED---   float b = ave_dphi2-ave_dphi*ave_dphi;
      //UNUSED---   float c = -(sum_detadphi/sumW2-ave_deta*ave_dphi);
      //UNUSED---   float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
      //UNUSED---   if (a+b+delta > 0) {
        //UNUSED---   axis1 = sqrt(0.5*(a+b+delta));
      //UNUSED---   }
      //UNUSED---   if (a+b-delta > 0) {  
      //UNUSED---   axis2 = sqrt(0.5*(a+b-delta));
      //UNUSED---   }
      //UNUSED---   if (c != 0) {
      //UNUSED---   tana = 0.5*(b-a+delta)/c;
      //UNUSED---   }
    //UNUSED---   }
    if (sumW_QC > 0) {
      ave_deta_QC = sum_deta_QC/sumW2_QC;
      ave_dphi_QC = sum_dphi_QC/sumW2_QC;
      ave_deta2_QC = sum_deta2_QC/sumW2_QC;
      ave_dphi2_QC = sum_dphi2_QC/sumW2_QC;

      ptD_QC = sqrt(sumW2_QC)/sumW_QC;

      float a_QC = ave_deta2_QC-ave_deta_QC*ave_deta_QC;
      float b_QC = ave_dphi2_QC-ave_dphi_QC*ave_dphi_QC;
      float c_QC = -(sum_detadphi_QC/sumW2_QC-ave_deta_QC*ave_dphi_QC);
      float delta_QC = sqrt(fabs((a_QC-b_QC)*(a_QC-b_QC)+4*c_QC*c_QC));

      if (a_QC+b_QC+delta_QC > 0) {
        axis1_QC = sqrt(0.5*(a_QC+b_QC+delta_QC));
      }
      if (a_QC+b_QC-delta_QC > 0) {
        axis2_QC = sqrt(0.5*(a_QC+b_QC-delta_QC));
      }
    }
    //UNUSED---   int nChg  = ijet->chargedMultiplicity();
    //UNUSED---   int nNeutral  = ijet->neutralMultiplicity();
   
	
	//------------variables defined compute QGL
	QGL=-1.0;
	int nPFCand_QC_ptCut=nChg_QC+nNeutral_ptCut;
		if(debug>1) std::cout<<"--computeQGL"<<endl;
	QGL=qglikeli_->computeQGLikelihood(ijet->pt(),*rho,ijet->eta(),nPFCand_QC_ptCut, ptD_QC ,axis1_QC, axis2_QC);
		if(debug>0) std::cout<<"--- "<<QGL <<endl;
		if(debug>1) std::cout<<"--- Pt"<<ijet->pt() <<endl;
		if(debug>1) std::cout<<"--- rho"<<*rho <<endl;
	values.push_back(QGL);
	
} //for - pat Jet

//common code
std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
edm::ValueMap<float>::Filler filler(*out);
if(!isPatJet_)
	filler.insert(pfjets,values.begin(),values.end());
else 
	filler.insert(patjets,values.begin(),values.end());
filler.fill();
iEvent.put(out);

}


// ------------ method called once each job just before starting event loop  ------------
void
QuarkGluonTagger2012::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
QuarkGluonTagger2012::endJob() {
}



//define this as a plug-in
DEFINE_FWK_MODULE(QuarkGluonTagger2012);
