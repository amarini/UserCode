#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

//#include "../interface/QGLikelihoodCalculator.h"




class QuarkGluonTagger2012 : public edm::EDProducer {
   public:
      explicit QuarkGluonTagger2012(const edm::ParameterSet&);
      ~QuarkGluonTagger2012();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data --------------------------
      edm::InputTag src_,srcRho_,srcRho2_;
      std::string jecService_;
      bool isPatJet_;
  //    QGLikelihoodCalculator *qglikeli_;
      const JetCorrector *JEC_;           
};

