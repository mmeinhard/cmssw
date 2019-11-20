// system include file
using namespace std;

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

class HTTSubjetFlavourTableProducer : public edm::stream::EDProducer<> {
    public:
        explicit HTTSubjetFlavourTableProducer(const edm::ParameterSet &iConfig) :
            name_(iConfig.getParameter<std::string>("name")),
            src_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src"))),
            cut_(iConfig.getParameter<std::string>("cut"), true),
            deltaR_(iConfig.getParameter<double>("deltaR")),
            subjetFlavourInfosToken_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("subjetFlavourInfos")))
        {
            produces<nanoaod::FlatTable>();
        }

        ~HTTSubjetFlavourTableProducer() override {};

        static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
            edm::ParameterSetDescription desc;
            desc.add<edm::InputTag>("src")->setComment("input Subjet collection");
            desc.add<edm::InputTag>("subjetFlavourInfos")->setComment("input flavour info collection");
            desc.add<std::string>("name")->setComment("name of the HTTSubjet FlatTable we are extending with flavour information");
            desc.add<std::string>("cut")->setComment("cut on input HTTSubjet collection");
            desc.add<double>("deltaR")->setComment("deltaR to match Subjets");
            descriptions.add("HTTSubjetFlavourTable", desc);
        }

    private:
        void produce(edm::Event&, edm::EventSetup const&) override ;

        std::string name_;
        edm::EDGetTokenT<std::vector<pat::Jet> > src_;
        const StringCutObjectSelector<pat::Jet> cut_;
        const double deltaR_;
        edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> subjetFlavourInfosToken_;

};

// ------------ method called to produce the data  ------------
void
HTTSubjetFlavourTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(src_, jets);
    
    edm::Handle<reco::JetFlavourInfoMatchingCollection> subjetFlavourInfos;
    iEvent.getByToken(subjetFlavourInfosToken_, subjetFlavourInfos);

    unsigned int ncand = 0;
    std::vector<int> partonFlavour;
    std::vector<uint8_t> hadronFlavour;
     
    for (const pat::Jet & jet : *jets) {
      if (!cut_(jet)) continue;
      ++ncand;
      bool matched = false;
      for (const reco::JetFlavourInfoMatching & subjetFlavourInfoMatching : *subjetFlavourInfos) {
        if (deltaR(jet.p4(), subjetFlavourInfoMatching.first->p4()) < deltaR_) {
          partonFlavour.push_back(subjetFlavourInfoMatching.second.getPartonFlavour());
          hadronFlavour.push_back(subjetFlavourInfoMatching.second.getHadronFlavour());
          matched = true;
          break;
        }
      }
      if (!matched) {
        partonFlavour.push_back(0);
        hadronFlavour.push_back(0);
      }
    }

    auto tab  = std::make_unique<nanoaod::FlatTable>(ncand, name_, false, true);
    tab->addColumn<int>("partonFlavour", partonFlavour, "flavour from parton matching", nanoaod::FlatTable::IntColumn);
    tab->addColumn<int>("hadronFlavour", hadronFlavour, "flavour from hadron ghost clustering", nanoaod::FlatTable::IntColumn);

    iEvent.put(std::move(tab));
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(HTTSubjetFlavourTableProducer);
