#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "RunSteps/ValidaThor/interface/ValHits.h"

class RunStepsRecHitsValidaThor : public edm::EDAnalyzer {

public:

    explicit RunStepsRecHitsValidaThor(const edm::ParameterSet&);
    ~RunStepsRecHitsValidaThor();
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

private:
    edm::InputTag src_;

};

RunStepsRecHitsValidaThor::RunStepsRecHitsValidaThor(const edm::ParameterSet& iConfig) {
    src_ = iConfig.getParameter<edm::InputTag>("src");

    std::cout << "------------------------------------------------------------" << std::endl
              << "-- Running RunSteps RunStepsRecHitsValidaThor v0.0" << std::endl
              << "------------------------------------------------------------" << std::endl;
}

RunStepsRecHitsValidaThor::~RunStepsRecHitsValidaThor() { }

void RunStepsRecHitsValidaThor::beginJob() { }

void RunStepsRecHitsValidaThor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    // Get the clusters
    edm::Handle< SiPixelClusterCollectionNew > clustersHandle;
    iEvent.getByLabel(src_, clustersHandle);

    // const edmNew::DetSetVector< SiPixelCluster >& clusters = *clustersHandle;

    // Get the links
    edm::Handle< edm::DetSetVector< PixelClusterSimLink > > clusterLinksHandle;
    iEvent.getByLabel(src_, clusterLinksHandle);

    const edm::DetSetVector< PixelClusterSimLink >& clusterLinks = *clusterLinksHandle;

    ValHitsCollection hitsCollection = ValHitsBuilder((edm::DetSetVector< PixelClusterSimLink >*) & clusterLinks);

    for (ValHitsCollection::const_iterator vhCollectionIter = hitsCollection.begin(); vhCollectionIter != hitsCollection.end(); ++vhCollectionIter) {

        // unsigned int detId = vhCollectionIter->first;
        ValHitsVector hitsVector = vhCollectionIter->second;


        for (ValHitsVector::const_iterator vhVectorIter = hitsVector.begin(); vhVectorIter != hitsVector.end(); ++vhVectorIter) {

            ValHit hit = *vhVectorIter;

            std::cout << hit.x << std::endl;

        }

    }

}

void RunStepsRecHitsValidaThor::endJob() { }

DEFINE_FWK_MODULE(RunStepsRecHitsValidaThor);
