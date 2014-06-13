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
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"

#include "RunSteps/ValidaThor/interface/ValHits.h"

#include <TH2F.h>

class RunStepsGoValidaThor : public edm::EDAnalyzer {

public:

    explicit RunStepsGoValidaThor(const edm::ParameterSet&);
    ~RunStepsGoValidaThor();
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

private:
    edm::InputTag rechits_, clusters_;
    bool useRecHits_;

    TH2D* xyhisto;
};

RunStepsGoValidaThor::RunStepsGoValidaThor(const edm::ParameterSet& iConfig) {
    rechits_ = iConfig.getParameter< edm::InputTag >("rechits");
    clusters_ = iConfig.getParameter< edm::InputTag >("clusters");
    useRecHits_ = iConfig.getParameter< bool >("useRecHits");

    std::cout << "------------------------------------------------------------" << std::endl
              << "-- Running RunSteps RunStepsGoValidaThor v0.0" << std::endl
              << "------------------------------------------------------------" << std::endl;

    // Use RecHits
    if (useRecHits_) std::cout << "INFO: Using RecHits" << std::endl;
    // Use Clusters
    else std::cout << "INFO: Using Clusters" << std::endl;


    // Make test histo
    edm::Service<TFileService> fs;
    xyhisto = fs->make<TH2D>("2D", "2D", 1000, 0., 0., 1000, 0., 0.);
}

RunStepsGoValidaThor::~RunStepsGoValidaThor() { }

void RunStepsGoValidaThor::beginJob() { }

void RunStepsGoValidaThor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    // Get the clusters
    edm::Handle< SiPixelClusterCollectionNew > clustersHandle;
    iEvent.getByLabel(clusters_, clustersHandle);
    // const edmNew::DetSetVector< SiPixelCluster >& clusters = *clustersHandle;

    // Get the cluster simlinks
    edm::Handle< edm::DetSetVector< PixelClusterSimLink > > clusterLinksHandle;
    iEvent.getByLabel(clusters_, clusterLinksHandle);
    const edm::DetSetVector< PixelClusterSimLink >& clusterLinks = *clusterLinksHandle;

    // Get the Geometry
    edm::ESHandle< TrackerGeometry > geomHandle;
    iSetup.get< TrackerDigiGeometryRecord >().get(geomHandle);
    const TrackerGeometry& tkGeom = *geomHandle;


    // Make selection on RecHits or Clusters
    ValHitsCollection hitsCollection;

    //Get the RecHits
    if (useRecHits_) {
        edm::Handle< SiPixelRecHitCollection > recHitsHandle;
        iEvent.getByLabel(rechits_, recHitsHandle);
        const edmNew::DetSetVector< SiPixelRecHit >& recHits = *recHitsHandle;

        hitsCollection = ValHitsBuilder((TrackerGeometry*) & tkGeom, (edm::DetSetVector< PixelClusterSimLink >*) & clusterLinks, (edmNew::DetSetVector< SiPixelRecHit >*) & recHits);
    }
    // Use Clusters
    else hitsCollection = ValHitsBuilder((TrackerGeometry*) & tkGeom, (edm::DetSetVector< PixelClusterSimLink >*) & clusterLinks);

    //////////////////////////////////////////////////////////////////////
    // Give hitsCollection to Validatator to play with the hits         //
    //////////////////////////////////////////////////////////////////////

    // Loop over the Hits
    for (ValHitsCollection::const_iterator vhCollectionIter = hitsCollection.begin(); vhCollectionIter != hitsCollection.end(); ++vhCollectionIter) {

        // unsigned int detId = vhCollectionIter->first;
        ValHitsVector hitsVector = vhCollectionIter->second;


        for (ValHitsVector::const_iterator vhVectorIter = hitsVector.begin(); vhVectorIter != hitsVector.end(); ++vhVectorIter) {

            ValHit hit = *vhVectorIter;

            /////////////////////////////////////////
            // Do things here to the hits          //
            /////////////////////////////////////////
            xyhisto->Fill(hit.globalPos.x(), hit.globalPos.y());
        }
    }

}

void RunStepsGoValidaThor::endJob() { }

DEFINE_FWK_MODULE(RunStepsGoValidaThor);
