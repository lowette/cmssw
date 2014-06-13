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
#include "RunSteps/ValidaThor/interface/ValidaThor.h"

#include <TH2F.h>

class RunStepsGoValidaThor : public edm::EDAnalyzer {

public:

    explicit RunStepsGoValidaThor(const edm::ParameterSet&);
    ~RunStepsGoValidaThor();
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

private:
    bool useRecHits_;
    ValidaThor* validaThor;

};

RunStepsGoValidaThor::RunStepsGoValidaThor(const edm::ParameterSet& iConfig) {
    useRecHits_ = iConfig.getParameter< bool >("useRecHits");

    std::cout << "------------------------------------------------------------" << std::endl
              << "-- Running RunSteps RunStepsGoValidaThor v0.0" << std::endl
              << "------------------------------------------------------------" << std::endl;

    // Use RecHits
    if (useRecHits_) std::cout << "INFO: Using RecHits" << std::endl;
    // Use Clusters
    else std::cout << "INFO: Using Clusters" << std::endl;

    validaThor = new ValidaThor();
}

RunStepsGoValidaThor::~RunStepsGoValidaThor() { }

void RunStepsGoValidaThor::beginJob() { }

void RunStepsGoValidaThor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    // Get the clusters
    edm::Handle< SiPixelClusterCollectionNew > clustersHandle;
    iEvent.getByLabel("siPixelClusters", clustersHandle);
    // const edmNew::DetSetVector< SiPixelCluster >& clusters = *clustersHandle;

    // Get the cluster simlinks
    edm::Handle< edm::DetSetVector< PixelClusterSimLink > > clusterLinksHandle;
    iEvent.getByLabel("siPixelClusters", clusterLinksHandle);
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
        iEvent.getByLabel("siPixelRecHits", recHitsHandle);
        const edmNew::DetSetVector< SiPixelRecHit >& recHits = *recHitsHandle;

        hitsCollection = ValHitsBuilder((TrackerGeometry*) & tkGeom, (edm::DetSetVector< PixelClusterSimLink >*) & clusterLinks, (edmNew::DetSetVector< SiPixelRecHit >*) & recHits);
    }
    // Use Clusters
    else hitsCollection = ValHitsBuilder((TrackerGeometry*) & tkGeom, (edm::DetSetVector< PixelClusterSimLink >*) & clusterLinks);

    // SimHit
    edm::Handle< edm::PSimHitContainer > simHits_BHandle;
    iEvent.getByLabel("g4SimHits", "TrackerHitsPixelBarrelLowTof", simHits_BHandle);
    const edm::PSimHitContainer& simHits_B = *simHits_BHandle;

    edm::Handle< edm::PSimHitContainer > simHits_EHandle;
    iEvent.getByLabel("g4SimHits", "TrackerHitsPixelEndcapLowTof", simHits_EHandle);
    const edm::PSimHitContainer& simHits_E = *simHits_EHandle;

    // SimTrack
    edm::Handle< edm::SimTrackContainer > simTracksHandle;
    iEvent.getByLabel("g4SimHits", simTracksHandle);
    const edm::SimTrackContainer& simTracks = *simTracksHandle;

    // SimVertex
    edm::Handle< edm::SimVertexContainer > simVerticesHandle;
    iEvent.getByLabel("g4SimHits", simVerticesHandle);
    const edm::SimVertexContainer& simVertices = *simVerticesHandle;

    // Validation module
    validaThor->analyze((ValHitsCollection*) & hitsCollection, (edm::PSimHitContainer*) & simHits_B, (edm::PSimHitContainer*) & simHits_E, (edm::SimTrackContainer*) & simTracks, (edm::SimVertexContainer*) & simVertices, (TrackerGeometry*) & tkGeom);
}

void RunStepsGoValidaThor::endJob() { }

DEFINE_FWK_MODULE(RunStepsGoValidaThor);
