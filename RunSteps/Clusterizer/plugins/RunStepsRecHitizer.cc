#include "RunSteps/Clusterizer/interface/RunStepsRecHitizer.h"
#include "RunSteps/Clusterizer/interface/PixelClusterizer.h"
#include "RunSteps/Clusterizer/interface/WeightedMeans2D.h"
#include "RunSteps/Clusterizer/interface/AdjacentHits.h"
#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"

#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
//#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonFwd.h"

//Needed for RecHits:
//#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
//#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
//#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include <vector>
#include <memory>
#include <string>
#include <iostream>

namespace cms {

    RunStepsRecHitizer::RunStepsRecHitizer(edm::ParameterSet const& conf) :
        conf_(conf),
        src_(conf.getParameter<edm::InputTag>("src")),
        algorithm_(conf.getParameter<std::string>("algorithm")),
        clusterSimLink_(conf.getParameter<bool>("clusterSimLink")) {

        // Names of the object that will be produced
        const std::string alias("siPixelClusters");

        // Objects that will be produced
        produces< SiPixelClusterCollectionNew >().setBranchAlias(alias);
        produces< edm::DetSetVector<PixelClusterSimLink> >().setBranchAlias(alias);

	//RecHit information: (from https://cmssdt.cern.ch/SDT/lxr/source/RecoLocalTracker/SiPixelRecHits/plugins/SiPixelRecHitConverter.cc)
        //tPixelCluster(consumes< edmNew::DetSetVector<SiPixelCluster> >( src_)) {
        //--- Declare to the EDM what kind of collections we will be making.
        produces<SiPixelRecHitCollection>();

        // Debug
        std::cout << "------------------------------------------------------------" << std::endl
                  << "-- Running RunSteps Clusterizer v0.5" << std::endl
		  << "       * 26/05/2014 RecHit info added (AO) " << std::endl
                  << "------------------------------------------------------------" << std::endl;

        // Set the algorithm to use
        if (algorithm_.compare("WeightedMeans2D") == 0) {
            std::cout << "Using the WeightedMeans2D algorithm" << std::endl;
            clusterizer_ = new WeightedMeans2D(conf);
        }
        else if (algorithm_.compare("AdjacentHits") == 0) {
            std::cout << "Using the AdjacentHits algorithm" << std::endl;
            clusterizer_ = new AdjacentHits(conf);
        }
        else {
            std::cout << "Using the default algorithm" << std::endl;
            clusterizer_ = new AdjacentHits(conf);
        }
    }

    RunStepsRecHitizer::~RunStepsRecHitizer() { }

    void RunStepsRecHitizer::beginJob() {
        edm::LogInfo("SiPixelClusterizer") << "[SiPixelClusterizer::beginJob]";
    }

    void RunStepsRecHitizer::produce(edm::Event & e, const edm::EventSetup & eventSetup) {
        // Get the Digis
        edm::Handle< edm::DetSetVector<PixelDigi> >  digis;
        e.getByLabel(src_, digis);

        // Get the simlinks for the Digis
        edm::Handle< edm::DetSetVector<PixelDigiSimLink> > pixelSimLinks;
        e.getByLabel(src_, pixelSimLinks);

        // Get the geometry
        edm::ESHandle<TrackerGeometry> geom;
        eventSetup.get<TrackerDigiGeometryRecord>().get(geom);

        // Global container for the clusters of each detector and the clusterLinks
        std::auto_ptr<SiPixelClusterCollectionNew> outputClusters( new SiPixelClusterCollectionNew() );
        std::vector<edm::DetSet<PixelClusterSimLink> > linksByDet;

	//Create empty output collection for the RecHits (from https://cmssdt.cern.ch/SDT/lxr/source/RecoLocalTracker/SiPixelRecHits/plugins/SiPixelRecHitConverter.cc)
	//std::auto_ptr<SiPixelRecHitCollectionNew> output(new SiPixelRecHitCollectionNew);
        
	// Step B*: create CPE
        //edm::ESHandle<PixelClusterParameterEstimator> hCPE;
        //std::string cpeName_ = conf_.getParameter<std::string>("CPE");
        //eventSetup.get<TkPixelCPERecord>().get(cpeName_,hCPE);
        //cpe_ = dynamic_cast< const PixelCPEBase* >(&(*hCPE));
	//cpe_ = &cpe;
    }
}
