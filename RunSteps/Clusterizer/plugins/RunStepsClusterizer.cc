#include "RunSteps/Clusterizer/interface/RunStepsClusterizer.h"
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

#include <vector>
#include <memory>
#include <string>
#include <iostream>

namespace cms {

    RunStepsClusterizer::RunStepsClusterizer(edm::ParameterSet const& conf) :
        conf_(conf),
        src_(conf.getParameter<edm::InputTag>("src")),
        algorithm_(conf.getParameter<std::string>("algorithm")),
        clusterSimLink_(conf.getParameter<bool>("clusterSimLink")) {

        // Names of the object that will be produced
        const std::string alias("siPixelClusters");

        // Objects that will be produced
        produces< SiPixelClusterCollectionNew >().setBranchAlias(alias);
        produces< edm::DetSetVector<PixelClusterSimLink> >().setBranchAlias(alias);

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

    RunStepsClusterizer::~RunStepsClusterizer() { }

    void RunStepsClusterizer::beginJob() {
        edm::LogInfo("SiPixelClusterizer") << "[SiPixelClusterizer::beginJob]";
    }

    void RunStepsClusterizer::produce(edm::Event & e, const edm::EventSetup & eventSetup) {
        // Get the Digis
        edm::Handle< edm::DetSetVector<PixelDigi> >  digis;
        e.getByLabel(src_, digis);

        // Get the simlinks for the Digis
        edm::Handle< edm::DetSetVector<PixelDigiSimLink> > pixelSimLinks;
        e.getByLabel(src_, pixelSimLinks);

        // Get the geometry
        edm::ESHandle<TrackerGeometry> geom;
        eventSetup.get<TrackerDigiGeometryRecord>().get(geom);

        // Global container for the clusters of each detector
        std::auto_ptr<SiPixelClusterCollectionNew> outputClusters( new SiPixelClusterCollectionNew() );

	   // Go over all the detectors
        for (edm::DetSetVector<PixelDigi>::const_iterator DSViter = digis->begin(); DSViter != digis->end(); ++DSViter) {
            DetId detIdObject(DSViter->detId());
            const GeomDetUnit* geoUnit = geom->idToDetUnit(detIdObject);
            const PixelGeomDetUnit * pixDet = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);
            if (!pixDet) assert(0);

            // Container for the clusters and clusterLinks that will be produced for this detector
            edmNew::DetSetVector<SiPixelCluster>::FastFiller clusters(*outputClusters, detIdObject);

            // Setup the clusterizer algorithm for this detector (see PixelClusterizer for more details)
            clusterizer_->setup(pixDet);

            // Pass the list of Digis to the main algorithm
            // This function will store the clusters in the previously created container
            clusterizer_->clusterizeDetUnit(*DSViter, pixelSimLinks, clusters);

            if (clusters.empty()) clusters.abort();
        }

        // Add the data to the output
        edm::OrphanHandle< edmNew::DetSetVector< SiPixelCluster > > clusterCollection = e.put(outputClusters);

        // Same for links

        if (clusterSimLink_) {
            std::vector< edm::DetSet< PixelClusterSimLink > > linksByDet;
            clusterizer_->makeLinks(clusterCollection, linksByDet);
            std::auto_ptr< edm::DetSetVector< PixelClusterSimLink > > outputLinks(new edm::DetSetVector< PixelClusterSimLink >(linksByDet));
            e.put(outputLinks);
        }
    }
}
