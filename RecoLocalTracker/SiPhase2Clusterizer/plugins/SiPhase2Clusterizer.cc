#include "RecoLocalTracker/SiPhase2Clusterizer/interface/SiPhase2Clusterizer.h"
#include "RecoLocalTracker/SiPhase2Clusterizer/interface/ClusterizerAlgorithm.h"
#include "RecoLocalTracker/SiPhase2Clusterizer/interface/PixelClusterSimLink.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"


#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <vector>
#include <memory>
#include <string>
#include <iostream>

namespace cms {

    SiPhase2Clusterizer::SiPhase2Clusterizer(edm::ParameterSet const& conf) :
        conf_(conf),
        src_(conf.getParameter< edm::InputTag >("src")),
        generateClusterSimLink_(conf.getParameter< bool >("clusterSimLink")) {

        // Names of the object that will be produced
        const std::string alias("siPixelClusters");

        // Objects that will be produced
        produces< SiPixelClusterCollectionNew >().setBranchAlias("siPixelClusters");

        // Optionally produce simlinks
        if (generateClusterSimLink_) produces< edm::DetSetVector<PixelClusterSimLink> >().setBranchAlias("siPixelClusters");

        // Debug
        // std::cout << "------------------------------------------------------------" << std::endl
                  // << "-- Running SiPhase2Clusterizer v1.0" << std::endl
                  // << "------------------------------------------------------------" << std::endl;

        // Set the algorithm to use
        clusterizer_ = new ClusterizerAlgorithm(conf);
    }

    SiPhase2Clusterizer::~SiPhase2Clusterizer() { }

    void SiPhase2Clusterizer::beginJob() {
        edm::LogInfo("SiPhase2Clusterizer") << "[SiPhase2Clusterizer::beginJob]";
    }

    void SiPhase2Clusterizer::produce(edm::Event & e, const edm::EventSetup & eventSetup) {

        // Get the Digis
        edm::Handle< edm::DetSetVector< PixelDigi > >  digis;
        e.getByLabel(src_, digis);

        // Get the simlinks for the Digis
        edm::Handle< edm::DetSetVector< PixelDigiSimLink > > pixelSimLinks;
        e.getByLabel(src_, pixelSimLinks);

        // Get the geometry
        edm::ESHandle< TrackerGeometry > geomHandle;
        eventSetup.get< TrackerDigiGeometryRecord >().get(geomHandle);
        const TrackerGeometry* tkGeom = &(*geomHandle);

        // Global container for the clusters of each modules
        std::auto_ptr< SiPixelClusterCollectionNew > outputClusters(new SiPixelClusterCollectionNew());

        // Go over all the modules
        for (edm::DetSetVector< PixelDigi >::const_iterator DSViter = digis->begin(); DSViter != digis->end(); ++DSViter) {

            // Geometry & detID
            DetId detId(DSViter->detId());
            const GeomDetUnit* geomDetUnit = tkGeom->idToDetUnit(detId);
            const PixelGeomDetUnit* pixDet = dynamic_cast< const PixelGeomDetUnit* >(geomDetUnit);
            if (!pixDet) assert(0);

            // Container for the clusters that will be produced for this modules
            edmNew::DetSetVector<SiPixelCluster>::FastFiller clusters(*outputClusters, DSViter->detId());

            // Setup the clusterizer algorithm for this detector (see ClusterizerAlgorithm for more details)
            clusterizer_->setup(pixDet);

            // Pass the list of Digis to the main algorithm
            // This function will store the clusters in the previously created container
            clusterizer_->clusterizeDetUnit(*DSViter, pixelSimLinks, clusters);

            if (clusters.empty()) clusters.abort();
        }

        // Add the data to the output
        edm::OrphanHandle< edmNew::DetSetVector<SiPixelCluster> > clusterCollection = e.put(outputClusters);

        // Do the same operations for the SimLinks if we have to generate them
        if (generateClusterSimLink_) {
            std::vector< edm::DetSet< PixelClusterSimLink > > linksByDet;
            clusterizer_->makeLinks(clusterCollection, linksByDet);
            std::auto_ptr< edm::DetSetVector< PixelClusterSimLink > > outputLinks(new edm::DetSetVector< PixelClusterSimLink >(linksByDet));
            e.put(outputLinks);
        }
    }
}

using cms::SiPhase2Clusterizer;
DEFINE_FWK_MODULE(SiPhase2Clusterizer);

