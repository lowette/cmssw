#include "RunSteps/Clusterizer/interface/RunStepsClusterizer.h"
#include "RunSteps/Clusterizer/interface/PixelClusterizer.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationService.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineService.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationForHLTService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include <memory>
#include <string>
#include <iostream>

namespace cms {

    RunStepsClusterizer::RunStepsClusterizer(edm::ParameterSet const& conf) :
        conf_(conf),
        clusterizer_(new PixelClusterizer(conf)),
        src_(conf.getParameter<edm::InputTag>("src")) {
        const std::string alias("siPixelClusters");
        produces< edm::DetSetVector<SiPixelCluster> >().setBranchAlias(alias);

        std::cout << "------------------------------------------------------------" << std::endl
                  << "-- Running RunSteps Clusterizer v0.2" << std::endl
                  << "------------------------------------------------------------" << std::endl;
    }

    RunStepsClusterizer::~RunStepsClusterizer() {
        delete clusterizer_;
    }

    void RunStepsClusterizer::beginJob() {
        edm::LogInfo("SiPixelClusterizer") << "[SiPixelClusterizer::beginJob]";
    }

    void RunStepsClusterizer::produce(edm::Event& e, const edm::EventSetup& eventSetup) {
        edm::Handle< edm::DetSetVector<PixelDigi> >  input;
        e.getByLabel(src_, input);

        edm::ESHandle<TrackerGeometry> geom;
        eventSetup.get<TrackerDigiGeometryRecord>().get(geom);

        // Global container for the clusters of each detector
        std::vector<edm::DetSet<SiPixelCluster> > clustersByDet;

        // Go over all the detectors
        for(edm::DetSetVector<PixelDigi>::const_iterator DSViter = input->begin(); DSViter != input->end(); ++DSViter) {
            DetId detIdObject(DSViter->detId());
            const GeomDetUnit* geoUnit = geom->idToDetUnit(detIdObject);
            const PixelGeomDetUnit * pixDet = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);
            if (!pixDet) assert(0);

            // Container for the clusters that will be produced for this detector
            edm::DetSet<SiPixelCluster> clusters(DSViter->detId());

            // Setup the clusterizer algorithm for this detector (see PixelClusterizer for more details)
            clusterizer_->setup(pixDet);

            // Pass the list of Digis to the main algorithm
            // This function will store the clusters in the previously created container
            clusterizer_->clusterizeDetUnit(*DSViter, pixDet, clusters.data);

            // Add the clusters for this detector to the global container
            clustersByDet.push_back(clusters);
        }

        std::auto_ptr< edm::DetSetVector<SiPixelCluster> > output(new edm::DetSetVector<SiPixelCluster>(clustersByDet));
        e.put(output);
    }
}
