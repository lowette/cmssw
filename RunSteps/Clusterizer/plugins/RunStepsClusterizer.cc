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

#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationService.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineService.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationForHLTService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//From HitEff.cc file
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "AnalysisDataFormats/SiStripClusterInfo/interface/SiStripClusterInfo.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "CalibTracker/Records/interface/SiStripQualityRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

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

	//get pixel Cluster Parameter Estimator    //Based on : https://github.com/cms-analysis/TrackingAnalysis-Cosmics/blob/master/src/HitEff.cc
	edm::ESHandle<PixelClusterParameterEstimator> pixelCPE;
	//es.get<TkPixelCPERecord>().get("PixelCPEfromTrackAngle", pixelCPE); 
	es.get<TkPixelCPERecord>().get("PixelCPETemplateReco", pixelCPE); 
	const PixelClusterParameterEstimator &pixelcpe(*pixelCPE);


        // Global container for the clusters of each detector and the clusterLinks
        std::auto_ptr<SiPixelClusterCollectionNew> outputClusters( new SiPixelClusterCollectionNew() );
        std::vector<edm::DetSet<PixelClusterSimLink> > linksByDet;

        // Go over all the detectors
        for (edm::DetSetVector<PixelDigi>::const_iterator DSViter = digis->begin(); DSViter != digis->end(); ++DSViter) {
            DetId detIdObject(DSViter->detId());
            const GeomDetUnit* geoUnit = geom->idToDetUnit(detIdObject);
            const PixelGeomDetUnit * pixDet = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);
            if (!pixDet) assert(0);

            // Container for the clusters and clusterLinks that will be produced for this detector
            edmNew::DetSetVector<SiPixelCluster>::FastFiller clusters(*outputClusters, detIdObject);

            edm::DetSet<PixelClusterSimLink> links(DSViter->detId());

            // Setup the clusterizer algorithm for this detector (see PixelClusterizer for more details)
            clusterizer_->setup(pixDet);

            // Pass the list of Digis to the main algorithm
            // This function will store the clusters in the previously created container
            clusterizer_->clusterizeDetUnit(*DSViter, pixelSimLinks, clusters, links.data);

            if (clusters.empty()) clusters.abort();

            // Add the clusters for this detector to the global container
            linksByDet.push_back(links);
        }

        // Add the data to the output
        e.put(outputClusters);

        // Same for links

        if (clusterSimLink_) {
            std::auto_ptr< edm::DetSetVector<PixelClusterSimLink> > outputLinks(new edm::DetSetVector<PixelClusterSimLink>(linksByDet));
            e.put(outputLinks);
        }
    }
}
