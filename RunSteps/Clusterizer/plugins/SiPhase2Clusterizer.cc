#include "RunSteps/Clusterizer/interface/SiPhase2Clusterizer.h"
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

#include <vector>
#include <memory>
#include <string>
#include <iostream>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace cms {

	SiPhase2Clusterizer::SiPhase2Clusterizer(edm::ParameterSet const& conf) : 
    		conf_(conf),
    		clusterizer_(new PixelClusterizer(conf)),
    		src_(conf.getParameter<edm::InputTag>( "src" )) {
		const std::string alias ("siPixelClusters"); 
    		produces< edm::DetSetVector<SiPixelCluster> >().setBranchAlias(alias); 
	}

  	SiPhase2Clusterizer::~SiPhase2Clusterizer() { 
    		delete clusterizer_;
  	}  

  	void SiPhase2Clusterizer::beginJob() {
    		edm::LogInfo("SiPixelClusterizer") << "[SiPixelClusterizer::beginJob]";
  	}
  
	void SiPhase2Clusterizer::produce(edm::Event& e, const edm::EventSetup& es) {
    		edm::Handle< edm::DetSetVector<PixelDigi> >  input;
    		e.getByLabel(src_, input);
    		
		edm::ESHandle<TrackerGeometry> geom;
    		es.get<TrackerDigiGeometryRecord>().get(geom);

    		
		std::vector<edm::DetSet<SiPixelCluster> > clustersByDet;
		
		for(edm::DetSetVector<PixelDigi>::const_iterator DSViter = input->begin(); DSViter != input->end(); DSViter++) {
      			DetId detIdObject(DSViter->detId());
      			const GeomDetUnit* geoUnit = geom->idToDetUnit(detIdObject);
      			const PixelGeomDetUnit * pixDet = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);
      			if (!pixDet) assert(0);
      			
			edm::DetSet<SiPixelCluster> clusters(DSViter->detId());
      			clusterizer_->clusterizeDetUnit(*DSViter, pixDet, clusters.data);
      			
			clustersByDet.push_back(clusters);
      		}

		std::auto_ptr< edm::DetSetVector<SiPixelCluster> > output(new edm::DetSetVector<SiPixelCluster>(clustersByDet));
		e.put(output);
  	}
} 
