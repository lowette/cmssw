#include "RecoLocalTracker/SiPixelClusterizer/interface/SiPixelClusterProducer.h"
#include "RecoLocalTracker/SiPixelClusterizer/interface/PixelThresholdClusterizer.h"

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

	SiPixelClusterProducer::SiPixelClusterProducer(edm::ParameterSet const& conf) : 
    		conf_(conf),
    		clusterizer_(new PixelThresholdClusterizer(conf)),
    		src_(conf.getParameter<edm::InputTag>( "src" )) {
    		produces<SiPixelClusterCollectionNew>(); 
	}

  	SiPixelClusterProducer::~SiPixelClusterProducer() { 
    		delete clusterizer_;
  	}  

  	void SiPixelClusterProducer::beginJob() {
    		edm::LogInfo("SiPixelClusterizer") << "[SiPixelClusterizer::beginJob]";
  	}
  
	void SiPixelClusterProducer::produce(edm::Event& e, const edm::EventSetup& es) {
    		edm::Handle< edm::DetSetVector<PixelDigi> >  input;
    		e.getByLabel(src_, input);
    		
		edm::ESHandle<TrackerGeometry> geom;
    		es.get<TrackerDigiGeometryRecord>().get(geom);

    		std::auto_ptr<SiPixelClusterCollectionNew> output(new SiPixelClusterCollectionNew());

    		run(*input, geom, *output);

    		e.put(output);
 	}
  	
	void SiPixelClusterProducer::run(const edm::DetSetVector<PixelDigi>  & input, edm::ESHandle<TrackerGeometry> & geom, edmNew::DetSetVector<SiPixelCluster> & output) {
    		for(edm::DetSetVector<PixelDigi>::const_iterator DSViter = input.begin(); DSViter != input.end(); DSViter++) {
      			DetId detIdObject(DSViter->detId());
      			const GeomDetUnit* geoUnit = geom->idToDetUnit(detIdObject);
      			const PixelGeomDetUnit * pixDet = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);
      			if (!pixDet) assert(0);
      			edmNew::DetSetVector<SiPixelCluster>::FastFiller spc(output, DSViter->detId());
      			clusterizer_->clusterizeDetUnit(*DSViter, pixDet, spc);
      			if (spc.empty()) spc.abort();
      		}
  	}
} 
