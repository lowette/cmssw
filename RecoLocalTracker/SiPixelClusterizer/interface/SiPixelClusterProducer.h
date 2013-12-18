#ifndef RecoLocalTracker_SiPixelClusterizer_SiPixelClusterProducer_h
#define RecoLocalTracker_SiPixelClusterizer_SiPixelClusterProducer_h

#include "RecoLocalTracker/SiPixelClusterizer/interface/PixelThresholdClusterizer.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace cms {

	class SiPixelClusterProducer : public edm::EDProducer {
  	
	public:
		explicit SiPixelClusterProducer(const edm::ParameterSet& conf);
    		virtual ~SiPixelClusterProducer();
    		void setupClusterizer();
    		virtual void beginJob( );
    		virtual void produce(edm::Event& e, const edm::EventSetup& c);

  	private:
    		edm::ParameterSet conf_;
    		PixelThresholdClusterizer* clusterizer_;
    		edm::InputTag src_;
  	};
}

#endif
