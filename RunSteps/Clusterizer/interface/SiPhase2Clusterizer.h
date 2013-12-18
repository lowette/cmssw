#ifndef RunSteps_Clusterizer_SiPhase2Clusterizer_h
#define RunSteps_Clusterizer_SiPhase2Clusterizer_h

#include "RunSteps/Clusterizer/interface/PixelClusterizer.h"
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

	class SiPhase2Clusterizer : public edm::EDProducer {
  	
	public:
		explicit SiPhase2Clusterizer(const edm::ParameterSet& conf);
    		virtual ~SiPhase2Clusterizer();
    		void setupClusterizer();
    		virtual void beginJob( );
    		virtual void produce(edm::Event& e, const edm::EventSetup& c);

  	private:
    		edm::ParameterSet conf_;
    		PixelClusterizer* clusterizer_;
    		edm::InputTag src_;
  	};
}

#endif
