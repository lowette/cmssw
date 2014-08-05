#ifndef RunSteps_Clusterizer_RunStepsClusterizer_h
#define RunSteps_Clusterizer_RunStepsClusterizer_h

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

#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include <string>

namespace cms {

	class RunStepsClusterizer : public edm::EDProducer {

	public:
		explicit RunStepsClusterizer(const edm::ParameterSet& conf);
    		virtual ~RunStepsClusterizer();
    		void setupClusterizer();
    		virtual void beginJob( );
    		virtual void produce(edm::Event& e, const edm::EventSetup& eventSetup);

  	private:
    		edm::ParameterSet conf_;
    		PixelClusterizer* clusterizer_;
            edm::InputTag src_;
    		std::string algorithm_;
            bool clusterSimLink_;
  	};
}

#endif
