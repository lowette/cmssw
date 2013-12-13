#ifndef RecoLocalTracker_SiPixelClusterizer_PixelThresholdClusterizer_H
#define RecoLocalTracker_SiPixelClusterizer_PixelThresholdClusterizer_H

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "RecoLocalTracker/SiPixelClusterizer/interface/SiPixelArrayBuffer.h"

#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationServiceBase.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>

class PixelGeomDetUnit;

class PixelThresholdClusterizer {

public:
	typedef edm::DetSet<PixelDigi>::const_iterator    DigiIterator;

	PixelThresholdClusterizer(edm::ParameterSet const& conf);
 	~PixelThresholdClusterizer();

  	void clusterizeDetUnit(const edm::DetSet<PixelDigi> & input, const PixelGeomDetUnit * pixDet, edmNew::DetSetVector<SiPixelCluster>::FastFiller& output);

private:
  	edm::ParameterSet conf_;
  	std::vector<SiPixelCluster> theClusters; 
  	SiPixelArrayBuffer theBuffer;
  	int  nrows_;
  	int  ncols_;
  	uint32_t detid_;
  	
	bool setup(const PixelGeomDetUnit * pixDet);
  	void copy_to_buffer(DigiIterator begin, DigiIterator end);
  	void clear_buffer(DigiIterator begin, DigiIterator end);
  	SiPixelCluster make_cluster(const SiPixelCluster::PixelPos& pix, edmNew::DetSetVector<SiPixelCluster>::FastFiller& output);
};

#endif
