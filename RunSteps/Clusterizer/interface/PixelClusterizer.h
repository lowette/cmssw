#ifndef RunSteps_Clusterizer_PixelClusterizer_H
#define RunSteps_Clusterizer_PixelClusterizer_H

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

#include "RunSteps/Clusterizer/interface/SiPixelArrayBuffer.h"

#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationServiceBase.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>

class PixelGeomDetUnit;

class PixelClusterizer {

public:
    virtual void clusterizeDetUnit(const edm::DetSet<PixelDigi> & input, const PixelGeomDetUnit* pixDet, std::vector<SiPixelCluster> & output) = 0;
    virtual bool setup(const PixelGeomDetUnit* pixDet) = 0;

};

#endif
