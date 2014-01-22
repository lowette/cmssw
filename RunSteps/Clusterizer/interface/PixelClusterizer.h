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
    typedef edm::DetSet<PixelDigi>::const_iterator DigiIterator;

    PixelClusterizer(edm::ParameterSet const& conf);
    ~PixelClusterizer();


    bool setup(const PixelGeomDetUnit* pixDet);
    void clusterizeDetUnit(const edm::DetSet<PixelDigi> & input, const PixelGeomDetUnit* pixDet, std::vector<SiPixelCluster> & output);

private:
    edm::ParameterSet conf_;
    std::vector<SiPixelCluster> theClusters;
    SiPixelArrayBuffer theBuffer;
    int  nrows_;
    int  ncols_;
    uint32_t detid_;

    int32_t DigiSearchWidthX_;
    int32_t DigiSearchWidthY_;
    int32_t MinimumWidthX_;
    int32_t MaximumWidthX_;
    int32_t MinimumWidthY_;
    int32_t MaximumWidthY_;
    int32_t MinimumNumberOfDigis_;
    int32_t MaximumNumberOfDigis_;
    bool SplitIfOverThresholds_;

    void copy_to_buffer(DigiIterator begin, DigiIterator end);
    void clear_buffer(DigiIterator begin, DigiIterator end);
    SiPixelCluster make_cluster(const SiPixelCluster::PixelPos& pix);
};

#endif
