#ifndef RunSteps_Clusterizer_ADJACENTHITS_H
#define RunSteps_Clusterizer_ADJACENTHITS_H

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"

#include "RunSteps/Clusterizer/interface/PixelClusterizer.h"
#include "RunSteps/Clusterizer/interface/SiPixelArrayBuffer.h"
#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"

#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationServiceBase.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>

class PixelClusterizer;

class PixelGeomDetUnit;

class PixelClusterSimLink;

class AdjacentHits : public PixelClusterizer {

public:
    typedef edm::DetSet<PixelDigi>::const_iterator DigiIterator;

    AdjacentHits(edm::ParameterSet const& conf);

    void setup(const PixelGeomDetUnit* pixDet);
    void clusterizeDetUnit(const edm::DetSet<PixelDigi> & pixelDigis, const edm::Handle< edm::DetSetVector< PixelDigiSimLink > > & pixelSimLinks, edmNew::DetSetVector<SiPixelCluster>::FastFiller & clusters);

private:
    void copy_to_buffer(DigiIterator begin, DigiIterator end);
    void clear_buffer(DigiIterator begin, DigiIterator end);

};

#endif
