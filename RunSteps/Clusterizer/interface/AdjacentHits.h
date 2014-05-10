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
    void clusterizeDetUnit(const edm::DetSet<PixelDigi> & pixelDigis, const edm::Handle< edm::DetSetVector< PixelDigiSimLink > > & pixelSimLinks, edmNew::DetSetVector<SiPixelCluster>::FastFiller & clusters, std::vector<PixelClusterSimLink> & links);

private:
    edm::ParameterSet conf_;
    SiPixelArrayBuffer hitArray;
    int nrows_;
    int ncols_;
    uint32_t detid_;

    void copy_to_buffer(DigiIterator begin, DigiIterator end);
    void clear_buffer(DigiIterator begin, DigiIterator end);
    unsigned int getSimTrackId(const edm::Handle< edm::DetSetVector< PixelDigiSimLink > > & pixelSimLinks, int channel);

    //

    struct AccretionCluster {
        static constexpr unsigned short MAXSIZE = 256;
        unsigned short adc[256];
        unsigned short x[256];
        unsigned short y[256];
        unsigned short xmin = 16000;
        unsigned short xmax = 0;
        unsigned short ymin = 16000;
        unsigned short ymax = 0;
        unsigned int isize = 0;
        unsigned int curr = 0;
        unsigned short top() const { return curr; }
        void pop() { ++curr; }
        bool empty() { return curr == isize; }
        bool add(SiPixelCluster::PixelPos const & p, unsigned short const iadc) {
            if (isize == MAXSIZE) return false;
            xmin = std::min(xmin, (unsigned short) p.row());
            xmax = std::max(xmax, (unsigned short) p.row());
            ymin = std::min(ymin, (unsigned short) p.col());
            ymax = std::max(ymax, (unsigned short) p.col());
            adc[isize] = iadc;
            x[isize] = p.row();
            y[isize++] = p.col();
            return true;
        }
        unsigned short size() { return isize; }
        unsigned short xsize() { return xmax - xmin + 1; }
        unsigned short ysize() { return ymax - ymin + 1; }
    };
};

#endif
