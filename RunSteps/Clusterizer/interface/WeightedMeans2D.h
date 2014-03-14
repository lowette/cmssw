#ifndef RunSteps_Clusterizer_2DWEIGHTEDMEANS_H
#define RunSteps_Clusterizer_2DWEIGHTEDMEANS_H

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

#include "RunSteps/Clusterizer/interface/PixelClusterizer.h"
#include "RunSteps/Clusterizer/interface/SiPixelArrayBuffer.h"

#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationServiceBase.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>

class PixelClusterizer;

class PixelGeomDetUnit;

class WeightedMeans2D : public PixelClusterizer {

public:
    typedef edm::DetSet<PixelDigi>::const_iterator DigiIterator;

    WeightedMeans2D(edm::ParameterSet const& conf);

    bool setup(const PixelGeomDetUnit* pixDet);
    void clusterizeDetUnit(const edm::DetSet<PixelDigi> & input, const PixelGeomDetUnit* pixDet, std::vector<SiPixelCluster> & output);

private:
    edm::ParameterSet conf_;
    std::vector<SiPixelCluster> theClusters;
    SiPixelArrayBuffer hitArray;
    SiPixelArrayBuffer weightArray;
    SiPixelArrayBuffer maskedArray;
    int nrows_;
    int ncols_;
    uint32_t detid_;

    SiPixelCluster make_cluster(int row, int col);
    void copy_to_buffer(DigiIterator begin, DigiIterator end);
    void clear_buffer(DigiIterator begin, DigiIterator end);
    unsigned int getLayerNumber(unsigned int & detid);

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
