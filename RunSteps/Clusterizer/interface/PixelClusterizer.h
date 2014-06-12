#ifndef RunSteps_Clusterizer_PixelClusterizer_H
#define RunSteps_Clusterizer_PixelClusterizer_H

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"

#include "RunSteps/Clusterizer/interface/SiPixelArrayBuffer.h"
#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"

#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationServiceBase.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>
#include <map>

class PixelGeomDetUnit;

class PixelClusterizer {

public:
    virtual void setup(const PixelGeomDetUnit* pixDet) = 0;
    virtual void clusterizeDetUnit(const edm::DetSet<PixelDigi> & pixelDigis, const edm::Handle< edm::DetSetVector< PixelDigiSimLink > > & pixelSimLinks, edmNew::DetSetVector<SiPixelCluster>::FastFiller & clusters) = 0;

    void makeLinks(edm::OrphanHandle< edmNew::DetSetVector<SiPixelCluster> > & clusters, std::vector<edm::DetSet<PixelClusterSimLink> > & linksByDet);

    unsigned int getSimTrackId(const edm::Handle< edm::DetSetVector< PixelDigiSimLink > > & pixelSimLinks, int channel);

public:
    edm::ParameterSet conf_;
    SiPixelArrayBuffer hitArray;
    int nrows_;
    int ncols_;
    uint32_t detid_;

    std::map< SiPixelCluster, std::vector< unsigned int > > tmpSimLinks;

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
