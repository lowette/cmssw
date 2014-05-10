#include "RunSteps/Clusterizer/interface/AdjacentHits.h"
#include "RunSteps/Clusterizer/interface/PixelClusterizer.h"
#include "RunSteps/Clusterizer/interface/SiPixelArrayBuffer.h"
#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelGainCalibrationOffline.h"

#include <stack>
#include <vector>
#include <iostream>

using namespace std;

AdjacentHits::AdjacentHits(edm::ParameterSet const& conf) :
    conf_(conf),
    nrows_(0),
    ncols_(0),
    detid_(0) {
    // Create a 2D matrix for this detector
    hitArray.setSize(nrows_, ncols_);
}

// Change the size of the 2D matrix for this detector unit
void AdjacentHits::setup(const PixelGeomDetUnit* pixDet) {
    const PixelTopology & topol = pixDet->specificTopology();
    nrows_ = topol.nrows();
    ncols_ = topol.ncolumns();
    if (nrows_ > hitArray.rows() || ncols_ > hitArray.columns()) hitArray.setSize(nrows_, ncols_);
}

// Go over the Digis and create clusters
void AdjacentHits::clusterizeDetUnit(const edm::DetSet<PixelDigi> & pixelDigis, const edm::Handle< edm::DetSetVector< PixelDigiSimLink > > & pixelSimLinks, edmNew::DetSetVector<SiPixelCluster>::FastFiller & clusters, std::vector<PixelClusterSimLink> & links) {

    // Get the det Id
    detid_ = pixelDigis.detId();

    // Fill the 2D matrix with the ADC values
    copy_to_buffer(pixelDigis.begin(), pixelDigis.end());

    // Loop over the Digis
    for (unsigned int row = 0; row < (unsigned int) nrows_; ++row) {
        for (unsigned int col = 0; col < (unsigned int)  ncols_; ++col) {

            // If the Digi is active
            if (hitArray(row, col)) {

                // Try to form a cluster

                // Create an entry for the pixelClusterLink
                PixelClusterSimLink link;

                // Add the simtrack of the Digi to the link
                link.addSimTrack(getSimTrackId(pixelSimLinks, PixelDigi::pixelToChannel(row, col)));

                // Set the value of this pixel to 0 as it is used to form a Digi
                hitArray.set(row, col, 0);

                // Create a temporary cluster (this allows to them easily form a "real" cluster with CMSSW data format)
                AccretionCluster acluster;

                // Create a pixel entry for the cluster
                SiPixelCluster::PixelPos firstpix(row, col);

                // Add the first pixel to the cluster
                acluster.add(firstpix, 255);

                // Go over all the pixels in the cluster
                while (!acluster.empty()) {

                    // Get the current pixel we are looking at
                    auto curInd = acluster.top();
                    acluster.pop();

                    // Look left and right
                    for (auto r = acluster.x[curInd] - 1; r <= acluster.x[curInd] + 1; ++r) {

                      // Look bottom and top
                        for (auto c = acluster.y[curInd] - 1; c <= acluster.y[curInd] + 1; ++c) {

                            // If the pixel is hit and has the same weight as the first pixel (probably from the same cluster)
                            if (hitArray(r, c)) {

                                // Add it to the cluster
                                SiPixelCluster::PixelPos newpix(r, c);
                                if (!acluster.add(newpix, 255)) break;

                                // And change its value
                                hitArray.set(newpix, 0);

                                // Add the simtrack of the Digi to the link
                                link.addSimTrack(getSimTrackId(pixelSimLinks, PixelDigi::pixelToChannel(row, col)));
                            }
                        }
                    }
                }

                // Form a "real" CMSSW cluster
                SiPixelCluster cluster(acluster.isize, acluster.adc, acluster.x, acluster.y, acluster.xmin, acluster.ymin);

                // Set the cluster of the link
                link.setCluster(cluster);

                // Add the cluster and the link to the list
                clusters.push_back(cluster);
                links.push_back(link);
            }
        }
    }

    // Reset the matrix
    clear_buffer(pixelDigis.begin(), pixelDigis.end());
}

void AdjacentHits::copy_to_buffer(DigiIterator begin, DigiIterator end) {
    // Copy the value of the Digis' ADC to the 2D matrix
    for (DigiIterator di = begin; di != end; ++di) hitArray.set(di->row(), di->column(), di->adc());
}

void AdjacentHits::clear_buffer(DigiIterator begin, DigiIterator end) {
    // Resets the matrix
    for (DigiIterator di = begin; di != end; ++di) hitArray.set(di->row(), di->column(), 0);
}

unsigned int AdjacentHits::getSimTrackId(const edm::Handle< edm::DetSetVector< PixelDigiSimLink > > & pixelSimLinks, int channel) {
    edm::DetSetVector<PixelDigiSimLink>::const_iterator isearch = pixelSimLinks->find(detid_);

    unsigned int simTrkId(0);
    if (isearch == pixelSimLinks->end()) return simTrkId;

    edm::DetSet<PixelDigiSimLink> link_detset = (*pixelSimLinks)[detid_];
    int iSimLink = 0;
    for (edm::DetSet<PixelDigiSimLink>::const_iterator it = link_detset.data.begin(); it != link_detset.data.end(); it++,iSimLink++) {
        if (channel == (int) it->channel()) {
            simTrkId = it->SimTrackId();
            break;
        }
    }
    return simTrkId;
}
