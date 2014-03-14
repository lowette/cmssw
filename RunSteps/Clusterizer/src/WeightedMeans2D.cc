#include "RunSteps/Clusterizer/interface/WeightedMeans2D.h"
#include "RunSteps/Clusterizer/interface/PixelClusterizer.h"
#include "RunSteps/Clusterizer/interface/SiPixelArrayBuffer.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelGainCalibrationOffline.h"

#include <stack>
#include <vector>
#include <iostream>

using namespace std;

WeightedMeans2D::WeightedMeans2D(edm::ParameterSet const& conf) :
    conf_(conf),
    nrows_(0),
    ncols_(0),
    detid_(0) {
    // Create a 2D matrix for this detector
    hitArray.setSize(nrows_, ncols_);
    // Create a 2D matrix for the weight of the pixels
    weightArray.setSize(nrows_, ncols_);
}

// Change the size of the 2D matrix for this detector unit
bool WeightedMeans2D::setup(const PixelGeomDetUnit* pixDet) {
    const PixelTopology & topol = pixDet->specificTopology();
    nrows_ = topol.nrows();
    ncols_ = topol.ncolumns();
    if (nrows_ > hitArray.rows() || ncols_ > hitArray.columns()) {
        hitArray.setSize(nrows_, ncols_);
        weightArray.setSize(nrows_, ncols_);
        maskedArray.setSize(nrows_, ncols_);
    }
    return true;
}

// Go over the Digis and create clusters
void WeightedMeans2D::clusterizeDetUnit(const edm::DetSet<PixelDigi> & input, const PixelGeomDetUnit* pixDet, std::vector<SiPixelCluster> & output) {
	if (!setup(pixDet)) return;

	detid_ = input.detId();
    // unsigned int layer = getLayerNumber(detid_);

	// Fill the 2D matrix with the ADC values
  	copy_to_buffer(input.begin(), input.end());

    // Compute the weight
    for (unsigned int row = 0; row < (unsigned int) nrows_; ++row) {
        for (unsigned int col = 0; col < (unsigned int)  ncols_; ++col) {
            int weight = 6 * hitArray(row, col)
                       + 3 * (hitArray(row + 1, col) + hitArray(row - 1, col) + hitArray(row, col + 1) + hitArray(row, col + 1))
                       + 1 * (hitArray(row + 1, col + 1) + hitArray(row - 1, col + 1) + hitArray(row - 1, col - 1) + hitArray(row + 1, col - 1));
            weightArray.set(row, col, weight);
            maskedArray.set(row, col, weight);
        }
    }

    // Apply a second mask
    for (unsigned int row = 0; row < (unsigned int) nrows_; ++row) {
        for (unsigned int col = 0; col < (unsigned int)  ncols_; ++col) {
            int weight = weightArray(row, col);

            if (weightArray(row - 1, col) < weight) maskedArray.set(row - 1, col, 0);
            if (weightArray(row + 1, col) < weight) maskedArray.set(row + 1, col, 0);
            if (weightArray(row, col - 1) < weight) maskedArray.set(row, col - 1, 0);
            if (weightArray(row, col + 1) < weight) maskedArray.set(row, col + 1, 0);
            if (weightArray(row - 1, col - 1) < weight) maskedArray.set(row - 1, col - 1, 0);
            if (weightArray(row - 1, col + 1) < weight) maskedArray.set(row - 1, col + 1, 0);
            if (weightArray(row + 1, col - 1) < weight) maskedArray.set(row + 1, col - 1, 0);
            if (weightArray(row + 1, col + 1) < weight) maskedArray.set(row + 1, col + 1, 0);
        }
    }

    // Loop over the Digis
    // cout << "--------------------" << endl;
    for (unsigned int row = 0; row < (unsigned int) nrows_; ++row) {
        for (unsigned int col = 0; col < (unsigned int)  ncols_; ++col) {
            // If the Digi is active
            if (maskedArray(row, col)) {
                // Try to form a cluster
                SiPixelCluster cluster = make_cluster(row, col);
                // Add the cluster to the list
                output.push_back(cluster);
                // cout << "- Pixel(" << layer << ", " << row << ", " << col << ")\t" << weightArray(row, col) << "\t" << maskedArray(row, col) << endl;
            }
        }
    }

    // Reset the matrix
	clear_buffer(input.begin(), input.end());
}

SiPixelCluster WeightedMeans2D::make_cluster(int row, int col) {
    // Get the weight of the pixel we are trying to form a cluster with
    int weight = maskedArray(row, col);

    // Set the value of this pixel to 0 as it is used to form a Digi
    maskedArray.set(row, col, 0);

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
                if (maskedArray(r, c) == weight) {
                    // Add it to the cluster
                    SiPixelCluster::PixelPos newpix(r, c);
                    if (!acluster.add(newpix, 255)) break;
                    // And change its value
                    maskedArray.set(newpix, 0);
                }
            }
        }
    }

    // Form a "real" CMSSW cluster
    return  SiPixelCluster(acluster.isize, acluster.adc, acluster.x, acluster.y, acluster.xmin, acluster.ymin);
}

void WeightedMeans2D::copy_to_buffer(DigiIterator begin, DigiIterator end) {
    // Copy the value of the Digis' ADC to the 2D matrix
    for (DigiIterator di = begin; di != end; ++di) hitArray.set(di->row(), di->column(), di->adc() / 255);
}

void WeightedMeans2D::clear_buffer(DigiIterator begin, DigiIterator end) {
    // Resets the matrix
    for (DigiIterator di = begin; di != end; ++di) {
        hitArray.set(di->row(), di->column(), 0);
        weightArray.set(di->row(), di->column(), 0);
        maskedArray.set(di->row(), di->column(), 0);
    }
}

unsigned int WeightedMeans2D::getLayerNumber(unsigned int & detid) {
    unsigned int layer = 999;
    DetId theDetId(detid);
    if (theDetId.det() == DetId::Tracker) {
        if (theDetId.subdetId() == PixelSubdetector::PixelBarrel) {
            PXBDetId pb_detId = PXBDetId(detid);
            layer = pb_detId.layer();
        }
        else if (theDetId.subdetId() == PixelSubdetector::PixelEndcap) {
            PXFDetId pf_detId = PXFDetId(detid);
            layer = 100 * pf_detId.side() + pf_detId.disk();
        }
        else std::cout << ">>> Invalid subdetId() = " << theDetId.subdetId() << std::endl;
    }
    return layer;
}
