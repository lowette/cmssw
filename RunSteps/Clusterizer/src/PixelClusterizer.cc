#include "RunSteps/Clusterizer/interface/PixelClusterizer.h"
#include "RunSteps/Clusterizer/interface/SiPixelArrayBuffer.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelGainCalibrationOffline.h"

#include <stack>
#include <vector>
#include <iostream>

using namespace std;

PixelClusterizer::PixelClusterizer(edm::ParameterSet const& conf) :
	conf_(conf),
	nrows_(0),
	ncols_(0),
	detid_(0) {
	// Create a 2D matrix for this detector
	theBuffer.setSize(nrows_, ncols_);
}

PixelClusterizer::~PixelClusterizer() {}

// Change the size of the 2D matrix for this detector unit
bool PixelClusterizer::setup(const PixelGeomDetUnit * pixDet) {
	const PixelTopology & topol = pixDet->specificTopology();
	nrows_ = topol.nrows();
	ncols_ = topol.ncolumns();
	if (nrows_ > theBuffer.rows() || ncols_ > theBuffer.columns()) theBuffer.setSize(nrows_, ncols_);
	return true;
}

// Go over the Digis and create clusters
void PixelClusterizer::clusterizeDetUnit(const edm::DetSet<PixelDigi> & input, const PixelGeomDetUnit * pixDet, std::vector<SiPixelCluster> & output) {
	if (!setup(pixDet)) return;

	detid_ = input.detId();

	// Fill the 2D matrix with the ADC values
  	copy_to_buffer(input.begin(), input.end());

  	// Go over the Digis
	for (DigiIterator it = input.begin(); it != input.end(); ++it) {
		// Check if the Digi can be used
		// Non hit pixels have a value of 0
		// Hit pixels have a value of 255
		// Hit pixels that have already been used in other clusters have a value of 127
      	if (theBuffer(it->row(), it->column()) == 255) {
      		// Use the Digi to form a cluster
			SiPixelCluster cluster = make_cluster(SiPixelCluster::PixelPos(it->row(), it->column()));
			// Add the cluster to the list (Here would come some cuts on the size, ...)
	  		output.push_back(cluster);
		}
    }
    // Reset the matrix
	clear_buffer(input.begin(), input.end());
}

void PixelClusterizer::copy_to_buffer(DigiIterator begin, DigiIterator end) {
	// Copy the value of the Digis' ADC to the 2D matrix
	for (DigiIterator di = begin; di != end; ++di) theBuffer.set_adc(di->row(), di->column(), di->adc());
}

void PixelClusterizer::clear_buffer(DigiIterator begin, DigiIterator end) {
	// Resets the matrix
	for (DigiIterator di = begin; di != end; ++di) theBuffer.set_adc(di->row(), di->column(), 0);
}

namespace {
	struct AccretionCluster {
		static constexpr unsigned short MAXSIZE = 256;
		unsigned short adc[256];
		unsigned short x[256];
		unsigned short y[256];
		unsigned short xmin = 16000;
		unsigned short ymin = 16000;
		unsigned int isize = 0;
		unsigned int curr = 0;
		unsigned short top() const { return curr; }
		void pop() { ++curr; }
		bool empty() { return curr == isize; }
		bool add(SiPixelCluster::PixelPos const & p, unsigned short const iadc) {
			if (isize == MAXSIZE) return false;
			xmin = std::min(xmin, (unsigned short) p.row());
			ymin = std::min(ymin, (unsigned short) p.col());
			adc[isize] = iadc;
			x[isize] = p.row();
			y[isize++] = p.col();
			return true;
		}
	};
}

SiPixelCluster PixelClusterizer::make_cluster( const SiPixelCluster::PixelPos& pix) {
	// Set the value of this pixel to 127 as it is used to form a Digi
	theBuffer.set_adc(pix, 127);

	// Create a temporary cluster (this allows to them easily form a "real" cluster with CMSSW data format)
  	AccretionCluster acluster;
  	// Add the pixel to the cluster
  	acluster.add(pix, 255);

  	// Go over all the pixels in the cluster
  	while (!acluster.empty()) {
  		// Get the current pixel we are looking at
      	auto curInd = acluster.top();
		acluster.pop();
		// Look left and right
      	for (auto r = acluster.x[curInd] - 1; r <= acluster.x[curInd] + 1; ++r) {
      		// Look bottom and top
	  		for (auto c = acluster.y[curInd] - 1; c <= acluster.y[curInd] + 1; ++c) {
	  			// If the pixel is hit and unused
	      		if (theBuffer(r, c) == 255) {
	      			// Add it to the cluster
		  			SiPixelCluster::PixelPos newpix(r, c);
		  			if (!acluster.add(newpix, theBuffer(r, c))) break;
		  			// And change its value
		  			theBuffer.set_adc(newpix, 127);
				}
	    	}
		}
    }

    // Form a "real" CMSSW cluster
	SiPixelCluster cluster(acluster.isize, acluster.adc, acluster.x, acluster.y, acluster.xmin, acluster.ymin);

	// Return it
	return cluster;
}
