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
	theBuffer.setSize(nrows_, ncols_);
}

PixelClusterizer::~PixelClusterizer() {}

bool PixelClusterizer::setup(const PixelGeomDetUnit * pixDet) {
	const PixelTopology & topol = pixDet->specificTopology();
	nrows_ = topol.nrows(); 
	ncols_ = topol.ncolumns();
	if (nrows_ > theBuffer.rows() || ncols_ > theBuffer.columns()) theBuffer.setSize(nrows_, ncols_);
	return true;   
}

void PixelClusterizer::clusterizeDetUnit(const edm::DetSet<PixelDigi> & input, const PixelGeomDetUnit * pixDet, std::vector<SiPixelCluster> & output) {
	if (!setup(pixDet)) return;
  	
	detid_ = input.detId();
  	copy_to_buffer(input.begin(), input.end());

	for (DigiIterator it = input.begin(); it != input.end(); ++it) {
      		if (it->adc() == 255) {
	  		SiPixelCluster cluster = make_cluster(SiPixelCluster::PixelPos(it->row(), it->column()));
		}
    	}
	clear_buffer(input.begin(), input.end());
}

void PixelClusterizer::copy_to_buffer(DigiIterator begin, DigiIterator end) {	
	for (DigiIterator di = begin; di != end; ++di) theBuffer.set_adc(di->row(), di->column(), di->adc());
}

void PixelClusterizer::clear_buffer(DigiIterator begin, DigiIterator end) {	
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
	theBuffer.set_adc(pix, 127);
  
  	AccretionCluster acluster;
  	acluster.add(pix, 255);
  
  	while (! acluster.empty()) {
      		auto curInd = acluster.top(); 
		acluster.pop();
      		for (auto r = acluster.x[curInd] - 1; r <= acluster.x[curInd] + 1; ++r) {
	  		for (auto c = acluster.y[curInd] - 1; c <= acluster.y[curInd] + 1; ++c) {
	      			if (theBuffer(r, c) == 255) {
		  			SiPixelCluster::PixelPos newpix(r, c);
		  			if (!acluster.add(newpix, theBuffer(r, c))) break;
		  			theBuffer.set_adc(newpix, 127);
				}

	    		}
		}
      
    	} 
  
	SiPixelCluster cluster(acluster.isize, acluster.adc, acluster.x, acluster.y, acluster.xmin, acluster.ymin);
 	
	return cluster;
}

