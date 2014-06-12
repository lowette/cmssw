#include "RunSteps/Clusterizer/interface/PixelClusterizer.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"

using namespace std;

void PixelClusterizer::makeLinks(edm::OrphanHandle< edmNew::DetSetVector<SiPixelCluster> > & clusters, std::vector<edm::DetSet<PixelClusterSimLink> > & linksByDet) {

    // Go over all the detectors
    for (edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter = clusters->begin(); DSViter != clusters->end(); ++DSViter) {
        DetId detIdObject(DSViter->detId());

        edm::DetSet<PixelClusterSimLink> links(DSViter->detId());

        for (edmNew::DetSet<SiPixelCluster>::const_iterator clustIt = DSViter->begin(); clustIt != DSViter->end(); clustIt++) {

            std::map< SiPixelCluster, std::vector< unsigned int > >::iterator it = tmpSimLinks.find(*clustIt);

            if (it != tmpSimLinks.end()) {
                std::vector< unsigned int > simTracks = it->second;

                PixelClusterSimLink simLink;
                simLink.setCluster(edmNew::makeRefTo(clusters, clustIt));
                simLink.setSimTracks(simTracks);
                links.data.push_back(simLink);
            }
        }

        linksByDet.push_back(links);
    }
}

unsigned int PixelClusterizer::getSimTrackId(const edm::Handle< edm::DetSetVector< PixelDigiSimLink > > & pixelSimLinks, int channel) {
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
