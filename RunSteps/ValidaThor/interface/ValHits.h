#ifndef RunSteps_ValidaThor_ValHit_H
#define RunSteps_ValidaThor_ValHit_H

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"

#include <vector>
#include <map>

typedef std::map< unsigned int, std::vector< struct ValHit > > ValHitsCollection;
typedef std::vector< struct ValHit > ValHitsVector;
typedef struct ValHit ValHit;

struct ValHit {
    double x, y, xx, xy, yy;
    std::vector< unsigned int > simTracks;
    edm::Ref< edmNew::DetSetVector< SiPixelCluster >, SiPixelCluster > cluster;
};

ValHitsCollection ValHitsBuilder(edm::DetSetVector< PixelClusterSimLink >* clusterLinks);
ValHitsCollection ValHitsBuilder(edm::DetSetVector< PixelClusterSimLink >* clusterLinks, edmNew::DetSetVector< SiPixelRecHit >* recHits);

#endif
