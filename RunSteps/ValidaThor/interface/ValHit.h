#ifndef RunSteps_ValidaThor_ValHit_H
#define RunSteps_ValidaThor_ValHit_H

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"

#include <vector>
#include <map>

struct HitSidor {
    double x, y;
    std::vector< unsigned int > simTracks;
};

class ValHit {

public:
    ValHit(edm::DetSetVector< PixelClusterSimLink >* clusterLinks);
    ValHit(edm::DetSetVector< PixelClusterSimLink >* clusterLinks, edmNew::DetSetVector< SiPixelRecHit >* recHits);

private:
    std::map< unsigned int, std::vector< struct HitSidor > > hits_;

};

#endif
