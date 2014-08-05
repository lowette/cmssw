#ifndef RUNSTEPS_PIXELCLUSTERSIMLINK_H_
#define RUNSTEPS_PIXELCLUSTERSIMLINK_H_

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Ref.h"

class PixelClusterSimLink {

public:
    PixelClusterSimLink() { };
    void setSimTracks(std::vector< unsigned int > simTrack) { simTracks_ = simTrack; };
    void setCluster(edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster>  cluster) { cluster_ = cluster; };

    std::vector< unsigned int > getSimTracks() { return simTracks_; };
    edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster>  getCluster() { return cluster_; };

    inline bool operator< ( const PixelClusterSimLink& other ) const { return true; }

private:
    edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster>  cluster_;
    std::vector< unsigned int > simTracks_;

};

#endif
