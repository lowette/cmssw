#ifndef RUNSTEPS_PIXELCLUSTERSIMLINK_H_
#define RUNSTEPS_PIXELCLUSTERSIMLINK_H_

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

class PixelClusterSimLink {

public:
    PixelClusterSimLink() { };
    void addSimTrack(unsigned int simTrack) { simTracks_.push_back(simTrack); };
    void setCluster(SiPixelCluster cluster) { cluster_ = cluster; };

    std::vector< unsigned int > getSimTracks() { return simTracks_; };
    SiPixelCluster getCluster() { return cluster_; };

    inline bool operator< ( const PixelClusterSimLink& other ) const { return true; }

private:
    SiPixelCluster cluster_;
    std::vector< unsigned int > simTracks_;

};

#endif
