#ifndef RUNSTEPS_CLUSTERIZER_CLASSES_H
#define RUNSTEPS_CLUSTERIZER_CLASSES_H

#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"

namespace {
    struct dictionary {
        edm::Wrapper<PixelClusterSimLink> PixelClusterLink1;
        edm::Wrapper< std::vector<PixelClusterSimLink>  > PixelClusterLink2;
        edm::Wrapper< edm::DetSet<PixelClusterSimLink> > PixelClusterLink3;
        edm::Wrapper< std::vector<edm::DetSet<PixelClusterSimLink> > > PixelClusterLink4;
        edm::Wrapper< edm::DetSetVector<PixelClusterSimLink> > PixelClusterLink5;
    };
}

#endif
