#include "RunSteps/ValidaThor/interface/ValHit.h"

#include <iostream>

ValHit::ValHit(edm::DetSetVector< PixelClusterSimLink >* clusterLinks) {

    edm::DetSetVector< PixelClusterSimLink >::const_iterator DSViter;
    edm::DetSet< PixelClusterSimLink >::const_iterator clusterLink;

    for (DSViter = clusterLinks->begin(); DSViter != clusterLinks->end(); ++DSViter) {

        // Loop over cluster links
        for (clusterLink = DSViter->data.begin(); clusterLink != DSViter->data.end(); ++clusterLink) {

            PixelClusterSimLink link = *clusterLink;

            // Get the cluster
            edm::Ref< edmNew::DetSetVector< SiPixelCluster >, SiPixelCluster > const& cluster = link.getCluster();

            // Get the parameters
            HitSidor newHit;
            newHit.x = cluster->x();
            newHit.y = cluster->y();
            newHit.simTracks = link.getSimTracks();


            // Add the Hit
            hits_[DSViter->detId()].push_back(newHit);

        }
    }
}

ValHit::ValHit(edm::DetSetVector< PixelClusterSimLink >* clusterLinks, edmNew::DetSetVector< SiPixelRecHit >* recHits) {

    edmNew::DetSetVector< SiPixelRecHit >::const_iterator DSViter;
    edmNew::DetSet< SiPixelRecHit >::const_iterator rechHitIter;

    edm::DetSet< PixelClusterSimLink >::const_iterator clusterLink;

    for (DSViter = recHits->begin(); DSViter != recHits->end(); ++DSViter) {

        // Get the cluster simlinks DetSet
        edm::DetSet< PixelClusterSimLink > clusters = (*clusterLinks)[DSViter->id()];

        // Loop over recHits
        for (rechHitIter = DSViter->begin(); rechHitIter != DSViter->end(); ++rechHitIter) {

            // Get the cluster
            edm::Ref< edmNew::DetSetVector< SiPixelCluster >, SiPixelCluster > const& cluster = rechHitIter->cluster();

            // Get the parameters
            HitSidor newHit;
            newHit.x = rechHitIter->localPosition().x();
            newHit.y = rechHitIter->localPosition().y();

            bool clusterFound = false;

            // Loop over the clusters
            for (clusterLink = clusters.begin(); clusterLink != clusters.end(); ++clusterLink) {
                PixelClusterSimLink link = *clusterLink;

                edm::Ref< edmNew::DetSetVector< SiPixelCluster >, SiPixelCluster > const& clusterFromLink = link.getCluster();

                // Compare the clusters
                if (cluster == clusterFromLink) {
                    newHit.simTracks = link.getSimTracks();
                    clusterFound = true;
                    break;
                }
            }

            // Add the Hit
            if (clusterFound) hits_[DSViter->detId()].push_back(newHit);

        }
    }
}
