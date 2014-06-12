#include "RunSteps/ValidaThor/interface/ValHits.h"

ValHitsCollection ValHitsBuilder(edm::DetSetVector< PixelClusterSimLink >* clusterLinks) {

    ValHitsCollection hits;

    edm::DetSetVector< PixelClusterSimLink >::const_iterator DSViter;
    edm::DetSet< PixelClusterSimLink >::const_iterator clusterLink;

    for (DSViter = clusterLinks->begin(); DSViter != clusterLinks->end(); ++DSViter) {

        // Loop over cluster links
        for (clusterLink = DSViter->data.begin(); clusterLink != DSViter->data.end(); ++clusterLink) {

            PixelClusterSimLink link = *clusterLink;

            // Get the cluster
            edm::Ref< edmNew::DetSetVector< SiPixelCluster >, SiPixelCluster > const& cluster = link.getCluster();

            // Get the parameters
            ValHit newHit;
            newHit.x = cluster->x();
            newHit.y = cluster->y();
            newHit.simTracks = link.getSimTracks();


            // Add the Hit
            hits[DSViter->detId()].push_back(newHit);

        }
    }

    return hits;
}

ValHitsCollection ValHitsBuilder(edm::DetSetVector< PixelClusterSimLink >* clusterLinks, edmNew::DetSetVector< SiPixelRecHit >* recHits) {

    ValHitsCollection hits;

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
            ValHit newHit;
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
            if (clusterFound) hits[DSViter->detId()].push_back(newHit);

        }
    }

    return hits;
}
