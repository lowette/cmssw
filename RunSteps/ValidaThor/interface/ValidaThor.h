#ifndef __RUNSTEPS_VALIDATHOR__
#define __RUNSTEPS_VALIDATHOR__

#include <memory>
#include <map>

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RunSteps/ValidaThor/interface/ValHits.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>
#include <THStack.h>
#include <TProfile.h>

int verbose = 2;
const int nTypes = 18;

class ValidaThor {

  typedef std::vector< std::pair< PSimHit , std::vector< ValHit > > > V_HIT_CLUSTERS;
  typedef std::map< int , V_HIT_CLUSTERS > M_TRK_HIT_CLUSTERS;

public:
    void analyze(ValHitsCollection* hitsCollection, edm::PSimHitContainer* simHits_B, edm::PSimHitContainer* simHits_E, edm::SimTrackContainer* simTracks, edm::SimVertexContainer* simVertices, TrackerGeometry* tkGeom);

private:
    TH2F* trackerLayout_;
    TH2F* trackerLayoutXY_;
    TH2F* trackerLayoutXYBar_;
    TH2F* trackerLayoutXYEC_;

    struct ClusterHistos {
        THStack* NumberOfClustersSource;
        TH1F* NumberOfClusterPixel;
        TH1F* NumberOfClusterStrip;

        TH1F* NumberOfClustersLink;

        TH1F* NumberOfMatchedHits[nTypes];
        TH1F* NumberOfMatchedClusters[nTypes];
        TH1F* hEfficiency[nTypes];
        TH1F* h_dx_Truth;
        TH1F* h_dy_Truth;

        THStack* ClustersSizeSource;
        TH1F* clusterSizePixel;
        TH1F* clusterSizeStrip;

        TH1F* clusterSize;
        TH1F* clusterSizeX;
        TH1F* clusterSizeY;

        TH1F* clusterShapeX;
        TH1F* clusterShapeY;
        TH2F* localPosXY;
        TH2F* globalPosXY;

        TH2F* localPosXYPixel;
        TH2F* localPosXYStrip;

        TH1F* digiType;
        TH2F* digiPosition;
    };

    std::map< unsigned int, ClusterHistos > layerHistoMap;

public:
    void createLayerHistograms(unsigned int iLayer);
    void createHistograms(unsigned int nLayer);
    unsigned int getLayerNumber(const TrackerGeometry* tkgeom, unsigned int& detid);
    unsigned int getLayerNumber(unsigned int& detid);
};

#endif
