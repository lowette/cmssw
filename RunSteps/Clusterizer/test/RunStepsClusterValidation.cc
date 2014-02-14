#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TProfile.h>

using namespace std;
using namespace edm;

class RunStepsClusterValidation : public EDAnalyzer {

public:

    explicit RunStepsClusterValidation(const ParameterSet&);
    ~RunStepsClusterValidation();
    virtual void beginJob();
    virtual void analyze(const Event&, const EventSetup&);
    virtual void endJob();

private:
    TH2F* trackerLayout_;
    TH2F* trackerLayoutXY_;
    TH2F* trackerLayoutXYBar_;
    TH2F* trackerLayoutXYEC_;

    struct ClusterHistos {
        TH1F* clusterSize_;
        TH1F* clusterSizeX_;
        TH1F* clusterSizeY_;
        TH1F* clusterShapeX_;
        TH1F* clusterShapeY_;
        TH2F* localPosXY_;
        TH2F* globalPosXY_;
    };

    std::map<unsigned int, ClusterHistos> layerHistoMap;

    InputTag src_;

public:
    void createLayerHistograms(unsigned int iLayer);
    void createHistograms(unsigned int nLayer);
    unsigned int getLayerNumber(const TrackerGeometry* tkgeom, unsigned int& detid);
    unsigned int getLayerNumber(unsigned int& detid);
};

RunStepsClusterValidation::RunStepsClusterValidation(const ParameterSet& iConfig) {
    src_ = iConfig.getParameter<InputTag>("src");

    std::cout << "------------------------------------------------------------" << std::endl
              << "-- Running RunSteps ClusterValidation v0.0" << std::endl
              << "------------------------------------------------------------" << std::endl;
}

RunStepsClusterValidation::~RunStepsClusterValidation() { }

void RunStepsClusterValidation::beginJob() {
    createHistograms(19);
}

void RunStepsClusterValidation::analyze(const Event& iEvent, const EventSetup& iSetup) {

    // Get the clusters
    Handle< DetSetVector<SiPixelCluster> > pixelClusters;
    iEvent.getByLabel(src_, pixelClusters);

    // Get the geometry
    ESHandle<TrackerGeometry> geomHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get(geomHandle);
    const TrackerGeometry*  tkGeom = &(*geomHandle);

    // Go over the detector units
    DetSetVector<SiPixelCluster>::const_iterator DSViter;
    for (DSViter = pixelClusters->begin(); DSViter != pixelClusters->end(); DSViter++) {
        // Get the detector unit's id
        unsigned int rawid = DSViter->detId();
        unsigned int layer = getLayerNumber(rawid);
        DetId detId(rawid);
        const GeomDetUnit* geomDetUnit = tkGeom->idToDetUnit(detId);

        if (!geomDetUnit) break;

        // Create histograms for the layer if they do not yet exist
        std::map<unsigned int, ClusterHistos>::iterator iPos = layerHistoMap.find(layer);
        if (iPos == layerHistoMap.end()) {
            createLayerHistograms(layer);
            iPos = layerHistoMap.find(layer);
        }

        // Go over the clusters in the detector unit
        DetSet<SiPixelCluster>::const_iterator cu;
        for (cu = DSViter->data.begin(); cu != DSViter->data.end(); ++cu) {
            // Get the cluster's size
            int size = cu->size();
            int sizeX = cu->sizeX();
            int sizeY = cu->sizeY();

            iPos->second.clusterSize_->Fill(size);
            iPos->second.clusterSizeX_->Fill(sizeX);
            iPos->second.clusterSizeY_->Fill(sizeY);
            // Get the cluster's local
            float x = cu->x();
            float y = cu->y();

            // Get the cluster's global position
            MeasurementPoint mp(x, y);
            LocalPoint lPos = geomDetUnit->topology().localPosition(mp);
            GlobalPoint gPos = geomDetUnit->surface().toGlobal(geomDetUnit->topology().localPosition(mp));

            // Fill the histograms
            iPos->second.localPosXY_->Fill(lPos.x(), lPos.y());
            iPos->second.globalPosXY_->Fill(gPos.x(), gPos.y());

            trackerLayout_->Fill(gPos.z(), gPos.perp());
            trackerLayoutXY_->Fill(gPos.y(), gPos.x());
            if (layer < 100) trackerLayoutXYBar_->Fill(gPos.y(), gPos.x());
            else trackerLayoutXYEC_->Fill(gPos.y(), gPos.x());

            // Get the pixels that form the Cluster
            const vector<SiPixelCluster::Pixel>& pixelsVec = cu->pixels();

            // Go over the pixels
            vector<SiPixelCluster::Pixel>::const_iterator pixelIt;
            for (pixelIt = pixelsVec.begin(); pixelIt != pixelsVec.end(); ++pixelIt) {
                //////////////////////////
                // NOT WORKING !!!!!!   //
                //////////////////////////
                iPos->second.clusterShapeX_->Fill(gPos.x() - pixelIt->x);
                iPos->second.clusterShapeY_->Fill(gPos.y() - pixelIt->y);
            }
        }
    }
}

void RunStepsClusterValidation::endJob() { }

// Create the histograms
void RunStepsClusterValidation::createLayerHistograms(unsigned int ival) {
    std::ostringstream fname1, fname2;

    Service<TFileService> fs;
    fs->file().cd("/");

    std::string tag;
    unsigned int id;
    if (ival < 100) {
        id = ival;
        fname1 << "Barrel";
        fname2 << "Layer_" << id;
        tag = "_layer_";
    } else {
        int side = ival / 100;
        id = ival - side*100;
        std::cout << " Creating histograms for Disc " << id << " with " << ival << std::endl;
        fname1 << "EndCap_Side_" << side;
        fname2 << "Disc_" << id;
        tag = "_disc_";
    }
    TFileDirectory td1 = fs->mkdir(fname1.str().c_str());
    TFileDirectory td = td1.mkdir(fname2.str().c_str());

    ClusterHistos local_histos;

    std::ostringstream htit1;
    htit1 << "ClusterSize" << tag.c_str() <<  id;
    local_histos.clusterSize_ = td.make<TH1F>(htit1.str().c_str(), htit1.str().c_str(), 1000, 0., 0.);

    std::ostringstream htit2;
    htit2 << "ClusterSizeX" << tag.c_str() <<  id;
    local_histos.clusterSizeX_ = td.make<TH1F>(htit2.str().c_str(), htit2.str().c_str(), 1000, 0., 0.);

    std::ostringstream htit3;
    htit3 << "ClusterSizeY" << tag.c_str() <<  id;
    local_histos.clusterSizeY_ = td.make<TH1F>(htit3.str().c_str(), htit3.str().c_str(), 1000, 0., 0.);

    std::ostringstream htit4;
    htit4 << "ClusterShapeX" << tag.c_str() <<  id;
    local_histos.clusterShapeX_ = td.make<TH1F>(htit4.str().c_str(), htit4.str().c_str(), 1000, 0., 0.);

    std::ostringstream htit5;
    htit5 << "ClusterShapeY" << tag.c_str() <<  id;
    local_histos.clusterShapeY_ = td.make<TH1F>(htit5.str().c_str(), htit5.str().c_str(), 1000, 0., 0.);

    std::ostringstream htit7;
    htit7 << "LocalPositionXY" << tag.c_str() <<  id;
    local_histos.localPosXY_ = td.make<TH2F>(htit7.str().c_str(), htit7.str().c_str(), 2000, 0., 0., 2000, 0., 0.);

    std::ostringstream htit8;
    htit8 << "GlobalPositionXY" << tag.c_str() <<  id;
    local_histos.globalPosXY_ = td.make<TH2F>(htit8.str().c_str(), htit8.str().c_str(), 2400, -120.0, 120.0, 2400, -120.0, 120.0);

    layerHistoMap.insert(std::make_pair(ival, local_histos));
    fs->file().cd("/");
}

void RunStepsClusterValidation::createHistograms(unsigned int nLayer) {
    Service<TFileService> fs;
    fs->file().cd("/");
    TFileDirectory td = fs->mkdir("Common");

    trackerLayout_ = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 1200, 0.0, 120.0);
    trackerLayoutXY_ = td.make<TH2F>("XVsY", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
    trackerLayoutXYBar_ = td.make<TH2F>("XVsYBar", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
    trackerLayoutXYEC_ = td.make<TH2F>("XVsYEC", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
}

unsigned int RunStepsClusterValidation::getLayerNumber(const TrackerGeometry* tkgeom, unsigned int& detid) {
    unsigned int layer = 999;
    DetId theDetId(detid);
    if (theDetId.subdetId() != 1) std::cout << ">>> Method1 : Det id " << theDetId.det() << " Subdet Id " << theDetId.subdetId() << std::endl;
    const PixelGeomDetUnit * theGeomDet = dynamic_cast<const PixelGeomDetUnit*> (tkgeom->idToDet(theDetId));

    const GeomDetUnit* it = tkgeom->idToDetUnit(DetId(theDetId));
    if (!it) std::cout << ">>> rawdetid " << detid << " GeomDetUnit " << it << " PixelGeomDetUnit " << theGeomDet << " DetId " << theDetId.det() << " Subdet Id " << theDetId.subdetId() << std::endl;
    if (it && it->type().isTracker()) {
        if (it->type().isBarrel()) {
            PXBDetId pb_detId = PXBDetId(detid);
            layer = pb_detId.layer();
        }
        else if (it->type().isEndcap()) {
            std::cout << " IAM HERE >>>>>>>>>>>>>>>>>>>>>>>>>>> " << std::endl;
            PXFDetId pf_detId = PXFDetId(detid);
            layer = 100*pf_detId.side() + pf_detId.disk();
        }
    }
    return layer;
}

unsigned int RunStepsClusterValidation::getLayerNumber(unsigned int & detid) {
    unsigned int layer = 999;
    DetId theDetId(detid);
    if (theDetId.det() == DetId::Tracker) {
        if (theDetId.subdetId() == PixelSubdetector::PixelBarrel) {
            PXBDetId pb_detId = PXBDetId(detid);
            layer = pb_detId.layer();
        }
        else if (theDetId.subdetId() == PixelSubdetector::PixelEndcap) {
            PXFDetId pf_detId = PXFDetId(detid);
            layer = 100 * pf_detId.side() + pf_detId.disk();
        }
        else std::cout << ">>> Invalid subdetId() = " << theDetId.subdetId() << std::endl;
    }
    return layer;
}

DEFINE_FWK_MODULE(RunStepsClusterValidation);
