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

class ClusterValidationTest : public edm::EDAnalyzer {

public:

    explicit ClusterValidationTest(const edm::ParameterSet&);
    ~ClusterValidationTest();
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

private:
    bool PRINT;

    TH2F* trackerLayout_;
    TH2F* trackerLayoutXY_;
    TH2F* trackerLayoutXYBar_;
    TH2F* trackerLayoutXYEC_;

    struct ClusterHistos {
        TH1F* clusterSize_;
        TH1F* clusterSizeX_;
        TH1F* clusterSizeY_;
        TH1F* numberOfPixels_;
        TH2F* localPosXY_;
        TH2F* globalPosXY_;
    };

    std::map<unsigned int, ClusterHistos> layerHistoMap;

    edm::InputTag src_;

public:
    void createLayerHistograms(unsigned int iLayer);
    void createHistograms(unsigned int nLayer);
    unsigned int getLayerNumber(const TrackerGeometry* tkgeom, unsigned int& detid);
    unsigned int getLayerNumber(unsigned int& detid);
};

ClusterValidationTest::ClusterValidationTest(const edm::ParameterSet& iConfig) {
    PRINT = iConfig.getUntrackedParameter<bool>("Verbosity",false);
    src_ = iConfig.getParameter<edm::InputTag>("src");
    if (PRINT) std::cout << ">>> Construct ClusterValidationTest " << std::endl;
}

ClusterValidationTest::~ClusterValidationTest() {
    if (PRINT) std::cout << ">>> Destroy ClusterValidationTest " << std::endl;
}

void ClusterValidationTest::beginJob() {
    using namespace edm;
    if (PRINT) std::cout << "Initialize ClusterValidationTest " << std::endl;
    createHistograms(19);
}

void ClusterValidationTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;

    //int event = iEvent.id().event();

    edm::Handle< edm::DetSetVector<SiPixelCluster> > pixelClusters;
    iEvent.getByLabel(src_, pixelClusters);

    edm::ESHandle<TrackerGeometry> geomHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get(geomHandle);
    const TrackerGeometry*  tkGeom = &(*geomHandle);

    edm::DetSetVector<SiPixelCluster>::const_iterator DSViter;
    for (DSViter = pixelClusters->begin(); DSViter != pixelClusters->end(); DSViter++) {
        unsigned int rawid = DSViter->detId();
        DetId detId(rawid);

        unsigned int layer = getLayerNumber(rawid);
        std::map<unsigned int, ClusterHistos>::iterator iPos = layerHistoMap.find(layer);
        if (iPos == layerHistoMap.end()) {
            createLayerHistograms(layer);
            iPos = layerHistoMap.find(layer);
        }
        const GeomDetUnit* geomDetUnit = tkGeom->idToDetUnit(detId);
        //const PixelGeomDetUnit* pixdet = (PixelGeomDetUnit*) geomDetUnit;

        edm::DetSet<SiPixelCluster>::const_iterator cu;
        for (cu = DSViter->data.begin(); cu != DSViter->data.end(); ++cu) {
            int size = cu->size();
            int sizeX = cu->sizeX();
            int sizeY = cu->sizeY();
            float x = cu->x();
            float y = cu->y();
            const vector<SiPixelCluster::Pixel>& pixelsVec = cu->pixels();
            float numPixels = pixelsVec.size();

            iPos->second.clusterSize_->Fill(size);
            iPos->second.clusterSizeX_->Fill(sizeX);
            iPos->second.clusterSizeY_->Fill(sizeY);
            iPos->second.numberOfPixels_->Fill(numPixels);

            MeasurementPoint mp(x, y);

            if (geomDetUnit) {
                LocalPoint lPos = geomDetUnit->topology().localPosition(mp);

                GlobalPoint pdPos = geomDetUnit->surface().toGlobal(geomDetUnit->topology().localPosition(mp));

                iPos->second.localPosXY_->Fill(lPos.x(), lPos.y());
                iPos->second.globalPosXY_->Fill(pdPos.x(), pdPos.y());

                trackerLayout_->Fill(pdPos.z(), pdPos.perp());
                trackerLayoutXY_->Fill(pdPos.y(), pdPos.x());
                if (layer < 100) trackerLayoutXYBar_->Fill(pdPos.y(), pdPos.x());
                else trackerLayoutXYEC_->Fill(pdPos.y(), pdPos.x());
            }
        }
    }
}

void ClusterValidationTest::endJob() { }

void ClusterValidationTest::createLayerHistograms(unsigned int ival) {
    std::ostringstream fname1, fname2;

    edm::Service<TFileService> fs;
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
    htit4 << "numberOfPixels" << tag.c_str() <<  id;
    local_histos.numberOfPixels_ = td.make<TH1F>(htit4.str().c_str(), htit4.str().c_str(), 1000, 0., 0.);

    std::ostringstream htit5;
    htit5 << "LocalPositionXY" << tag.c_str() <<  id;
    local_histos.localPosXY_ = td.make<TH2F>(htit5.str().c_str(), htit5.str().c_str(), 2000, 0., 0., 2000, 0., 0.);

    std::ostringstream htit6;
    htit6 << "GlobalPositionXY" << tag.c_str() <<  id;
    local_histos.globalPosXY_ = td.make<TH2F>(htit6.str().c_str(), htit6.str().c_str(), 2400, -120.0, 120.0, 2400, -120.0, 120.0);

    layerHistoMap.insert(std::make_pair(ival, local_histos));
    fs->file().cd("/");
}

void ClusterValidationTest::createHistograms(unsigned int nLayer) {
    edm::Service<TFileService> fs;
    fs->file().cd("/");
    TFileDirectory td = fs->mkdir("Common");

    trackerLayout_ = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 1200, 0.0, 120.0);
    trackerLayoutXY_ = td.make<TH2F>("XVsY", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
    trackerLayoutXYBar_ = td.make<TH2F>("XVsYBar", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
    trackerLayoutXYEC_ = td.make<TH2F>("XVsYEC", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
}

unsigned int ClusterValidationTest::getLayerNumber(const TrackerGeometry* tkgeom, unsigned int& detid) {
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

unsigned int ClusterValidationTest::getLayerNumber(unsigned int& detid) {
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

DEFINE_FWK_MODULE(ClusterValidationTest);

