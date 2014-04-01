#include <memory>

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>
#include <THStack.h>
#include <TProfile.h>

int verbose=2;

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
        THStack* NumberOfClustersSource;
        TH1F* NumberOfClusterPixel;
        TH1F* NumberOfClusterStrip;

        TH1F* NumberOfClustersLink;

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
              << "-- Running RunSteps ClusterValidation v0.2" << std::endl
              << "------------------------------------------------------------" << std::endl;
}

RunStepsClusterValidation::~RunStepsClusterValidation() { }

void RunStepsClusterValidation::beginJob() {
    createHistograms(19);
}

void RunStepsClusterValidation::analyze(const Event& iEvent, const EventSetup& iSetup) {

    // Simulated information
    // SimHit
    Handle<PSimHitContainer> simHits;
    iEvent.getByLabel("g4SimHits","TrackerHitsPixelBarrelLowTof" ,simHits);

    // SimTrack
    Handle<SimTrackContainer> simTracks;
    iEvent.getByLabel("g4SimHits",simTracks);

    // SimVertex
    Handle<SimVertexContainer> simVertices;
    iEvent.getByLabel("g4SimHits", simVertices);

    // Get the clusters
    Handle< DetSetVector<SiPixelCluster> > pixelClusters;
    iEvent.getByLabel(src_, pixelClusters);

    // Get the links
    Handle< DetSetVector<PixelClusterSimLink> > clusterLinks;
    iEvent.getByLabel(src_, clusterLinks);

    // Get the geometry
    ESHandle<TrackerGeometry> geomHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get(geomHandle);
    const TrackerGeometry* tkGeom = &(*geomHandle);

    // Go over the detector units
    DetSetVector<SiPixelCluster>::const_iterator DSViter;
    for (DSViter = pixelClusters->begin(); DSViter != pixelClusters->end(); DSViter++) {
        // Clusters
        unsigned int nClusters = 0;

        // Get the detector unit's id
        unsigned int rawid = DSViter->detId();
        unsigned int layer = getLayerNumber(rawid);
        DetId detId(rawid);

        // Get the geometry of the tracker
        const GeomDetUnit* geomDetUnit = tkGeom->idToDetUnit(detId);
        const PixelGeomDetUnit* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(geomDetUnit);
        const PixelTopology & topol = theGeomDet->specificTopology();

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

            iPos->second.clusterSize->Fill(size);
            iPos->second.clusterSizeX->Fill(sizeX);
            iPos->second.clusterSizeY->Fill(sizeY);

            // Get the cluster's local
            float x = cu->x();
            float y = cu->y();

            // Get the cluster's global position
            MeasurementPoint mp(x, y);
            LocalPoint lPos = geomDetUnit->topology().localPosition(mp);
            GlobalPoint gPos = geomDetUnit->surface().toGlobal(geomDetUnit->topology().localPosition(mp));

            // Fill the histograms
            iPos->second.localPosXY->Fill(lPos.x(), lPos.y());
            iPos->second.globalPosXY->Fill(gPos.x(), gPos.y());

            trackerLayout_->Fill(gPos.z(), gPos.perp());
            trackerLayoutXY_->Fill(gPos.y(), gPos.x());
            if (layer < 100) trackerLayoutXYBar_->Fill(gPos.y(), gPos.x());
            else trackerLayoutXYEC_->Fill(gPos.y(), gPos.x());

            // Pixel module
            if (topol.ncolumns() == 32) {
                iPos->second.localPosXYPixel->Fill(lPos.x(), lPos.y());
                iPos->second.clusterSizePixel->Fill(size);
            }
            // Strip module
            else if (topol.ncolumns() == 2) {
                iPos->second.localPosXYStrip->Fill(lPos.x(), lPos.y());
                iPos->second.clusterSizeStrip->Fill(size);
            }

            // Get the pixels that form the Cluster
            const vector<SiPixelCluster::Pixel>& pixelsVec = cu->pixels();

            // Go over the pixels
            for (vector<SiPixelCluster::Pixel>::const_iterator pixelIt = pixelsVec.begin(); pixelIt != pixelsVec.end(); ++pixelIt) {
                SiPixelCluster::Pixel PDigi = (SiPixelCluster::Pixel) *pixelIt;

                iPos->second.digiPosition->Fill(PDigi.x, PDigi.y);

                //////////////////////////
                // NOT WORKING !!!!!!   //
                //////////////////////////
                iPos->second.clusterShapeX->Fill(lPos.x() - PDigi.x);
                iPos->second.clusterShapeY->Fill(lPos.y() - PDigi.y);
            }

            ++nClusters;
        }

        // Pixel module
        if (topol.ncolumns() == 32) iPos->second.NumberOfClusterPixel->Fill(nClusters);
        // Strip module
        else if (topol.ncolumns() == 2) iPos->second.NumberOfClusterStrip->Fill(nClusters);
    }

    // Go over the links
    DetSetVector<PixelClusterSimLink>::const_iterator DSViterLinks;
    for (DSViterLinks = clusterLinks->begin(); DSViterLinks != clusterLinks->end(); DSViterLinks++) {

        // Get the detector unit's id
        unsigned int rawid = DSViterLinks->detId();
        unsigned int layer = getLayerNumber(rawid);

        // Create histograms for the layer if they do not yet exist
        std::map<unsigned int, ClusterHistos>::iterator iPos = layerHistoMap.find(layer);
        if (iPos == layerHistoMap.end()) {
            createLayerHistograms(layer);
            iPos = layerHistoMap.find(layer);
        }

        unsigned int nLinks = 0;

        // Go over the links in the detector unit
	std::vector< unsigned int > simTrackID;
        DetSet<PixelClusterSimLink>::const_iterator link;
	PixelClusterSimLink li;
	uint trkID=-1;
	uint sizeLink=0;
	bool combinatoric=false;

	if(verbose>1) cout << "### SimLinks ###" << endl;

        for (link = DSViterLinks->data.begin(); link != DSViterLinks->data.end(); ++link) {
            // Use the link
	    combinatoric=false;
            nLinks++;
	    li = *link;
	    simTrackID    = li.getSimTracks();
	    sizeLink      = simTrackID.size();
	    if(verbose>1) cout << sizeLink << " SimTracks | " ;

	    for(uint i=0 ; i<sizeLink ; i++) {	      
	      if(verbose>1) cout << simTrackID[i] << " | " ;
	      if(i==0) trkID=simTrackID[i];
	      else if(simTrackID[i]!=trkID) combinatoric=true; 
	    }
	    if(combinatoric && verbose>1) cout << " COMBINATORIC !!! ";
	    cout << endl;

	    // Find SimHit corresponding to SimTrack matched to cluster
	    // I should build a map <SimTrack ID , vector<SimHit> (all layers) >
	    /*
	    if(!combinatoric) {
	      for(

		  }
	    */

        }

        iPos->second.NumberOfClustersLink->Fill(nLinks);

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
    }
    else {
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

    std::ostringstream histoName;

    histoName.str("");
    histoName << "Number_of_Clusters_Pixel" << tag.c_str() <<  id;
    local_histos.NumberOfClusterPixel = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 21, 0., 20.);
    histoName.str("");
    histoName << "Number_of_Clusters_Strip" << tag.c_str() <<  id;
    local_histos.NumberOfClusterStrip = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 21, 0., 20.);

    local_histos.NumberOfClusterPixel->SetFillColor(kBlue);
    local_histos.NumberOfClusterStrip->SetFillColor(kRed);

    histoName.str("");
    histoName << "Number_of_Clusters" << tag.c_str() <<  id;
    local_histos.NumberOfClustersSource = td.make<THStack>(histoName.str().c_str(), histoName.str().c_str());
    local_histos.NumberOfClustersSource->Add(local_histos.NumberOfClusterPixel);
    local_histos.NumberOfClustersSource->Add(local_histos.NumberOfClusterStrip);


    histoName.str("");
    histoName << "Number_of_Clusters_Link" << tag.c_str() <<  id;
    local_histos.NumberOfClustersLink = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 21, 0., 20.);


    histoName.str("");
    histoName << "ClusterSize" << tag.c_str() <<  id;
    local_histos.clusterSize = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 50, 0., 50.);
    histoName.str("");
    histoName << "ClusterSize_Pixel" << tag.c_str() <<  id;
    local_histos.clusterSizePixel = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 50, 0., 50.);
    histoName.str("");
    histoName << "ClusterSize_Strip" << tag.c_str() <<  id;
    local_histos.clusterSizeStrip = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 50, 0., 50.);

    local_histos.clusterSizePixel->SetFillColor(kBlue);
    local_histos.clusterSizeStrip->SetFillColor(kRed);

    histoName.str("");
    histoName << "Clusters_Size_Source" << tag.c_str() <<  id;
    local_histos.ClustersSizeSource = td.make<THStack>(histoName.str().c_str(), histoName.str().c_str());
    local_histos.ClustersSizeSource->Add(local_histos.clusterSizePixel);
    local_histos.ClustersSizeSource->Add(local_histos.clusterSizeStrip);

    histoName.str("");
    histoName << "ClusterSizeX" << tag.c_str() <<  id;
    local_histos.clusterSizeX = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 1000, 0., 0.);
    histoName.str("");
    histoName << "ClusterSizeY" << tag.c_str() <<  id;
    local_histos.clusterSizeY = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 1000, 0., 0.);

    histoName.str("");
    histoName << "ClusterShapeX" << tag.c_str() <<  id;
    local_histos.clusterShapeX = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 1000, 0., 0.);
    histoName.str("");
    histoName << "ClusterShapeY" << tag.c_str() <<  id;
    local_histos.clusterShapeY = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 1000, 0., 0.);

    histoName.str("");
    histoName << "LocalPositionXY" << tag.c_str() <<  id;
    local_histos.localPosXY = td.make<TH2F>(histoName.str().c_str(), histoName.str().c_str(), 2000, 0., 0., 2000, 0., 0.);

    histoName.str("");
    histoName << "GlobalPositionXY" << tag.c_str() <<  id;
    local_histos.globalPosXY = td.make<TH2F>(histoName.str().c_str(), histoName.str().c_str(), 2400, -120.0, 120.0, 2400, -120.0, 120.0);
    histoName.str("");
    histoName << "LocalPositionXY_Pixel" << tag.c_str() <<  id;
    local_histos.localPosXYPixel = td.make<TH2F>(histoName.str().c_str(), histoName.str().c_str(), 2000, 0., 0., 2000, 0., 0.);
    histoName.str("");
    histoName << "LocalPositionXY_Strip" << tag.c_str() <<  id;
    local_histos.localPosXYStrip = td.make<TH2F>(histoName.str().c_str(), histoName.str().c_str(), 2000, 0., 0., 2000, 0., 0.);

    histoName.str("");
    histoName << "Digi_type_" << tag.c_str() <<  id;
    local_histos.digiType = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 2, 0, 1);
    histoName.str("");
    histoName << "Digi_position_" << tag.c_str() <<  id;
    local_histos.digiPosition = td.make<TH2F>(histoName.str().c_str(), histoName.str().c_str(), 2000, 0., 0., 2000, 0., 0.);

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
