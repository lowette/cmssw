#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
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
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
#include "DataFormats/DetId/interface/DetId.h"

//RecHit Stuff
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"

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

class RunStepsDigiValidation_RecHits : public EDAnalyzer {

public:
    explicit RunStepsDigiValidation_RecHits(const ParameterSet&);
    ~RunStepsDigiValidation_RecHits();
    virtual void beginJob();
    virtual void analyze(const Event&, const EventSetup&);
    virtual void endJob();

private:
    InputTag src_;
    InputTag simG4_;

public:
    //void createLayerHistograms(unsigned int iLayer);
    //void createHistograms(unsigned int nLayer);
    //unsigned int getSimTrackId(Handle<DetSetVector<PixelDigiSimLink> >&, DetId& detId, unsigned int& channel);
    //int matchedSimTrack(Handle<SimTrackContainer>& SimTk, unsigned int simTrkId);
    void initializeVariables();
    //unsigned int getMaxPosition(vector<float>& charge_vec);
    //unsigned int getLayerNumber(const TrackerGeometry* tkgeom, unsigned int& detid);
    //unsigned int getLayerNumber(unsigned int& detid);
    //int isPrimary(const SimTrack& simTrk, Handle<PSimHitContainer>& simHits);
    //int isPrimary(const SimTrack& simTrk, const PSimHit& simHit);
    //void fillMatchedSimTrackHistos(DigiHistos& digiHistos, const SimTrack& simTk, int ptype, unsigned int layer);

};

RunStepsDigiValidation_RecHits::RunStepsDigiValidation_RecHits(const ParameterSet& iConfig) {
    src_ =  iConfig.getParameter<InputTag>("src");
    simG4_ = iConfig.getParameter<InputTag>("simG4");

    cout << "------------------------------------------------------------" << endl
         << "-- Running RunSteps DigiValidation v0.1" << endl
         << "------------------------------------------------------------" << endl;
}
RunStepsDigiValidation_RecHits::~RunStepsDigiValidation_RecHits() { }

void RunStepsDigiValidation_RecHits::beginJob() {
    //createHistograms(19);
}

void RunStepsDigiValidation_RecHits::analyze(const Event& iEvent, const EventSetup& iSetup) {
    //int run       = iEvent.id().run();
    //int event     = iEvent.id().event();
    //int lumiBlock = iEvent.luminosityBlock();
    //int bx        = iEvent.bunchCrossing();
    //int orbit     = iEvent.orbitNumber();

    // Get digis
    Handle< DetSetVector<PixelDigi> > pixelDigis;
    iEvent.getByLabel(src_, pixelDigis);

    // Get simlink data

    Handle< DetSetVector<PixelDigiSimLink> > pixelSimLinks;
    iEvent.getByLabel(src_,   pixelSimLinks);

    // Get event setup (to get global transformation)
    ESHandle<TrackerGeometry> geomHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get( geomHandle );
    const TrackerGeometry*  tkGeom = &(*geomHandle);

    // Get PSimHits
    Handle<PSimHitContainer> simHits;
    iEvent.getByLabel("g4SimHits","TrackerHitsPixelBarrelLowTof" ,simHits);

    Handle<SimTrackContainer> simTracks;
    iEvent.getByLabel("g4SimHits",simTracks);

    // SimVertex
    Handle<SimVertexContainer> simVertices;
    iEvent.getByLabel("g4SimHits", simVertices);

    // Get RecHits:
    Handle<SiPixelRecHitCollection> siPixelRecHit;
    iEvent.getByLabel("siPixelRecHits", siPixelRecHit);

    initializeVariables();

    vector<int> processTypes;

	//Loop over rechits:
	for(SiPixelRecHitCollection::const_iterator recHitItr = siPixelRecHit->begin(); recHitItr != siPixelRecHit->end(); ++ recHitItr){
	    const SiPixelRecHit* pixelHit = &(*recHit);

	    int RecHitPos = pixelHit.pos();
	    std::cout << " Position of RecHit : " << RecHitPos << std::endl;
	}		

    // Loop over Sim Tracks and Fill relevant histograms
//    int nTracks = 0;

  /*  for (SimTrackContainer::const_iterator simTrkItr = simTracks->begin(); simTrkItr != simTracks->end(); ++simTrkItr) {
        int type = isPrimary((*simTrkItr), simHits);
        processTypes.push_back(type);

        // remove neutrinos
        if (simTrkItr->charge() == 0) continue;

        nTracks++;
        float simTk_pt =  simTrkItr->momentum().pt();
        float simTk_eta = simTrkItr->momentum().eta();
        float simTk_phi = simTrkItr->momentum().phi();

        simTrackPt_->Fill(simTk_pt);
        simTrackEta_->Fill(simTk_eta);
        simTrackPhi_->Fill(simTk_phi);
        if (fabs(simTk_eta) < 1.0) simTrackPtBar_->Fill(simTk_pt);
        else if (fabs(simTk_eta) > 1.6) simTrackPtEC_->Fill(simTk_pt);

        if (type == 1) {
            simTrackPtP_->Fill(simTk_pt);
            if (fabs(simTk_eta) < 1.0) simTrackPtPBar_->Fill(simTk_pt);
            else if (fabs(simTk_eta) > 1.6) simTrackPtPEC_->Fill(simTk_pt);
            simTrackEtaP_->Fill(simTk_eta);
            simTrackPhiP_->Fill(simTk_phi);
        }
        else if (type == 0) {
            simTrackPtS_->Fill(simTk_pt);
            if (fabs(simTk_eta) < 1.0) simTrackPtSBar_->Fill(simTk_pt);
            else if (fabs(simTk_eta) > 1.6) simTrackPtSEC_->Fill(simTk_pt);
            simTrackEtaS_->Fill(simTk_eta);
            simTrackPhiS_->Fill(simTk_phi);
        }
    }

    nSimTracks_->Fill(nTracks++);
*/
    // Loop Over Digis and Fill Histograms
  /*  for (DetSetVector<PixelDigi>::const_iterator DSViter = pixelDigis->begin(); DSViter != pixelDigis->end(); ++DSViter) {
        // Detector id
        unsigned int rawid = DSViter->id;
        DetId detId(rawid);
        unsigned int layer = getLayerNumber(rawid);

        // Tracker geometry
        const GeomDetUnit* geomDetUnit = tkGeom->idToDetUnit(detId);
        const PixelGeomDetUnit* pixdet = dynamic_cast<const PixelGeomDetUnit*>(geomDetUnit);
        const PixelTopology & topol = pixdet->specificTopology();

        // Get the histogram for this layer
        map<unsigned int, DigiHistos>::iterator iPos = layerHistoMap.find(layer);

        if (iPos == layerHistoMap.end()) {
            createLayerHistograms(layer);
            iPos = layerHistoMap.find(layer);
        }

        int col_last = -1;
        int row_last = -1;
        int nDigiP = 0;
        int nDigiS = 0;
        int nDigiA = 0;

        vector<MyCluster> cluster_vec;

        MyCluster cluster;
        cluster.charge = 0.0;
        cluster.width = 0;
        cluster.trkType = false;
        cluster.trkPt = -999.0;
        cluster.trkEta = -999.0;
        cluster.strip_charges.clear();

        // Loop over the Digis
        for (DetSet<PixelDigi>::const_iterator di = DSViter->data.begin(); di != DSViter->data.end(); ++di) {
            int adc = di->adc();    // charge, modifued to unsiged short
            int col = di->column(); // column
            int row = di->row();    // row

            unsigned int channel = PixelChannelIdentifier::pixelToChannel(row, col);
            unsigned int simTkId = getSimTrackId(pixelSimLinks, detId, channel);

            // Get the Digi position
            MeasurementPoint mp(row + 0.5, col + 0.5 );

            if (geomDetUnit) {
                LocalPoint lPos = geomDetUnit->topology().localPosition(mp);
                GlobalPoint pdPos = geomDetUnit->surface().toGlobal(geomDetUnit->topology().localPosition(mp)) ;

                // Fill some histograms
                iPos->second.LocalPosition->Fill(lPos.x(), lPos.y());
                iPos->second.YposVsXpos->Fill(pdPos.y(), pdPos.x());
                iPos->second.RVsZpos->Fill(pdPos.z(), pdPos.perp());
                trackerLayout_->Fill(pdPos.z(), pdPos.perp());
                trackerLayoutXY_->Fill(pdPos.y(), pdPos.x());

                if (layer < 100) trackerLayoutXYBar_->Fill(pdPos.y(), pdPos.x());
                else trackerLayoutXYEC_->Fill(pdPos.y(), pdPos.x());

                // // Pixel module
                // if (topol.ncolumns() == 32) iPos->second.LocalPositionPixel->Fill(nClusterP + nClusterS);
                // // Strip module
                // else if (topol.ncolumns() == 2) iPos->second.LocalPositionStrip->Fill(nClusterP + nClusterS);
            }

            // Do the simtrack matching
            int iSimTrk = matchedSimTrack(simTracks, simTkId);
            int primaryTrk = -1;

            if (iSimTrk != -1) {
                primaryTrk = processTypes[iSimTrk];
                iPos->second.simTkIndx.insert(iSimTrk);
            }

            // Fill the charge of the Digi (not really useful as it is 0 or 255 for us)
            iPos->second.DigiCharge->Fill(adc);

            if (primaryTrk == 1) {
                iPos->second.DigiChargePrimary->Fill(adc);
                nDigiP++;
            }
            else if (primaryTrk == 0){
                iPos->second.DigiChargeSecondary->Fill(adc);
                nDigiS++;
            }

            nDigiA++;

            // Form a cluster of continuous Digis in one row
            if (col_last == -1) {
                cluster.charge = adc;
                cluster.width  = 1;
                cluster.trkType = primaryTrk;
                cluster.trkPt = (*simTracks)[iSimTrk].momentum().pt();
                cluster.trkEta = (*simTracks)[iSimTrk].momentum().eta();
                cluster.strip_charges.clear();
                cluster.strip_charges.push_back(adc);
            }
            else {
                if (abs(col - col_last) == 1 && row == row_last) {
                    cluster.charge += adc;
                    cluster.width++;
                    cluster.strip_charges.push_back(adc);
                }
                else {
                    cluster_vec.push_back(cluster);
                    cluster.charge = adc;
                    cluster.width = 1;
                    cluster.trkType = primaryTrk;
                    cluster.trkPt = (*simTracks)[iSimTrk].momentum().pt();
                    cluster.trkEta = (*simTracks)[iSimTrk].momentum().eta();
                    cluster.strip_charges.clear();
                    cluster.strip_charges.push_back(adc);
                }
            }

            col_last = col;
            row_last = row;
        }

        cluster_vec.push_back(cluster);

        // Fill some histograms about Digis
        iPos->second.NumberOfDigisPrimary->Fill(nDigiP);
        iPos->second.NumberOfDigisSecondary->Fill(nDigiS);
        // iPos->second.NumberOfDigis->Fill(nDigiP + nDigiS);
        iPos->second.NumberOfDigis->Fill(nDigiA);
        iPos->second.totNDigis += nDigiP + nDigiS;

        // Pixel module
        if (topol.ncolumns() == 32) iPos->second.NumberOfDigisPixel->Fill(nDigiA);
        // Strip module
        else if (topol.ncolumns() == 2) iPos->second.NumberOfDigisStrip->Fill(nDigiA);

        int nClusterP = 0;
        int nClusterS = 0;
        int nClusterA = 0;

        for (vector<MyCluster>::iterator ic = cluster_vec.begin(); ic != cluster_vec.end(); ++ic) {
            // Cluster parameters
            float cl_charge = ic->charge;
            int cl_width = ic->width;
            int cl_type = ic->trkType;
            float trk_pt = ic->trkPt;
            float trk_eta = ic->trkEta;
            vector<float> str_charges = ic->strip_charges;
            unsigned int max_pos = getMaxPosition(str_charges);

            // Compute the cluster shape
            if (max_pos != 999) {
                for (unsigned int ival = 0; ival < str_charges.size(); ++ival) {
                    int pos = ival - max_pos;
                    iPos->second.ClusterShape->Fill(pos, str_charges[ival]);
                    if (cl_type == 1) iPos->second.ClusterShapePrimary->Fill(pos, str_charges[ival]);
                    else if (cl_type == 0) iPos->second.ClusterShapeSecondart->Fill(pos, str_charges[ival]);
                }
            }

            // Fill some plots
            iPos->second.ClusterCharge->Fill(cl_charge);
            iPos->second.ClusterWidth->Fill(cl_width);
            iPos->second.ClusterWidthVsSimTrkPt->Fill(trk_pt, cl_width);
            iPos->second.ClusterWidthVsSimTrkEta->Fill(trk_eta, cl_width);

            if (cl_type == 1) {
                iPos->second.ClusterChargePrimary->Fill(cl_charge);
                iPos->second.ClusterWidthPrimary->Fill(cl_width);
                iPos->second.ClusterWidthVsSimTrkPtPrimary->Fill(trk_pt, cl_width);
                iPos->second.ClusterWidthVsSimTrkEtaPrimary->Fill(trk_eta, cl_width);
                nClusterP++;
            }
            else if (cl_type == 0) {
                iPos->second.ClusterChargeSecondary->Fill(cl_charge);
                iPos->second.ClusterWidthSecondary->Fill(cl_width);
                iPos->second.ClusterWidthVsSimTrkPtSecondary->Fill(trk_pt, cl_width);
                iPos->second.ClusterWidthVsSimTrkEtaSecondary->Fill(trk_eta, cl_width);
                nClusterS++;
            }

            nClusterA++;
        }

        // Fill some plots
        iPos->second.NumberOfClustersPrimary->Fill(nClusterP);
        iPos->second.NumberOfClustersSecondary->Fill(nClusterS);
        // iPos->second.NumberOfClusters->Fill(nClusterP + nClusterS);
        iPos->second.NumberOfClusters->Fill(nClusterA);
        iPos->second.totNClusters += nClusterP + nClusterS;

        // Pixel module
        if (topol.ncolumns() == 32) iPos->second.NumberOfClustersPixel->Fill(nClusterA);
        // Strip module
        else if (topol.ncolumns() == 2) iPos->second.NumberOfClustersStrip->Fill(nClusterA);
    }*/

    // Fill Layer Level Histograms
  /*  for (map<unsigned int, DigiHistos>::iterator iPos  = layerHistoMap.begin(); iPos != layerHistoMap.end(); ++iPos) {
        DigiHistos local_histos = iPos->second;
        local_histos.TotalNumberOfClusters->Fill(local_histos.totNClusters);
        local_histos.TotalNumberOfDigis->Fill(local_histos.totNDigis);
        local_histos.NumberOfSimHits->Fill(local_histos.totSimHits);

        local_histos.NumberOfMatchedSimHits->Fill(local_histos.totMatchedSimHits);
        local_histos.NumberOfMatchedSimHitsPrimary->Fill(local_histos.totMatchedSimHitsPrimary);
        local_histos.NumberOfMatchedSimHitsSecondary->Fill(local_histos.totMatchedSimHitsSecondary);

        if (local_histos.totSimHits) {
            float eff;
            eff = local_histos.totMatchedSimHits * 1. / local_histos.totSimHits;
            local_histos.DigiEfficiency->Fill(eff);
            eff = local_histos.totMatchedSimHitsPrimary * 1. / local_histos.totSimHitsPrimary;
            local_histos.DigiEfficiencyPrimary->Fill(eff);
            eff = local_histos.totMatchedSimHitsSecondary * 1. / local_histos.totSimHitsSecondary;
            local_histos.DigiEfficiencySecondary->Fill(eff);
        }

        for (set<int>::iterator ii = local_histos.simTkIndx.begin(); ii != local_histos.simTkIndx.end(); ++ii) {
            unsigned int index = (*ii);
            int pid = processTypes[index];
            fillMatchedSimTrackHistos(local_histos, (*simTracks.product())[index], pid, iPos->first);
        }
    }*/
}

void RunStepsDigiValidation_RecHits::endJob() { }

/*void RunStepsDigiValidation::createLayerHistograms(unsigned int ival) {
    ostringstream fname1, fname2;

    Service<TFileService> fs;
    fs->file().cd("/");

    string tag;
    unsigned int id;
    if (ival < 100) {
        id = ival;
        fname1 << "Barrel";
        fname2 << "Layer_" << id;
        tag = "_layer_";
    }
    else {
        int side = ival / 100;
        id = ival - side * 100;
        cout << " Creating histograms for Disc " << id << " with " << ival << endl;
        fname1 << "EndCap_Side_" << side;
        fname2 << "Disc_" << id;
        tag = "_disc_";
    }

    TFileDirectory td1 = fs->mkdir(fname1.str().c_str());
    TFileDirectory td = td1.mkdir(fname2.str().c_str());

    DigiHistos local_histos;
    ostringstream htit1;
    htit1 << "NumberOfDigis" << tag.c_str() << id;
    local_histos.NumberOfDigis = td.make<TH1F>(htit1.str().c_str(), htit1.str().c_str(), 51, -0.5, 50.5);
    htit1.str("");
    htit1 << "NumberOfDigisPrimary" << tag.c_str() << id;
    local_histos.NumberOfDigisPrimary = td.make<TH1F>(htit1.str().c_str(), htit1.str().c_str(), 51, -0.5, 50.5);
    htit1.str("");
    htit1 << "NumberOfDigisSecondary" << tag.c_str() << id;
    local_histos.NumberOfDigisSecondary = td.make<TH1F>(htit1.str().c_str(), htit1.str().c_str(), 51, -0.5, 50.5);
    htit1.str("");
    htit1 << "NumberOfDigis_Pixel" << tag.c_str() << id;
    local_histos.NumberOfDigisPixel = td.make<TH1F>(htit1.str().c_str(), htit1.str().c_str(), 51, -0.5, 50.5);
    htit1.str("");
    htit1 << "NumberOfDigis_Strip" << tag.c_str() << id;
    local_histos.NumberOfDigisStrip = td.make<TH1F>(htit1.str().c_str(), htit1.str().c_str(), 51, -0.5, 50.5);


    ostringstream htit2;
    htit2 << "DigiCharge" << tag.c_str() << id;
    local_histos.DigiCharge = td.make<TH1F>(htit2.str().c_str(), htit2.str().c_str(), 261, -0.5, 260.5);
    htit2.str("");
    htit2 << "DigiChargePrimary" << tag.c_str() << id;
    local_histos.DigiChargePrimary = td.make<TH1F>(htit2.str().c_str(), htit2.str().c_str(), 261, -0.5, 260.5);
    htit2.str("");
    htit2 << "DigiChargeSecondary" << tag.c_str() << id;
    local_histos.DigiChargeSecondary = td.make<TH1F>(htit2.str().c_str(), htit2.str().c_str(), 261, -0.5, 260.5);

    ostringstream htit3;
    htit3 << "NumberOfClusters" << tag.c_str() << id;
    local_histos.NumberOfClusters = td.make<TH1F>(htit3.str().c_str(), htit3.str().c_str(), 51, -0.5, 50.5);
    htit3.str("");
    htit3 << "NumberOfClustersPrimary" << tag.c_str() << id;
    local_histos.NumberOfClustersPrimary = td.make<TH1F>(htit3.str().c_str(), htit3.str().c_str(), 51, -0.5, 50.5);
    htit3.str("");
    htit3 << "NumberOfClustersSecondary" << tag.c_str() << id;
    local_histos.NumberOfClustersSecondary = td.make<TH1F>(htit3.str().c_str(), htit3.str().c_str(), 51, -0.5, 50.5);
    htit3.str("");
    htit3 << "NumberOfClusters_Pixel" << tag.c_str() << id;
    local_histos.NumberOfClustersPixel = td.make<TH1F>(htit3.str().c_str(), htit3.str().c_str(), 51, -0.5, 50.5);
    htit3.str("");
    htit3 << "NumberOfClusters_Strip" << tag.c_str() << id;
    local_histos.NumberOfClustersStrip = td.make<TH1F>(htit3.str().c_str(), htit3.str().c_str(), 51, -0.5, 50.5);

    ostringstream htit4;
    htit4 << "ClusterCharge" << tag.c_str() << id;
    local_histos.ClusterCharge = td.make<TH1F>(htit4.str().c_str(), htit4.str().c_str(), 1041, -0.5, 1040.5);
    htit4.str("");
    htit4 << "ClusterChargePrimary" << tag.c_str() << id;
    local_histos.ClusterChargePrimary = td.make<TH1F>(htit4.str().c_str(), htit4.str().c_str(), 1041, -0.5,1040.5);
    htit4.str("");
    htit4 << "ClusterChargeSecondary" << tag.c_str() << id;
    local_histos.ClusterChargeSecondary = td.make<TH1F>(htit4.str().c_str(), htit4.str().c_str(), 1041, -0.5, 1040.5);

    ostringstream htit5;
    htit5 << "ClusterWidth" << tag.c_str() << id;
    local_histos.ClusterWidth = td.make<TH1F>(htit5.str().c_str(), htit5.str().c_str(), 16, -0.5, 15.5);
    htit5.str("");
    htit5 << "ClusterWidthPrimary" << tag.c_str() << id;
    local_histos.ClusterWidthPrimary = td.make<TH1F>(htit5.str().c_str(), htit5.str().c_str(), 16, -0.5, 15.5);
    htit5.str("");
    htit5 << "ClusterWidthSecondary" << tag.c_str() << id;
    local_histos.ClusterWidthSecondary = td.make<TH1F>(htit5.str().c_str(), htit5.str().c_str(), 16, -0.5, 15.5);

    ostringstream htit6;
    htit6 << "TotalNumberOfDigis" << tag.c_str() << id;
    local_histos.TotalNumberOfDigis = td.make<TH1F>(htit6.str().c_str(), htit6.str().c_str(), 100, -0.5, 100.5);

    ostringstream htit7;
    htit7 << "TotalNumberOfClusters" << tag.c_str() << id;
    local_histos.TotalNumberOfClusters = td.make<TH1F>(htit7.str().c_str(), htit7.str().c_str(), 100, -0.5, 100.5);

    ostringstream htit8;
    htit8 << "ClusterShape" << tag.c_str() << id;
    local_histos.ClusterShape = td.make<TH1F>(htit8.str().c_str(), htit8.str().c_str(), 21, -20.5, 20.5);
    htit8.str("");
    htit8<< "ClusterShapePrimary" << tag.c_str() << id;
    local_histos.ClusterShapePrimary = td.make<TH1F>(htit8.str().c_str(), htit8.str().c_str(), 21, -20.5, 20.5);
    htit8.str("");
    htit8 << "ClusterShapeSecondart" << tag.c_str() << id;
    local_histos.ClusterShapeSecondart = td.make<TH1F>(htit8.str().c_str(), htit8.str().c_str(), 21, -20.5, 20.5);

    ostringstream htit9;
    htit9 << "NumberOfSimHits" << tag.c_str() << id;
    local_histos.NumberOfSimHits = td.make<TH1F>(htit9.str().c_str(), htit9.str().c_str(), 201, -0.5, 200.5);
    htit9.str("");
    htit9 << "NumberOfMatchedSimHits" << tag.c_str() << id;
    local_histos.NumberOfMatchedSimHits = td.make<TH1F>(htit9.str().c_str(), htit9.str().c_str(), 201, -0.5, 200.5);
    htit9.str("");
    htit9 << "NumberOfMatchedSimHitsPrimary" << tag.c_str() << id;
    local_histos.NumberOfMatchedSimHitsPrimary = td.make<TH1F>(htit9.str().c_str(), htit9.str().c_str(), 201, -0.5, 200.5);
    htit9.str("");
    htit9 << "NumberOfMatchedSimHitsSecondary" << tag.c_str() << id;
    local_histos.NumberOfMatchedSimHitsSecondary = td.make<TH1F>(htit9.str().c_str(), htit9.str().c_str(), 201, -0.5, 200.5);

    ostringstream htit10;
    htit10 << "DigiEfficiency" << tag.c_str() << id;
    local_histos.DigiEfficiency = td.make<TH1F>(htit10.str().c_str(), htit10.str().c_str(), 55, -0.05, 1.05);
    htit10.str("");
    htit10 << "DigiEfficiencyPrimary" << tag.c_str() << id;
    local_histos.DigiEfficiencyPrimary = td.make<TH1F>(htit10.str().c_str(), htit10.str().c_str(), 55, -0.05, 1.05);
    htit10.str("");
    htit10 << "DigiEfficiencySecondary" << tag.c_str() << id;
    local_histos.DigiEfficiencySecondary = td.make<TH1F>(htit10.str().c_str(), htit10.str().c_str(), 55, -0.05, 1.05);

    ostringstream htit11;
    htit11 << "YposVsXpos" << tag.c_str() << id;
    local_histos.YposVsXpos = td.make<TH2F>(htit11.str().c_str(), htit11.str().c_str(), 240, -120.0, 120.0, 240, -120.0, 120.0);

    ostringstream htit12;
    htit12 << "RVsZpos" << tag.c_str() << id;
    local_histos.RVsZpos = td.make<TH2F>(htit12.str().c_str(), htit12.str().c_str(), 600, -300.0, 300.0, 120, 0.0, 120.0);

    ostringstream htit13;
    htit13 << "DigiChargeMatched" << tag.c_str() << id;
    local_histos.DigiChargeMatched = td.make<TH1F>(htit13.str().c_str(), htit13.str().c_str(), 261, -0.5, 260.5);

    ostringstream htit14;
    htit14 << "ClusterWidthVsSimTrkPt" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkPt = td.make<TProfile>(htit14.str().c_str(),htit14.str().c_str(),56, -0.5, 55.5,-0.5,15.5);
    htit14.str("");
    htit14 << "ClusterWidthVsSimTrkPtPrimary" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkPtPrimary = td.make<TProfile>(htit14.str().c_str(),htit14.str().c_str(),56, -0.5, 55.5,-0.5,15.5);
    htit14.str("");
    htit14 << "ClusterWidthVsSimTkrPtS" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkPtSecondary = td.make<TProfile>(htit14.str().c_str(),htit14.str().c_str(),56, -0.5, 55.5,-0.5,15.5);

    ostringstream htit15;
    htit15 << "ClusterWidthVsSimTrkEta" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkEta = td.make<TProfile>(htit15.str().c_str(),htit15.str().c_str(),50, -2.5, 2.5,-0.5,15.5);
    htit15.str("");
    htit15 << "ClusterWidthVsSimTrkEtaPrimary" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkEtaPrimary = td.make<TProfile>(htit15.str().c_str(),htit15.str().c_str(),50, -2.5, 2.5,-0.5,15.5);
    htit15.str("");
    htit15 << "ClusterWidthVsSimTkrEtaS" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkEtaSecondary = td.make<TProfile>(htit15.str().c_str(),htit15.str().c_str(),50, -2.5, 2.5,-0.5,15.5);

    ostringstream htit16;
    htit16 << "MatchedSimTrackPt" << tag.c_str() << id;
    local_histos.matchedSimTrackPt_  = td.make<TH1F>(htit16.str().c_str(),htit16.str().c_str(),101,-0.5,100.5);
    htit16.str("");
    htit16 << "MatchedSimTrackPtP" << tag.c_str() << id;
    local_histos.matchedSimTrackPtPrimary_  = td.make<TH1F>(htit16.str().c_str(),htit16.str().c_str(),101,-0.5,100.5);
    htit16.str("");
    htit16 << "MatchedSimTrackPtS" << tag.c_str() << id;
    local_histos.matchedSimTrackPtSecondary_  = td.make<TH1F>(htit16.str().c_str(),htit16.str().c_str(),101,-0.5,100.5);

    ostringstream htit17;
    htit17 << "MatchedSimTrackEta" << tag.c_str() << id;
    local_histos.matchedSimTrackEta_  = td.make<TH1F>(htit17.str().c_str(),  htit17.str().c_str(), 50, -2.5, 2.5);
    htit17.str("");
    htit17 << "MatchedSimTrackEtaP" << tag.c_str() << id;
    local_histos.matchedSimTrackEtaPrimary_  = td.make<TH1F>(htit17.str().c_str(),  htit17.str().c_str(), 50, -2.5, 2.5);
    htit17.str("");
    htit17 << "MatchedSimTrackEtaS" << tag.c_str() << id;
    local_histos.matchedSimTrackEtaSecondary_  = td.make<TH1F>(htit17.str().c_str(),  htit17.str().c_str(), 50, -2.5, 2.5);

    ostringstream htit18;
    htit18 << "MatchedSimTrackPhi" << tag.c_str() << id;
    local_histos.matchedSimTrackPhi_  = td.make<TH1F>(htit18.str().c_str(),  htit18.str().c_str(), 160, -3.2, 3.2);
    htit18.str("");
    htit18 << "MatchedSimTrackPhiP" << tag.c_str() << id;
    local_histos.matchedSimTrackPhiPrimary_  = td.make<TH1F>(htit18.str().c_str(),  htit18.str().c_str(), 160, -3.2, 3.2);
    htit18.str("");
    htit18 << "MatchedSimTrackPhiS" << tag.c_str() << id;
    local_histos.matchedSimTrackPhiSecondary_  = td.make<TH1F>(htit18.str().c_str(),  htit18.str().c_str(), 160, -3.2, 3.2);

    ostringstream htit19;
    htit19 << "LocalPosition" << tag.c_str() << id;
    local_histos.LocalPosition = td.make<TH2F>(htit19.str().c_str(),htit19.str().c_str(),10000, -5, 5 , 10000, -5 ,5);

    layerHistoMap.insert( make_pair(ival, local_histos));

    fs->file().cd("/");

    layerHistoMap[ival].totNDigis = 0;
    layerHistoMap[ival].totNClusters = 0;
    layerHistoMap[ival].totSimHits = 0;
    layerHistoMap[ival].totSimHitsPrimary = 0;
    layerHistoMap[ival].totSimHitsSecondary = 0;
    layerHistoMap[ival].totMatchedSimHits = 0;
    layerHistoMap[ival].totMatchedSimHitsPrimary = 0;
    layerHistoMap[ival].totMatchedSimHitsSecondary = 0;

    layerHistoMap[ival].simTkIndx.clear();

}/*
/*
void RunStepsDigiValidation::createHistograms(unsigned int nLayer) {
    Service<TFileService> fs;
    fs->file().cd("/");
    TFileDirectory td = fs->mkdir("Common");

    nSimTracks_  = td.make<TH1F>("nSimTracks", "Number of Sim Tracks" , 201, -0.5, 200.5);
    simTrackPt_  = td.make<TH1F>("SimTrackPt", "Pt of Sim Tracks", 101, -0.5, 100.5);
    simTrackPtBar_  = td.make<TH1F>("SimTrackPtBar", "Pt of Sim Tracks( eta < 1.0)", 101, -0.5, 100.5);
    simTrackPtEC_  = td.make<TH1F>("SimTrackPtEC", "Pt of Sim Tracks( eta > 1.6)", 101, -0.5, 100.5);
    simTrackEta_ =  td.make<TH1F>("SimTrackEta", "Eta of Sim Tracks", 50, -2.5, 2.5);
    simTrackPhi_ =  td.make<TH1F>("SimTrackPhi", "Phi of Sim Tracks", 160, -3.2, 3.2);

    simTrackPtP_  = td.make<TH1F>("SimTrackPtP", "Pt of Primary Sim Tracks", 101, -0.5, 100.5);
    simTrackPtPBar_  = td.make<TH1F>("SimTrackPtPBar", "Pt of Primary Sim Tracks( eta < 1.0)", 101, -0.5, 100.5);
    simTrackPtPEC_  = td.make<TH1F>("SimTrackPtPEC", "Pt of Primary Sim Tracks( eta > 1.6)", 101, -0.5, 100.5);
    simTrackEtaP_ =  td.make<TH1F>("SimTrackEtaP", "Eta of Primary Sim Tracks", 50, -2.5, 2.5);
    simTrackPhiP_ =  td.make<TH1F>("SimTrackPhiP", "Phi of Primary Sim Tracks", 160, -3.2, 3.2);

    simTrackPtS_  = td.make<TH1F>("SimTrackPtS", "Pt of Secondary Sim Tracks", 101, -0.5, 100.5);
    simTrackPtSBar_  = td.make<TH1F>("SimTrackPtSBar", "Pt of Secondary Sim Tracks( eta < 1.0)", 101, -0.5, 100.5);
    simTrackPtSEC_  = td.make<TH1F>("SimTrackPtSEC", "Pt of Secondary Sim Tracks( eta > 1.6)", 101, -0.5, 100.5);
    simTrackEtaS_ =  td.make<TH1F>("SimTrackEtaS", "Eta of Secondary Sim Tracks", 50, -2.5, 2.5);
    simTrackPhiS_ =  td.make<TH1F>("SimTrackPhiS", "Phi of Secondary Sim Tracks", 160, -3.2, 3.2);

    trackerLayout_ = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 1200, 0.0, 120.0);
    trackerLayoutXY_ = td.make<TH2F>("XVsY", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
    trackerLayoutXYBar_ = td.make<TH2F>("XVsYBar", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
    trackerLayoutXYEC_ = td.make<TH2F>("XVsYEC", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
}

unsigned int RunStepsDigiValidation::getSimTrackId(Handle<DetSetVector<PixelDigiSimLink> >& pixelSimLinks, DetId& detId, unsigned int& channel) {
    DetSetVector<PixelDigiSimLink>::const_iterator isearch = pixelSimLinks->find(detId);

    unsigned int simTrkId(0);
    if (isearch == pixelSimLinks->end()) return simTrkId;

    DetSet<PixelDigiSimLink> link_detset = (*pixelSimLinks)[detId];
    int iSimLink = 0;
    for (DetSet<PixelDigiSimLink>::const_iterator it = link_detset.data.begin(); it != link_detset.data.end(); it++,iSimLink++) {
        if (channel == it->channel()) {
            simTrkId = it->SimTrackId();
            break;
        }
    }
    return simTrkId;
}

int RunStepsDigiValidation::matchedSimTrack(Handle<SimTrackContainer>& SimTk, unsigned int simTrkId) {
    SimTrackContainer sim_tracks = (*SimTk.product());
    for(unsigned int it = 0; it < sim_tracks.size(); it++) {
        if (sim_tracks[it].trackId() == simTrkId) return it;
    }
    return -1;
}

void RunStepsDigiValidation::initializeVariables() {
    for (map<unsigned int, DigiHistos>::iterator iPos = layerHistoMap.begin(); iPos != layerHistoMap.end(); iPos++) {
        iPos->second.totNDigis    = 0;
        iPos->second.totNClusters = 0;
        iPos->second.totSimHits    = 0;
        iPos->second.totSimHitsPrimary   = 0;
        iPos->second.totSimHitsSecondary   = 0;
        iPos->second.totMatchedSimHits = 0;
        iPos->second.totMatchedSimHitsPrimary = 0;
        iPos->second.totMatchedSimHitsSecondary = 0;
        iPos->second.simTkIndx.clear();
    }
}

unsigned int RunStepsDigiValidation::getLayerNumber(const TrackerGeometry* tkgeom,unsigned int& detid) {
    unsigned int layer = 999;
    DetId theDetId(detid);
    if (theDetId.subdetId() != 1) cout << ">>> Method1 : Det id " << theDetId.det() << " Subdet Id " << theDetId.subdetId() << endl;
    const PixelGeomDetUnit * theGeomDet = dynamic_cast<const PixelGeomDetUnit*> ( tkgeom->idToDet(theDetId) );

    const GeomDetUnit* it = tkgeom->idToDetUnit(DetId(theDetId));
    if (!it) cout << ">>> rawdetid " << detid << " GeomDetUnit " << it << " PixelGeomDetUnit " << theGeomDet << " DetId " << theDetId.det() << " Subdet Id " << theDetId.subdetId() << endl;

    if (it && it->type().isTracker()) {
        if (it->type().isBarrel()) {
            PXBDetId pb_detId = PXBDetId(detid);
            layer = pb_detId.layer();
        }
        else if (it->type().isEndcap()) {
            cout << " IAM HERE >>>>>>>>>>>>>>>>>>>>>>>>>>> " << endl;
            PXFDetId pf_detId = PXFDetId(detid);
            layer = 100 * pf_detId.side() + pf_detId.disk();
        }
    }
    return layer;
}

unsigned int RunStepsDigiValidation::getLayerNumber(unsigned int& detid) {
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
        else cout << ">>> Invalid subdetId() = " << theDetId.subdetId() << endl;
    }
    return layer;
}

unsigned int RunStepsDigiValidation::getMaxPosition(vector<float>& charge_vec) {
    unsigned int ipos = 999;
    float max_val = 0.0;
    for (unsigned int ival = 0; ival < charge_vec.size(); ++ival) {
        if (charge_vec[ival] > max_val) {
            max_val = charge_vec[ival];
            ipos = ival;
        }
    }
    return ipos;
}

int RunStepsDigiValidation::isPrimary(const SimTrack& simTrk, Handle<PSimHitContainer>& simHits) {
    int result = -1;
    unsigned int trkId = simTrk.trackId();
    if (trkId > 0) {
        int vtxIndx = simTrk.vertIndex();
        for (PSimHitContainer::const_iterator iHit = simHits->begin(); iHit != simHits->end(); ++iHit) {
            if (trkId == iHit->trackId()) {
                int ptype = iHit->processType();
                if ((vtxIndx == 0) && (ptype == 2 || ptype == 7 || ptype == 9 || ptype == 11 || ptype == 15)) result = 1;
                else result = 0;
                break;
            }
        }
    }
    return result;
}

int RunStepsDigiValidation::isPrimary(const SimTrack& simTrk, const PSimHit& simHit) {
    int result = -1;
    unsigned int trkId = simTrk.trackId();
    if (trkId > 0) {
        int vtxIndx = simTrk.vertIndex();
        int ptype = simHit.processType();
        if ((vtxIndx == 0) && (ptype == 2 || ptype == 7 || ptype == 9 || ptype == 11 || ptype == 15)) result = 1;
        else result = 0;
    }
    return result;
}

void RunStepsDigiValidation::fillMatchedSimTrackHistos(DigiHistos& digiHistos, const SimTrack& simTk, int ptype, unsigned int layer){
    float pt =  simTk.momentum().pt();
    float eta = simTk.momentum().eta();
    float phi = simTk.momentum().phi();

    if (layer < 100 && fabs(eta) < 1.0) {
        digiHistos.matchedSimTrackPt_->Fill(pt);
        if (ptype == 1) digiHistos.matchedSimTrackPtPrimary_->Fill(pt);
        else if (ptype == 0) digiHistos.matchedSimTrackPtSecondary_->Fill(pt);
    }
    else if (layer > 100 && fabs(eta) > 1.6) {
        digiHistos.matchedSimTrackPt_->Fill(pt);
        if (ptype == 1) digiHistos.matchedSimTrackPtPrimary_->Fill(pt);
        else if (ptype == 0) digiHistos.matchedSimTrackPtSecondary_->Fill(pt);
    }
    if (ptype == 1) {
        digiHistos.matchedSimTrackEtaPrimary_->Fill(eta);
        digiHistos.matchedSimTrackPhiPrimary_->Fill(phi);
    }
    else if (ptype == 0) {
        digiHistos.matchedSimTrackEtaSecondary_->Fill(eta);
        digiHistos.matchedSimTrackPhiSecondary_->Fill(phi);
    }
    digiHistos.matchedSimTrackEta_->Fill(eta);
    digiHistos.matchedSimTrackPhi_->Fill(phi);
}*/

DEFINE_FWK_MODULE(RunStepsDigiValidation_RecHits);
