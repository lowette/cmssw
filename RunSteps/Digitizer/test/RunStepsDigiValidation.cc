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
#include <THStack.h>
#include <TProfile.h>

int verbose=2;
const int nTypes=18;

using namespace std;
using namespace edm;

class RunStepsDigiValidation : public EDAnalyzer {

  typedef vector< pair< PSimHit , vector<PixelDigi> > > V_HIT_DIGI;
  typedef map< int , V_HIT_DIGI > M_TRK_HIT_DIGI;

public:
    explicit RunStepsDigiValidation(const ParameterSet&);
    ~RunStepsDigiValidation();
    virtual void beginJob();
    virtual void analyze(const Event&, const EventSetup&);
    virtual void endJob();

private:
    TH1F* nSimTracks_;

    TH1F* simTrackPt_;
    TH1F* simTrackPtBar_;
    TH1F* simTrackPtEC_;
    TH1F* simTrackEta_;
    TH1F* simTrackPhi_;

    TH1F* simTrackPtP_;
    TH1F* simTrackPtPBar_;
    TH1F* simTrackPtPEC_;
    TH1F* simTrackEtaP_;
    TH1F* simTrackPhiP_;

    TH1F* simTrackPtS_;
    TH1F* simTrackPtSBar_;
    TH1F* simTrackPtSEC_;
    TH1F* simTrackEtaS_;
    TH1F* simTrackPhiS_;

    TH2F* trackerLayout_;
    TH2F* trackerLayoutXY_;
    TH2F* trackerLayoutXYBar_;
    TH2F* trackerLayoutXYEC_;

    struct DigiHistos {
        // THStack* NumberOfDigisType;
        // THStack* NumberOfDigisSource;
        TH1F* NumberOfDigisPrimary; // Secondary
        TH1F* NumberOfDigisSecondary; // Primary
        // TH1F* NumberOfDigisPixel; // Pixel
        // TH1F* NumberOfDigisStrip; // Strip

      // Truth Matching
        TH1F* NumberOfMatchedHits[nTypes];
        TH1F* NumberOfMatchedDigis[nTypes];
        TH1F* hEfficiency[nTypes];
        TH1F* h_dx_Truth;
        TH1F* h_dy_Truth;

        TH1F* DigiCharge;
        // TH1F* DigiChargePrimary;
        // TH1F* DigiChargeSecondary;

        TH1F* DigiChargeMatched;

        TH1F* NumberOfClusters;
        // TH1F* NumberOfClustersPrimary; // Primary
        // TH1F* NumberOfClustersSecondary; // Secondary
        // TH1F* NumberOfClustersPixel; // Pixel
        // TH1F* NumberOfClustersStrip; // Strip

        TH1F* ClusterCharge;
        // TH1F* ClusterChargePrimary;
        // TH1F* ClusterChargeSecondary;

        TH1F* ClusterWidth;
        // TH1F* ClusterWidthPrimary;
        // TH1F* ClusterWidthSecondary;

        TH1F* TotalNumberOfDigis;
        TH1F* TotalNumberOfClusters;

        TH1F* ClusterShape;
        TH1F* ClusterShapePrimary;
        TH1F* ClusterShapeSecondart;

        TH1F* NumberOfSimHits;
        TH1F* NumberOfMatchedSimHits;
        TH1F* NumberOfMatchedSimHitsPrimary;
        TH1F* NumberOfMatchedSimHitsSecondary;

        TH1F* DigiEfficiency;
        TH1F* DigiEfficiencyPrimary;
        TH1F* DigiEfficiencySecondary;

        TH2F* YposVsXpos;
        TH2F* RVsZpos;
        TH2F* LocalPosition;

        TProfile* ClusterWidthVsSimTrkPt;
        TProfile* ClusterWidthVsSimTrkPtPrimary;
        TProfile* ClusterWidthVsSimTrkPtSecondary;

        TProfile* ClusterWidthVsSimTrkEta;
        TProfile* ClusterWidthVsSimTrkEtaPrimary;
        TProfile* ClusterWidthVsSimTrkEtaSecondary;

        TH1F* matchedSimTrackPt_;
        TH1F* matchedSimTrackEta_;
        TH1F* matchedSimTrackPhi_;

        TH1F* matchedSimTrackPtPrimary_;
        TH1F* matchedSimTrackEtaPrimary_;
        TH1F* matchedSimTrackPhiPrimary_;

        TH1F* matchedSimTrackPtSecondary_;
        TH1F* matchedSimTrackEtaSecondary_;
        TH1F* matchedSimTrackPhiSecondary_;

        int totNDigis;
        int totNClusters;

        int totSimHits;
        int totSimHitsPrimary;
        int totSimHitsSecondary;

        int totMatchedSimHits;
        int totMatchedSimHitsPrimary;
        int totMatchedSimHitsSecondary;

        set<int> simTkIndx;
    };

    struct MyCluster {
        float charge;
        int width;
        bool trkType;
        float trkPt;
        float trkEta;
        vector<float> strip_charges;
    };

    map<unsigned int, DigiHistos> layerHistoMap;

    InputTag src_;
    InputTag simG4_;

public:
    void createLayerHistograms(unsigned int iLayer);
    void createHistograms(unsigned int nLayer);
    unsigned int getSimTrackId(Handle<DetSetVector<PixelDigiSimLink> >&, DetId& detId, unsigned int& channel);
    int matchedSimTrack(Handle<SimTrackContainer>& SimTk, unsigned int simTrkId);
    void initializeVariables();
    unsigned int getMaxPosition(vector<float>& charge_vec);
    unsigned int getLayerNumber(const TrackerGeometry* tkgeom, unsigned int& detid);
    unsigned int getLayerNumber(unsigned int& detid);
    int isPrimary(const SimTrack& simTrk, Handle<PSimHitContainer>& simHits);
    int isPrimary(const SimTrack& simTrk, const PSimHit& simHit);
    void fillMatchedSimTrackHistos(DigiHistos& digiHistos, const SimTrack& simTk, int ptype, unsigned int layer);

};

RunStepsDigiValidation::RunStepsDigiValidation(const ParameterSet& iConfig) {
    src_ =  iConfig.getParameter<InputTag>("src");
    simG4_ = iConfig.getParameter<InputTag>("simG4");

    cout << "------------------------------------------------------------" << endl
         << "-- Running RunSteps DigiValidation v0.1" << endl
         << "------------------------------------------------------------" << endl;
}
RunStepsDigiValidation::~RunStepsDigiValidation() { }

void RunStepsDigiValidation::beginJob() {
    createHistograms(19);
}

void RunStepsDigiValidation::analyze(const Event& iEvent, const EventSetup& iSetup) {
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
    if(verbose>1) cout << "- Getting SimHits" << endl << "-- TrackerHitsPixelBarrelLowTof" << endl;
    Handle<PSimHitContainer> simHits_B;
    iEvent.getByLabel("g4SimHits","TrackerHitsPixelBarrelLowTof" ,simHits_B);

    if(verbose>1) cout << "-- TrackerHitsPixelEndcapLowTof" << endl;
    Handle<PSimHitContainer> simHits_E;
    iEvent.getByLabel("g4SimHits","TrackerHitsPixelEndcapLowTof" ,simHits_E);

    // SimTracks
    if(verbose>1) cout << "- Getting SimTracks" << endl;
    Handle<SimTrackContainer> simTracks;
    iEvent.getByLabel("g4SimHits",simTracks);

    // SimVertices
    if(verbose>1) cout << "- Getting SimVertices" << endl;
    Handle<SimVertexContainer> simVertices;
    iEvent.getByLabel("g4SimHits", simVertices);

    initializeVariables();

    vector<int> processTypes;

    ////////////////////////////////
    // MAP SIM HITS TO SIM TRACKS //
    ////////////////////////////////

    // all hits ; type-2 hits ; primary hits ; secondary hits
    int    nMatchedHits[nTypes]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //
    Local3DPoint pos_hit;
    double x_hit=0,y_hit=0,z_hit=0,/*x_cl=0,y_cl=0,*/dx=0,dy=0;
    bool found_hits=false;
    bool fill_dtruth=false;
    //
    // cluster and hit informations
    //unsigned int trkID=-1;
    //unsigned int sizeLink=0;
    //unsigned int rawid=0;
    //unsigned int layer=0;
    unsigned int simh_detid=0;
    unsigned int simh_layer=0;
    int          simh_type=0;
    //unsigned int nLinks = 0;
    //bool combinatoric=false;

    vector<PixelDigi> matched_digis;
    V_HIT_DIGI matched_hits;
    M_TRK_HIT_DIGI map_hits;

    // Fill the map
    int nHits=0;
    for (PSimHitContainer::const_iterator iHit = simHits_B->begin(); iHit != simHits_B->end(); ++iHit) {
      map_hits[iHit->trackId()].push_back( make_pair(*iHit , matched_digis) );
      nHits++ ;
    }
    for (PSimHitContainer::const_iterator iHit = simHits_E->begin(); iHit != simHits_E->end(); ++iHit) {
      map_hits[iHit->trackId()].push_back( make_pair(*iHit , matched_digis) );
      nHits++ ;
    }

    if(verbose>1) cout << endl << "- Number of SimHits in the event : " << nHits << endl;

    // Loop over Sim Tracks and Fill relevant histograms
    int nTracks = 0;

    if(verbose>1) cout << "- SimTracks : " ;

    for (SimTrackContainer::const_iterator simTrkItr = simTracks->begin(); simTrkItr != simTracks->end(); ++simTrkItr) {

      if(verbose>1) cout << simTrkItr->trackId() << " | " ;

        int type = isPrimary((*simTrkItr), simHits_B);
        if(type==-1) type = isPrimary((*simTrkItr), simHits_E);

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

    if(verbose>1) cout << endl;

    nSimTracks_->Fill(nTracks++);

    // Loop Over Digis and Fill Histograms
    if(verbose>1) cout << "- Start looping over DetSetVector<PixelDigi>" << endl;

    for (DetSetVector<PixelDigi>::const_iterator DSViter = pixelDigis->begin(); DSViter != pixelDigis->end(); ++DSViter) {

        // Detector id
        unsigned int rawid = DSViter->id;
        DetId detId(rawid);
        unsigned int layer = getLayerNumber(rawid);

        // Tracker geometry
        const GeomDetUnit* geomDetUnit = tkGeom->idToDetUnit(detId);
        // const PixelGeomDetUnit* pixdet = dynamic_cast<const PixelGeomDetUnit*>(geomDetUnit);
        // const PixelTopology & topol = pixdet->specificTopology();

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

        // Loop over the links in the detector unit
	if(verbose>1) cout << endl << endl << "-- DetId=" << rawid << endl;

        // Loop over the Digis
        for (DetSet<PixelDigi>::const_iterator di = DSViter->data.begin(); di != DSViter->data.end(); ++di) {
            int adc = di->adc();    // charge, modifued to unsiged short
            int col = di->column(); // column
            int row = di->row();    // row

            unsigned int channel = PixelChannelIdentifier::pixelToChannel(row, col);
            unsigned int simTkId = getSimTrackId(pixelSimLinks, detId, channel);

	    if(verbose>1) cout << endl << "--- Digi :"
			       << " col="     << col      
			       << " row="     << row 
			       << " ch="      << channel 
			       << " simTkId=" << simTkId
			       << endl;

            // Get the Digi position
            MeasurementPoint mp(row + 0.5, col + 0.5 );
	    LocalPoint lPos;
	    GlobalPoint pdPos;

            if (geomDetUnit) {
                lPos  = geomDetUnit->topology().localPosition(mp);
                pdPos = geomDetUnit->surface().toGlobal(geomDetUnit->topology().localPosition(mp));

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
                // iPos->second.DigiChargePrimary->Fill(adc);
                nDigiP++;
            }
            else if (primaryTrk == 0) {
                // iPos->second.DigiChargeSecondary->Fill(adc);
                nDigiS++;
            }
	    else cout << "--- primaryTrk=" << primaryTrk 
		      << " | simTkId="     << simTkId
		      << " | iSimTrk="     << iSimTrk
		      << endl;

            nDigiA++;

	    // Get matched hits
	    if(verbose>1) cout << "--- Getting hits matched to SimTrack (id " << simTkId << ")" << endl;
	    matched_hits = map_hits[simTkId];

	    if(verbose>1) {
	      cout << "     number of hits matched to the SimTrack = " << matched_hits.size() ;
	      if(matched_hits.size()!=0) cout << " ids(" ;
	      
	      // printout list of SimHits matched to the SimTrack
	      for (unsigned int iH=0 ; iH<matched_hits.size() ; iH++) {
		cout << matched_hits[iH].first.detUnitId() ;
		if(iH<matched_hits.size()-1) cout << "," ;
		else cout << ")" ;
	      }
	      cout << endl;
	    }
	    
	    // digi matching quantities
	    for(int iM=0 ; iM<nTypes ; iM++)
	      nMatchedHits[iM] = 0;

	    // Loop over matched SimHits
	    if(verbose>2) cout << "--- start looping over matched hits" << endl;
	    for (unsigned int iH=0 ; iH<matched_hits.size() ; iH++) {
	      
	      if(verbose>2) cout << "---- iteration #" << iH << endl;
	      
	      // Consider only SimHits with same DetID as current cluster
	      simh_detid = matched_hits[iH].first.detUnitId();
	      if(simh_detid!=rawid) continue;
	      else found_hits=true;
	      
	      // Map current digi to current SimHit
	      map_hits[ simTkId ][ iH ].second.push_back(*di);

	      simh_layer = getLayerNumber( simh_detid );
	      simh_type  = matched_hits[iH].first.processType();
	      pos_hit    = matched_hits[iH].first.localPosition();
	      x_hit      = pos_hit.x();
	      y_hit      = pos_hit.y();
	      z_hit      = pos_hit.z();

	      if(simh_type>=0 && simh_type<17) nMatchedHits[simh_type]++ ;
	      nMatchedHits[17]++ ;
	      
	      if(simh_type==2) {
		dx = x_hit - lPos.x();
		dy = y_hit - lPos.y();
		if(fill_dtruth==true) fill_dtruth=false; // eliminates cases with several type-2 hits
		fill_dtruth=true; // toggle filling of the histo only when a type-2 hit is found
	      }
	      
	      //if(simh_type == 2 || simh_type == 7 || simh_type == 9 || simh_type == 11 || simh_type == 15)
	      
	      if(verbose>1) cout << "----- SimHit #" << iH
				 << " type="    << simh_type
		//<< " s_id="    << simh_detid
				 << " s_lay="   << simh_layer 
				 << " c_lay="   << layer
				 << " s("   << x_hit    << " , " << y_hit    << " , " << z_hit    << ")"
		//<< " c_g(" << gPos.x() << " , " << gPos.y() << " , " << gPos.z() << ")"
				 << endl;
	      
	    } // end loop over matched SimHits

	    // Number of matched hits (per type)
	    if(verbose>2) cout << "--- Filling NumberOfMatchedHits histograms" << endl;
	    for(int iM=0 ; iM<nTypes ; iM++)
	      iPos->second.NumberOfMatchedHits[iM]-> Fill(nMatchedHits[iM]);
	    
	    // Position resolution
	    if(fill_dtruth) {
	      if(verbose>2) cout << "--- Filling dx,dy histograms" << endl;
	      iPos->second.h_dx_Truth->Fill(dx);
	      iPos->second.h_dy_Truth->Fill(dy);
	    }
	    
	    if(!found_hits && verbose>1) cout << "----- FOUND NO MATCHED HITS" << endl;

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
        // iPos->second.NumberOfDigis->Fill(nDigiA);
        iPos->second.totNDigis += nDigiP + nDigiS;

        // Pixel module
        // if (topol.ncolumns() == 32) iPos->second.NumberOfDigisPixel->Fill(nDigiA);
        // Strip module
        // else if (topol.ncolumns() == 2) iPos->second.NumberOfDigisStrip->Fill(nDigiA);

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
                // iPos->second.ClusterChargePrimary->Fill(cl_charge);
                // iPos->second.ClusterWidthPrimary->Fill(cl_width);
                iPos->second.ClusterWidthVsSimTrkPtPrimary->Fill(trk_pt, cl_width);
                iPos->second.ClusterWidthVsSimTrkEtaPrimary->Fill(trk_eta, cl_width);
                nClusterP++;
            }
            else if (cl_type == 0) {
                // iPos->second.ClusterChargeSecondary->Fill(cl_charge);
                // iPos->second.ClusterWidthSecondary->Fill(cl_width);
                iPos->second.ClusterWidthVsSimTrkPtSecondary->Fill(trk_pt, cl_width);
                iPos->second.ClusterWidthVsSimTrkEtaSecondary->Fill(trk_eta, cl_width);
                nClusterS++;
            }

            nClusterA++;
        }

        // Fill some plots
        // iPos->second.NumberOfClustersPrimary->Fill(nClusterP);
        // iPos->second.NumberOfClustersSecondary->Fill(nClusterS);
        // iPos->second.NumberOfClusters->Fill(nClusterP + nClusterS);
        iPos->second.NumberOfClusters->Fill(nClusterA);
        iPos->second.totNClusters += nClusterP + nClusterS;

        // Pixel module
        // if (topol.ncolumns() == 32) iPos->second.NumberOfClustersPixel->Fill(nClusterA);
        // Strip module
        // else if (topol.ncolumns() == 2) iPos->second.NumberOfClustersStrip->Fill(nClusterA);
    }

    // Fill Layer Level Histograms
    for (map<unsigned int, DigiHistos>::iterator iPos  = layerHistoMap.begin(); iPos != layerHistoMap.end(); ++iPos) {
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
    }

    ////////////////////////////////////
    // COMPUTE CLUSTERIZER EFFICIENCY //
    ////////////////////////////////////
    
    if(verbose>1) cout << "- Enter efficiency computation" << endl;

    // Iterate over the map of hits & clusters
    M_TRK_HIT_DIGI::const_iterator iMapHits;

    // Counters
    int nTrackHits=0;
    int countHit=0;
    int nMatchedDigis=0;
    int nTotalHits=0;
    int nMatchHits=0;
    float efficiency=0;

    // Hit informations
    unsigned int theHit_id=0;   
    unsigned int theHit_layer=0;
    unsigned int theHit_type=0;

    // Prepare the map of counters for efficiency
    std::map<unsigned int, DigiHistos>::iterator iPos;
    map< unsigned int , vector<vector<int>> > map_effi;
    map< unsigned int , vector<vector<int>> >::const_iterator iMapEffi;
    vector<int> init_counter;
    for(int iM=0 ; iM<2 ; iM++) init_counter.push_back(0);

    // Loop over the entries in the map of hits & clusters
    if(verbose>1) cout << "- loop over map of hits & clusters (size=" << map_hits.size() << ")" << endl;

    for( iMapHits=map_hits.begin() ; iMapHits!=map_hits.end() ; iMapHits++) {

      if(verbose>1) cout << "-- SimTrack ID=" << iMapHits->first << endl;
      nTrackHits = (iMapHits->second).size(); 
      countHit  += nTrackHits;

      // Loop over the hits matched to the current SimTrack ID
      for(int iH=0 ; iH<nTrackHits ; iH++) {

	// Current SimHit
	if(verbose>1) cout << "--- SimHit #" << iH << endl;
	PSimHit theHit( ((iMapHits->second)[iH]).first );
	theHit_id    = theHit.detUnitId();
	theHit_layer = getLayerNumber(theHit_id);
	theHit_type  = theHit.processType();
	//if(verbose>1) cout << theHit_id << theHit_layer << theHit_type << endl;

	// Clusters matched to the SimHit
	matched_digis = ((iMapHits->second)[iH]).second;
	nMatchedDigis = matched_digis.size();

	// Find layer in map of histograms
	if(verbose>2) cout << "--- Find layer=" << theHit_layer << " in map of histograms" << endl;
	iPos = layerHistoMap.find(theHit_layer);
        if (iPos == layerHistoMap.end()) {
	  if(verbose>2) cout << "---- add layer in the map" << endl;
	  createLayerHistograms(theHit_layer);
	  iPos = layerHistoMap.find(theHit_layer);
        }

	// Fill Histograms
	(iPos->second.NumberOfMatchedDigis[17])->Fill( nMatchedDigis );
	if(theHit_type<17)
	  (iPos->second.NumberOfMatchedDigis[theHit_type])->Fill( nMatchedDigis );

	if( nMatchedDigis==0 ) {
	  if(verbose>1) cout << "---- No Digi Matched" << endl;
	}
	else {
	  if(verbose>1) cout << "---- Yes Digi Matched = " << nMatchedDigis << endl;
	}
	
	if( map_effi.find(theHit_layer)==map_effi.end() )
	  for(int iT=0 ; iT<nTypes ; iT++) {
	    map_effi[theHit_layer].push_back(init_counter);
	    if(verbose>2) cout << "----- type #" << iT << " layer=" << theHit_layer << " map size=" << map_effi.size() << endl;
	  }
	(map_effi[theHit_layer][theHit_type][0])++ ; // total number of hits of this type in this layer
	if(nMatchedDigis>0) (map_effi[theHit_layer][theHit_type][1])++ ; // number of hits matched to >=1 digi(s)

	(map_effi[theHit_layer][17][0])++ ; // total number of hits of this type in this layer
	if(nMatchedDigis>0) (map_effi[theHit_layer][17][1])++ ; // number of hits matched to >=1 digi(s)
      }

    }

    // Fill histograms from the map_effi
    if(verbose>1) cout << "- fill [per layer] effi histo from effi map (size=" << map_effi.size() << ")" << endl;
    for( iMapEffi=map_effi.begin() ; iMapEffi!=map_effi.end() ; iMapEffi++ ) {
      
      iPos = layerHistoMap.find(iMapEffi->first);
      if(verbose>1) cout << "-- layer=" << iMapEffi->first << endl;

      for(int iT=0 ; iT<nTypes ; iT++) {
	nTotalHits = iMapEffi->second[iT][0];
	nMatchHits = iMapEffi->second[iT][1]; // number of hits in a given layer that have at least 1 digi matched
	efficiency = nTotalHits!=0 ? float(nMatchHits)/float(nTotalHits) : -1 ;
	if(efficiency>=0) (iPos->second.hEfficiency[iT])->Fill( efficiency );
	if(verbose>1) cout << "--- type #"   << iT 
			   << " nTotalHits=" << nTotalHits 
			   << " nMatchHits=" << nMatchHits 
			   << " efficiency=" << efficiency
			   << endl;
      }
    }

    // Check if all event's SimHits are mapped
    if( countHit != nHits )
      if(verbose>1) cout << "---- Missing hits in the efficiency computation : " 
			 << countHit << " != " << nHits << endl;

}

void RunStepsDigiValidation::endJob() { }

void RunStepsDigiValidation::createLayerHistograms(unsigned int ival) {
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
    ostringstream histoName;

    histoName << "NumberOfDigisPrimary" << tag.c_str() << id;
    local_histos.NumberOfDigisPrimary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 51, -0.5, 50.5);
    histoName.str("");
    histoName << "NumberOfDigisSecondary" << tag.c_str() << id;
    local_histos.NumberOfDigisSecondary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 51, -0.5, 50.5);
    // histoName.str("");
    // histoName << "NumberOfDigis_Pixel" << tag.c_str() << id;
    // local_histos.NumberOfDigisPixel = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 51, -0.5, 50.5);
    // histoName.str("");
    // histoName << "NumberOfDigis_Strip" << tag.c_str() << id;
    // local_histos.NumberOfDigisStrip = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 51, -0.5, 50.5);

    local_histos.NumberOfDigisPrimary->SetFillColor(kBlue);
    local_histos.NumberOfDigisSecondary->SetFillColor(kRed);

    // local_histos.NumberOfDigisPixel->SetFillColor(kBlue);
    // local_histos.NumberOfDigisStrip->SetFillColor(kRed);

    // histoName.str("");
    // histoName << "NumberOfDigisType" << tag.c_str() << id;
    // local_histos.NumberOfDigisType = td.make<THStack>(histoName.str().c_str(), histoName.str().c_str());
    // local_histos.NumberOfDigisType->Add(local_histos.NumberOfDigisPrimary);
    // local_histos.NumberOfDigisType->Add(local_histos.NumberOfDigisSecondary);

    // histoName.str("");
    // histoName << "NumberOfDigisSource" << tag.c_str() << id;
    // local_histos.NumberOfDigisSource = td.make<THStack>(histoName.str().c_str(), histoName.str().c_str());
    // local_histos.NumberOfDigisSource->Add(local_histos.NumberOfDigisPixel);
    // local_histos.NumberOfDigisSource->Add(local_histos.NumberOfDigisStrip);

    histoName.str("");
    histoName << "DigiCharge" << tag.c_str() << id;
    local_histos.DigiCharge = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 261, -0.5, 260.5);
    // histoName.str("");
    // histoName << "DigiChargePrimary" << tag.c_str() << id;
    // local_histos.DigiChargePrimary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 261, -0.5, 260.5);
    // histoName.str("");
    // histoName << "DigiChargeSecondary" << tag.c_str() << id;
    // local_histos.DigiChargeSecondary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 261, -0.5, 260.5);

    histoName.str("");
    histoName << "NumberOfClusters" << tag.c_str() << id;
    local_histos.NumberOfClusters = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 51, -0.5, 50.5);
    // histoName.str("");
    // histoName << "NumberOfClustersPrimary" << tag.c_str() << id;
    // local_histos.NumberOfClustersPrimary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 51, -0.5, 50.5);
    // histoName.str("");
    // histoName << "NumberOfClustersSecondary" << tag.c_str() << id;
    // local_histos.NumberOfClustersSecondary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 51, -0.5, 50.5);
    // histoName.str("");
    // histoName << "NumberOfClusters_Pixel" << tag.c_str() << id;
    // local_histos.NumberOfClustersPixel = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 51, -0.5, 50.5);
    // histoName.str("");
    // histoName << "NumberOfClusters_Strip" << tag.c_str() << id;
    // local_histos.NumberOfClustersStrip = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 51, -0.5, 50.5);

    histoName.str("");
    histoName << "ClusterCharge" << tag.c_str() << id;
    local_histos.ClusterCharge = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 1041, -0.5, 1040.5);
    // histoName.str("");
    // histoName << "ClusterChargePrimary" << tag.c_str() << id;
    // local_histos.ClusterChargePrimary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 1041, -0.5,1040.5);
    // histoName.str("");
    // histoName << "ClusterChargeSecondary" << tag.c_str() << id;
    // local_histos.ClusterChargeSecondary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 1041, -0.5, 1040.5);

    histoName.str("");
    histoName << "ClusterWidth" << tag.c_str() << id;
    local_histos.ClusterWidth = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 16, -0.5, 15.5);
    // histoName.str("");
    // histoName << "ClusterWidthPrimary" << tag.c_str() << id;
    // local_histos.ClusterWidthPrimary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 16, -0.5, 15.5);
    // histoName.str("");
    // histoName << "ClusterWidthSecondary" << tag.c_str() << id;
    // local_histos.ClusterWidthSecondary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 16, -0.5, 15.5);

    histoName.str("");
    histoName << "TotalNumberOfDigis" << tag.c_str() << id;
    local_histos.TotalNumberOfDigis = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 100, -0.5, 100.5);

    histoName.str("");
    histoName << "TotalNumberOfClusters" << tag.c_str() << id;
    local_histos.TotalNumberOfClusters = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 100, -0.5, 100.5);

    histoName.str("");
    histoName << "ClusterShape" << tag.c_str() << id;
    local_histos.ClusterShape = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 21, -20.5, 20.5);
    histoName.str("");
    histoName<< "ClusterShapePrimary" << tag.c_str() << id;
    local_histos.ClusterShapePrimary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 21, -20.5, 20.5);
    histoName.str("");
    histoName << "ClusterShapeSecondart" << tag.c_str() << id;
    local_histos.ClusterShapeSecondart = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 21, -20.5, 20.5);

    histoName.str("");
    histoName << "NumberOfSimHits" << tag.c_str() << id;
    local_histos.NumberOfSimHits = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 201, -0.5, 200.5);
    histoName.str("");
    histoName << "NumberOfMatchedSimHits" << tag.c_str() << id;
    local_histos.NumberOfMatchedSimHits = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 201, -0.5, 200.5);
    histoName.str("");
    histoName << "NumberOfMatchedSimHitsPrimary" << tag.c_str() << id;
    local_histos.NumberOfMatchedSimHitsPrimary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 201, -0.5, 200.5);
    histoName.str("");
    histoName << "NumberOfMatchedSimHitsSecondary" << tag.c_str() << id;
    local_histos.NumberOfMatchedSimHitsSecondary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 201, -0.5, 200.5);

    histoName.str("");
    histoName << "DigiEfficiency" << tag.c_str() << id;
    local_histos.DigiEfficiency = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 55, -0.05, 1.05);
    histoName.str("");
    histoName << "DigiEfficiencyPrimary" << tag.c_str() << id;
    local_histos.DigiEfficiencyPrimary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 55, -0.05, 1.05);
    histoName.str("");
    histoName << "DigiEfficiencySecondary" << tag.c_str() << id;
    local_histos.DigiEfficiencySecondary = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 55, -0.05, 1.05);

    histoName.str("");
    histoName << "YposVsXpos" << tag.c_str() << id;
    local_histos.YposVsXpos = td.make<TH2F>(histoName.str().c_str(), histoName.str().c_str(), 240, -120.0, 120.0, 240, -120.0, 120.0);

    histoName.str("");
    histoName << "RVsZpos" << tag.c_str() << id;
    local_histos.RVsZpos = td.make<TH2F>(histoName.str().c_str(), histoName.str().c_str(), 600, -300.0, 300.0, 120, 0.0, 120.0);

    histoName.str("");
    histoName << "DigiChargeMatched" << tag.c_str() << id;
    local_histos.DigiChargeMatched = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 261, -0.5, 260.5);

    histoName.str("");
    histoName << "ClusterWidthVsSimTrkPt" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkPt = td.make<TProfile>(histoName.str().c_str(),histoName.str().c_str(),56, -0.5, 55.5,-0.5,15.5);
    histoName.str("");
    histoName << "ClusterWidthVsSimTrkPtPrimary" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkPtPrimary = td.make<TProfile>(histoName.str().c_str(),histoName.str().c_str(),56, -0.5, 55.5,-0.5,15.5);
    histoName.str("");
    histoName << "ClusterWidthVsSimTkrPtS" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkPtSecondary = td.make<TProfile>(histoName.str().c_str(),histoName.str().c_str(),56, -0.5, 55.5,-0.5,15.5);

    histoName.str("");
    histoName << "ClusterWidthVsSimTrkEta" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkEta = td.make<TProfile>(histoName.str().c_str(),histoName.str().c_str(),50, -2.5, 2.5,-0.5,15.5);
    histoName.str("");
    histoName << "ClusterWidthVsSimTrkEtaPrimary" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkEtaPrimary = td.make<TProfile>(histoName.str().c_str(),histoName.str().c_str(),50, -2.5, 2.5,-0.5,15.5);
    histoName.str("");
    histoName << "ClusterWidthVsSimTkrEtaS" << tag.c_str() << id;
    local_histos.ClusterWidthVsSimTrkEtaSecondary = td.make<TProfile>(histoName.str().c_str(),histoName.str().c_str(),50, -2.5, 2.5,-0.5,15.5);

    histoName.str("");
    histoName << "MatchedSimTrackPt" << tag.c_str() << id;
    local_histos.matchedSimTrackPt_  = td.make<TH1F>(histoName.str().c_str(),histoName.str().c_str(),101,-0.5,100.5);
    histoName.str("");
    histoName << "MatchedSimTrackPtP" << tag.c_str() << id;
    local_histos.matchedSimTrackPtPrimary_  = td.make<TH1F>(histoName.str().c_str(),histoName.str().c_str(),101,-0.5,100.5);
    histoName.str("");
    histoName << "MatchedSimTrackPtS" << tag.c_str() << id;
    local_histos.matchedSimTrackPtSecondary_  = td.make<TH1F>(histoName.str().c_str(),histoName.str().c_str(),101,-0.5,100.5);

    histoName.str("");
    histoName << "MatchedSimTrackEta" << tag.c_str() << id;
    local_histos.matchedSimTrackEta_  = td.make<TH1F>(histoName.str().c_str(),  histoName.str().c_str(), 50, -2.5, 2.5);
    histoName.str("");
    histoName << "MatchedSimTrackEtaP" << tag.c_str() << id;
    local_histos.matchedSimTrackEtaPrimary_  = td.make<TH1F>(histoName.str().c_str(),  histoName.str().c_str(), 50, -2.5, 2.5);
    histoName.str("");
    histoName << "MatchedSimTrackEtaS" << tag.c_str() << id;
    local_histos.matchedSimTrackEtaSecondary_  = td.make<TH1F>(histoName.str().c_str(),  histoName.str().c_str(), 50, -2.5, 2.5);

    histoName.str("");
    histoName << "MatchedSimTrackPhi" << tag.c_str() << id;
    local_histos.matchedSimTrackPhi_  = td.make<TH1F>(histoName.str().c_str(),  histoName.str().c_str(), 160, -3.2, 3.2);
    histoName.str("");
    histoName << "MatchedSimTrackPhiP" << tag.c_str() << id;
    local_histos.matchedSimTrackPhiPrimary_  = td.make<TH1F>(histoName.str().c_str(),  histoName.str().c_str(), 160, -3.2, 3.2);
    histoName.str("");
    histoName << "MatchedSimTrackPhiS" << tag.c_str() << id;
    local_histos.matchedSimTrackPhiSecondary_  = td.make<TH1F>(histoName.str().c_str(),  histoName.str().c_str(), 160, -3.2, 3.2);

    histoName.str("");
    histoName << "LocalPosition" << tag.c_str() << id;
    local_histos.LocalPosition = td.make<TH2F>(histoName.str().c_str(),histoName.str().c_str(),10000, -5, 5 , 10000, -5 ,5);


    // Truth Matching 

    string name_types[nTypes] = {"Undefined","Unknown","Primary","Hadronic",
				 "Decay","Compton","Annihilation","EIoni",
				 "HIoni","MuIoni","Photon","MuPairProd",
				 "Conversions","EBrem","SynchrotronRadiation",
				 "MuBrem","MuNucl","AllTypes"};

    for(int iM=0 ; iM<nTypes ; iM++) {
      histoName.str("");
      histoName << "NumberOfMatchedHits_" << name_types[iM] << tag.c_str() <<  id;
      local_histos.NumberOfMatchedHits[iM] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 21, 0., 20.);

      histoName.str("");
      histoName << "NumberOfMatchedDigis_" << name_types[iM] << tag.c_str() <<  id;
      local_histos.NumberOfMatchedDigis[iM] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 21, 0., 20.);

      histoName.str("");
      histoName << "Efficiency_" << name_types[iM] << tag.c_str() <<  id;
      local_histos.hEfficiency[iM] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 110, 0., 1.1);
    }

    histoName.str("");
    histoName << "DeltaX_simhit_digi" << tag.c_str() <<  id;
    local_histos.h_dx_Truth = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 1000, 0., 0.);

    histoName.str("");
    histoName << "DeltaY_simhit_digi" << tag.c_str() <<  id;
    local_histos.h_dy_Truth = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 1000, 0., 0.);


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
}

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
}

DEFINE_FWK_MODULE(RunStepsDigiValidation);
