#include <memory>
#include <map>

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "RunSteps/Clusterizer/interface/PixelClusterSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

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
const int nMatch=4;


using namespace std;
using namespace edm;

class RunStepsClusterValidation : public EDAnalyzer {

  typedef vector< pair< PSimHit , vector<SiPixelCluster> > > V_HIT_CLUSTERS;
  typedef map< int , V_HIT_CLUSTERS >                        M_TRK_HIT_CLUSTERS;

public:

    explicit RunStepsClusterValidation(const ParameterSet&);
    ~RunStepsClusterValidation();
    virtual void beginJob();
    virtual void analyze(const Event&, const EventSetup&);
    virtual void endJob();
    int isPrimary(const SimTrack& simTrk, Handle<PSimHitContainer>& simHits);
    int isPrimary(const SimTrack& simTrk, const PSimHit& simHit);

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

        TH1F* NumberOfMatchedHits[4];
        TH1F* h_dx_Truth;
        TH1F* h_dy_Truth;

        TH1F* hEfficiency[17];

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

  verbose=2;

    // Simulated information
    // SimHit
    Handle<PSimHitContainer> simHits_B;
    iEvent.getByLabel("g4SimHits","TrackerHitsPixelBarrelLowTof" ,simHits_B);

    Handle<PSimHitContainer> simHits_E;
    iEvent.getByLabel("g4SimHits","TrackerHitsPixelEndcapLowTof" ,simHits_E);

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

    ////////////////////////////////
    // MAP SIM HITS TO SIM TRACKS //
    ////////////////////////////////
    vector<SiPixelCluster> matched_clusters;
    V_HIT_CLUSTERS         matched_hits;
    M_TRK_HIT_CLUSTERS     map_hits;

    // Fill the map
    int nHits=0;
    for (PSimHitContainer::const_iterator iHit = simHits_B->begin(); iHit != simHits_B->end(); ++iHit) {
      map_hits[iHit->trackId()].push_back( make_pair(*iHit , matched_clusters) ) ;
      nHits++ ;
    }
    for (PSimHitContainer::const_iterator iHit = simHits_E->begin(); iHit != simHits_E->end(); ++iHit) {
      map_hits[iHit->trackId()].push_back( make_pair(*iHit , matched_clusters) ) ;
      nHits++ ;
    }

    if(verbose>1) cout << endl << "-- Number of SimHits in the event : " << nHits << endl;
    
    //////////////////////////////////
    // LOOP OVER CLUSTER COLLECTION //
    //////////////////////////////////

    // Loop over the detector units

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

        // Loop over the clusters in the detector unit
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

            // Loop over the pixels
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

    /////////////////////////////
    // LOOP OVER CLUSTER LINKS //
    /////////////////////////////

    DetSetVector<PixelClusterSimLink>::const_iterator DSViterLinks;
    DetSet<PixelClusterSimLink>::const_iterator iterLinks;
    std::vector< unsigned int > simTrackID;
    PixelClusterSimLink link;
    SiPixelCluster cluster;

    // cluster and hit informations
    unsigned int trkID=-1;
    unsigned int sizeLink=0;
    unsigned int rawid=0;
    unsigned int layer=0;
    unsigned int simh_detid=0;
    unsigned int simh_layer=0;
    int          simh_type=0;
    unsigned int nLinks = 0;
    bool combinatoric=false;

    // matching quantities
    //
    // all hits ; type-2 hits ; primary hits ; secondary hits
    int    nMatchedHits[nMatch]={0,0,0,0}; 
    //
    Local3DPoint pos_hit;
    double x_hit=0,y_hit=0,z_hit=0,x_cl=0,y_cl=0,dx=0,dy=0;
    bool found_hits=false;
    bool fill_dtruth=false;

    if(verbose>1) cout << endl << "-- Enter loop over links" << endl;

    for (DSViterLinks = clusterLinks->begin(); DSViterLinks != clusterLinks->end(); DSViterLinks++) {

        trkID=-1;
	sizeLink=0;
	combinatoric=false;

	// Get the detector unit's id
        rawid = DSViterLinks->detId();
        layer = getLayerNumber(rawid);
        DetId detId(rawid);

        // Get the geometry of the tracker
        const GeomDetUnit* geomDetUnit = tkGeom->idToDetUnit(detId);

        if (!geomDetUnit) break;

        // Create histograms for the layer if they do not yet exist
        std::map<unsigned int, ClusterHistos>::iterator iPos = layerHistoMap.find(layer);
        if (iPos == layerHistoMap.end()) {
            createLayerHistograms(layer);
            iPos = layerHistoMap.find(layer);
        }

        // Loop over the links in the detector unit
	if(verbose>1) cout << endl << endl << "--- DetId=" << rawid << endl;

        for (iterLinks = DSViterLinks->data.begin(); iterLinks != DSViterLinks->data.end(); ++iterLinks) {

            // Link informations
	    combinatoric=false;
            nLinks++;
	    link = *iterLinks;
	    simTrackID = link.getSimTracks();
	    sizeLink   = simTrackID.size();

	    // cluster matching quantities
	    for(int iM=0 ; iM<nMatch ; iM++)
	      nMatchedHits[iM] = 0;

	    // Cluster informations
	    cluster = link.getCluster();
	    x_cl    = cluster.x();
	    y_cl    = cluster.y();
            MeasurementPoint mp(x_cl, y_cl);
            LocalPoint lPos  = geomDetUnit->topology().localPosition(mp);
            //GlobalPoint gPos = geomDetUnit->surface().toGlobal(geomDetUnit->topology().localPosition(mp));

	    if(verbose>1) cout << endl << "---- Cluster size="  << cluster.size()    << " | " << sizeLink << " SimTracks: ids=(" ;

	    for(unsigned int i=0 ; i<sizeLink ; i++) {	      
	      if(verbose>1) {
		cout << simTrackID[i] ;
		if(i<sizeLink-1) cout << "," ;
	      }
	      if(i==0) trkID=simTrackID[i];
	      else if(simTrackID[i]!=trkID) combinatoric=true; 
	    }

	    if(verbose>1) {
	      if(combinatoric) cout << ") COMBINATORIC !!! ";
	      cout << ")" << endl << "     cluster local position = (" << lPos.x() << " , " << lPos.y() << " , " << "n/a"    << ")"
		   << "   layer=" << layer
		   << endl;
	    }

	    // Get matched SimHits from the map
	    if(!combinatoric) {
	      
	      matched_hits = map_hits[trkID];

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

	      found_hits=false;
	      fill_dtruth=false;

	      // Loop over matched SimHits

	      for (unsigned int iH=0 ; iH<matched_hits.size() ; iH++) {
		
		// Consider only SimHits with same DetID as current cluster
		simh_detid = matched_hits[iH].first.detUnitId();
		if(simh_detid!=rawid) continue;
		else found_hits=true;

		// Map current cluster to current SimHit for efficiency study
		map_hits[ trkID ][ iH ].second.push_back( cluster );
		
		simh_layer = getLayerNumber( simh_detid );
		simh_type  = matched_hits[iH].first.processType();
		pos_hit    = matched_hits[iH].first.localPosition();
		x_hit      = pos_hit.x();
		y_hit      = pos_hit.y();
		z_hit      = pos_hit.z();

		nMatchedHits[0]++ ;
		if(simh_type==2) {
		  nMatchedHits[1]++ ;
		  dx = x_hit - lPos.x();
		  dy = y_hit - lPos.y();
		  if(fill_dtruth==true) fill_dtruth=false; // eliminates cases with several type-2 hits
		  fill_dtruth=true; // toggle filling of the histo only when a type-2 hit is found
		}

		if(simh_type == 2 || simh_type == 7 || simh_type == 9 || simh_type == 11 || simh_type == 15)
		  nMatchedHits[2]++ ;
		else 
		  nMatchedHits[3]++ ;
		
		if(verbose>1) cout << "----- SimHit #" << iH
				   << " type="    << simh_type
		  //<< " s_id="    << simh_detid
				   << " s_lay="   << simh_layer 
				   << " c_lay="   << layer
				   << " s("   << x_hit    << " , " << y_hit    << " , " << z_hit    << ")"
		  //<< " c_g(" << gPos.x() << " , " << gPos.y() << " , " << gPos.z() << ")"
				   << endl;
		
	      } // end loop over matched SimHits

	      if(!found_hits && verbose>1) cout << "----- FOUND NO MATCHED HITS" << endl;

	    } // endif !combinatoric

	    // Number of matched hits (per type)
	    for(int iM=0 ; iM<nMatch ; iM++)
	      iPos->second.NumberOfMatchedHits[iM]-> Fill(nMatchedHits[iM]);

	    // Position resolution
	    if(fill_dtruth) {
	      iPos->second.h_dx_Truth->Fill(dx);
	      iPos->second.h_dy_Truth->Fill(dy);
	    }

        } // end loop over links within a single DetID

        iPos->second.NumberOfClustersLink         -> Fill(nLinks);

    } // end loop over all links


    ////////////////////////////////////
    // COMPUTE CLUSTERIZER EFFICIENCY //
    ////////////////////////////////////
    
    if(verbose>1) cout << "- Enter efficiency computation" << endl;

    // Iterate over the map of hits & clusters
    M_TRK_HIT_CLUSTERS::const_iterator iMapHits;

    // Counters
    int nTrackHits=0;
    int countHit=0;
    int nMatchedClusters=0;
    int nTotalHits=0;
    int nMatchHits=0;
    float efficiency=0;

    // Hit informations
    unsigned int theHit_id=0;   
    unsigned int theHit_layer=0;
    unsigned int theHit_type=0;

    // Prepare the map of counters for efficiency
    std::map<unsigned int, ClusterHistos>::iterator iPos;
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
	matched_clusters = ((iMapHits->second)[iH]).second;
	nMatchedClusters = matched_clusters.size();

	if( nMatchedClusters==0 ) {
	  if(verbose>1) cout << "---- No Cluster Matched" << endl;
	}
	else {
	  if(verbose>1) cout << "---- Yes Cluster Matched = " << nMatchedClusters << endl;
	}
	
	if( map_effi.find(theHit_layer)==map_effi.end() )
	  for(int iT=0 ; iT<17 ; iT++) {
	    map_effi[theHit_layer].push_back(init_counter);
	    if(verbose>2) cout << "----- type #" << iT << " layer=" << theHit_layer << " map size=" << map_effi.size() << endl;
	  }
	(map_effi[theHit_layer][theHit_type][0])++ ; // total number of hits of this type in this layer
	if(nMatchedClusters>0) (map_effi[theHit_layer][theHit_type][1])++ ; // number of hits matched to >=1 cluster(s)
      }

    }

    // Fill histograms from the map_effi
    if(verbose>1) cout << "- fill [per layer] effi histo from effi map (size=" << map_effi.size() << ")" << endl;
    for( iMapEffi=map_effi.begin() ; iMapEffi!=map_effi.end() ; iMapEffi++ ) {
      
      iPos = layerHistoMap.find(iMapEffi->first);
      if(verbose>1) cout << "-- layer=" << iMapEffi->first << endl;

      for(int iT=0 ; iT<17 ; iT++) {
	nTotalHits = iMapEffi->second[iT][0];
	nMatchHits = iMapEffi->second[iT][1];
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

    // Truth Matching 
    string name_match[nMatch] = {"AllType", "Type2", "Primary", "Secondary"};
    for(int iM=0 ; iM<nMatch ; iM++) {
      histoName.str("");
      histoName << "NumberOfMatchedHits" << name_match[iM] << tag.c_str() <<  id;
      local_histos.NumberOfMatchedHits[iM] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 21, 0., 20.);
    }

    histoName.str("");
    histoName << "DeltaX_simhit_cluster" << tag.c_str() <<  id;
    local_histos.h_dx_Truth = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 1000, 0., 0.);

    histoName.str("");
    histoName << "DeltaY_simhit_cluster" << tag.c_str() <<  id;
    local_histos.h_dy_Truth = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 1000, 0., 0.);

    string name_types[17] = {"Undefined","Unknown","Primary","Hadronic",
			     "Decay","Compton","Annihilation","EIoni",
			     "HIoni","MuIoni","Photon","MuPairProd",
			     "Conversions","EBrem","SynchrotronRadiation",
			     "MuBrem","MuNucl"};

    for(int iM=0 ; iM<17 ; iM++) {
      histoName.str("");
      histoName << "Efficiency" << name_types[iM] << tag.c_str() <<  id;
      local_histos.hEfficiency[iM] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 110, 0., 1.1);
    }

    // Cluster topology

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

int RunStepsClusterValidation::isPrimary(const SimTrack& simTrk, Handle<PSimHitContainer>& simHits) {
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

int RunStepsClusterValidation::isPrimary(const SimTrack& simTrk, const PSimHit& simHit) {
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

DEFINE_FWK_MODULE(RunStepsClusterValidation);
