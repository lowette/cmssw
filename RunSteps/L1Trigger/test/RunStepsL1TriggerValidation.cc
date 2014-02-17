#include <map>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include <TH1D.h>
#include <TH2D.h>

using namespace std;
using namespace edm;

class RunStepsL1TriggerValidation : public EDAnalyzer {

public:

    explicit RunStepsL1TriggerValidation(const ParameterSet& iConfig);
    virtual ~RunStepsL1TriggerValidation();
    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const Event& iEvent, const EventSetup& iSetup);

private:

    // TrackingParticle and TrackingVertex
    TH2D* hSimVtx_XY;
    TH2D* hSimVtx_RZ;

    TH1D* hTPart_Pt;
    TH1D* hTPart_Eta_Pt10;
    TH1D* hTPart_Phi_Pt10;

    // Global positions of TTClusters
    TH2D* hCluster_Barrel_XY;
    TH2D* hCluster_Barrel_XY_Zoom;
    TH2D* hCluster_Endcap_Fw_XY;
    TH2D* hCluster_Endcap_Bw_XY;
    TH2D* hCluster_RZ;
    TH2D* hCluster_Endcap_Fw_RZ_Zoom;
    TH2D* hCluster_Endcap_Bw_RZ_Zoom;

    TH1D* hCluster_IMem_Barrel;
    TH1D* hCluster_IMem_Endcap;
    TH1D* hCluster_OMem_Barrel;
    TH1D* hCluster_OMem_Endcap;

    TH1D* hCluster_Gen_Barrel;
    TH1D* hCluster_Unkn_Barrel;
    TH1D* hCluster_Comb_Barrel;
    TH1D* hCluster_Gen_Endcap;
    TH1D* hCluster_Unkn_Endcap;
    TH1D* hCluster_Comb_Endcap;

    TH1D* hCluster_Gen_Eta;
    TH1D* hCluster_Unkn_Eta;
    TH1D* hCluster_Comb_Eta;

    TH2D* hCluster_PID;
    TH2D* hCluster_W;

    TH1D* hTPart_Eta_INormalization;
    TH1D* hTPart_Eta_ICW_1;
    TH1D* hTPart_Eta_ICW_2;
    TH1D* hTPart_Eta_ICW_3;
    TH1D* hTPart_Eta_ONormalization;
    TH1D* hTPart_Eta_OCW_1;
    TH1D* hTPart_Eta_OCW_2;
    TH1D* hTPart_Eta_OCW_3;

    // Global positions of TTStubs
    TH2D* hStub_Barrel_XY;
    TH2D* hStub_Barrel_XY_Zoom;
    TH2D* hStub_Endcap_Fw_XY;
    TH2D* hStub_Endcap_Bw_XY;
    TH2D* hStub_RZ;
    TH2D* hStub_Endcap_Fw_RZ_Zoom;
    TH2D* hStub_Endcap_Bw_RZ_Zoom;

    TH1D* hStub_Barrel;
    TH1D* hStub_Endcap;

    TH1D* hStub_Gen_Barrel;
    TH1D* hStub_Unkn_Barrel;
    TH1D* hStub_Comb_Barrel;
    TH1D* hStub_Gen_Endcap;
    TH1D* hStub_Unkn_Endcap;
    TH1D* hStub_Comb_Endcap;

    TH1D* hStub_Gen_Eta;
    TH1D* hStub_Unkn_Eta;
    TH1D* hStub_Comb_Eta;

    TH1D* hStub_PID;
    TH2D* hStub_Barrel_W;
    TH2D* hStub_Barrel_O;
    TH2D* hStub_Endcap_W;
    TH2D* hStub_Endcap_O;

    // Stub finding coverage
    TH1D* hTPart_Eta_Pt10_Normalization;
    TH1D* hTPart_Eta_Pt10_NumPS;
    TH1D* hTPart_Eta_Pt10_Num2S;

    // Denominator for Stub Prod Eff
    map< unsigned int, TH1D* > mapCluLayer_hTPart_Pt;
    map< unsigned int, TH1D* > mapCluLayer_hTPart_Eta_Pt10;
    map< unsigned int, TH1D* > mapCluLayer_hTPart_Phi_Pt10;
    map< unsigned int, TH1D* > mapCluDisk_hTPart_Pt;
    map< unsigned int, TH1D* > mapCluDisk_hTPart_Eta_Pt10;
    map< unsigned int, TH1D* > mapCluDisk_hTPart_Phi_Pt10;

    // Cluster size
    map< unsigned int, TH1D* > mapCluSizeBarrel;
    map< unsigned int, TH1D* > mapCluSizeBarrel_Pixel;
    map< unsigned int, TH1D* > mapCluSizeBarrel_Strip;
    map< unsigned int, TH1D* > mapCluSizeEndcap;
    map< unsigned int, TH1D* > mapCluSizeEndcap_Pixel;
    map< unsigned int, TH1D* > mapCluSizeEndcap_Strip;

    // Local Positions
    map< unsigned int, TH2D* > mapLocalPositionBarrel;
    map< unsigned int, TH2D* > mapLocalPositionBarrel_Pixel;
    map< unsigned int, TH2D* > mapLocalPositionBarrel_Strip;
    map< unsigned int, TH2D* > mapLocalPositionEndcap;
    map< unsigned int, TH2D* > mapLocalPositionEndcap_Pixel;
    map< unsigned int, TH2D* > mapLocalPositionEndcap_Strip;

    // Numerator for Stub Prod Eff
    map< unsigned int, TH1D* > mapStubLayer_hTPart_Pt;
    map< unsigned int, TH1D* > mapStubLayer_hTPart_Eta_Pt10;
    map< unsigned int, TH1D* > mapStubLayer_hTPart_Phi_Pt10;
    map< unsigned int, TH1D* > mapStubDisk_hTPart_Pt;
    map< unsigned int, TH1D* > mapStubDisk_hTPart_Eta_Pt10;
    map< unsigned int, TH1D* > mapStubDisk_hTPart_Phi_Pt10;

    // Comparison of Stubs to TrackingParticles
    map< unsigned int, TH2D* > mapStubLayer_hStub_InvPt_TPart_InvPt;
    map< unsigned int, TH2D* > mapStubLayer_hStub_Pt_TPart_Pt;
    map< unsigned int, TH2D* > mapStubLayer_hStub_Eta_TPart_Eta;
    map< unsigned int, TH2D* > mapStubLayer_hStub_Phi_TPart_Phi;
    map< unsigned int, TH2D* > mapStubDisk_hStub_InvPt_TPart_InvPt;
    map< unsigned int, TH2D* > mapStubDisk_hStub_Pt_TPart_Pt;
    map< unsigned int, TH2D* > mapStubDisk_hStub_Eta_TPart_Eta;
    map< unsigned int, TH2D* > mapStubDisk_hStub_Phi_TPart_Phi;

    // Residuals
    map< unsigned int, TH2D* > mapStubLayer_hStub_InvPtRes_TPart_Eta;
    map< unsigned int, TH2D* > mapStubLayer_hStub_PtRes_TPart_Eta;
    map< unsigned int, TH2D* > mapStubLayer_hStub_EtaRes_TPart_Eta;
    map< unsigned int, TH2D* > mapStubLayer_hStub_PhiRes_TPart_Eta;
    map< unsigned int, TH2D* > mapStubDisk_hStub_InvPtRes_TPart_Eta;
    map< unsigned int, TH2D* > mapStubDisk_hStub_PtRes_TPart_Eta;
    map< unsigned int, TH2D* > mapStubDisk_hStub_EtaRes_TPart_Eta;
    map< unsigned int, TH2D* > mapStubDisk_hStub_PhiRes_TPart_Eta;

    // Stub Width vs Pt
    map< unsigned int, TH2D* > mapStubLayer_hStub_W_TPart_Pt;
    map< unsigned int, TH2D* > mapStubLayer_hStub_W_TPart_InvPt;
    map< unsigned int, TH2D* > mapStubDisk_hStub_W_TPart_Pt;
    map< unsigned int, TH2D* > mapStubDisk_hStub_W_TPart_InvPt;

    // Containers of parameters passed by python configuration file
    ParameterSet config;

    bool testedGeometry;
    bool DebugMode;
};

RunStepsL1TriggerValidation::RunStepsL1TriggerValidation(ParameterSet const& iConfig) :
    config(iConfig) {
    DebugMode = iConfig.getParameter< bool >("DebugMode");
    cout << "------------------------------------------------------------" << endl
         << "-- Running RunSteps L1TriggerValidation v0.1" << endl
        << "------------------------------------------------------------" << endl;
}

RunStepsL1TriggerValidation::~RunStepsL1TriggerValidation() { }

void RunStepsL1TriggerValidation::analyze(const Event& iEvent, const EventSetup& iSetup) {
    // Geometry handles etc
    ESHandle< TrackerGeometry > GeometryHandle;
    ESHandle< StackedTrackerGeometry > StackedGeometryHandle;
    const StackedTrackerGeometry* theStackedGeometry;
    StackedTrackerGeometry::StackContainerIterator StackedTrackerIterator;

    // Geometry setup
    // Set pointers to Geometry
    iSetup.get< TrackerDigiGeometryRecord >().get(GeometryHandle);

    // Set pointers to Stacked Modules
    iSetup.get< StackedTrackerGeometryRecord >().get(StackedGeometryHandle);
    theStackedGeometry = StackedGeometryHandle.product(); // Note this is different from the "global" geometry

    // Magnetic Field
    ESHandle< MagneticField > magneticFieldHandle;
    iSetup.get< IdealMagneticFieldRecord >().get(magneticFieldHandle);
    const MagneticField* theMagneticField = magneticFieldHandle.product();
    double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();

    // TrackingParticles
    Handle< vector< TrackingParticle > > TrackingParticleHandle;
    iEvent.getByLabel("mix", "MergedTrackTruth", TrackingParticleHandle);
    Handle< vector< TrackingVertex > > TrackingVertexHandle;
    iEvent.getByLabel("mix", "MergedTrackTruth", TrackingVertexHandle);

    // Track Trigger
    Handle< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > > > PixelDigiTTClusterHandle;
    Handle< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > > > PixelDigiTTClusterInclusiveHandle;
    Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > PixelDigiTTStubHandle;
    // NOTE: the InputTag for the "Accepted" clusters is different from the "Inclusive" one
    iEvent.getByLabel("TTStubsFromPixelDigis", "ClusterAccepted", PixelDigiTTClusterHandle);
    iEvent.getByLabel("TTClustersFromPixelDigis", "ClusterInclusive", PixelDigiTTClusterInclusiveHandle);
    iEvent.getByLabel("TTStubsFromPixelDigis", "StubAccepted", PixelDigiTTStubHandle);

    // Track Trigger MC Truth
    Handle< TTClusterAssociationMap< Ref_PixelDigi_ > > MCTruthTTClusterHandle;
    Handle< TTStubAssociationMap< Ref_PixelDigi_ > > MCTruthTTStubHandle;
    iEvent.getByLabel("TTClusterAssociatorFromPixelDigis", "ClusterAccepted", MCTruthTTClusterHandle);
    iEvent.getByLabel("TTStubAssociatorFromPixelDigis", "StubAccepted", MCTruthTTStubHandle);

    //////////////////////////////
    // COLLECT STUB INFORMATION //
    //////////////////////////////

    // Eta coverage
    // Go on only if there are TrackingParticles
    if (TrackingParticleHandle->size() > 0) {
        // Loop over TrackingParticles
        unsigned int tpCnt = 0;
        for (vector< TrackingParticle >::const_iterator iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
            // Make the pointer
            Ptr< TrackingParticle > tempTPPtr(TrackingParticleHandle, tpCnt++);

            // Search the cluster MC map
            vector< Ref< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >, TTCluster< Ref_PixelDigi_ > > > theseClusters = MCTruthTTClusterHandle->findTTClusterRefs(tempTPPtr);

            if (theseClusters.size() > 0) {
                bool normIClu = false;
                bool normOClu = false;

                // Loop over the Clusters
                for (unsigned int jc = 0; jc < theseClusters.size(); jc++) {
                    // Check if it is good
                    bool genuineClu = MCTruthTTClusterHandle->isGenuine(theseClusters.at(jc));
                    if (!genuineClu) continue;

                    unsigned int stackMember = theseClusters.at(jc)->getStackMember();
                    unsigned int clusterWidth = theseClusters.at(jc)->findWidth();

                    if (stackMember == 0) {
                        if (normIClu == false) {
                            hTPart_Eta_INormalization->Fill(fabs(tempTPPtr->momentum().eta()));
                            normIClu = true;
                        }
                        if (clusterWidth == 1) hTPart_Eta_ICW_1->Fill(fabs(tempTPPtr->momentum().eta()));
                        else if (clusterWidth == 2) hTPart_Eta_ICW_2->Fill(fabs(tempTPPtr->momentum().eta()));
                        else hTPart_Eta_ICW_3->Fill(fabs(tempTPPtr->momentum().eta()));
                    }
                    else if (stackMember == 1) {
                        if (normOClu == false) {
                            hTPart_Eta_ONormalization->Fill(fabs(tempTPPtr->momentum().eta()));
                            normOClu = true;
                        }

                        if (clusterWidth == 1) hTPart_Eta_OCW_1->Fill(fabs(tempTPPtr->momentum().eta()));
                        else if (clusterWidth == 2) hTPart_Eta_OCW_2->Fill(fabs(tempTPPtr->momentum().eta()));
                        else hTPart_Eta_OCW_3->Fill(fabs(tempTPPtr->momentum().eta()));
                    }
                } // End of loop over clusters
            }

            // Search the stub MC truth map
            vector< Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > theseStubs = MCTruthTTStubHandle->findTTStubRefs(tempTPPtr);

            if (tempTPPtr->p4().pt() <= 10) continue;

            if (theseStubs.size() > 0) {
                bool normStub = false;

                // Loop over the Stubs
                for (unsigned int js = 0; js < theseStubs.size(); js++) {
                    // Check if it is good
                    bool genuineStub = MCTruthTTStubHandle->isGenuine(theseStubs.at(js));
                    if (!genuineStub) continue;

                    if (normStub == false) {
                        hTPart_Eta_Pt10_Normalization->Fill(fabs(tempTPPtr->momentum().eta()));
                        normStub = true;
                    }

                    // Classify the stub
                    StackedTrackerDetId stDetId(theseStubs.at(js)->getDetId());
                    // Check if there are PS modules in seed or candidate
                    const GeomDetUnit* det0 = theStackedGeometry->idToDetUnit(stDetId, 0);
                    const GeomDetUnit* det1 = theStackedGeometry->idToDetUnit(stDetId, 1);
                    // Find pixel pitch and topology related information
                    const PixelGeomDetUnit* pix0 = dynamic_cast< const PixelGeomDetUnit* >(det0);
                    const PixelGeomDetUnit* pix1 = dynamic_cast< const PixelGeomDetUnit* >(det1);
                    const PixelTopology* top0 = dynamic_cast< const PixelTopology* >(&(pix0->specificTopology()));
                    const PixelTopology* top1 = dynamic_cast< const PixelTopology* >(&(pix1->specificTopology()));
                    int cols0 = top0->ncolumns();
                    int cols1 = top1->ncolumns();
                    int ratio = cols0/cols1; // This assumes the ratio is integer!

                    // 2S Modules
                    if (ratio == 1) hTPart_Eta_Pt10_Num2S->Fill(fabs(tempTPPtr->momentum().eta()));
                    // PS
                    else hTPart_Eta_Pt10_NumPS->Fill(fabs(tempTPPtr->momentum().eta()));
                } // End of loop over the Stubs generated by this TrackingParticle
            }
        } // End of loop over TrackingParticles
    }

    // Maps to store TrackingParticle information
    map< unsigned int, vector< Ptr< TrackingParticle > > > tpPerLayer;
    map< unsigned int, vector< Ptr< TrackingParticle > > > tpPerDisk;

    ///////////////////////
    //                   //
    // Interesting STUFF //
    //                   //
    ///////////////////////

    // Loop over the input Clusters
    typename edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >::const_iterator inputIter;
    typename edmNew::DetSet< TTCluster< Ref_PixelDigi_ > >::const_iterator contentIter;

    for (inputIter = PixelDigiTTClusterHandle->begin(); inputIter != PixelDigiTTClusterHandle->end(); ++inputIter) {
        for (contentIter = inputIter->begin(); contentIter != inputIter->end(); ++contentIter) {
            // Make the reference to be put in the map
            Ref< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >, TTCluster< Ref_PixelDigi_ > > tempCluRef = edmNew::makeRefTo(PixelDigiTTClusterHandle, contentIter);

            StackedTrackerDetId detIdClu(tempCluRef->getDetId());
            unsigned int memberClu = tempCluRef->getStackMember();
            bool genuineClu = MCTruthTTClusterHandle->isGenuine(tempCluRef);
            bool combinClu = MCTruthTTClusterHandle->isCombinatoric(tempCluRef);
            //bool unknownClu = MCTruthTTClusterHandle->isUnknown(tempCluRef);
            int partClu = 999999999;

            if (genuineClu) {
                Ptr< TrackingParticle > thisTP = MCTruthTTClusterHandle->findTrackingParticlePtr(tempCluRef);
                partClu = thisTP->pdgId();
            }

            unsigned int widClu = tempCluRef->findWidth();
            GlobalPoint posClu  = theStackedGeometry->findAverageGlobalPosition(&(*tempCluRef));
            MeasurementPoint locClu = tempCluRef->findAverageLocalCoordinates();

            hCluster_RZ->Fill(posClu.z(), posClu.perp());

            if (detIdClu.isBarrel()) {
                if (memberClu == 0) hCluster_IMem_Barrel->Fill(detIdClu.iLayer());
                else hCluster_OMem_Barrel->Fill(detIdClu.iLayer());

                if (genuineClu) hCluster_Gen_Barrel->Fill(detIdClu.iLayer());
                else if (combinClu) hCluster_Comb_Barrel->Fill(detIdClu.iLayer());
                else hCluster_Unkn_Barrel->Fill(detIdClu.iLayer());

                hCluster_Barrel_XY->Fill(posClu.x(), posClu.y());
                hCluster_Barrel_XY_Zoom->Fill(posClu.x(), posClu.y());

                mapLocalPositionBarrel[detIdClu.iLayer()]->Fill(locClu.x(), locClu.y());
                mapLocalPositionBarrel_Pixel[detIdClu.iLayer()]->Fill(locClu.x(), locClu.y());
                mapLocalPositionBarrel_Strip[detIdClu.iLayer()]->Fill(locClu.x(), locClu.y());

                mapCluSizeBarrel[detIdClu.iLayer()]->Fill(widClu);
                mapCluSizeBarrel_Pixel[detIdClu.iLayer()]->Fill(widClu);
                mapCluSizeBarrel_Strip[detIdClu.iLayer()]->Fill(widClu);
            }
            else if (detIdClu.isEndcap()) {
                if (memberClu == 0) hCluster_IMem_Endcap->Fill(detIdClu.iDisk());
                else hCluster_OMem_Endcap->Fill(detIdClu.iDisk());

                if (genuineClu) hCluster_Gen_Endcap->Fill(detIdClu.iDisk());
                else if (combinClu) hCluster_Comb_Endcap->Fill(detIdClu.iDisk());
                else hCluster_Unkn_Endcap->Fill(detIdClu.iDisk());

                if (posClu.z() > 0) {
                    hCluster_Endcap_Fw_XY->Fill(posClu.x(), posClu.y());
                    hCluster_Endcap_Fw_RZ_Zoom->Fill(posClu.z(), posClu.perp());
                }
                else {
                    hCluster_Endcap_Bw_XY->Fill(posClu.x(), posClu.y());
                    hCluster_Endcap_Bw_RZ_Zoom->Fill(posClu.z(), posClu.perp());
                }

                mapLocalPositionEndcap[detIdClu.iDisk()]->Fill(locClu.x(), locClu.y());
                mapLocalPositionEndcap_Pixel[detIdClu.iDisk()]->Fill(locClu.x(), locClu.y());
                mapLocalPositionEndcap_Strip[detIdClu.iDisk()]->Fill(locClu.x(), locClu.y());

                mapCluSizeEndcap[detIdClu.iDisk()]->Fill(widClu);
                mapCluSizeEndcap_Pixel[detIdClu.iDisk()]->Fill(widClu);
                mapCluSizeEndcap_Strip[detIdClu.iDisk()]->Fill(widClu);
            }

            // Another way of looking at MC truth
            if (genuineClu) hCluster_Gen_Eta->Fill(fabs(posClu.eta()));
            else if (combinClu) hCluster_Comb_Eta->Fill(fabs(posClu.eta()));
            else hCluster_Unkn_Eta->Fill(fabs(posClu.eta()));

            hCluster_PID->Fill(partClu, memberClu);
            hCluster_W->Fill(widClu, memberClu);

            // Store Track information in maps, skip if the Cluster is not good
            if (!genuineClu && !combinClu) continue;

            vector< Ptr< TrackingParticle > > theseTPs = MCTruthTTClusterHandle->findTrackingParticlePtrs(tempCluRef);

            for (unsigned int i = 0; i < theseTPs.size(); i++) {
                Ptr< TrackingParticle > tpPtr = theseTPs.at(i);

                if (tpPtr.isNull()) continue;

                // Get the corresponding vertex and reject the track
                // if its vertex is outside the beampipe
                if (tpPtr->vertex().rho() >= 2.) continue;

                if (detIdClu.isBarrel()) {
                    if (tpPerLayer.find(detIdClu.iLayer()) == tpPerLayer.end()) {
                        vector< Ptr< TrackingParticle > > tempVec;
                        tpPerLayer.insert(make_pair(detIdClu.iLayer(), tempVec));
                    }
                    tpPerLayer[detIdClu.iLayer()].push_back(tpPtr);
                }
                else if (detIdClu.isEndcap()) {
                    if (tpPerDisk.find(detIdClu.iDisk()) == tpPerDisk.end()) {
                        vector< Ptr< TrackingParticle > > tempVec;
                        tpPerDisk.insert(make_pair(detIdClu.iDisk(), tempVec));
                    }
                    tpPerDisk[detIdClu.iDisk()].push_back(tpPtr);
                }
            }
        }
    } // End of Loop over TTClusters

    // Clean the maps for TrackingParticles and fill histograms
    map< unsigned int, vector< Ptr< TrackingParticle > > >::iterator iterTPPerLayer;
    map< unsigned int, vector< Ptr< TrackingParticle > > >::iterator iterTPPerDisk;

    for (iterTPPerLayer = tpPerLayer.begin(); iterTPPerLayer != tpPerLayer.end(); ++iterTPPerLayer) {
        // Remove duplicates, if any
        vector< Ptr< TrackingParticle > > tempVec = iterTPPerLayer->second;
        sort(tempVec.begin(), tempVec.end());
        tempVec.erase(unique(tempVec.begin(), tempVec.end()), tempVec.end());

        // Loop over the TrackingParticles in this piece of the map
        for (unsigned int i = 0; i < tempVec.size(); i++) {
            if (tempVec.at(i).isNull()) continue;
            TrackingParticle thisTP = *(tempVec.at(i));
            mapCluLayer_hTPart_Pt[iterTPPerLayer->first]->Fill(thisTP.p4().pt());

            if (thisTP.p4().pt() > 10.) {
                mapCluLayer_hTPart_Eta_Pt10[iterTPPerLayer->first]->Fill(thisTP.momentum().eta());
                mapCluLayer_hTPart_Phi_Pt10[iterTPPerLayer->first]->Fill(thisTP.momentum().phi() > M_PI ? thisTP.momentum().phi() - 2. * M_PI : thisTP.momentum().phi());
            }
        }
    }

    for (iterTPPerDisk = tpPerDisk.begin(); iterTPPerDisk != tpPerDisk.end(); ++iterTPPerDisk) {
        // Remove duplicates, if any
        vector< Ptr< TrackingParticle > > tempVec = iterTPPerDisk->second;
        sort(tempVec.begin(), tempVec.end());
        tempVec.erase(unique(tempVec.begin(), tempVec.end()), tempVec.end());

        // Loop over the TrackingParticles in this piece of the map
        for (unsigned int i = 0; i < tempVec.size(); i++) {
            if (tempVec.at(i).isNull()) continue;
            TrackingParticle thisTP = *(tempVec.at(i));
            mapCluDisk_hTPart_Pt[iterTPPerDisk->first]->Fill(thisTP.p4().pt());

            if (thisTP.p4().pt() > 10.) {
                mapCluDisk_hTPart_Eta_Pt10[iterTPPerDisk->first]->Fill(thisTP.momentum().eta());
                mapCluDisk_hTPart_Phi_Pt10[iterTPPerDisk->first]->Fill(thisTP.momentum().phi() > M_PI ? thisTP.momentum().phi() - 2. * M_PI : thisTP.momentum().phi());
            }
        }
    }

    //////////////////////////////
    // COLLECT STUB INFORMATION //
    //////////////////////////////

    // Maps to store TrackingParticle information
    map< unsigned int, vector< Ptr< TrackingParticle > > > tpPerStubLayer;
    map< unsigned int, vector< Ptr< TrackingParticle > > > tpPerStubDisk;

    // Loop over the input Stubs
    typename edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::const_iterator otherInputIter;
    typename edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::const_iterator otherContentIter;

    for (otherInputIter = PixelDigiTTStubHandle->begin(); otherInputIter != PixelDigiTTStubHandle->end(); ++otherInputIter) {
        for (otherContentIter = otherInputIter->begin(); otherContentIter != otherInputIter->end(); ++otherContentIter) {
            // Make the reference to be put in the map
            Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > tempStubRef = edmNew::makeRefTo(PixelDigiTTStubHandle, otherContentIter);

            StackedTrackerDetId detIdStub(tempStubRef->getDetId());

            bool genuineStub = MCTruthTTStubHandle->isGenuine(tempStubRef);
            bool combinStub = MCTruthTTStubHandle->isCombinatoric(tempStubRef);
            //bool unknownStub = MCTruthTTStubHandle->isUnknown(tempStubRef);
            int partStub = 999999999;

            if (genuineStub) {
                Ptr< TrackingParticle > thisTP = MCTruthTTStubHandle->findTrackingParticlePtr(tempStubRef);
                partStub = thisTP->pdgId();
            }
            double displStub = tempStubRef->getTriggerDisplacement();
            double offsetStub = tempStubRef->getTriggerOffset();
            GlobalPoint posStub = theStackedGeometry->findGlobalPosition(&(*tempStubRef));

            hStub_RZ->Fill(posStub.z(), posStub.perp());

            if (detIdStub.isBarrel()) {
                hStub_Barrel->Fill(detIdStub.iLayer());

                if (genuineStub) hStub_Gen_Barrel->Fill(detIdStub.iLayer());
                else if (combinStub) hStub_Comb_Barrel->Fill(detIdStub.iLayer());
                else hStub_Unkn_Barrel->Fill(detIdStub.iLayer());

                hStub_Barrel_XY->Fill(posStub.x(), posStub.y());
                hStub_Barrel_XY_Zoom->Fill(posStub.x(), posStub.y());
            }
            else if (detIdStub.isEndcap()) {
                hStub_Endcap->Fill(detIdStub.iDisk());

                if (genuineStub) hStub_Gen_Endcap->Fill(detIdStub.iDisk());
                else if (combinStub) hStub_Comb_Endcap->Fill(detIdStub.iDisk());
                else hStub_Unkn_Endcap->Fill(detIdStub.iDisk());

                if (posStub.z() > 0) {
                    hStub_Endcap_Fw_XY->Fill(posStub.x(), posStub.y());
                    hStub_Endcap_Fw_RZ_Zoom->Fill(posStub.z(), posStub.perp());
                }
                else {
                    hStub_Endcap_Bw_XY->Fill(posStub.x(), posStub.y());
                    hStub_Endcap_Bw_RZ_Zoom->Fill(posStub.z(), posStub.perp());
                }
            }

            // Another way of looking at MC truth
            if (genuineStub) hStub_Gen_Eta->Fill(fabs(posStub.eta()));
            else if (combinStub) hStub_Comb_Eta->Fill(fabs(posStub.eta()));
            else hStub_Unkn_Eta->Fill(fabs(posStub.eta()));

            hStub_PID->Fill(partStub);

            // Store Track information in maps, skip if the Cluster is not good
            if (!genuineStub) continue;

            Ptr< TrackingParticle > tpPtr = MCTruthTTStubHandle->findTrackingParticlePtr(tempStubRef);

            // Get the corresponding vertex and reject the track
            // if its vertex is outside the beampipe
            if (tpPtr->vertex().rho() >= 2.0) continue;

            if (detIdStub.isBarrel()) {
                if (tpPerStubLayer.find(detIdStub.iLayer()) == tpPerStubLayer.end()) {
                    vector< Ptr< TrackingParticle > > tempVec;
                    tpPerStubLayer.insert(make_pair(detIdStub.iLayer(), tempVec));
                }

                tpPerStubLayer[detIdStub.iLayer()].push_back(tpPtr);

                hStub_Barrel_W->Fill(detIdStub.iLayer(), displStub - offsetStub);
                hStub_Barrel_O->Fill(detIdStub.iLayer(), offsetStub);
            }
            else if (detIdStub.isEndcap()) {
                if (tpPerStubDisk.find(detIdStub.iDisk()) == tpPerStubDisk.end()) {
                    vector< Ptr< TrackingParticle > > tempVec;
                    tpPerStubDisk.insert(make_pair(detIdStub.iDisk(), tempVec));
                }

                tpPerStubDisk[detIdStub.iDisk()].push_back(tpPtr);

                hStub_Endcap_W->Fill(detIdStub.iDisk(), displStub - offsetStub);
                hStub_Endcap_O->Fill(detIdStub.iDisk(), offsetStub);
            }

            // Compare to TrackingParticle

            if (tpPtr.isNull()) continue; // This prevents to fill the vector if the TrackingParticle is not found

            TrackingParticle thisTP = *tpPtr;

            double simPt = thisTP.p4().pt();
            double simEta = thisTP.momentum().eta();
            double simPhi = thisTP.momentum().phi();
            double recPt = theStackedGeometry->findRoughPt(mMagneticFieldStrength, &(*tempStubRef));
            double recEta = theStackedGeometry->findGlobalDirection(&(*tempStubRef)).eta();
            double recPhi = theStackedGeometry->findGlobalDirection(&(*tempStubRef)).phi();

            if (simPhi > M_PI) simPhi -= 2*M_PI;
            if (recPhi > M_PI) recPhi -= 2*M_PI;

            if (detIdStub.isBarrel()) {
                mapStubLayer_hStub_InvPt_TPart_InvPt[detIdStub.iLayer()]->Fill(1./simPt, 1./recPt);
                mapStubLayer_hStub_Pt_TPart_Pt[detIdStub.iLayer()]->Fill(simPt, recPt);
                mapStubLayer_hStub_Eta_TPart_Eta[detIdStub.iLayer()]->Fill(simEta, recEta);
                mapStubLayer_hStub_Phi_TPart_Phi[detIdStub.iLayer()]->Fill(simPhi, recPhi);

                mapStubLayer_hStub_InvPtRes_TPart_Eta[detIdStub.iLayer()]->Fill(simEta, 1./recPt - 1./simPt);
                mapStubLayer_hStub_PtRes_TPart_Eta[detIdStub.iLayer()]->Fill(simEta, recPt - simPt);
                mapStubLayer_hStub_EtaRes_TPart_Eta[detIdStub.iLayer()]->Fill(simEta, recEta - simEta);
                mapStubLayer_hStub_PhiRes_TPart_Eta[detIdStub.iLayer()]->Fill(simEta, recPhi - simPhi);

                mapStubLayer_hStub_W_TPart_Pt[detIdStub.iLayer()]->Fill(simPt, displStub - offsetStub);
                mapStubLayer_hStub_W_TPart_InvPt[detIdStub.iLayer()]->Fill(1./simPt, displStub - offsetStub);
            }
            else if (detIdStub.isEndcap()) {
                mapStubDisk_hStub_InvPt_TPart_InvPt[detIdStub.iDisk()]->Fill(1./simPt, 1./recPt);
                mapStubDisk_hStub_Pt_TPart_Pt[detIdStub.iDisk()]->Fill(simPt, recPt);
                mapStubDisk_hStub_Eta_TPart_Eta[detIdStub.iDisk()]->Fill(simEta, recEta);
                mapStubDisk_hStub_Phi_TPart_Phi[detIdStub.iDisk()]->Fill(simPhi, recPhi);

                mapStubDisk_hStub_InvPtRes_TPart_Eta[detIdStub.iDisk()]->Fill(simEta, 1./recPt - 1./simPt);
                mapStubDisk_hStub_PtRes_TPart_Eta[detIdStub.iDisk()]->Fill(simEta, recPt - simPt);
                mapStubDisk_hStub_EtaRes_TPart_Eta[detIdStub.iDisk()]->Fill(simEta, recEta - simEta);
                mapStubDisk_hStub_PhiRes_TPart_Eta[detIdStub.iDisk()]->Fill(simEta, recPhi - simPhi);

                mapStubDisk_hStub_W_TPart_Pt[detIdStub.iDisk()]->Fill(simPt, displStub - offsetStub);
                mapStubDisk_hStub_W_TPart_InvPt[detIdStub.iDisk()]->Fill(1./simPt, displStub - offsetStub);
            }
        }
    } // End of loop over TTStubs

    // Clean the maps for TrackingParticles and fill histograms
    map< unsigned int, vector< Ptr< TrackingParticle > > >::iterator iterTPPerStubLayer;
    map< unsigned int, vector< Ptr< TrackingParticle > > >::iterator iterTPPerStubDisk;

    for (iterTPPerStubLayer = tpPerStubLayer.begin(); iterTPPerStubLayer != tpPerStubLayer.end(); ++iterTPPerStubLayer) {
        // Remove duplicates, if any
        vector< Ptr< TrackingParticle > > tempVec = iterTPPerStubLayer->second;
        sort(tempVec.begin(), tempVec.end());
        tempVec.erase(unique(tempVec.begin(), tempVec.end()), tempVec.end());

        // Loop over the TrackingParticles in this piece of the map
        for (unsigned int i = 0; i < tempVec.size(); i++) {
            if (tempVec.at(i).isNull()) continue;
            TrackingParticle thisTP = *(tempVec.at(i));
            mapStubLayer_hTPart_Pt[iterTPPerStubLayer->first]->Fill(thisTP.p4().pt());
            if (thisTP.p4().pt() > 10.) {
                mapStubLayer_hTPart_Eta_Pt10[iterTPPerStubLayer->first]->Fill(thisTP.momentum().eta());
                mapStubLayer_hTPart_Phi_Pt10[iterTPPerStubLayer->first]->Fill(thisTP.momentum().phi() > M_PI ? thisTP.momentum().phi() - 2. * M_PI : thisTP.momentum().phi());
            }
        }
    }

    for (iterTPPerStubDisk = tpPerStubDisk.begin(); iterTPPerStubDisk != tpPerStubDisk.end(); ++iterTPPerStubDisk) {
        // Remove duplicates, if any
        vector< Ptr< TrackingParticle > > tempVec = iterTPPerStubDisk->second;
        sort(tempVec.begin(), tempVec.end());
        tempVec.erase(unique(tempVec.begin(), tempVec.end()), tempVec.end());

        // Loop over the TrackingParticles in this piece of the map
        for (unsigned int i = 0; i < tempVec.size(); i++) {
            if (tempVec.at(i).isNull()) continue;
            TrackingParticle thisTP = *(tempVec.at(i));
            mapStubDisk_hTPart_Pt[iterTPPerStubDisk->first]->Fill(thisTP.p4().pt());
            if (thisTP.p4().pt() > 10.) {
                mapStubDisk_hTPart_Eta_Pt10[iterTPPerStubDisk->first]->Fill(thisTP.momentum().eta());
                mapStubDisk_hTPart_Phi_Pt10[iterTPPerStubDisk->first]->Fill(thisTP.momentum().phi() > M_PI ? thisTP.momentum().phi() - 2. * M_PI : thisTP.momentum().phi());
            }
        }
    }

    ////////////////////////////
    // SPECTRUM OF SIM TRACKS //
    // WITHIN PRIMARY VERTEX  //
    // CONSTRAINTS            //
    ///////////////////////////

    // Go on only if there are TrackingParticles
    if (TrackingParticleHandle->size() != 0) {
        // Loop over TrackingParticles
        for (vector<TrackingParticle>::const_iterator iterTrackingParticles = TrackingParticleHandle->begin(); iterTrackingParticles != TrackingParticleHandle->end(); ++iterTrackingParticles) {
            // Get the corresponding vertex
            // Assume perfectly round beamspot
            // Correct and get the correct TrackingParticle Vertex position wrt beam center
            if (iterTrackingParticles->vertex().rho() >= 2) continue;

            // First of all, check beamspot and correction
            hSimVtx_XY->Fill(iterTrackingParticles->vertex().x(), iterTrackingParticles->vertex().y());
            hSimVtx_RZ->Fill(iterTrackingParticles->vertex().z(), iterTrackingParticles->vertex().rho());

            // Here we have only tracks form primary vertices
            // Check Pt spectrum and pseudorapidity for over-threshold tracks
            hTPart_Pt->Fill(iterTrackingParticles->p4().pt());

            if (iterTrackingParticles->p4().pt() > 10.) {
                hTPart_Eta_Pt10->Fill(iterTrackingParticles->momentum().eta());
                hTPart_Phi_Pt10->Fill(iterTrackingParticles->momentum().phi() > M_PI ? iterTrackingParticles->momentum().phi() - 2. * M_PI : iterTrackingParticles->momentum().phi());
            }
        } // End of Loop over TrackingParticles
    } // End of if (TrackingParticleHandle->size() != 0)

} // End of analyze()

void RunStepsL1TriggerValidation::endJob() {
    cerr << " RunStepsL1TriggerValidation::endJob" << endl;
}

void RunStepsL1TriggerValidation::beginJob() {
    testedGeometry = false;

    ostringstream histoName;
    ostringstream histoTitle;

    // Things to be done before entering the event Loop
    cerr << " RunStepsL1TriggerValidation::beginJob" << endl;

    // Book histograms etc
    Service<TFileService> fs;

    // Prepare for LogXY Plots
    int NumBins = 200;
    double MinPt = 0.;
    double MaxPt = 100.;

    double* BinVec = new double[NumBins+1];
    for (int iBin = 0; iBin < NumBins + 1; iBin++) {
        double temp = pow(10, (- NumBins + iBin) / (MaxPt - MinPt));
        BinVec[iBin] = temp;
    }

    fs->file().cd("/");
    TFileDirectory tdCommonDir = fs->mkdir("Common");

    // TrackingParticle and TrackingVertex
    hSimVtx_XY = tdCommonDir.make<TH2D>("hSimVtx_XY", "SimVtx y vs. x", 200, -0.4, 0.4, 200, -0.4, 0.4);
    hSimVtx_RZ = tdCommonDir.make<TH2D>("hSimVtx_RZ", "SimVtx #rho vs. z", 200, -50, 50, 200, 0, 0.4);
    hSimVtx_XY->Sumw2();
    hSimVtx_RZ->Sumw2();

    hTPart_Pt = tdCommonDir.make<TH1D>("hTPart_Pt", "TPart p_{T}", 200, 0, 50);
    hTPart_Eta_Pt10 = tdCommonDir.make<TH1D>("hTPart_Eta_Pt10", "TPart #eta (p_{T} > 10 GeV/c)", 180, -M_PI, M_PI);
    hTPart_Phi_Pt10 = tdCommonDir.make<TH1D>("hTPart_Phi_Pt10", "TPart #phi (p_{T} > 10 GeV/c)", 180, -M_PI, M_PI);
    hTPart_Pt->Sumw2();
    hTPart_Eta_Pt10->Sumw2();
    hTPart_Phi_Pt10->Sumw2();

    // Global position of TTCluster
    hCluster_Barrel_XY = tdCommonDir.make<TH2D>("hCluster_Barrel_XY", "TTCluster Barrel y vs. x", 960, -120, 120, 960, -120, 120);
    hCluster_Barrel_XY_Zoom = tdCommonDir.make<TH2D>("hCluster_Barrel_XY_Zoom", "TTCluster Barrel y vs. x", 960, 30, 60, 960, -15, 15);
    hCluster_Endcap_Fw_XY = tdCommonDir.make<TH2D>("hCluster_Endcap_Fw_XY", "TTCluster Forward Endcap y vs. x", 960, -120, 120, 960, -120, 120);
    hCluster_Endcap_Bw_XY = tdCommonDir.make<TH2D>("hCluster_Endcap_Bw_XY", "TTCluster Backward Endcap y vs. x", 960, -120, 120, 960, -120, 120);
    hCluster_RZ = tdCommonDir.make<TH2D>("hCluster_RZ", "TTCluster #rho vs. z", 900, -300, 300, 480, 0, 120);
    hCluster_Endcap_Fw_RZ_Zoom = tdCommonDir.make<TH2D>("hCluster_Endcap_Fw_RZ_Zoom", "TTCluster Forward Endcap #rho vs. z", 960, 140, 170, 960, 30, 60);
    hCluster_Endcap_Bw_RZ_Zoom = tdCommonDir.make<TH2D>("hCluster_Endcap_Bw_RZ_Zoom", "TTCluster Backward Endcap #rho vs. z", 960, -170, -140, 960, 70, 100);
    hCluster_Barrel_XY->Sumw2();
    hCluster_Barrel_XY_Zoom->Sumw2();
    hCluster_Endcap_Fw_XY->Sumw2();
    hCluster_Endcap_Bw_XY->Sumw2();
    hCluster_RZ->Sumw2();
    hCluster_Endcap_Fw_RZ_Zoom->Sumw2();
    hCluster_Endcap_Bw_RZ_Zoom->Sumw2();

    hCluster_IMem_Barrel = tdCommonDir.make<TH1D>("hCluster_IMem_Barrel", "Inner TTCluster Stack", 12, -0.5, 11.5);
    hCluster_IMem_Endcap = tdCommonDir.make<TH1D>("hCluster_IMem_Endcap", "Inner TTCluster Stack", 12, -0.5, 11.5);
    hCluster_OMem_Barrel = tdCommonDir.make<TH1D>("hCluster_OMem_Barrel", "Outer TTCluster Stack", 12, -0.5, 11.5);
    hCluster_OMem_Endcap = tdCommonDir.make<TH1D>("hCluster_OMem_Endcap", "Outer TTCluster Stack", 12, -0.5, 11.5);
    hCluster_IMem_Barrel->Sumw2();
    hCluster_IMem_Endcap->Sumw2();
    hCluster_OMem_Barrel->Sumw2();
    hCluster_OMem_Endcap->Sumw2();

    hCluster_Gen_Barrel = tdCommonDir.make<TH1D>("hCluster_Gen_Barrel", "Genuine TTCluster Stack", 12, -0.5, 11.5);
    hCluster_Unkn_Barrel = tdCommonDir.make<TH1D>("hCluster_Unkn_Barrel", "Unknown TTCluster Stack", 12, -0.5, 11.5);
    hCluster_Comb_Barrel = tdCommonDir.make<TH1D>("hCluster_Comb_Barrel", "Combinatorial TTCluster Stack", 12, -0.5, 11.5);
    hCluster_Gen_Endcap = tdCommonDir.make<TH1D>("hCluster_Gen_Endcap", "Genuine TTCluster Stack", 12, -0.5, 11.5);
    hCluster_Unkn_Endcap = tdCommonDir.make<TH1D>("hCluster_Unkn_Endcap", "Unknown TTCluster Stack", 12, -0.5, 11.5);
    hCluster_Comb_Endcap = tdCommonDir.make<TH1D>("hCluster_Comb_Endcap", "Combinatorial TTCluster Stack", 12, -0.5, 11.5);
    hCluster_Gen_Barrel->Sumw2();
    hCluster_Unkn_Barrel->Sumw2();
    hCluster_Comb_Barrel->Sumw2();
    hCluster_Gen_Endcap->Sumw2();
    hCluster_Unkn_Endcap->Sumw2();
    hCluster_Comb_Endcap->Sumw2();

    hCluster_Gen_Eta = tdCommonDir.make<TH1D>("hCluster_Gen_Eta", "Genuine TTCluster #eta", 90, 0, M_PI);
    hCluster_Unkn_Eta = tdCommonDir.make<TH1D>("hCluster_Unkn_Eta", "Unknown TTCluster #eta", 90, 0, M_PI);
    hCluster_Comb_Eta = tdCommonDir.make<TH1D>("hCluster_Comb_Eta", "Combinatorial TTCluster #eta", 90, 0, M_PI);
    hCluster_Gen_Eta->Sumw2();
    hCluster_Unkn_Eta->Sumw2();
    hCluster_Comb_Eta->Sumw2();

    hCluster_PID = tdCommonDir.make<TH2D>("hCluster_PID", "TTCluster PID (Member)", 501, -250.5, 250.5, 2, -0.5, 1.5);
    hCluster_W = tdCommonDir.make<TH2D>("hCluster_W", "TTCluster Width (Member)", 10, -0.5, 9.5, 2, -0.5, 1.5);
    hCluster_PID->Sumw2();
    hCluster_W->Sumw2();

    hTPart_Eta_INormalization = tdCommonDir.make<TH1D>("hTPart_Eta_INormalization", "TParticles vs. TPart #eta", 90, 0, M_PI);
    hTPart_Eta_ICW_1 = tdCommonDir.make<TH1D>("hTPart_Eta_ICW_1", "CW 1 vs. TPart #eta", 90, 0, M_PI);
    hTPart_Eta_ICW_2 = tdCommonDir.make<TH1D>("hTPart_Eta_ICW_2", "CW 2 vs. TPart #eta", 90, 0, M_PI);
    hTPart_Eta_ICW_3 = tdCommonDir.make<TH1D>("hTPart_Eta_ICW_3", "CW 3 or more vs. TPart #eta", 90, 0, M_PI);
    hTPart_Eta_INormalization->Sumw2();
    hTPart_Eta_ICW_1->Sumw2();
    hTPart_Eta_ICW_2->Sumw2();
    hTPart_Eta_ICW_3->Sumw2();

    hTPart_Eta_ONormalization = tdCommonDir.make<TH1D>("hTPart_Eta_ONormalization", "TParticles vs. TPart #eta", 90, 0, M_PI);
    hTPart_Eta_OCW_1 = tdCommonDir.make<TH1D>("hTPart_Eta_OCW_1", "CW 1 vs. TPart #eta", 90, 0, M_PI);
    hTPart_Eta_OCW_2 = tdCommonDir.make<TH1D>("hTPart_Eta_OCW_2", "CW 2 vs. TPart #eta", 90, 0, M_PI);
    hTPart_Eta_OCW_3 = tdCommonDir.make<TH1D>("hTPart_Eta_OCW_3", "CW 3 or more vs. TPart #eta", 90, 0, M_PI);
    hTPart_Eta_ONormalization->Sumw2();
    hTPart_Eta_OCW_1->Sumw2();
    hTPart_Eta_OCW_2->Sumw2();
    hTPart_Eta_OCW_3->Sumw2();

    // Global position of TTStub
    hStub_Barrel_XY = tdCommonDir.make<TH2D>("hStub_Barrel_XY", "TTStub Barrel y vs. x", 960, -120, 120, 960, -120, 120);
    hStub_Barrel_XY_Zoom = tdCommonDir.make<TH2D>("hStub_Barrel_XY_Zoom", "TTStub Barrel y vs. x", 960, 30, 60, 960, -15, 15);
    hStub_Endcap_Fw_XY = tdCommonDir.make<TH2D>("hStub_Endcap_Fw_XY", "TTStub Forward Endcap y vs. x", 960, -120, 120, 960, -120, 120);
    hStub_Endcap_Bw_XY = tdCommonDir.make<TH2D>("hStub_Endcap_Bw_XY", "TTStub Backward Endcap y vs. x", 960, -120, 120, 960, -120, 120);
    hStub_RZ = tdCommonDir.make<TH2D>("hStub_RZ", "TTStub #rho vs. z", 900, -300, 300, 480, 0, 120);
    hStub_Endcap_Fw_RZ_Zoom = tdCommonDir.make<TH2D>("hStub_Endcap_Fw_RZ_Zoom", "TTStub Forward Endcap #rho vs. z", 960, 140, 170, 960, 30, 60);
    hStub_Endcap_Bw_RZ_Zoom = tdCommonDir.make<TH2D>("hStub_Endcap_Bw_RZ_Zoom", "TTStub Backward Endcap #rho vs. z", 960, -170, -140, 960, 70, 100);
    hStub_Barrel_XY->Sumw2();
    hStub_Barrel_XY_Zoom->Sumw2();
    hStub_Endcap_Fw_XY->Sumw2();
    hStub_Endcap_Bw_XY->Sumw2();
    hStub_RZ->Sumw2();
    hStub_Endcap_Fw_RZ_Zoom->Sumw2();
    hStub_Endcap_Bw_RZ_Zoom->Sumw2();

    hStub_Barrel = tdCommonDir.make<TH1D>("hStub_Barrel", "TTStub Stack", 12, -0.5, 11.5);
    hStub_Endcap = tdCommonDir.make<TH1D>("hStub_Endcap", "TTStub Stack", 12, -0.5, 11.5);
    hStub_Barrel->Sumw2();
    hStub_Endcap->Sumw2();

    hStub_Gen_Barrel = tdCommonDir.make<TH1D>("hStub_Gen_Barrel", "Genuine TTStub Stack", 12, -0.5, 11.5);
    hStub_Unkn_Barrel = tdCommonDir.make<TH1D>("hStub_Unkn_Barrel", "Unknown  TTStub Stack", 12, -0.5, 11.5);
    hStub_Comb_Barrel = tdCommonDir.make<TH1D>("hStub_Comb_Barrel", "Combinatorial TTStub Stack", 12, -0.5, 11.5);
    hStub_Gen_Endcap  = tdCommonDir.make<TH1D>("hStub_Gen_Endcap", "Genuine TTStub Stack", 12, -0.5, 11.5);
    hStub_Unkn_Endcap = tdCommonDir.make<TH1D>("hStub_Unkn_Endcap", "Unknown  TTStub Stack", 12, -0.5, 11.5);
    hStub_Comb_Endcap = tdCommonDir.make<TH1D>("hStub_Comb_Endcap", "Combinatorial TTStub Stack", 12, -0.5, 11.5);
    hStub_Gen_Barrel->Sumw2();
    hStub_Unkn_Barrel->Sumw2();
    hStub_Comb_Barrel->Sumw2();
    hStub_Gen_Endcap->Sumw2();
    hStub_Unkn_Endcap->Sumw2();
    hStub_Comb_Endcap->Sumw2();

    hStub_Gen_Eta = tdCommonDir.make<TH1D>("hStub_Gen_Eta", "Genuine TTStub #eta", 90, 0, M_PI);
    hStub_Unkn_Eta = tdCommonDir.make<TH1D>("hStub_Unkn_Eta", "Unknown TTStub #eta", 90, 0, M_PI);
    hStub_Comb_Eta = tdCommonDir.make<TH1D>("hStub_Comb_Eta", "Combinatorial TTStub #eta", 90, 0, M_PI);
    hStub_Gen_Eta->Sumw2();
    hStub_Unkn_Eta->Sumw2();
    hStub_Comb_Eta->Sumw2();

    hStub_PID = tdCommonDir.make<TH1D>("hStub_PID", "TTStub PID", 501, -250.5, 250.5);
    hStub_Barrel_W = tdCommonDir.make<TH2D>("hStub_Barrel_W", "TTStub Post-Corr Displacement (Layer)", 12, -0.5, 11.5, 43, -10.75, 10.75);
    hStub_Barrel_O = tdCommonDir.make<TH2D>("hStub_Barrel_O", "TTStub Offset (Layer)", 12, -0.5, 11.5, 43, -10.75, 10.75);
    hStub_Endcap_W = tdCommonDir.make<TH2D>("hStub_Endcap_W", "TTStub Post-Corr Displacement (Layer)", 12, -0.5, 11.5, 43, -10.75, 10.75);
    hStub_Endcap_O = tdCommonDir.make<TH2D>("hStub_Endcap_O", "TTStub Offset (Layer)", 12, -0.5, 11.5, 43, -10.75, 10.75);

    hStub_PID->Sumw2();
    hStub_Barrel_W->Sumw2();
    hStub_Barrel_O->Sumw2();
    hStub_Endcap_W->Sumw2();
    hStub_Endcap_O->Sumw2();

    hTPart_Eta_Pt10_Normalization = tdCommonDir.make<TH1D>("hTPart_Eta_Pt10_Normalization", "TParticles vs. TPart #eta", 90, 0, M_PI);
    hTPart_Eta_Pt10_NumPS = tdCommonDir.make<TH1D>("hTPart_Eta_Pt10_NumPS", "PS Stubs vs. TPart #eta", 90, 0, M_PI);
    hTPart_Eta_Pt10_Num2S = tdCommonDir.make<TH1D>("hTPart_Eta_Pt10_Num2S", "2S Stubs vs. TPart #eta", 90, 0, M_PI);
    hTPart_Eta_Pt10_Normalization->Sumw2();
    hTPart_Eta_Pt10_NumPS->Sumw2();
    hTPart_Eta_Pt10_Num2S->Sumw2();

    // Stub Production Efficiency and comparison to TrackingParticle
    for (unsigned int stackIdx = 0; stackIdx < 12; stackIdx++) {

        ////////// BARREL

        fs->file().cd("/");
        TFileDirectory tdBarrelDir = fs->mkdir("Barrel");

        ostringstream directoryName;
        directoryName << "Layer_" << stackIdx;
        TFileDirectory tdBarrerLayerDir = tdBarrelDir.mkdir(directoryName.str().c_str());

        // Local position
        histoName.str("");
        histoName << "Local_Pos_XY_Barrel" << stackIdx;
        histoTitle.str("");
        histoTitle << "Local_Pos_XY_Barrel" << stackIdx;
        mapLocalPositionBarrel[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 2000, 0., 0., 2000, 0., 0.);

        histoName.str("");
        histoName << "Local_Pos_XY_Barrel_Pixel" << stackIdx;
        histoTitle.str("");
        histoTitle << "Local_Pos_XY_Barrel_Pixel " << stackIdx;
        mapLocalPositionBarrel_Pixel[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 2000, 0., 0., 2000, 0., 0.);

        histoName.str("");
        histoName << "Local_Pos_XY_Barrel_Strip" << stackIdx;
        histoTitle.str("");
        histoTitle << "Local_Pos_XY_Barrel_Strip " << stackIdx;
        mapLocalPositionBarrel_Strip[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 2000, 0., 0., 2000, 0., 0.);

        // Cluster size
        histoName.str("");
        histoName << "ClusterSize_Barrel" << stackIdx;
        histoTitle.str("");
        histoTitle << "ClusterSize Barrel " << stackIdx;
        mapCluSizeBarrel[stackIdx] = tdBarrerLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 50, 0, 50);

        histoName.str("");
        histoName << "ClusterSize_Barrel_Pixel" << stackIdx;
        histoTitle.str("");
        histoTitle << "ClusterSize Barrel Pixel " << stackIdx;
        mapCluSizeBarrel_Pixel[stackIdx] = tdBarrerLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 50, 0, 50);

        histoName.str("");
        histoName << "ClusterSize_Barrel_Strip" << stackIdx;
        histoTitle.str("");
        histoTitle << "ClusterSize Barrel Strip " << stackIdx;
        mapCluSizeBarrel_Strip[stackIdx] = tdBarrerLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 50, 0, 50);

        // Denominators
        histoName.str("");
        histoName << "hTPart_Pt_Clu_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart p_{T}, Cluster, Barrel Stack " << stackIdx;
        mapCluLayer_hTPart_Pt[stackIdx] = tdBarrerLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 200, 0, 50);

        histoName.str("");
        histoName << "hTPart_Eta_Pt10_Clu_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart #eta (p_{T} > 10 GeV/c), Cluster, Barrel Stack " << stackIdx;
        mapCluLayer_hTPart_Eta_Pt10[stackIdx] = tdBarrerLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI);

        histoName.str("");
        histoName << "hTPart_Phi_Pt10_Clu_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart #phi (p_{T} > 10 GeV/c), Cluster, Barrel Stack " << stackIdx;
        mapCluLayer_hTPart_Phi_Pt10[stackIdx] = tdBarrerLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI);

        mapCluLayer_hTPart_Pt[stackIdx]->Sumw2();
        mapCluLayer_hTPart_Eta_Pt10[stackIdx]->Sumw2();
        mapCluLayer_hTPart_Phi_Pt10[stackIdx]->Sumw2();

        // Numerators GeV/c
        histoName.str("");
        histoName << "hTPart_Pt_Stub_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart p_{T}, Stub, Barrel Stack " << stackIdx;
        mapStubLayer_hTPart_Pt[stackIdx] = tdBarrerLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 200, 0, 50);

        histoName.str("");
        histoName << "hTPart_Eta_Pt10_Stub_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart #eta (p_{T} > 10 GeV/c), Stub, Barrel Stack " << stackIdx;
        mapStubLayer_hTPart_Eta_Pt10[stackIdx] = tdBarrerLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI);

        histoName.str("");
        histoName << "hTPart_Phi_Pt10_Stub_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart #phi (p_{T} > 10 GeV/c), Stub, Barrel Stack " << stackIdx;
        mapStubLayer_hTPart_Phi_Pt10[stackIdx] = tdBarrerLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI);

        mapStubLayer_hTPart_Pt[stackIdx]->Sumw2();
        mapStubLayer_hTPart_Eta_Pt10[stackIdx]->Sumw2();
        mapStubLayer_hTPart_Phi_Pt10[stackIdx]->Sumw2();

        // Comparison to TrackingParticle
        histoName.str("");
        histoName << "hStub_InvPt_TPart_InvPt_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub p_{T}^{-1} vs. TPart p_{T}^{-1}, Barrel Stack " << stackIdx;
        mapStubLayer_hStub_InvPt_TPart_InvPt[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 200, 0.0, 0.8, 200, 0.0, 0.8);
        mapStubLayer_hStub_InvPt_TPart_InvPt[stackIdx]->GetXaxis()->Set(NumBins, BinVec);
        mapStubLayer_hStub_InvPt_TPart_InvPt[stackIdx]->GetYaxis()->Set(NumBins, BinVec);
        mapStubLayer_hStub_InvPt_TPart_InvPt[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_Pt_TPart_Pt_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub p_{T} vs. TPart p_{T}, Barrel Stack " << stackIdx;
        mapStubLayer_hStub_Pt_TPart_Pt[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 100, 0, 50, 100, 0, 50);
        mapStubLayer_hStub_Pt_TPart_Pt[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_Eta_TPart_Eta_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub #eta vs. TPart #eta, Barrel Stack " << stackIdx;
        mapStubLayer_hStub_Eta_TPart_Eta[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 180, -M_PI, M_PI);
        mapStubLayer_hStub_Eta_TPart_Eta[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_Phi_TPart_Phi_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub #phi vs. TPart #phi, Barrel Stack " << stackIdx;
        mapStubLayer_hStub_Phi_TPart_Phi[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 180, -M_PI, M_PI);
        mapStubLayer_hStub_Phi_TPart_Phi[stackIdx]->Sumw2();

        // Residuals
        histoName.str("");
        histoName << "hStub_InvPtRes_TPart_Eta_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub p_{T}^{-1} - TPart p_{T}^{-1} vs. TPart #eta, Barrel Stack " << stackIdx;
        mapStubLayer_hStub_InvPtRes_TPart_Eta[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 100, -2.0, 2.0);
        mapStubLayer_hStub_InvPtRes_TPart_Eta[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_PtRes_TPart_Eta_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub p_{T} - TPart p_{T} vs. TPart #eta, Barrel Stack " << stackIdx;
        mapStubLayer_hStub_PtRes_TPart_Eta[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 100, -40, 40);
        mapStubLayer_hStub_PtRes_TPart_Eta[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_EtaRes_TPart_Eta_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub #eta - TPart #eta vs. TPart #eta, Barrel Stack " << stackIdx;
        mapStubLayer_hStub_EtaRes_TPart_Eta[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 100, -2, 2);
        mapStubLayer_hStub_EtaRes_TPart_Eta[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_PhiRes_TPart_Eta_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub #phi - TPart #phi vs. TPart #eta, Barrel Stack " << stackIdx;
        mapStubLayer_hStub_PhiRes_TPart_Eta[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 100, -0.5, 0.5);
        mapStubLayer_hStub_PhiRes_TPart_Eta[stackIdx]->Sumw2();

        // Stub Width vs. Pt
        histoName.str("");
        histoName << "hStub_W_TPart_Pt_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub Width vs. TPart p_{T}, Barrel Stack " << stackIdx;
        mapStubLayer_hStub_W_TPart_Pt[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 200, 0, 50, 41, -10.25, 10.25);
        mapStubLayer_hStub_W_TPart_Pt[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_W_TPart_InvPt_L" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub Width vs. TPart p_{T}^{-1}, Barrel Stack " << stackIdx;
        mapStubLayer_hStub_W_TPart_InvPt[stackIdx] = tdBarrerLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 200, 0, 0.8, 41, -10.25, 10.25);
        mapStubLayer_hStub_W_TPart_InvPt[stackIdx]->GetXaxis()->Set(NumBins, BinVec);
        mapStubLayer_hStub_W_TPart_InvPt[stackIdx]->Sumw2();

        ////////// ENDCAP

        fs->file().cd("/");
        TFileDirectory tdEndcapDir = fs->mkdir("Endcap");

        ostringstream directoryName2;
        directoryName2 << "Layer_" << stackIdx;
        TFileDirectory tdEndcapLayerDir = tdEndcapDir.mkdir(directoryName2.str().c_str());

        // Local position
        histoName.str("");
        histoName << "Local_Pos_XY_Endcap" << stackIdx;
        histoTitle.str("");
        histoTitle << "Local_Pos_XY_Endcap " << stackIdx;
        mapLocalPositionEndcap[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 2000, 0., 0., 2000, 0., 0.);

        histoName.str("");
        histoName << "Local_Pos_XY_Endcap_Pixel" << stackIdx;
        histoTitle.str("");
        histoTitle << "Local_Pos_XY_Endcap_Pixel " << stackIdx;
        mapLocalPositionEndcap_Pixel[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 2000, 0., 0., 2000, 0., 0.);

        histoName.str("");
        histoName << "Local_Pos_XY_Endcap_Strip" << stackIdx;
        histoTitle.str("");
        histoTitle << "Local_Pos_XY_Endcap_Strip " << stackIdx;
        mapLocalPositionEndcap_Strip[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 2000, 0., 0., 2000, 0., 0.);

        // Cluster size
        histoName.str("");
        histoName << "ClusterSize_Endcap" << stackIdx;
        histoTitle.str("");
        histoTitle << "ClusterSize Endcap " << stackIdx;
        mapCluSizeEndcap[stackIdx] = tdEndcapLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 50, 0, 50);

        histoName.str("");
        histoName << "ClusterSize_Endcap_Pixel" << stackIdx;
        histoTitle.str("");
        histoTitle << "ClusterSize Endcap Pixel " << stackIdx;
        mapCluSizeEndcap_Pixel[stackIdx] = tdEndcapLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 50, 0, 50);

        histoName.str("");
        histoName << "ClusterSize_Endcap_Strip" << stackIdx;
        histoTitle.str("");
        histoTitle << "ClusterSize Endcap Strip " << stackIdx;
        mapCluSizeEndcap_Strip[stackIdx] = tdEndcapLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 50, 0, 50);

        // Denominators
        histoName.str("");
        histoName << "hTPart_Pt_Clu_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart p_{T}, Cluster, Endcap Stack " << stackIdx;
        mapCluDisk_hTPart_Pt[stackIdx] = tdEndcapLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 200, 0, 50);

        histoName.str("");
        histoName << "hTPart_Eta_Pt10_Clu_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart #eta (p_{T} > 10 GeV/c), Cluster, Endcap Stack " << stackIdx;
        mapCluDisk_hTPart_Eta_Pt10[stackIdx] = tdEndcapLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI);

        histoName.str("");
        histoName << "hTPart_Phi_Pt10_Clu_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart #phi (p_{T} > 10 GeV/c), Cluster, Endcap Stack " << stackIdx;
        mapCluDisk_hTPart_Phi_Pt10[stackIdx] = tdEndcapLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI);

        mapCluDisk_hTPart_Pt[stackIdx]->Sumw2();
        mapCluDisk_hTPart_Eta_Pt10[stackIdx]->Sumw2();
        mapCluDisk_hTPart_Phi_Pt10[stackIdx]->Sumw2();

        // Numerators GeV/c
        histoName.str("");
        histoName << "hTPart_Pt_Stub_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart p_{T}, Stub, Endcap Stack " << stackIdx;
        mapStubDisk_hTPart_Pt[stackIdx] = tdEndcapLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 200, 0, 50);

        histoName.str("");
        histoName << "hTPart_Eta_Pt10_Stub_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart #eta (p_{T} > 10 GeV/c), Stub, Endcap Stack " << stackIdx;
        mapStubDisk_hTPart_Eta_Pt10[stackIdx] = tdEndcapLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI);

        histoName.str("");
        histoName << "hTPart_Phi_Pt10_Stub_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "TPart #phi (p_{T} > 10 GeV/c), Stub, Endcap Stack " << stackIdx;
        mapStubDisk_hTPart_Phi_Pt10[stackIdx] = tdEndcapLayerDir.make<TH1D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI);

        mapStubDisk_hTPart_Pt[stackIdx]->Sumw2();
        mapStubDisk_hTPart_Eta_Pt10[stackIdx]->Sumw2();
        mapStubDisk_hTPart_Phi_Pt10[stackIdx]->Sumw2();

        // Comparison to TrackingParticle
        histoName.str("");
        histoName << "hStub_InvPt_TPart_InvPt_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub p_{T}^{-1} vs. TPart p_{T}^{-1}, Endcap Stack " << stackIdx;
        mapStubDisk_hStub_InvPt_TPart_InvPt[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 200, 0.0, 0.8, 200, 0.0, 0.8);
        mapStubDisk_hStub_InvPt_TPart_InvPt[stackIdx]->GetXaxis()->Set(NumBins, BinVec);
        mapStubDisk_hStub_InvPt_TPart_InvPt[stackIdx]->GetYaxis()->Set(NumBins, BinVec);
        mapStubDisk_hStub_InvPt_TPart_InvPt[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_Pt_TPart_Pt_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub p_{T} vs. TPart p_{T}, Endcap Stack " << stackIdx;
        mapStubDisk_hStub_Pt_TPart_Pt[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 100, 0, 50, 100, 0, 50);
        mapStubDisk_hStub_Pt_TPart_Pt[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_Eta_TPart_Eta_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub #eta vs. TPart #eta, Endcap Stack " << stackIdx;
        mapStubDisk_hStub_Eta_TPart_Eta[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 180, -M_PI, M_PI);
        mapStubDisk_hStub_Eta_TPart_Eta[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_Phi_TPart_Phi_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub #phi vs. TPart #phi, Endcap Stack " << stackIdx;
        mapStubDisk_hStub_Phi_TPart_Phi[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 180, -M_PI, M_PI);
        mapStubDisk_hStub_Phi_TPart_Phi[stackIdx]->Sumw2();

        // Residuals
        histoName.str("");
        histoName << "hStub_InvPtRes_TPart_Eta_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub p_{T}^{-1} - TPart p_{T}^{-1} vs. TPart #eta, Endcap Stack " << stackIdx;
        mapStubDisk_hStub_InvPtRes_TPart_Eta[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 100, -2.0, 2.0);
        mapStubDisk_hStub_InvPtRes_TPart_Eta[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_PtRes_TPart_Eta_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub p_{T} - TPart p_{T} vs. TPart #eta, Endcap Stack " << stackIdx;
        mapStubDisk_hStub_PtRes_TPart_Eta[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 100, -40, 40);
        mapStubDisk_hStub_PtRes_TPart_Eta[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_EtaRes_TPart_Eta_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub #eta - TPart #eta vs. TPart #eta, Endcap Stack " << stackIdx;
        mapStubDisk_hStub_EtaRes_TPart_Eta[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 100, -2, 2);
        mapStubDisk_hStub_EtaRes_TPart_Eta[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_PhiRes_TPart_Eta_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub #phi - TPart #phi vs. TPart #eta, Endcap Stack " << stackIdx;
        mapStubDisk_hStub_PhiRes_TPart_Eta[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 180, -M_PI, M_PI, 100, -0.5, 0.5);
        mapStubDisk_hStub_PhiRes_TPart_Eta[stackIdx]->Sumw2();

        // Stub Width vs. Pt
        histoName.str("");
        histoName << "hStub_W_TPart_Pt_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub Width vs. TPart p_{T}, Endcap Stack " << stackIdx;
        mapStubDisk_hStub_W_TPart_Pt[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 200, 0, 50, 41, -10.25, 10.25);
        mapStubDisk_hStub_W_TPart_Pt[stackIdx]->Sumw2();

        histoName.str("");
        histoName << "hStub_W_TPart_InvPt_D" << stackIdx;
        histoTitle.str("");
        histoTitle << "Stub Width vs. TPart p_{T}^{-1}, Endcap Stack " << stackIdx;
        mapStubDisk_hStub_W_TPart_InvPt[stackIdx] = tdEndcapLayerDir.make<TH2D>(histoName.str().c_str(), histoTitle.str().c_str(), 200, 0, 0.8, 41, -10.25, 10.25);
        mapStubDisk_hStub_W_TPart_InvPt[stackIdx]->GetXaxis()->Set(NumBins, BinVec);
        mapStubDisk_hStub_W_TPart_InvPt[stackIdx]->Sumw2();
    }

    fs->file().cd("/");
}

DEFINE_FWK_MODULE(RunStepsL1TriggerValidation);
