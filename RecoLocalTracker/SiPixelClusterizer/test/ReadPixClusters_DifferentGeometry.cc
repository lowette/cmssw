// File: ReadPixClusters.cc
// Description: T0 test the pixel clusters. 
// Author: Danek Kotlinski 
// Creation Date:  Initial version. 3/06
// Modify to work with CMSSW354, 11/03/10 d.k.
//--------------------------------------------
#include <memory>
#include <string>
#include <iostream>

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
//#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerDetUnit.h"
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// For ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

#define HISTOS

using namespace std;

//=============================================================================

class ReadPixClusters_DifferentGeometry : public edm::EDAnalyzer {
public:
  
  explicit ReadPixClusters_DifferentGeometry(const edm::ParameterSet& conf);  
  virtual ~ReadPixClusters_DifferentGeometry();
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  virtual void beginRun(const edm::EventSetup& iSetup);
  virtual void beginJob();
  virtual void endJob();
  
private:
  edm::ParameterSet conf_;
  edm::InputTag src_;
  bool printLocal;
  int countEvents, countAllEvents;
  double sumClusters;

  //TFile* hFile;
};

/////////////////////////////////////////////////////////////////
// Contructor, empty.
ReadPixClusters_DifferentGeometry::ReadPixClusters_DifferentGeometry(edm::ParameterSet const& conf) 
  : conf_(conf), src_(conf.getParameter<edm::InputTag>( "src" )) { 
  printLocal = conf.getUntrackedParameter<bool>("Verbosity",true);
  //src_ =  conf.getParameter<edm::InputTag>( "src" );
  cout<<" Construct "<<printLocal<<endl;

}
// Virtual destructor needed.
ReadPixClusters_DifferentGeometry::~ReadPixClusters_DifferentGeometry() { }  

// ------------ method called at the begining   ------------
void ReadPixClusters_DifferentGeometry::beginRun(const edm::EventSetup& iSetup) {
  cout << "beginRun -  PixelClusterTest " <<printLocal<<endl;
}

// ------------ method called at the begining   ------------
void ReadPixClusters_DifferentGeometry::beginJob() {
  cout << "Initialize PixelClusterTest " <<printLocal<<endl;

  // NEW way to use root (from 2.0.0?)
  edm::Service<TFileService> fs;

  //=====================================================================

  countEvents=0;
  countAllEvents=0;
  sumClusters=0.;

}
// ------------ method called to at the end of the job  ------------
void ReadPixClusters_DifferentGeometry::endJob(){
  sumClusters = sumClusters/float(countEvents);
  cout << " End PixelClusTest, events all/with hits=  " << countAllEvents<<"/"<<countEvents<<" "<<sumClusters<<" "<<printLocal<<endl;

}
//////////////////////////////////////////////////////////////////
// Functions that gets called by framework every event
void ReadPixClusters_DifferentGeometry::analyze(const edm::Event& e, const edm::EventSetup& es) {
  using namespace edm;

  /// Geometry handles etc
  //edm::ESHandle< TrackerGeometry >                GeometryHandle;
  edm::ESHandle< StackedTrackerGeometry >         StackedGeometryHandle;
  const StackedTrackerGeometry*                   theStackedGeometry;
  //StackedTrackerGeometry::StackContainerIterator  StackedTrackerIterator;
  // Get event setup 
  //edm::ESHandle<TrackerGeometry> geom;
  //es.get<TrackerDigiGeometryRecord>().get( geom );
  //const TrackerGeometry& theTracker(*geom);


  /// Geometry setup
  /// Set pointers to Geometry
  //es.get< TrackerDigiGeometryRecord >().get(GeometryHandle);
  /// Set pointers to Stacked Modules
  es.get< StackedTrackerGeometryRecord >().get(StackedGeometryHandle);
  theStackedGeometry = StackedGeometryHandle.product(); /// Note this is different 
                                                        /// from the "global" geometry

  /// Magnetic Field
  edm::ESHandle< MagneticField > magneticFieldHandle;
  es.get< IdealMagneticFieldRecord >().get(magneticFieldHandle);
  const MagneticField* theMagneticField = magneticFieldHandle.product();
  double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();


  countAllEvents++;
  int run       = e.id().run();
  int event     = e.id().event();

  int lumiBlock = e.luminosityBlock();
  int bx        = e.bunchCrossing();
  int orbit     = e.orbitNumber();

  // Get Cluster Collection from InputTag
  edm::Handle< edmNew::DetSetVector<SiPixelCluster> > clusters;
  e.getByLabel( src_ , clusters);

  const edmNew::DetSetVector<SiPixelCluster>& ClusterCollection = *clusters;     
  int sizeOfClusterCollection = ClusterCollection.size();
<<<<<<< HEAD
=======

  std::cout << " ******************************************************************* " << std::endl;
  std::cout << " **** Looking at Different geometry test                         *** " << std::endl;
  std::cout << " ******************************************************************** " << std::endl;
>>>>>>> 4c5493c1fe530d2b2193c01b3381749f50c361f1
  
  if(printLocal) cout<<"run "<<run<<" event "<<event<<" bx "<<bx<<" lumi "<<lumiBlock<<" orbit "<<orbit<<" "<<sizeOfClusterCollection<<endl;  

  // Select events with pixels
  if(sizeOfClusterCollection<1) 
    std::cout << " Event encountered with pixel dets!" << std::endl;
    return; // skip events with  pixel dets
  //if(sizeOfClusterCollection<4) return; // skip events with few pixel dets

  countEvents++;
  /*
  int numberOfDetUnits = 0;
  int numberOfClusters = 0;
  int numberOfPixels = 0;
  int numberOfDetUnits1 = 0;
  int numOfClustersPerDet1=0;        
  int numOfClustersPerLay1=0;        
  int numberOfDetUnits2 = 0;
  int numOfClustersPerDet2=0;        
  int numOfClustersPerLay2=0;        
  int numberOfDetUnits3 = 0;
  int numOfClustersPerDet3=0;        
  int numOfClustersPerLay3=0;        

  int numOfPixPerLay1=0;     
  int numOfPixPerLay2=0;     
  int numOfPixPerLay3=0;     

  //int numOfPixPerDet1=0;   //Unused Value
  //int numOfPixPerDet2=0;  //Unused Value
  //int numOfPixPerDet3=0;  //Unused value
      
  //int numOfPixPerLink11=0;  //Unused value
  //int numOfPixPerLink12=0;  //Unused value
  //int numOfPixPerLink21=0;  //Unused value
  //int numOfPixPerLink22=0;  //Unused value
  //SK:unused  int numOfPixPerLink3=0;  

  int maxClusPerDet=0;
  //int maxPixPerDet=0;  //Unused value
  unsigned int maxPixPerClu=0;

  int numOfClustersPerDisk1=0;  
  int numOfClustersPerDisk2=0;  
  int numOfClustersPerDisk3=0;  
  int numOfClustersPerDisk4=0;  
  //int numOfPixPerDisk1=0;  //Unused value
  //int numOfPixPerDisk2=0;  //Unused value
  //int numOfPixPerDisk3=0;  //Unused value
  //int numOfPixPerDisk4=0;  //Unused value
        
  float aveCharge1 = 0., aveCharge2 = 0., aveCharge3 = 0., aveCharge4 = 0., aveCharge5 = 0.;

  //static int module1[416][160] = {{0}};  //Unused value
  //static int module2[416][160] = {{0}};  //Unused value
  //static int module3[416][160] = {{0}};  //Unused value
  */

  //L1TrackTriggerCode:
  /// Track Trigger
  edm::Handle< L1TkCluster_PixelDigi_Collection > PixelDigiL1TkClusterHandle;
  edm::Handle< L1TkStub_PixelDigi_Collection >    PixelDigiL1TkStubHandle;
  edm::Handle< L1TkStub_PixelDigi_Collection >    PixelDigiL1TkFailedStubHandle;
  iEvent.getByLabel( "L1TkClustersFromPixelDigis",             PixelDigiL1TkClusterHandle );
  iEvent.getByLabel( "L1TkStubsFromPixelDigis", "StubsPass",   PixelDigiL1TkStubHandle );
  iEvent.getByLabel( "L1TkStubsFromPixelDigis", "StubsFail",   PixelDigiL1TkFailedStubHandle );

  /// Go on only if there are L1TkCluster from PixelDigis
  if ( PixelDigiL1TkClusterHandle->size() > 0 )
  {
    std::cout << " Number of L1TkClusters from PixelDigis : " << PixelDigiL1TkClusterHandle->size() << std::endl;
    /// Loop over L1TkClusters
    L1TkCluster_PixelDigi_Collection::const_iterator iterL1TkCluster;
    for ( iterL1TkCluster = PixelDigiL1TkClusterHandle->begin();
          iterL1TkCluster != PixelDigiL1TkClusterHandle->end();
          ++iterL1TkCluster )
    {
  
  // get vector of detunit ids
  //--- Loop over detunits.   (Question: Isn't this a loop over the SiPixelClusters?? ~ Line 601 in L1TrackTrigger code??)
  //edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=ClusterCollection.begin();
  //for ( ; DSViter != ClusterCollection.end() ; DSViter++) {
    //bool valid = false;  //Unused Value
    unsigned int detid = iterL1TkCluster->getDetId();
    // Det id
    StackedTrackerDetId detId = DetId(detid);       // Get the Detid object
    unsigned int detType=detId.det(); // det type, pixel=1
    unsigned int subid=detId.subdetId(); //subdetector type, barrel=1
 
    if(printLocal)
      cout<<"Det: "<<detId.rawId()<<" "<<detId.null()<<" "<<detType<<" "<<subid<<endl;    

    if(detType!=1) continue; // look only at pixels
    ++numberOfDetUnits;
  
    //const GeomDetUnit * genericDet = geom->idToDet(detId);
    //const PixelGeomDetUnit * pixDet = 
    //dynamic_cast<const PixelGeomDetUnit*>(genericDet);
    
    //L1TrackTrigger geom-detector code:
    StackedTrackerDetUnit* detUnitIt = *StackedTrackerIterator;
    StackedTrackerDetId stackDetId = detUnitIt->Id();
    assert(detUnitIt == theStackedGeometry->idToStack(stackDetId));

    uint3a2_t rawId0 = theStackedGeometry->idToDet(stackDetId, 0)->geographicalId().rawId();
    uint32_t rawId1 = theStackedGeometry->idToDet(stackDetId, 1)->geographicalId().rawId();
    std::cout << " rawId0 : " << rawId0 << " , rawId1 : " << rawId1 << std::endl;

    // Get the geom-detector
    //const PixelGeomDetUnit * theGeomDet =dynamic_cast<const PixelGeomDetUnit*> (theStackedGeometry.idToDet(detId) );
    //double detZ = theGeomDet->surface().position().z();
    //double detR = theGeomDet->surface().position().perp();

    //const BoundPlane& plane = theGeomDet->surface(); //for transf.
    
    //double detThick = theGeomDet->specificSurface().bounds().thickness();
    //int cols = theGeomDet->specificTopology().ncolumns();
    //int rows = theGeomDet->specificTopology().nrows();
    
    //*****************************//
    //   Specific Pixel Topology   //
    //*****************************//
    
    //const PixelTopology * topol = &(theGeomDet->specificTopology());

    /*
    // barrel ids
    unsigned int layerC=0;
    unsigned int ladderC=0;
    unsigned int zindex=0;
    int shell  = 0; // shell id // Shell { mO = 1, mI = 2 , pO =3 , pI =4 };
    int sector = 0; // 1-8
    int ladder = 0; // 1-22
    int layer  = 0; // 1-3
    int module = 0; // 1-4
    bool half  = false; // 

    // Endcap ids
    unsigned int disk=0; //1,2,3
    unsigned int blade=0; //1-24
    unsigned int zindexF=0; //
    unsigned int side=0; //size=1 for -z, 2 for +z
    unsigned int panel=0; //panel=1

    edmNew::DetSet<SiPixelCluster>::const_iterator clustIt;

    // Subdet id, pix barrel=1, forward=2
    if(subid==2) {  // forward

      PXFDetId pdetId = PXFDetId(detid);       
      disk=pdetId.disk(); //1,2,3
      blade=pdetId.blade(); //1-24
      zindexF=pdetId.module(); //
      side=pdetId.side(); //size=1 for -z, 2 for +z
      panel=pdetId.panel(); //panel=1
      
      if(printLocal) cout<<" forward det, disk "<<disk<<", blade "<<blade<<", module "<<zindexF<<", side "<<side<<", panel "<<panel<<" pos = "<<detZ<<" "<<detR<<endl;

    } 
    else if (subid==1) {  // barrel
     
      PXBDetId pdetId = PXBDetId(detid);
      //unsigned int detTypeP=pdetId.det();
      //unsigned int subidP=pdetId.subdetId();
      // Barell layer = 1,2,3
      layerC=pdetId.layer();
      // Barrel ladder id 1-20,32,44.
      ladderC=pdetId.ladder();
      // Barrel Z-index=1,8
      zindex=pdetId.module();

      // Convert to online 
      PixelBarrelName pbn(pdetId);
      // Shell { mO = 1, mI = 2 , pO =3 , pI =4 };
      PixelBarrelName::Shell sh = pbn.shell(); //enum
      sector = pbn.sectorName();
      ladder = pbn.ladderName();
      layer  = pbn.layerName();
      module = pbn.moduleName();
      half  = pbn.isHalfModule();
      shell = int(sh);
      // change the module sign for z<0
      if(shell==1 || shell==2) module = -module;
      // change ladeer sign for Outer )x<0)
      if(shell==1 || shell==3) ladder = -ladder;
      
      if(printLocal) { 
	cout<<" Barrel layer, ladder, module "<<layerC<<" "<<ladderC<<" "<<zindex<<" "<<sh<<"("<<shell<<") "<<sector<<" "<<layer<<" "<<ladder<<" "<<module<<" "<<half<< endl;
	//cout<<" Barrel det, thick "<<detThick<<" "<<" layer, ladder, module "<<layer<<" "<<ladder<<" "<<zindex<<endl;
	//cout<<" col/row, pitch "<<cols<<" "<<rows<<" "<<pitchX<<" "<<pitchY<<endl;
      }      
      
    } // if subid
    */
    
    //if(printLocal) {
    //  cout<<"List clusters : "<<endl;
    //   cout<<"Num Charge Size SizeX SizeY X Y Xmin Xmax Ymin Ymax Edge"<<endl;
    //}

    // Loop over clusters
    //for (clustIt = DSViter->begin(); clustIt != DSViter->end(); clustIt++)  {
    //   sumClusters++;
    //   numberOfClusters++;
    //   float ch = float(clustIt->charge())/1000.; // convert ke to electrons
    //   int size = clustIt->size();
    //   int sizeX = clustIt->sizeX(); //x=row=rfi, 
    //   int sizeY = clustIt->sizeY(); //y=col=z_global
    //   float x = clustIt->x(); // cluster position as float (int+0.5)
    //   float y = clustIt->y(); // analog average
    //   // Returns int index of the cluster min/max  
    //   int minPixelRow = clustIt->minPixelRow(); //x
    //   int maxPixelRow = clustIt->maxPixelRow();
    //   int minPixelCol = clustIt->minPixelCol(); //y
    //   int maxPixelCol = clustIt->maxPixelCol();

    //   //unsigned int geoId = clustIt->geographicalId(); // always 0?!

    //   // edge method moved to topologu class
    //   bool edgeHitX = (topol->isItEdgePixelInX(minPixelRow)) || 
    // 	(topol->isItEdgePixelInX(maxPixelRow)); 
    //   bool edgeHitY = (topol->isItEdgePixelInY(minPixelCol)) || 
    // 	(topol->isItEdgePixelInY(maxPixelCol)); 

    //   bool edgeHitX2 = false; // edge method moved 
    //   bool edgeHitY2 = false; // to topologu class
            
    //   if(printLocal) 
    // 	cout<<numberOfClusters<<" "<<ch<<" "<<size<<" "<<sizeX<<" "<<sizeY<<" "<<x<<" "<<y<<" "<<minPixelRow<<" "<<maxPixelRow<<" "<<minPixelCol<<" "<<maxPixelCol<<" "<<edgeHitX<<" "<<edgeHitY<<endl;

    //   // Get the pixels in the Cluster
    //   const vector<SiPixelCluster::Pixel>& pixelsVec = clustIt->pixels();
    //   if(printLocal) cout<<" Pixels in this cluster "<<endl;
    //   bool bigInX=false, bigInY=false;
    //   // Look at pixels in this cluster. ADC is calibrated, in electrons
    //   bool edgeInX = false; // edge method moved 
    //   bool edgeInY = false; // to topologu class
    //   //SK:unused      bool cluBigInX = false; // does this clu include a big pixel
    //   //SK:unused      bool cluBigInY = false; // does this clu include a big pixel
    //   //int noisy = 0;

    //   if(pixelsVec.size()>maxPixPerClu) maxPixPerClu = pixelsVec.size();
 
    //   for (unsigned int i = 0;  i < pixelsVec.size(); ++i) { // loop over pixels
    // 	numberOfPixels++;
    // 	float pixx = pixelsVec[i].x; // index as float=iteger, row index
    // 	float pixy = pixelsVec[i].y; // same, col index
    // 	float adc = (float(pixelsVec[i].adc)/1000.);
    // 	//int chan = PixelChannelIdentifier::pixelToChannel(int(pixx),int(pixy));
    // 	//bool binInX = (RectangularPixelTopology::isItBigPixelInX(int(pixx)));
    // 	//bool bigInY = (RectangularPixelTopology::isItBigPixelInY(int(pixy)));
    // 	//int roc = rocId(int(pixy),int(pixx));  // column, row
	
    // 	edgeInX = topol->isItEdgePixelInX(int(pixx));
    // 	edgeInY = topol->isItEdgePixelInY(int(pixy));
	
    // 	if(printLocal) cout<<i<<" "<<pixx<<" "<<pixy<<" "<<adc<<" "<<bigInX<<" "<<bigInY<<" "<<edgeInX<<" "<<edgeInY<<endl;
	
    // 	if(edgeInX) edgeHitX2=true;
    // 	if(edgeInY) edgeHitY2=true; 
    // 	//SK:unused	if(bigInX) cluBigInX=true;
    // 	//SK:unused	if(bigInY) cluBigInY=true;

    //   } // pixel loop      

    //   if(edgeHitX != edgeHitX2) 
    // 	cout<<" wrong egdeX "<<edgeHitX<<" "<<edgeHitX2<<endl;
    //   if(edgeHitY != edgeHitY2) 
    // 	cout<<" wrong egdeY "<<edgeHitY<<" "<<edgeHitY2<<endl;

    // }
    // clusters 
    
    // if(numOfClustersPerDet1>maxClusPerDet) maxClusPerDet = numOfClustersPerDet1;
    // if(numOfClustersPerDet2>maxClusPerDet) maxClusPerDet = numOfClustersPerDet2;
    // if(numOfClustersPerDet3>maxClusPerDet) maxClusPerDet = numOfClustersPerDet3;

    // if(printLocal) {
    //   if(layer==1) 
    // 	cout<<"Lay1: number of clusters per det = "<<numOfClustersPerDet1<<endl;
    //   else if(layer==2) 
    // 	cout<<"Lay2: number of clusters per det = "<<numOfClustersPerDet1<<endl;
    //   else if(layer==3) 
    // 	cout<<"Lay3: number of clusters per det = "<<numOfClustersPerDet1<<endl;
    // } // end if printLocal
    
  } // detunits loop

  // if( printLocal ) {
  //   cout<<"run "<<run<<" event "<<event<<" bx "<<bx<<" lumi "<<lumiBlock<<" orbit "<<orbit<<" num "<<countEvents<<endl;   
  //   cout<<"Num of pix "<<numberOfPixels<<" num of clus "<<numberOfClusters<<" num of dets "<<sizeOfClusterCollection<<" max clus per det "<<maxClusPerDet<<" max pix per clu "<<maxPixPerClu<<" count "<<countEvents<<endl;
  //   cout<<"Number of clusters per      Lay1,2,3: "<<numOfClustersPerLay1<<" "<<numOfClustersPerLay2<<" "<<numOfClustersPerLay3<<endl;
  //   cout<<"Number of pixels per        Lay1,2,3: "<<numOfPixPerLay1<<" "<<numOfPixPerLay2<<" "<<numOfPixPerLay3<<endl;
  //   cout<<"Number of dets with clus in Lay1,2,3: "<<numberOfDetUnits1<<" "<<numberOfDetUnits2<<" "<<numberOfDetUnits3<<endl;
  //   cout<<"Number of clus in disks 1,2,3,4     : "<<numOfClustersPerDisk1<<" "<<numOfClustersPerDisk2<<" "<<numOfClustersPerDisk3<<" "<<numOfClustersPerDisk4<<endl;
  //   aveCharge1 /= float(numOfClustersPerLay1);
  //   aveCharge2 /= float(numOfClustersPerLay2);
  //   aveCharge3 /= float(numOfClustersPerLay3);
  //   aveCharge4 /= float(numOfClustersPerDisk2 + numOfClustersPerDisk3);
  //   aveCharge5 /= float(numOfClustersPerDisk1 + numOfClustersPerDisk4);
  //   cout<<" Average charge "<<aveCharge1<<" "<<aveCharge2<<" "<<aveCharge3<<" "<<aveCharge4<<" "<<aveCharge5<<endl;
  // }
    
} // end 

//define this as a plug-in
DEFINE_FWK_MODULE(ReadPixClusters_DifferentGeometry);
