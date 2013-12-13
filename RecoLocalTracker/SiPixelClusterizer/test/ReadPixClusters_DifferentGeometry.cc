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
  edm::ESHandle< StackedTrackerGeometry >         StackedGeometryHandle;
  //const StackedTrackerGeometry*                   theStackedGeometry;

  /// Geometry setup
  /// Set pointers to Geometry
  es.get< StackedTrackerGeometryRecord >().get(StackedGeometryHandle);
  //theStackedGeometry = StackedGeometryHandle.product(); /// Note this is different from the "global" geometry

  /// Magnetic Field
  //edm::ESHandle< MagneticField > magneticFieldHandle;
  //es.get< IdealMagneticFieldRecord >().get(magneticFieldHandle);
  //const MagneticField* theMagneticField = magneticFieldHandle.product();
  //double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();

  countAllEvents++;
  int run       = e.id().run();
  int event     = e.id().event();

  int lumiBlock = e.luminosityBlock();
  int bx        = e.bunchCrossing();
  int orbit     = e.orbitNumber();

  std::cout << " ******************************************************************* " << std::endl;
  std::cout << " **** Looking at Different geometry test                         *** " << std::endl;
  std::cout << " ******************************************************************** " << std::endl;
  
  cout<<"run "<<run<<" event "<<event<<" bx "<<bx<<" lumi "<<lumiBlock<<" orbit "<<orbit<<" magnetic field strenth " << endl;  
    
} // end 

//define this as a plug-in
DEFINE_FWK_MODULE(ReadPixClusters_DifferentGeometry);
