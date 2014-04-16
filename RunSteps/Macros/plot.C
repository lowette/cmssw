#include <iostream>

#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"

using namespace std;

vector<TString> getHistoNames(TString);

int plot(TString level="Digis", TString path="../Output/", TString pathOut="../Plots/")
{

  TString file  = "" ;
  TString suffix= "" ;

  if(level.Contains("Digi")) {
    suffix = "Digis" ;
    file   = "DigiValidation.root" ;
  }
  else if(level.Contains("Clus")) {
    suffix = "Clusters" ;
    file   = "ClusterValidation.root" ;
  }

  TString nameF = path+"/"+file;

  TFile* f = new TFile(nameF, "READ");
  if(f->IsZombie()) {
    cout << "ERROR !!! FILE " << nameF << " IS ZOMBIE !!!" << endl;
    return 2;
  }

  // Names of Parts of Tracker
  const u_int nPart=3;
  TString namePart[nPart]={"Barrel", "EndCap_Side_1", "EndCap_Side_2"};

  // Names of Layers/Discs
  vector<TString> nameLayer[nPart][2];
  for(u_int iL=1 ; iL<=10 ; iL++) {
    nameLayer[0][0].push_back( TString(Form( "Layer_%d" , iL )) );
    nameLayer[0][1].push_back( TString(Form( "layer_%d" , iL )) );
  }
  for(u_int iL=1 ; iL<=8 ; iL++) {
    nameLayer[1][0].push_back( TString(Form( "Disc_%d" , iL )) );
    nameLayer[1][1].push_back( TString(Form( "disc_%d" , iL )) );
    nameLayer[2][0].push_back( TString(Form( "Disc_%d" , iL )) );
    nameLayer[2][1].push_back( TString(Form( "disc_%d" , iL )) );
  }

  // Names of histograms
  vector<TString> myNames = getHistoNames(suffix);

  // Prepare reading from file
  TDirectory* myDir;
  TString nameDir="";

  vector<TH1F> histos;
  TH1F* hTemp=0;
  TString nameHisto="";

  // Loop over input histograms
  for(u_int iPart=0 ; iPart<nPart ; iPart++) {
    for(u_int iL=0 ; iL<nameLayer[iPart][0].size() ; iL++) {

      nameDir   = "analysis/"+namePart[iPart]+"/"+nameLayer[iPart][0][iL] ;

      myDir = f->GetDirectory(nameDir);
      if(!myDir) {
	cout << "ERROR : no such TDirectory named '" << nameDir << "'" << endl;
	continue;
      }

      for(u_int iH=0 ; iH<myNames.size() ; iH++) {
	nameHisto = myNames[iH]+"_"+nameLayer[iPart][1][iL];
	
	myDir->GetObject(nameHisto , hTemp );

	if(!hTemp) {
	  cout << "ERROR : no such plot named '" << nameHisto << "'" << endl;
	  continue;
	}

	//histos.push_back( *hTemp );
	TCanvas c("c","c",5,30,800,600);
	hTemp->Draw();
	c.Print(pathOut+"/"+myNames[iH]+"_"+namePart[iPart]+"_"+nameLayer[iPart][1][iL]+".pdf");
      }

    }
  }

  cout << "THE END" << endl;
  return 0;
}

vector<TString> getHistoNames(TString suffix)
{
  vector<TString> myNames;

  const u_int nH=5;
  TString nameHistos[nH]={"NumberOfMatchedHits", "NumberOfMatched"+suffix, "Efficiency", 
			  "DeltaX_simhit_cluster", "DeltaY_simhit_cluster"};

  //const u_int nSuffix1=4;
  const u_int nSuffix2=18;
  //TString suffix1[nSuffix1] = {"AllType", "Primary", "Secondary", "Type2"};
  TString suffix2[nSuffix2] = {"Undefined","Unknown","Primary","Hadronic",
			       "Decay","Compton","Annihilation","EIoni",
			       "HIoni","MuIoni","Photon","MuPairProd",
			       "Conversions","EBrem","SynchrotronRadiation",
			       "MuBrem","MuNucl","AllTypes"};

  for(u_int iH=0 ; iH<nH ; iH++) {
    
    //if(iH==0)      
    //for(u_int iS=0 ; iS<nSuffix1 ; iS++)
    //myNames.push_back(nameHistos[iH]+"_"+suffix1[iS]);
	  
    //else if(iH==1)      
    
    if(iH<=2)
      for(u_int iS=0 ; iS<nSuffix2 ; iS++)
	myNames.push_back(nameHistos[iH]+"_"+suffix2[iS]);

    else
      myNames.push_back(nameHistos[iH]);
  }

  return myNames;
}
