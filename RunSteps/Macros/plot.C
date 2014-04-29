#include <iostream>

#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TROOT.h"

#include "tdrstyle.h"

using namespace std;

int verbose=3;

vector<TString> getHistoNames(TString);

int plot(TString level="Digis", TString path="../Output/", TString pathOut="../Plots/", bool doDraw=true)
{

  // Style
  //setTDRStyle();
  //gStyle->SetOptTitle(1);
  gStyle->SetTitleSize(0.05, "XYZ");
  //gStyle->SetTitleYOffset(1.2);

  // Inputs and outputs
  if(verbose>1) cout << "- inputs and outputs" << endl;

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

  // Parts of Tracker
  if(verbose>1) cout << "- parts of tracker" << endl;
  
  const u_int nPart=3;
  const u_int nLayers[nPart]={10,8,8};
  const u_int iStartLay[nPart]={4,3,3};  // to plot values outside pixels
  TString namePart[nPart]={"Barrel", "EndCap_Side_1", "EndCap_Side_2"};

  // Efficiencies, Responses, Resolutions
  if(verbose>1) cout << "- values" << endl;

  const u_int nVar=5; // effi,respX,resolX,respY,resolY (nErr=nominal/error)
  const u_int nErr=2; // nominal/error

  TString g_name[  nVar]={"Efficiency", "ResponseX", "ResolutionX", "ResponseY", "ResolutionY"};
  TString g_titleY[nVar]={"Efficiency", "Response in x [cm]", "Resolution in x [cm]",
			  "Response in y [cm]", "Resolution in y [cm]"};

  vector<double> values[nVar][nErr][nPart];
  
  for(u_int iV=0 ; iV<nVar ; iV++)
    for(u_int iE=0 ; iE<nErr ; iE++)
      for(u_int iP=0 ; iP<nPart ; iP++)
	for(u_int iL=0 ; iL<nLayers[iP] ; iL++)
	  values[iV][iE][iP].push_back(0);

  TGraphErrors g_effi[nVar][nPart];

  // Names of Layers/Discs
  if(verbose>1) cout << "- names of layers/discs" << endl;

  const u_int nCase=2;
  TString nameTypeLayer[nPart][nCase] = { {"Layer","layer"} , {"Disc","disc"} , {"Disc","disc"} };
  vector<TString> nameLayer[nPart][nCase];
  
  // Indices of layers/discs
  vector<double*> idxLayers[nErr];
  double idxTB[ nErr][nLayers[0]];
  double idxTE1[nErr][nLayers[1]];
  double idxTE2[nErr][nLayers[2]];

  // Loop over TB, TE1, TE2
  for(u_int iP=0 ; iP<nPart ; iP++) {

    // Put the array-pointer into the global vector
    for(u_int iE=0 ; iE<nErr ; iE++) {
      idxLayers[iE].push_back( idxTB[ iE] );
      idxLayers[iE].push_back( idxTE1[iE] );
      idxLayers[iE].push_back( idxTE2[iE] );
    }

    // Loop over all layers in part #iP, define their names with/without case
    for(u_int iL=0 ; iL<nLayers[iP] ; iL++) {
      for(u_int iC=0 ; iC<nCase ; iC++) {
	nameLayer[iP][iC].push_back( nameTypeLayer[iP][iC]+"_"+TString(Form( "%d" , iL+1 )) );
      }

      // Define the index for layer #iL (in part #iP)
      idxLayers[0][iP][iL] = iL+1;
      idxLayers[1][iP][iL] = 0.5;
    }
  }


  ////////////////////
  // GET HISTOGRAMS //
  ////////////////////

  // Names of histograms
  vector<TString> myNames = getHistoNames(suffix);

  // Prepare reading from file
  TDirectory* myDir;
  TString nameDir="";

  //vector<TH1F> histos;
  TH1F* hTemp=0;
  TString nameHisto="";

  // Loop over input histograms
  if(verbose>1) cout << "- loop over input histograms" << endl;

  double tempMax=0;
  double g_max[nVar][nPart];

  for(u_int iP=0 ; iP<nPart ; iP++) {

    for(u_int iV=0 ; iV<nVar ; iV++) g_max[iV][iP]=0;

    for(u_int iL=0 ; iL<nLayers[iP] ; iL++) {

      nameDir   = "analysis/"+namePart[iP]+"/"+nameLayer[iP][0][iL] ;
      if(verbose>2) cout << "-- nameDir=" << nameDir << endl;

      myDir = f->GetDirectory(nameDir);
      if(!myDir) {
	cout << "ERROR : no such TDirectory named '" << nameDir << "'" << endl;
	continue;
      }

      if(verbose>2) cout << "-- loop over histograms (myNames.size()=" << myNames.size() << ")" << endl;

      for(u_int iH=0 ; iH<myNames.size() ; iH++) {

	nameHisto = myNames[iH]+"_"+nameLayer[iP][1][iL];
	if(verbose>2) cout << "--- histo : " << nameHisto << endl;
	
	myDir->GetObject(nameHisto , hTemp );

	if(!hTemp) {
	  cout << "ERROR : no such plot named '" << nameHisto << "'" << endl;
	  continue;
	}

	// get values
	if(verbose>2) cout << "--- nEntries=" << hTemp->GetEntries() << endl
			   << "--- getting values" << endl;
	
	// effi,respX,resolX,respY,resolY
	if( nameHisto.Contains("Efficiency") && nameHisto.Contains("Primary") ) {
	  values[0][0][iP][iL] = hTemp->GetMean();
	  values[0][1][iP][iL] = hTemp->GetMeanError();
	}
	  
	if(nameHisto.Contains("DeltaX")) {
	  if(verbose>2) cout << "--- enter DeltaX" << endl;

	  values[1][0][iP][iL] = hTemp->GetMean();
	  values[1][1][iP][iL] = hTemp->GetMeanError();
	  values[2][0][iP][iL] = hTemp->GetRMS();
	  values[2][1][iP][iL] = hTemp->GetRMSError();
	}
	  
	if(nameHisto.Contains("DeltaY")) {
	  values[3][0][iP][iL] = hTemp->GetMean();
	  values[3][1][iP][iL] = hTemp->GetMeanError();
	  values[4][0][iP][iL] = hTemp->GetRMS();
	  values[4][1][iP][iL] = hTemp->GetRMSError();
	}

	if(verbose>2) cout << "--- got values" << endl;	

	//histos.push_back( *hTemp );
	TCanvas c("c","c",5,30,800,600);
	if(doDraw) {
	  hTemp->Draw();
	  c.Print(pathOut+"/"+suffix+"/"+myNames[iH]+"_"+namePart[iP]+"_"+nameLayer[iP][1][iL]+".pdf");
	}

      }

    }
  }

  // Build first the graphs and find maxima
  if(verbose>1) cout << "- build first the graphs and find maxima" << endl;

  // Determine maxima per variable and part of the Tracker
  for(u_int iP=0 ; iP<nPart ; iP++) {
    for(u_int iL=0 ; iL<nLayers[iP] ; iL++) {
	for(u_int iV=0 ; iV<nVar ; iV++) {
	  
	  tempMax = values[iV][0][iP][iL];

	  if(verbose>2)
	    cout << "--- "    << g_name[iV] 
		 << " "       << namePart[iP] 
		 << " layer " << iL 
		 << " : val=" << tempMax << endl;
	  
	  if( tempMax > g_max[iV][iP] ) g_max[iV][iP] = tempMax;
	}
    }
  }

  double g_maximum[nVar];

  for(u_int iV=0 ; iV<nVar ; iV++) {

    g_maximum[iV]=0;

    for(u_int iP=0 ; iP<nPart ; iP++) {

      // make arrays on the fly to fill TGraphErrors constructor arguments
      double g_x[nLayers[iP]], g_y[nLayers[iP]], g_err_x[nLayers[iP]], g_err_y[nLayers[iP]];
      for(u_int iL=0 ; iL<nLayers[iP] ; iL++) {
	g_x[iL]     = idxLayers[0][iP][iL];
	g_err_x[iL] = idxLayers[1][iP][iL];
	g_y[iL]     = values[iV][0][iP][iL];
	g_err_y[iL] = values[iV][1][iP][iL];
      }

      // effi, respX,resolX,err_respX,err_resolX, respY,resolY,err_respY,err_resolY;
      g_effi[iV][iP] = TGraphErrors(    nLayers[iP]-iStartLay[iP] , 
					(iStartLay[iP]+g_x),
					(iStartLay[iP]+g_y),
					(iStartLay[iP]+g_err_x),
					(iStartLay[iP]+g_err_y)
				       );

      if( g_max[iV][iP]>g_maximum[iV] ) g_maximum[iV] = g_max[iV][iP];

    }
  }


  // Draw graphs : effi(TB,TE1,TE2), response(TB,TE1,TE2), resolution(TB,TE1,TE2)
  if(verbose>1) cout << "- draw graphs" << endl;

  for(u_int iV=0 ; iV<nVar ; iV++) {

    for(u_int iP=0 ; iP<nPart ; iP++) {

      TCanvas c_graph("cg","cg",10,10,800,600);
      
      if(verbose>2) cout << "--- " << g_name[iV] 
			 << " "    << namePart[iP] 
			 << " maxima : global=" << g_maximum[iV] 
			 << " local=" << g_max[iV][iP]
			 << endl;

      //if(iV>0) g_effi[iV][iP].SetMaximum( 1.05*g_maximum[iV] );
      if(iV>0) g_effi[iV][iP].SetMaximum( 1.05*g_max[iV][iP] );

      //g_effi[iV][iP].SetMarkerStyle( kOpenSquare );      
      //g_effi[iV][iP].SetMarkerSize( 2 );      
      //g_effi[iV][iP].SetMarkerColor( kAzure-9 );      
      g_effi[iV][iP].SetTitle( suffix+" "+g_name[iV]+" "+namePart[iP] );
      g_effi[iV][iP].GetXaxis()->SetTitle( nameTypeLayer[iP][0] );
      g_effi[iV][iP].GetYaxis()->SetTitle( g_titleY[iV] );

      if(verbose>2) cout << "--- graph mean = " << g_effi[iV][iP].GetMean() << endl;

      g_effi[iV][iP].Draw("AP");

      c_graph.Print(pathOut+"/"+suffix+"/SummaryPlots/"+g_name[iV]+"_"+suffix+"_"+namePart[iP]+".pdf");
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
