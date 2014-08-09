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
#include "TF1.h"

#include "tdrstyle.h"

using namespace std;

int verbose=3;

vector<TString> getHistoNames(TString,TString,TString);

int plot(TString level="Digis", TString file="../Output/DigiValidation_Phase2.root", TString pathOut="../Plots/", bool doDraw=true)
{

  // Style
  //setTDRStyle();
  //gStyle->SetOptTitle(1);
  gStyle->SetTitleSize(0.05, "XYZ");
  //gStyle->SetTitleYOffset(1.2);

  // Inputs and outputs
  if(verbose>1) cout << "- inputs and outputs" << endl;

  //TString file  = "" ;
  TString suffix= "" ;
  TString suffix0= "" ;
  TString reco  = "" ;

  if(level.Contains("Digi")) {
    suffix0 = suffix = "Digis" ;
    reco   = "digi" ;
    //file   = "DigiValidation.root" ;
  }
  else if(level.Contains("Clus")) {
    suffix0 = suffix = "Clusters" ;
    reco   = "cluster" ;
    //file   = "ClusterValidation.root" ;
  }
  else if(level.Contains("Rec")) {
    suffix0= "Clusters" ;
    suffix = "RecHits" ;
    reco   = "cluster" ;
    //file   = "RecHitValidation.root" ;
  }

  TFile* f = new TFile(file, "READ");
  if(f->IsZombie()) {
    cout << "ERROR !!! INPUT FILE " << file << " IS ZOMBIE !!!" << endl;
    return 2;
  }

  TString nameF_o = pathOut+"/"+suffix+"/SummaryPlots/"+suffix+"Summary.root";
  TFile* fout = new TFile(nameF_o, "RECREATE");
  if(fout->IsZombie()) {
    cout << "ERROR !!! OUTPUT FILE " << nameF_o << " IS ZOMBIE !!!" << endl;
    return 2;
  }

  // Parts of Tracker
  if(verbose>1) cout << "- parts of tracker" << endl;

  const UInt_t nFit=2;
  const UInt_t nPart=3;
  const UInt_t nLayers[nPart]={10,8,8};
  const UInt_t iStartLay[nPart]={4,3,3};  // to plot values outside pixels
  TString namePart[nPart]={"Barrel", "EndCap_Side_1", "EndCap_Side_2"};
  TString nameFit[nFit]={"histo","fit"};

  // Efficiencies, Responses, Resolutions
  if(verbose>1) cout << "- values" << endl;

  const UInt_t nVar=9; // effi,respX,resolX,respY,resolY (nErr=nominal/error)
  const UInt_t nErr=2; // nominal/error

  TString g_name[  nVar]={"Efficiency", "ResponseX", "ResolutionX", "ResponseY", "ResolutionY", 
			  "RH_ResponseX", "RH_ResponseY", "RH_ResolutionX", "RH_ResolutionY"};

  TString g_titleY[nVar]={"Efficiency", 
			  "Response in x [cm]", "Resolution in x [cm]",
			  "Response in y [cm]", "Resolution in y [cm]", 
			  "RecHit Response in x [cm]", "RecHit Resolution in x [cm]",
			  "RecHit Response in y [cm]", "RecHit Resolution in y [cm]"};

  const UInt_t nPS=3;
  TString name_PS[nPS] = {"_AllMod", "_PixelMod", "_StripMod"};

  UInt_t nPS_L=3;
  if(suffix.Contains("Digi")) {
    nPS_L=1;
    for(UInt_t iPS=0 ; iPS<nPS ; iPS++)
      name_PS[iPS] = "";  
  }

  vector<double> values[nVar][nErr][nPart][nFit][nPS_L];
  
  for(UInt_t iV=0 ; iV<nVar ; iV++)
    for(UInt_t iE=0 ; iE<nErr ; iE++)
      for(UInt_t iP=0 ; iP<nPart ; iP++)
	for(UInt_t iF=0 ; iF<nFit ; iF++)
	  for(UInt_t iL=0 ; iL<nLayers[iP] ; iL++)
	    for(UInt_t iPS=0 ; iPS<nPS_L ; iPS++)
	      values[iV][iE][iP][iF][iPS].push_back(0);

  TGraphErrors g_effi[nVar][nPart][nFit][nPS_L];

  // Names of Layers/Discs
  if(verbose>1) cout << "- names of layers/discs" << endl;

  const UInt_t nCase=2;
  TString nameTypeLayer[nPart][nCase] = { {"Layer","layer"} , {"Disc","disc"} , {"Disc","disc"} };
  vector<TString> nameLayer[nPart][nCase];
  
  // Indices of layers/discs
  vector<double*> idxLayers[nErr];
  double idxTB[ nErr][nLayers[0]];
  double idxTE1[nErr][nLayers[1]];
  double idxTE2[nErr][nLayers[2]];

  // Loop over TB, TE1, TE2
  for(UInt_t iP=0 ; iP<nPart ; iP++) {

    // Put the array-pointer into the global vector
    for(UInt_t iE=0 ; iE<nErr ; iE++) {
      idxLayers[iE].push_back( idxTB[ iE] );
      idxLayers[iE].push_back( idxTE1[iE] );
      idxLayers[iE].push_back( idxTE2[iE] );
    }

    // Loop over all layers in part #iP, define their names with/without case
    for(UInt_t iL=0 ; iL<nLayers[iP] ; iL++) {
      for(UInt_t iC=0 ; iC<nCase ; iC++) {
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
  vector<TString> myNames = getHistoNames(suffix0, reco, suffix);

  // Prepare reading from file
  TDirectory* myDir;
  TString nameDir="";

  //vector<TH1F> histos;
  TH1F* hTemp=0;
  TString nameHisto="";
  TF1 *fTemp  = 0;
  
  double mean     = 0;
  double err_mean = 0;
  double rms      = 0;
  double err_rms  = 0;
  
  double fMean      = 0;
  double fSigma     = 0;
  double err_fMean  = 0;
  double err_fSigma = 0;

  // Loop over input histograms
  if(verbose>1) cout << "- loop over input histograms" << endl;

  double tempMax=0;
  double g_max[nVar][nPart][nFit][nPS_L];

  for(UInt_t iP=0 ; iP<nPart ; iP++) {

    for(UInt_t iV=0 ; iV<nVar ; iV++) 
      for(UInt_t iF=0 ; iF<nFit ; iF++) 
	for(UInt_t iPS=0 ; iPS<nPS_L ; iPS++)
	  g_max[iV][iP][iF][iPS]=0;

    for(UInt_t iL=0 ; iL<nLayers[iP] ; iL++) {
      
      nameDir   = "analysis/"+namePart[iP]+"/"+nameLayer[iP][0][iL] ;
      if(verbose>2) cout << "-- nameDir=" << nameDir << endl;

      myDir = f->GetDirectory(nameDir);
      if(!myDir) {
	cout << "ERROR : no such TDirectory named '" << nameDir << "'" << endl;
	continue;
      }

      if(verbose>2) cout << "-- loop over histograms (myNames.size()=" << myNames.size() << ")" << endl;

      for(UInt_t iH=0 ; iH<myNames.size() ; iH++) {

	nameHisto = myNames[iH]+"_"+nameLayer[iP][1][iL];
	if(verbose>2) cout << "--- histo : " << nameHisto << endl;
	
	myDir->GetObject(nameHisto , hTemp );

	if(!hTemp) {
	  cout << "ERROR : no such plot named '" << nameHisto << "'" << endl;
	  continue;
	}

	UInt_t idxPS=0;
	if(nameHisto.Contains("AllMod"))   idxPS=0;
	if(nameHisto.Contains("PixelMod")) idxPS=1;
	if(nameHisto.Contains("StripMod")) idxPS=2;

	// fit the plots
	if( (nameHisto.Contains("Efficiency") && nameHisto.Contains("Primary")) 
	    || nameHisto.Contains("DeltaX") || nameHisto.Contains("DeltaY") ) {

	  hTemp->Fit("gaus");
	  fTemp  = hTemp->GetFunction("gaus");

	  //continue; // nadir test

	  mean     = hTemp->GetMean();
	  err_mean = hTemp->GetMeanError();
	  rms      = hTemp->GetRMS();
	  err_rms  = hTemp->GetRMSError();

	  if(fTemp) {
	    fMean      = fTemp->GetParameter(1);
	    fSigma     = fTemp->GetParameter(2);
	    err_fMean  = fTemp->GetParError(1);
	    err_fSigma = fTemp->GetParError(2);
	  }
	  else {
	    fMean = fSigma = err_fMean = err_fSigma = 0;
	    if(verbose>1) cout << "!!! FIT FAILED !!!" << endl;
	  }
	}
	else mean = err_mean = rms = err_rms = fMean = fSigma = err_fMean = err_fSigma = 0;

	// get values
	if(verbose>2) cout << "--- nEntries=" << hTemp->GetEntries() << endl
			   << "--- getting values" << endl;
	
	// effi,respX,resolX,respY,resolY
	if( nameHisto.Contains("Efficiency") && nameHisto.Contains("Primary") ) {
	  if(verbose>2) cout << "--- enter Efficiency Primary" << endl;

	  values[0][0][iP][0][iL][idxPS] = mean;
	  values[0][1][iP][0][iL][idxPS] = err_mean;
	  values[0][0][iP][1][iL][idxPS] = fMean;
	  values[0][1][iP][1][iL][idxPS] = err_fMean;
	}
	  
	if(nameHisto.Contains("DeltaX")) {
	  if(nameHisto.Contains("simhit")) {
	    if(verbose>2) cout << "--- enter DeltaX" << endl;
	    
	    values[1][0][iP][0][iL][idxPS] = mean;
	    values[1][1][iP][0][iL][idxPS] = err_mean;
	    values[2][0][iP][0][iL][idxPS] = rms;
	    values[2][1][iP][0][iL][idxPS] = err_rms;
	    
	    values[1][0][iP][1][iL][idxPS] = fMean;
	    values[1][1][iP][1][iL][idxPS] = err_fMean;
	    values[2][0][iP][1][iL][idxPS] = fSigma;
	    values[2][1][iP][1][iL][idxPS] = err_fSigma;
	  }
	  else if(nameHisto.Contains("RecHit")) {
	    values[5][0][iP][0][iL][idxPS] = mean;
	    values[5][1][iP][0][iL][idxPS] = err_mean;
	    values[6][0][iP][0][iL][idxPS] = rms;
	    values[6][1][iP][0][iL][idxPS] = err_rms;
	    
	    values[5][0][iP][1][iL][idxPS] = fMean;
	    values[5][1][iP][1][iL][idxPS] = err_fMean;
	    values[6][0][iP][1][iL][idxPS] = fSigma;
	    values[6][1][iP][1][iL][idxPS] = err_fSigma;
	  }
	}
	if(nameHisto.Contains("DeltaY")) {
	  if(nameHisto.Contains("simhit")) {
	    if(verbose>2) cout << "--- enter DeltaY" << endl;
	    
	    values[3][0][iP][0][iL][idxPS] = mean;
	    values[3][1][iP][0][iL][idxPS] = err_mean;
	    values[4][0][iP][0][iL][idxPS] = rms;
	    values[4][1][iP][0][iL][idxPS] = err_rms;
	    
	    values[3][0][iP][1][iL][idxPS] = fMean;
	    values[3][1][iP][1][iL][idxPS] = err_fMean;
	    values[4][0][iP][1][iL][idxPS] = fSigma;
	    values[4][1][iP][1][iL][idxPS] = err_fSigma;
	  }
	  else if(nameHisto.Contains("RecHit")) {
	    values[7][0][iP][0][iL][idxPS] = mean;
	    values[7][1][iP][0][iL][idxPS] = err_mean;
	    values[8][0][iP][0][iL][idxPS] = rms;
	    values[8][1][iP][0][iL][idxPS] = err_rms;
	    
	    values[7][0][iP][1][iL][idxPS] = fMean;
	    values[7][1][iP][1][iL][idxPS] = err_fMean;
	    values[8][0][iP][1][iL][idxPS] = fSigma;
	    values[8][1][iP][1][iL][idxPS] = err_fSigma;
	  }
	}
	
	if(verbose>2) cout << "--- got values" << endl;	

	//histos.push_back( *hTemp );
	TCanvas c("c","c",5,30,800,600);
	if(doDraw) {
	  hTemp->Draw();
	  c.Print(pathOut+"/"+suffix+"/"+myNames[iH]+"_"+namePart[iP]+"_"+nameLayer[iP][1][iL]+".pdf");
	  if(iH>2) c.Print(pathOut+"/"+suffix+"/"+myNames[iH]+"_"+namePart[iP]+"_"+nameLayer[iP][1][iL]+".png");
	}

      }

    }
  }

  //return 42; // nadir test

  // Build first the graphs and find maxima
  if(verbose>1) cout << "- build first the graphs and find maxima" << endl;

  // Determine maxima per variable and part of the Tracker
  for(UInt_t iF=0 ; iF<nFit ; iF++) {
    for(UInt_t iP=0 ; iP<nPart ; iP++) {
      for(UInt_t iL=0 ; iL<nLayers[iP] ; iL++) {
	for(UInt_t iV=0 ; iV<nVar ; iV++) {
	  for(UInt_t iPS=0 ; iPS<nPS_L ; iPS++) {
	  
	    tempMax = values[iV][0][iP][iF][iL][iPS];
	  
	    if(verbose>2)
	      cout << "--- "    << g_name[iV] 
		   << " "       << namePart[iP] 
		   << " layer " << iL 
		   << " : val=" << tempMax << endl;
	  
	    if( tempMax > g_max[iV][iP][iF][iPS] ) g_max[iV][iP][iF][iPS] = tempMax;
	  }
	}
      }
    }
  }

  if(verbose>1) cout << "Build the graphs" << endl;

  double g_maximum[nVar][nFit][nPS_L];

  for(UInt_t iPS=0 ; iPS<nPS_L ; iPS++) {
    for(UInt_t iF=0 ; iF<nFit ; iF++) {
      for(UInt_t iV=0 ; iV<nVar ; iV++) {

	g_maximum[iV][iF][iPS]=0;

	for(UInt_t iP=0 ; iP<nPart ; iP++) {

	  if(verbose>1) cout << "iF=" << iF << " iV=" << iV << " iP=" << iP << endl;

	  // make arrays on the fly to fill TGraphErrors constructor arguments
	  double g_x[nLayers[iP]], g_y[nLayers[iP]], g_err_x[nLayers[iP]], g_err_y[nLayers[iP]];
	  for(UInt_t iL=0 ; iL<nLayers[iP] ; iL++) {
	    g_x[iL]     = idxLayers[0][iP][iL];
	    g_err_x[iL] = idxLayers[1][iP][iL];
	    g_y[iL]     = values[iV][0][iP][iF][iL][iPS]; 
	    g_err_y[iL] = values[iV][1][iP][iF][iL][iPS]; 
	  }
	
	  // effi, respX,resolX,err_respX,err_resolX, respY,resolY,err_respY,err_resolY;
	  g_effi[iV][iP][iF][iPS] = TGraphErrors(    nLayers[iP]-iStartLay[iP] , 
						     (iStartLay[iP]+g_x),
						     (iStartLay[iP]+g_y),
						     (iStartLay[iP]+g_err_x),
						     (iStartLay[iP]+g_err_y)
						     );
	
	  if( g_max[iV][iP][iF][iPS]>g_maximum[iV][iF][iPS] ) g_maximum[iV][iF][iPS] = g_max[iV][iP][iF][iPS];
	
	}
      }
    }
  }

  // Draw graphs : effi(TB,TE1,TE2), response(TB,TE1,TE2), resolution(TB,TE1,TE2)
  if(verbose>1) cout << "- draw graphs" << endl;

  fout->cd();

  for(UInt_t iV=0 ; iV<nVar ; iV++) {

    for(UInt_t iPS=0 ; iPS<nPS_L ; iPS++) {

      for(UInt_t iP=0 ; iP<nPart ; iP++) {

	for(UInt_t iF=0 ; iF<nFit ; iF++) {

	  TCanvas c_graph("cg","cg",10,10,800,600);
      
	  if(verbose>2) cout << "--- " << g_name[iV] 
			     << " "    << namePart[iP] 
			     << " maxima : global=" << g_maximum[iV][iF][iPS]
			     << " local=" << g_max[iV][iP][iF][iPS]
			     << endl;
	  
	  //if(iV>0) g_effi[iV][iP][iF].SetMaximum( 1.05*g_maximum[iV] );
	  if(iV>0) g_effi[iV][iP][iF][iPS].SetMaximum( 1.05*g_max[iV][iP][iF][iPS] );

	  //g_effi[iV][iP][iF].SetMarkerStyle( kOpenSquare );      
	  //g_effi[iV][iP][iF].SetMarkerSize( 2 );      
	  //g_effi[iV][iP][iF].SetMarkerColor( kAzure-9 );      
	  g_effi[iV][iP][iF][iPS].SetTitle( suffix+" "+g_name[iV]+" "+namePart[iP] );
	  g_effi[iV][iP][iF][iPS].GetXaxis()->SetTitle( nameTypeLayer[iP][0] );
	  g_effi[iV][iP][iF][iPS].GetYaxis()->SetTitle( g_titleY[iV] );

	  if(verbose>2) cout << "--- graph mean = " << g_effi[iV][iP][iF][iPS].GetMean() << endl;

	  g_effi[iV][iP][iF][iPS].Draw("AP");
	  g_effi[iV][iP][iF][iPS].Write( "g_"+suffix+"_"+g_name[iV]+"_"+name_PS[iPS]+"_"+namePart[iP]+"_"+nameFit[iF] );

	  c_graph.Print(pathOut+"/"+suffix+"/SummaryPlots/"+
			g_name[iV]+"_"+suffix+"_"+name_PS[iPS]+"_"+namePart[iP]+"_"+nameFit[iF]
			+".pdf");
	}
      }
    }
  }

  fout->Write();

  f->Close();
  fout->Close();

  cout << "THE END" << endl;
  return 0;
}

vector<TString> getHistoNames(TString suffix0, TString reco, TString suffix)
{

  vector<TString> myNames;

  const UInt_t nPS=3;
  TString name_PS[nPS] = {"_AllMod", "_PixelMod", "_StripMod"};

  UInt_t nPS_loop=3;
  if(suffix.Contains("Digi")) {
    nPS_loop=1;
    for(UInt_t iPS=0 ; iPS<nPS ; iPS++)
      name_PS[iPS] = "";  
  }

  const UInt_t nH=5;
  TString commonNames[nH]={"NumberOfMatchedHits", "NumberOfMatched"+suffix0, "Efficiency", 
			  "DeltaX_simhit_"+reco, "DeltaY_simhit_"+reco};
  TString rechitNames[2]={"DeltaX_Cluster_RecHit_", "DeltaY_Cluster_RecHit_"};

  vector<TString> nameHistos;

  for(UInt_t iH=0 ; iH<nH ; iH++) {
    for(UInt_t iPS=0 ; iPS<nPS_loop ; iPS++)
      nameHistos.push_back(commonNames[iH]+name_PS[iPS]);
  }

  if(suffix.Contains("Rec"))
    for(UInt_t iH=0 ; iH<2 ; iH++)  nameHistos.push_back(rechitNames[iH]);

  //const UInt_t nSuffix1=4;
  const UInt_t nSuffix2=18;
  //TString suffix1[nSuffix1] = {"AllType", "Primary", "Secondary", "Type2"};
  TString suffix2[nSuffix2] = {"Undefined","Unknown","Primary","Hadronic",
			       "Decay","Compton","Annihilation","EIoni",
			       "HIoni","MuIoni","Photon","MuPairProd",
			       "Conversions","EBrem","SynchrotronRadiation",
			       "MuBrem","MuNucl","AllTypes"};

  for(UInt_t iH=0 ; iH<nameHistos.size() ; iH++) {
    
    //if(iH==0)      
    //for(UInt_t iS=0 ; iS<nSuffix1 ; iS++)
    //myNames.push_back(nameHistos[iH]+"_"+suffix1[iS]);
	  
    //else if(iH==1)      
    
    if(iH<=2*nPS_loop)
      for(UInt_t iS=0 ; iS<nSuffix2 ; iS++)
	myNames.push_back(nameHistos[iH]+"_"+suffix2[iS]);

    else
      myNames.push_back(nameHistos[iH]);
  }

  return myNames;
}
