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

int compare()
{
  
  gStyle->SetTitleSize(0.05, "XYZ");

  const int nT=2;
  const int nD=2;
  const int nC=2;
  const int nP=3;
  const int nV=5;
  const int nCase=2;

  TFile* file[ nT][nD][nC];
  TGraph graph[nV][nP][nT][nD][nC];
  TGraph* gTemp=0;

  TString type[ nT]={"Digis","Clusters"};
  TString algoD[nD]={"DigiRun1","DigiPhase2"};
  TString algoC[nC]={"AdjacentHits","WeightedMean"};

  TString g_name[  nV]={"Efficiency", "ResponseX", "ResolutionX", "ResponseY", "ResolutionY"};
  TString g_titleY[nV]={"Efficiency", "Response in x [cm]", "Resolution in x [cm]",
			"Response in y [cm]", "Resolution in y [cm]"};

  TString namePart[nP]={"Barrel", "EndCap_Side_1", "EndCap_Side_2"};
  TString nameTypeLayer[nP][nCase] = { {"Layer","layer"} , {"Disc","disc"} , {"Disc","disc"} };

  int mycolor[nD][nC]={ {kBlue,kViolet+2} , {kViolet-1,kRed} };

  TString path="../Plots/";
  TString filename="";

  // Graph maxima
  double tempMax=0;
  double g_max[nV][nP];
  for(int iV=0 ; iV<nV ; iV++)
    for(int iP=0 ; iP<nP ; iP++)
      g_max[iV][iP]=0;

  for(int iT=0 ; iT<nT ; iT++) {
    for(int iD=0 ; iD<nD ; iD++) {
      for(int iC=0 ; iC<nC ; iC++) {

	filename = path+"/"+algoD[iD]+"_"+algoC[iC]+"/"+type[iT]+"/SummaryPlots/"+type[iT]+"Summary.root" ;

	file[iT][iD][iC] = new TFile(filename , "READ");

	if(file[iT][iD][iC]->IsZombie()) {
	  cout << "ERROR !!! OUTPUT FILE " << filename << " IS ZOMBIE !!!" << endl;
	  return 2;
	}

	// Get all kinds of graphs from file (iT,iD,iC)
	for(int iV=0 ; iV<nV ; iV++) {
	  for(int iP=0 ; iP<nP ; iP++) {
	    file[ iT][iD][iC]->GetObject( "g_"+type[iT]+"_"+g_name[iV]+"_"+namePart[iP] , gTemp );
	    graph[iV][iP][iT][iD][iC] = (*gTemp);
	    tempMax = graph[iV][iP][iT][iD][iC].GetMaximum();
	    if(tempMax > g_max[iV][iP]) g_max[iV][iP] = tempMax;
	  }
	}
	// got the graphs from file (iT,iD,iC)
      }
    }
  }  

  // Produce the summary plots for digis, then for clusters
  for(int iT=0 ; iT<nT ; iT++) {

    // Loop over all kinds of plot (iV,iP)=(variable,part of the tracker)
    for(int iV=0 ; iV<nV ; iV++) {
      for(int iP=0 ; iP<nP ; iP++) {

	TCanvas c_graph("cg","cg",10,10,800,600);

	// Loop over all digitizers and clusterizers to compare on the same plot
	for(int iD=0 ; iD<nD ; iD++) {
	  for(int iC=0 ; iC<nC ; iC++) {

	    if(iV>0) graph[iV][iP][iT][iD][iC].SetMaximum( 1.05*g_max[iV][iP] );
	    
	    graph[iV][iP][iT][iD][iC].SetTitle( type[iT]+" "+g_name[iV]+" "+namePart[iP] );
	    graph[iV][iP][iT][iD][iC].GetXaxis()->SetTitle( nameTypeLayer[iP][0] );
	    graph[iV][iP][iT][iD][iC].GetYaxis()->SetTitle( g_titleY[iV] );

	    graph[iV][iP][iT][iD][iC].SetMarkerStyle(kPlus);
	    graph[iV][iP][iT][iD][iC].SetMarkerSize(1.2);
	    graph[iV][iP][iT][iD][iC].SetMarkerColor(mycolor[iD][iC]);

	    if(iD==0 && iC==0) graph[iV][iP][iT][iD][iC].Draw("AP");
	    else               graph[iV][iP][iT][iD][iC].Draw("PSAME");

	  }
	}
	
	// Plot the summary plot for type iT, variable iV, tracker part iP
	c_graph.Print(path+"/ComparisonPlots/graph_"+type[iT]+"_"+g_name[iV]+"_"+namePart[iP]+".pdf");
	
      }
    }
  }
  
  // Close all input files
  /*
  for(int iT=0 ; iT<nT ; iT++)
    for(int iD=0 ; iD<nD ; iD++)
      for(int iC=0 ; iC<nC ; iC++)
	file[iT][iD][iC]->Close();
  */

  return 0;
}
