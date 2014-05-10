#include <iostream>
#include <limits>

#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLegend.h"

#include "tdrstyle.h"

using namespace std;

int verbose=3;

double getMinimum(TGraph, double, double);
double getMaximum(TGraph, double, double);

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
  TString algoD_r[nD]={"Run1","Phase2"};
  TString algoC[nC]={"AdjacentHits","WeightedMean"};
  TString algoC_r[nC]={"AH","WM"};

  TString g_name[  nV]={"Efficiency", "ResponseX", "ResolutionX", "ResponseY", "ResolutionY"};
  TString g_titleY[nV]={"Efficiency", "Response in x [cm]", "Resolution in x [cm]",
			"Response in y [cm]", "Resolution in y [cm]"};

  TString namePart[nP]={"Barrel", "EndCap_Side_1", "EndCap_Side_2"};
  TString nameTypeLayer[nP][nCase] = { {"Layer","layer"} , {"Disc","disc"} , {"Disc","disc"} };
  double g_xmin[nP]={5,4,4};
  double g_xmax[nP]={10,8,8};

  int mycolor[nD][nC]={ {kBlue,kViolet} , {kPink,kRed} };
  int mystyle[nD][nC]={ {kOpenSquare,kFullCircle} , {kOpenTriangleUp,kFullDiamond} };

  TString path="../Plots/";
  TString filename="";

  // Graph maxima
  double tempMax=0, tempMin=0;
  double g_max[nV][nP][nT], g_min[nV][nP][nT];
  for(int iV=0 ; iV<nV ; iV++)
    for(int iP=0 ; iP<nP ; iP++)
      for(int iT=0 ; iT<nT ; iT++) {
	g_max[iV][iP][iT]=0;
	g_min[iV][iP][iT]=9999999999;
      }

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

	    tempMax = getMaximum( graph[iV][iP][iT][iD][iC] , g_xmin[iP] , g_xmax[iP]);
	    tempMin = getMinimum( graph[iV][iP][iT][iD][iC] , g_xmin[iP] , g_xmax[iP] );

	    if(tempMax > g_max[iV][iP][iT]) g_max[iV][iP][iT] = tempMax;
	    if(tempMin < g_min[iV][iP][iT]) g_min[iV][iP][iT] = tempMin;

	    if(verbose>1)// && tempMax==0) 
	      cout << algoD[iD]+"_"+algoC[iC]
		   << " : GRAPH = "  << "g_"+type[iT]+"_"+g_name[iV]+"_"+namePart[iP] 
		   << " ; MAX = " << tempMax << endl;

	    if(verbose>1)// && tempMin==0) 
	      cout << " ; MIN = " << tempMin << endl;

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

	TLegend leg(0.80 , 0.75 , 0.975 , 0.95);
	leg.SetMargin(0.15);
	leg.SetLineColor(1);
	leg.SetTextColor(1);
	leg.SetTextFont(42);
	leg.SetTextSize(0.03);
	leg.SetShadowColor(kWhite);
	leg.SetFillColor(kWhite);

	// Loop over all digitizers and clusterizers to compare on the same plot
	for(int iD=0 ; iD<nD ; iD++) {
	  for(int iC=0 ; iC<nC ; iC++) {

	    if(iV>0) {
	      graph[iV][iP][iT][iD][iC].SetMaximum( 1.05*g_max[iV][iP][iT] );
	      graph[iV][iP][iT][iD][iC].SetMinimum( g_min[iV][iP][iT] );
	    }

	    if(verbose>1) cout << "PLOT : "   << type[iT]+" "+g_name[iV]+" "+namePart[iP] 
			       << " ; MAX : " << 1.05*g_max[iV][iP][iT]          
			       << endl;
	    
	    // Labels
	    graph[iV][iP][iT][iD][iC].SetTitle( type[iT]+" "+g_name[iV]+" "+namePart[iP] );
	    graph[iV][iP][iT][iD][iC].GetXaxis()->SetTitle( nameTypeLayer[iP][0] );
	    graph[iV][iP][iT][iD][iC].GetYaxis()->SetTitle( g_titleY[iV] );

	    // Style
	    int linestyle;
	    if(iC==0) linestyle=7;
	    else      linestyle=6;
	    graph[iV][iP][iT][iD][iC].SetLineStyle(linestyle);
	    graph[iV][iP][iT][iD][iC].SetLineColor(mycolor[iD][iC]);
	    graph[iV][iP][iT][iD][iC].SetMarkerColor(mycolor[iD][iC]);
	    graph[iV][iP][iT][iD][iC].SetMarkerStyle(mystyle[iD][iC]);
	    graph[iV][iP][iT][iD][iC].SetMarkerSize(1.2);

	    // Plot
	    if(iD==0 && iC==0) graph[iV][iP][iT][iD][iC].Draw("APL");
	    else               graph[iV][iP][iT][iD][iC].Draw("PLSAME");

	    // Legend
	    leg.AddEntry( &(graph[iV][iP][iT][iD][iC]) , algoD_r[iD]+" "+algoC_r[iC] , "P");

	  }
	}

	// Plot the summary plot for type iT, variable iV, tracker part iP
	leg.Draw();
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

double getMinimum(TGraph g, double xMin, double xMax)
{

  double minimum=DBL_MAX;

  int N = g.GetN();
  double* X = g.GetX();
  double* Y = g.GetY();

  for(int i=0 ; i<N ; i++) {
    if(X[i]<xMin || X[i]>xMax) continue;
    
    if(Y[i] < minimum) minimum = Y[i];
  }

  return minimum;
}

double getMaximum(TGraph g, double xMin, double xMax)
{

  double maximum=-DBL_MAX;
  
  int N = g.GetN();
  double* X = g.GetX();
  double* Y = g.GetY();

  for(int i=0 ; i<N ; i++) {
    if(X[i]<xMin || X[i]>xMax) continue;
    
    if(Y[i] > maximum) maximum = Y[i];
  }

  return maximum;
}
