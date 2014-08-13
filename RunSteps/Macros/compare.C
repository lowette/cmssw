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

  const UInt_t nT=2;
  const UInt_t nD=2;
  const UInt_t nC=2;
  const UInt_t nP=3;
  const UInt_t nV=5;
  const UInt_t nCase=2;
  const UInt_t nPS=3;

  TFile* file[ nT][nD][nC];
  TGraph graph[nV][nP][nT][nD][nC][nPS];
  TGraph* gTemp=0;
  
  // Module type
  vector<TString> namePS[nT];
  namePS[0].push_back("");
  namePS[1].push_back("_AllMod");
  namePS[1].push_back("_PixelMod");
  namePS[1].push_back("_StripMod");

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
  double g_max[nV][nP][nT][nPS], g_min[nV][nP][nT][nPS];
  for(UInt_t iV=0 ; iV<nV ; iV++)
    for(UInt_t iP=0 ; iP<nP ; iP++)
      for(UInt_t iT=0 ; iT<nT ; iT++)
	for(UInt_t iPS=0 ; iPS<namePS[iT].size() ; iPS++) {
	  g_max[iV][iP][iT][iPS]=0;
	  g_min[iV][iP][iT][iPS]=9999999999;
      }

  for(UInt_t iT=0 ; iT<nT ; iT++) {
    for(UInt_t iD=0 ; iD<nD ; iD++) {
      for(UInt_t iC=0 ; iC<nC ; iC++) {

	filename = path+"/"+algoD[iD]+"_"+algoC[iC]+"/"+type[iT]+"/SummaryPlots/"+type[iT]+"Summary.root" ;

	file[iT][iD][iC] = new TFile(filename , "READ");

	if(file[iT][iD][iC]->IsZombie()) {
	  cout << "ERROR !!! OUTPUT FILE " << filename << " IS ZOMBIE !!!" << endl;
	  return 2;
	}

	// Get all kinds of graphs from file (iT,iD,iC)
	for(UInt_t iV=0 ; iV<nV ; iV++) {
	  for(UInt_t iP=0 ; iP<nP ; iP++) {
	    for(UInt_t iPS=0 ; iPS<namePS[iT].size() ; iPS++) {

	      file[ iT][iD][iC]->GetObject( "g_"+type[iT]+"_"+g_name[iV]+"_"+namePS[iT][iPS]+"_"+namePart[iP] , gTemp );
	      graph[iV][iP][iT][iD][iC][iPS] = (*gTemp);

	      tempMax = getMaximum( graph[iV][iP][iT][iD][iC][iPS] , g_xmin[iP] , g_xmax[iP]);
	      tempMin = getMinimum( graph[iV][iP][iT][iD][iC][iPS] , g_xmin[iP] , g_xmax[iP] );

	      if(tempMax > g_max[iV][iP][iT][iPS]) g_max[iV][iP][iT][iPS] = tempMax;
	      if(tempMin < g_min[iV][iP][iT][iPS]) g_min[iV][iP][iT][iPS] = tempMin;

	      if(verbose>1)// && tempMax==0) 
		cout << algoD[iD]+"_"+algoC[iC]
		     << " : GRAPH = "  << "g_"+type[iT]+"_"+g_name[iV]+"_"+namePart[iP] 
		     << " ; MAX = " << tempMax << endl;

	      if(verbose>1)// && tempMin==0) 
		cout << " ; MIN = " << tempMin << endl;

	    }
	  }
	}
	// got the graphs from file (iT,iD,iC)
      }
    }
  }  

  // Produce the summary plots for digis, then for clusters
  for(UInt_t iT=0 ; iT<nT ; iT++) {

    // Loop over all kinds of plot (iV,iP)=(variable,part of the tracker)
    for(UInt_t iV=0 ; iV<nV ; iV++) {
      for(UInt_t iP=0 ; iP<nP ; iP++) {
	for(UInt_t iPS=0 ; iPS<namePS[iT].size() ; iPS++) {

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
	  for(UInt_t iD=0 ; iD<nD ; iD++) {
	    for(UInt_t iC=0 ; iC<nC ; iC++) {

	      if(iV>0) {
		graph[iV][iP][iT][iD][iC][iPS].SetMaximum( 1.05*g_max[iV][iP][iT][iPS] );
		graph[iV][iP][iT][iD][iC][iPS].SetMinimum( g_min[iV][iP][iT][iPS] );
	      }

	      if(verbose>1) cout << "PLOT : "   << type[iT]+" "+g_name[iV]+" "+namePart[iP] 
				 << " ; MAX : " << 1.05*g_max[iV][iP][iT][iPS]          
				 << endl;
	    
	      // Labels
	      graph[iV][iP][iT][iD][iC][iPS].SetTitle( type[iT]+" "+g_name[iV]+" "+namePart[iP] );
	      graph[iV][iP][iT][iD][iC][iPS].GetXaxis()->SetTitle( nameTypeLayer[iP][0] );
	      graph[iV][iP][iT][iD][iC][iPS].GetYaxis()->SetTitle( g_titleY[iV] );

	      // Style
	      int linestyle;
	      if(iC==0) linestyle=7;
	      else      linestyle=6;
	      graph[iV][iP][iT][iD][iC][iPS].SetLineStyle(linestyle);
	      graph[iV][iP][iT][iD][iC][iPS].SetLineColor(mycolor[iD][iC]);
	      graph[iV][iP][iT][iD][iC][iPS].SetMarkerColor(mycolor[iD][iC]);
	      graph[iV][iP][iT][iD][iC][iPS].SetMarkerStyle(mystyle[iD][iC]);
	      graph[iV][iP][iT][iD][iC][iPS].SetMarkerSize(1.2);

	      // Plot
	      if(iD==0 && iC==0) graph[iV][iP][iT][iD][iC][iPS].Draw("APL");
	      else               graph[iV][iP][iT][iD][iC][iPS].Draw("PLSAME");

	      // Legend
	      leg.AddEntry( &(graph[iV][iP][iT][iD][iC][iPS]) , algoD_r[iD]+" "+algoC_r[iC] , "P");

	    }
	  }

	  // Plot the summary plot for type iT, variable iV, tracker part iP
	  leg.Draw();
	  c_graph.Print(path+"/ComparisonPlots/graph_"+type[iT]+"_"+g_name[iV]+"_"+namePart[iP]+namePS[iT][iPS]+".pdf");
	
	} // end loop over iPS module types
      } // end loop over iP parts
    } // end loop over iV variables
  } // end loop over object type
  
  // Close all input files
  /*
  for(UInt_t iT=0 ; iT<nT ; iT++)
    for(UInt_t iD=0 ; iD<nD ; iD++)
      for(UInt_t iC=0 ; iC<nC ; iC++)
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
