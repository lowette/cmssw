
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TText.h"


void rechitresiduals() {

  TFile * f = TFile::Open("rechits_validation.root");

  TH1F * h1p = (TH1F *) f->FindObjectAny("Cluster_Size_Pixel_layer_1");
  TH1F * h2p = (TH1F *) f->FindObjectAny("Cluster_Size_Pixel_layer_2");
  TH1F * h3p = (TH1F *) f->FindObjectAny("Cluster_Size_Pixel_layer_3");
  h1p->Add(h2p);
  h1p->Add(h3p);
  TH1F * h1s = (TH1F *) f->FindObjectAny("Cluster_Size_Strip_layer_1");
  TH1F * h2s = (TH1F *) f->FindObjectAny("Cluster_Size_Strip_layer_2");
  TH1F * h3s = (TH1F *) f->FindObjectAny("Cluster_Size_Strip_layer_3");
  h1s->Add(h2s);
  h1s->Add(h3s);
  TH1F * h4s = (TH1F *) f->FindObjectAny("Cluster_Size_Strip_layer_4");
  TH1F * h5s = (TH1F *) f->FindObjectAny("Cluster_Size_Strip_layer_5");
  TH1F * h6s = (TH1F *) f->FindObjectAny("Cluster_Size_Strip_layer_6");
  h4s->Add(h5s);
  h4s->Add(h6s);

  h1p->GetXaxis()->SetRange(1,10);
  h1s->GetXaxis()->SetRange(1,10);
  h4s->GetXaxis()->SetRange(1,10);

  TCanvas * c0 = new TCanvas();
  TH1 * hn1 = h1p->DrawNormalized();
  hn1->SetMaximum(2.0);
  hn1->SetXTitle("Cluster size");
  hn1->SetYTitle("Arbitrary units");
  h1p->SetLineColor(kRed+1);
  h1s->SetLineColor(kOrange-3);
  h4s->SetLineColor(kBlack);
  h1p->SetLineWidth(2);
  h1s->SetLineWidth(2);
  h4s->SetLineWidth(3);
  h1p->SetLineStyle(1);
  h1s->SetLineStyle(2);
  h4s->SetLineStyle(3);

  h1p->DrawNormalized("same");
  h1s->DrawNormalized("same");
  h4s->DrawNormalized("same");

  TLegend * l0 = new TLegend(0.63,0.72,0.91,0.91);
  l0->SetBorderSize(0);
  l0->SetFillStyle(0);
  l0->AddEntry(h1p,"Macro-pixels","L");
  l0->AddEntry(h1s,"PS strips","L");
  l0->AddEntry(h4s,"2S strips","L");
  l0->Draw();

  TText * t1 = new TText(.18,.96,"CMS");
  TText * t2 = new TText(.26,.96,"Phase2 Simulation");
  t1->SetNDC();
  t2->SetNDC();
  t1->SetTextFont(62);
  t2->SetTextFont(52);
  t1->SetTextSize(0.035);
  t2->SetTextSize(0.033);
  t1->Draw();
  t2->Draw();

  c0->SetLogy();
  c0->SaveAs("clustersize_p_vs_s_vs_2s.png");
  c0->SaveAs("clustersize_p_vs_s_vs_2s.pdf");

  h1p = (TH1F *) f->FindObjectAny("Delta_X_Pixel_layer_1");
  h2p = (TH1F *) f->FindObjectAny("Delta_X_Pixel_layer_2");
  h3p = (TH1F *) f->FindObjectAny("Delta_X_Pixel_layer_3");
  h1p->Add(h2p);
  h1p->Add(h3p);
  h1s = (TH1F *) f->FindObjectAny("Delta_X_Strip_layer_1");
  h2s = (TH1F *) f->FindObjectAny("Delta_X_Strip_layer_2");
  h3s = (TH1F *) f->FindObjectAny("Delta_X_Strip_layer_3");
  h1s->Add(h2s);
  h1s->Add(h3s);
  h4s = (TH1F *) f->FindObjectAny("Delta_X_Strip_layer_4");
  h5s = (TH1F *) f->FindObjectAny("Delta_X_Strip_layer_5");
  h6s = (TH1F *) f->FindObjectAny("Delta_X_Strip_layer_6");
  h4s->Add(h5s);
  h4s->Add(h6s);

  TCanvas * c1 = new TCanvas();
  hn1 = h1p->DrawNormalized();
  hn1->SetMaximum(0.1);
  hn1->SetXTitle("#Deltax(cluster,simhit) [cm]");
  hn1->SetYTitle("Arbitrary units");
  h1p->SetLineColor(kRed+1);
  h1s->SetLineColor(kOrange-3);
  h4s->SetLineColor(kBlack);
  h1p->SetLineWidth(2);
  h1s->SetLineWidth(2);
  h4s->SetLineWidth(3);
  h1p->SetLineStyle(1);
  h1s->SetLineStyle(2);
  h4s->SetLineStyle(3);

  h1p->DrawNormalized("same");
  h1p->Fit("gaus","N","",-0.006,0.006);
  h1s->DrawNormalized("same");
  h1s->Fit("gaus","N","",-0.006,0.006);
  h4s->DrawNormalized("same");
  h4s->Fit("gaus","N","",-0.006,0.006);

  TLegend * l1 = new TLegend(0.63,0.72,0.91,0.91);
  l1->SetBorderSize(0);
  l1->SetFillStyle(0);
  l1->AddEntry(h1p,"Macro-pixels","L");
  l1->AddEntry(h1s,"PS strips","L");
  l1->AddEntry(h4s,"2S strips","L");
  l1->Draw();

  t1->Draw();
  t2->Draw();

  c1->SaveAs("rechitresiduals_p_vs_s_vs_2s.png");
  c1->SaveAs("rechitresiduals_p_vs_s_vs_2s.pdf");

  TH1F * h1p1 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_1_ClS_1"); h1p1->Rebin();
  TH1F * h2p1 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_2_ClS_1"); h2p1->Rebin();
  TH1F * h3p1 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_3_ClS_1"); h3p1->Rebin();
  h1p1->Add(h2p1);
  h1p1->Add(h3p1);
  TH1F * h1s1 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_1_ClS_1"); h1s1->Rebin();
  TH1F * h2s1 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_2_ClS_1"); h2s1->Rebin();
  TH1F * h3s1 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_3_ClS_1"); h3s1->Rebin();
  h1p1->Add(h1s1);
  h1p1->Add(h2s1);
  h1p1->Add(h3s1);
  TH1F * h4s1 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_4_ClS_1"); h4s1->Rebin();
  TH1F * h5s1 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_5_ClS_1"); h5s1->Rebin();
  TH1F * h6s1 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_6_ClS_1"); h6s1->Rebin();
  h1p1->Add(h4s1);
  h1p1->Add(h5s1);
  h1p1->Add(h6s1);

  TH1F * h1p2 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_1_ClS_2"); h1p2->Rebin();
  TH1F * h2p2 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_2_ClS_2"); h2p2->Rebin();
  TH1F * h3p2 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_3_ClS_2"); h3p2->Rebin();
  h1p2->Add(h2p2);
  h1p2->Add(h3p2);
  TH1F * h1s2 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_1_ClS_2"); h1s2->Rebin();
  TH1F * h2s2 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_2_ClS_2"); h2s2->Rebin();
  TH1F * h3s2 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_3_ClS_2"); h3s2->Rebin();
  h1p2->Add(h1s2);
  h1p2->Add(h2s2);
  h1p2->Add(h3s2);
  TH1F * h4s2 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_4_ClS_2"); h4s2->Rebin();
  TH1F * h5s2 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_5_ClS_2"); h5s2->Rebin();
  TH1F * h6s2 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_6_ClS_2"); h6s2->Rebin();
  h1p2->Add(h4s2);
  h1p2->Add(h5s2);
  h1p2->Add(h6s2);

  TH1F * h1p3 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_1_ClS_3"); h1p3->Rebin();
  TH1F * h2p3 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_2_ClS_3"); h2p3->Rebin();
  TH1F * h3p3 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_3_ClS_3"); h3p3->Rebin();
  h1p3->Add(h2p3);
  h1p3->Add(h3p3);
  TH1F * h1s3 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_1_ClS_3"); h1s3->Rebin();
  TH1F * h2s3 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_2_ClS_3"); h2s3->Rebin();
  TH1F * h3s3 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_3_ClS_3"); h3s3->Rebin();
  h1p3->Add(h1s3);
  h1p3->Add(h2s3);
  h1p3->Add(h3s3);
  TH1F * h4s3 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_4_ClS_3"); h4s3->Rebin();
  TH1F * h5s3 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_5_ClS_3"); h5s3->Rebin();
  TH1F * h6s3 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_6_ClS_3"); h6s3->Rebin();
  h1p3->Add(h4s3);
  h1p3->Add(h5s3);
  h1p3->Add(h6s3);

  TH1F * h1p4 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_1_ClS_4"); h1p4->Rebin();
  TH1F * h2p4 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_2_ClS_4"); h2p4->Rebin();
  TH1F * h3p4 = (TH1F *) f->FindObjectAny("Pull_X_Pixel_layer_3_ClS_4"); h3p4->Rebin();
  h1p4->Add(h2p4);
  h1p4->Add(h3p4);
  TH1F * h1s4 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_1_ClS_4"); h1s4->Rebin();
  TH1F * h2s4 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_2_ClS_4"); h2s4->Rebin();
  TH1F * h3s4 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_3_ClS_4"); h3s4->Rebin();
  h1p4->Add(h1s4);
  h1p4->Add(h2s4);
  h1p4->Add(h3s4);
  TH1F * h4s4 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_4_ClS_4"); h4s4->Rebin();
  TH1F * h5s4 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_5_ClS_4"); h5s4->Rebin();
  TH1F * h6s4 = (TH1F *) f->FindObjectAny("Pull_X_Strip_layer_6_ClS_4"); h6s4->Rebin();
  h1p4->Add(h4s4);
  h1p4->Add(h5s4);
  h1p4->Add(h6s4);

  TCanvas * c2 = new TCanvas();
  TH1 * hn2 = h1p1->DrawNormalized();
  hn2->SetMaximum(0.13);
  hn2->SetXTitle("Residual pull");
  hn2->SetYTitle("Arbitrary units");
  h1p1->SetLineColor(kRed+1);
  h1p2->SetLineColor(kOrange-3);
  h1p3->SetLineColor(kBlack);
  h1p4->SetLineColor(kViolet-1);
  h1p1->SetLineWidth(2);
  h1p2->SetLineWidth(2);
  h1p3->SetLineWidth(3);
  h1p4->SetLineWidth(2);
  h1p1->SetLineStyle(1);
  h1p2->SetLineStyle(2);
  h1p3->SetLineStyle(3);
  h1p4->SetLineStyle(4);

  h1p1->DrawNormalized("same");
  h1p1->Fit("gaus","N","",-2,2);
  h1p2->DrawNormalized("same");
  h1p2->Fit("gaus","N","",-2,2);
  h1p3->DrawNormalized("same");
  h1p3->Fit("gaus","N","",-2,2);
  h1p4->DrawNormalized("same");
  h1p4->Fit("gaus","N","",-2,2);

  TLegend * l2 = new TLegend(0.6,0.67,0.98,0.91);
  l2->SetBorderSize(0);
  l2->SetFillStyle(0);
  l2->AddEntry(h1p1,"Cluster size 1","L");
  l2->AddEntry(h1p2,"Cluster size 2","L");
  l2->AddEntry(h1p3,"Cluster size 3","L");
  l2->AddEntry(h1p4,"Cluster size 4","L");
  l2->Draw();

  t1->Draw();
  t2->Draw();

  c2->SaveAs("rechitpull_clustersize.png");
  c2->SaveAs("rechitpull_clustersize.pdf");

}
