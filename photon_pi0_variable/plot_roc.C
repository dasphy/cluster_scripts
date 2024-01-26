#include <iostream>
#include <string>
using namespace std;

void plot_roc() {

  TFile  *f1 = new TFile("../5_20_GeV/bdt_ROC_TGraph.root");
  TGraph *g1 = (TGraph*)f1->Get("bdt_ROC_curve");
  TFile  *f2 = new TFile("../20_40_GeV/bdt_ROC_TGraph.root");
  TGraph *g2 = (TGraph*)f2->Get("bdt_ROC_curve");
  TFile  *f3 = new TFile("../40_60_GeV/bdt_ROC_TGraph.root");
  TGraph *g3 = (TGraph*)f3->Get("bdt_ROC_curve");
  TFile  *f4 = new TFile("../60_80_GeV/bdt_ROC_TGraph.root");
  TGraph *g4 = (TGraph*)f4->Get("bdt_ROC_curve");
  TFile  *f5 = new TFile("../80_100_GeV/bdt_ROC_TGraph.root");
  TGraph *g5 = (TGraph*)f5->Get("bdt_ROC_curve");

  TCanvas *mc = new TCanvas("mc", "", 0, 0, 1000, 600);
  mc->SetGrid();
  //mc->SetLogy();

  g1->SetLineColor(kBlack); g1->SetLineWidth(4);
  g2->SetLineColor(kRed); g2->SetLineWidth(4);
  g3->SetLineColor(kBlue); g3->SetLineWidth(4);
  g4->SetLineColor(kGreen+2); g4->SetLineWidth(4);
  g5->SetLineColor(kOrange+2); g5->SetLineWidth(4);

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("BDT ROC Curve");
  mg->Add(g1);
  mg->Add(g2);
  mg->Add(g3);
  mg->Add(g4);
  mg->Add(g5);
  mg->GetXaxis()->SetRangeUser(.75, 1.);
  mg->GetXaxis()->SetTitle("Signal efficiency");
  mg->GetYaxis()->SetTitle("Background Rejection (1-efficiency)");

  mg->SetMinimum(.95);
  mg->Draw("AC");

  TLegend *leg = new TLegend(0.1783567,0.226087,0.3687375,0.426087);
  leg->SetLineColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->AddEntry(g1, "5-20 GeV", "l");
  leg->AddEntry(g2, "20-40 GeV", "l");
  leg->AddEntry(g3, "40-60 GeV", "l");
  leg->AddEntry(g4, "60-80 GeV", "l");
  leg->AddEntry(g5, "80-100 GeV", "l");
  leg->SetTextSize(0.035);
  leg->Draw();

  float E[] = {12.5, 30, 50, 70, 90};
  float bkg_rej[] = {1., 0.998845, 0.996878, 0.978121, 0.948538};
  
  TGraph *gr = new TGraph();
  for (int i = 0; i < 5; i ++) {
    gr->SetPoint(i, E[i], bkg_rej[i]);
  }

  TCanvas *c2 = new TCanvas("c2", "", 0, 0, 1000, 600);
  c2->SetGrid();

  gr->SetTitle("Background rejection at signal eff 97%");
  gr->SetMarkerColor(kBlue+1);
  gr->SetMarkerStyle(34);
  gr->SetMarkerSize(2.5);
  gr->GetXaxis()->SetTitle("Truth Particle Energy [GeV]");
  gr->GetYaxis()->SetTitle("Bkg rejection at signal eff 97%");
  gr->SetMinimum(0.91);
  gr->SetMaximum(1.01);
  gr->Draw("AP");

}
