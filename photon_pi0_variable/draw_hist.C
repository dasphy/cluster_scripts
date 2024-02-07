#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TObjArray.h>

void draw_hist() {
    TFile *file1 = new TFile("./240207/gamma_variable.root");
    TFile *file2 = new TFile("./240207/pi0_variable.root");

    TTree *tree1 = dynamic_cast<TTree*>(file1->Get("variable"));
    TTree *tree2 = dynamic_cast<TTree*>(file2->Get("variable"));

    // get all leaves
    TObjArray *leaves1 = tree1->GetListOfLeaves();
    TObjArray *leaves2 = tree2->GetListOfLeaves();
    //cout << "leaves1->GetEntries() is : " << leaves1->GetEntries() << endl;
    //cout << "leaves2->GetEntries() is : " << leaves2->GetEntries() << endl;

    // 0     E_cluster
    // 1-12  E_frac
    // 13-24 w_theta
    // 25-36 w_module
    // 37-48 theta_diff
    // 49-60 E_ratio
    // 61-72 delta_E

    int nbin = 200;
    float xmin = 4000;
    float xmax = 102000;

    TString _bkg = "_bkg";

    for (Int_t i = 0; i < leaves1->GetEntries(); ++ i) {
        TLeaf *leaf1 = dynamic_cast<TLeaf*>(leaves1->At(i));
        TLeaf *leaf2 = dynamic_cast<TLeaf*>(leaves2->At(i));
        // get var names
        TString variableName1 = leaf1->GetName();
        TString variableName2 = leaf2->GetName() + _bkg;

        // create histo
        TH1F *histogram1 = new TH1F(variableName1.Data(), Form("%s", variableName1.Data()), nbin, xmin, xmax);
        TH1F *histogram2 = new TH1F(variableName2.Data(), Form("%s", variableName2.Data()), nbin, xmin, xmax);
        // fill histo
        tree1->Project(variableName1.Data(), variableName1.Data());
        tree2->Project(variableName2.Data(), variableName1.Data());  // use "variableName1.Data" in the second para

        // draw logY
        TCanvas *canvas = new TCanvas("canvas", "", 1000, 650);
        canvas->SetLogy();
        histogram1->SetLineColor(kBlue);
        histogram2->SetLineColor(kRed);
        histogram1->SetLineWidth(2);
        histogram2->SetLineWidth(2);
        histogram1->SetStats(0);
        histogram2->SetStats(0);
        histogram1->GetXaxis()->SetTitle(variableName1.Data());
        histogram1->Draw();
        histogram2->Draw("sames");
        TLegend *leg = new TLegend(0.6725,0.7248,0.7807,0.8704);
        leg->SetLineColor(kWhite);
        leg->SetFillColor(kWhite);
        leg->AddEntry(histogram1, "photon", "l");
        leg->AddEntry(histogram2, "pi 0", "l");
        leg->SetTextSize(0.035);
        leg->Draw();
        canvas->SaveAs(Form("%s_logY.pdf", variableName1.Data()));

        // draw linear
        TCanvas *canvas2 = new TCanvas("canvas2", "", 1000, 650);
        histogram1->GetYaxis()->SetRangeUser(0, 1.35 * histogram1->GetMaximum());
        histogram1->Draw();
        histogram2->Draw("sames");
        leg->Draw();
        canvas2->SaveAs(Form("%s.pdf", variableName1.Data()));

        delete canvas;
        delete canvas2;

    }
}
