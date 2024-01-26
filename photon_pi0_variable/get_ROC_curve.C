#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>

void DrawROCCurve(TString inputFileName, TString treeName, TString methodName, TString outputFileName) {
    // Open the ROOT file containing the output of the TMVA training
    TFile inputFile(inputFileName);

    // Get the test tree
    TTree* testTree = dynamic_cast<TTree*>(inputFile.Get(treeName));
    if (!testTree) {
        std::cerr << "Test tree not found!" << std::endl;
        return;
    }

    // Get the BDT response from the test tree
    Float_t bdtResponse;
    testTree->SetBranchAddress(methodName, &bdtResponse);

    // Create arrays to store true positive rate (TPR) and false positive rate (FPR)
    Float_t tpr[1000], fpr[1000];
    Int_t nPoints = 1000;

    // Loop over different BDT response thresholds
    for (Int_t i = 0; i < nPoints; ++i) {
        //if (nPoints%100==0)  std::cout << "Processing number of points : " << nPoints << std::endl;
        Float_t threshold = static_cast<Float_t>(i) / nPoints;

        // Initialize counters
        Int_t truePositive = 0, falsePositive = 0, trueNegative = 0, falseNegative = 0;

        // Loop over events in the test tree
        Long64_t nEntries = testTree->GetEntries();
        for (Long64_t j = 0; j < nEntries; ++ j) {
            testTree->GetEntry(j);

            // Classify events based on the BDT response and threshold
            Bool_t isSignal     = (bdtResponse > threshold);
            Bool_t isTrueSignal = (testTree->GetBranch("classID")->GetLeaf("classID")->GetValue() == 0);

            // Update counters
            if (isTrueSignal && isSignal) truePositive++;
            else if (!isTrueSignal && isSignal) falsePositive++;
            else if (isTrueSignal && !isSignal) falseNegative++;
            else trueNegative++;

        }

        // Calculate TPR and FPR for the current threshold
        tpr[i] = static_cast<Float_t>(truePositive) / (truePositive + falseNegative);
        fpr[i] = 1. - static_cast<Float_t>(falsePositive) / (falsePositive + trueNegative);
        //tpr[i] = static_cast<Float_t>(truePositive) / (truePositive + falseNegative);
        //fpr[i] = static_cast<Float_t>(falsePositive) / (falsePositive + trueNegative);

        if (tpr[i]>0.965&&tpr[i]<0.975)  std::cout << "when sig eff is 0.97 bkg rej is : " << fpr[i] << std::endl;

    }

    // Create a TGraph for the ROC curve
    TGraph* graph = new TGraph(nPoints, tpr, fpr);

    // Create a canvas and draw the ROC curve
    TCanvas* canvas = new TCanvas("canvas", "ROC Curve", 1000, 600);
    //graph->SetLineColor(kBlue);
    //graph->SetLineWidth(2);
    graph->SetMarkerStyle(8);
    graph->SetMarkerColor(kRed);
    graph->SetTitle("BDT ROC Curve");
    graph->GetXaxis()->SetTitle("Signal efficiency");
    graph->GetYaxis()->SetTitle("Background rejection (1-efficiency)");
    graph->Draw("AP*");

    // Save the ROC curve as a PDF file
    canvas->SaveAs(outputFileName);

    TFile outputFile("roc_test.root", "RECREATE");
    //TFile outputFile("bdt_ROC_TGraph.root", "RECREATE");
    graph->Write("bdt_ROC_curve");
    outputFile.Close();

    // Clean up
    delete canvas;
    delete graph;
}

int get_ROC_curve() {
    // Provide the input ROOT file name, test tree name, BDT method name, and output PDF file name
    TString inputFileName = "./TMVA_BDT.root";
    TString treeName = "dataset/TestTree";
    TString methodName = "BDT";
    TString outputFileName = "./roc_curve.png";

    // Draw the ROC curve
    DrawROCCurve(inputFileName, treeName, methodName, outputFileName);

    return 0;
}

