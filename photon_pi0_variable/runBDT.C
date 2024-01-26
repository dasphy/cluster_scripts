#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TFile.h"
#include "TTree.h"

void bdt() {
    // Create a factory object
    TMVA::Factory factory("TMVA_BDT", TFile::Open("TMVA_BDT.root", "RECREATE"), "AnalysisType=Classification");
    // Create a DataLoader
    TMVA::DataLoader loader("dataset");

    //// Add variables to the DataLoader
    loader.AddVariable("E_cluster", 'F');
    loader.AddVariable("E_frac_layer0", 'F');
    loader.AddVariable("E_frac_layer1", 'F');
    loader.AddVariable("E_frac_layer2", 'F');
    loader.AddVariable("E_frac_layer3", 'F');
    loader.AddVariable("E_frac_layer4", 'F');
    loader.AddVariable("E_frac_layer5", 'F');
    loader.AddVariable("E_frac_layer6", 'F');
    loader.AddVariable("E_frac_layer7", 'F');
    loader.AddVariable("E_frac_layer8", 'F');
    loader.AddVariable("E_frac_layer9", 'F');
    loader.AddVariable("E_frac_layer10", 'F');
    loader.AddVariable("E_frac_layer11", 'F');
    loader.AddVariable("theta_diff_layer0", 'F');
    loader.AddVariable("theta_diff_layer1", 'F');
    loader.AddVariable("theta_diff_layer2", 'F');
    loader.AddVariable("theta_diff_layer3", 'F');
    loader.AddVariable("theta_diff_layer4", 'F');
    loader.AddVariable("theta_diff_layer5", 'F');
    loader.AddVariable("theta_diff_layer6", 'F');
    loader.AddVariable("theta_diff_layer7", 'F');
    loader.AddVariable("theta_diff_layer8", 'F');
    loader.AddVariable("theta_diff_layer9", 'F');
    loader.AddVariable("theta_diff_layer10", 'F');
    loader.AddVariable("theta_diff_layer11", 'F');
    loader.AddVariable("E_ratio_layer0", 'F');
    loader.AddVariable("E_ratio_layer1", 'F');
    loader.AddVariable("E_ratio_layer2", 'F');
    loader.AddVariable("E_ratio_layer3", 'F');
    loader.AddVariable("E_ratio_layer4", 'F');
    loader.AddVariable("E_ratio_layer5", 'F');
    loader.AddVariable("E_ratio_layer6", 'F');
    loader.AddVariable("E_ratio_layer7", 'F');
    loader.AddVariable("E_ratio_layer8", 'F');
    loader.AddVariable("E_ratio_layer9", 'F');
    loader.AddVariable("E_ratio_layer10", 'F');
    loader.AddVariable("E_ratio_layer11", 'F');
    loader.AddVariable("delta_E_layer0", 'F');
    loader.AddVariable("delta_E_layer1", 'F');
    loader.AddVariable("delta_E_layer2", 'F');
    loader.AddVariable("delta_E_layer3", 'F');
    loader.AddVariable("delta_E_layer4", 'F');
    loader.AddVariable("delta_E_layer5", 'F');
    loader.AddVariable("delta_E_layer6", 'F');
    loader.AddVariable("delta_E_layer7", 'F');
    loader.AddVariable("delta_E_layer8", 'F');
    loader.AddVariable("delta_E_layer9", 'F');
    loader.AddVariable("delta_E_layer10", 'F');
    loader.AddVariable("delta_E_layer11", 'F');
    loader.AddVariable("w_theta_layer0", 'F');
    loader.AddVariable("w_theta_layer1", 'F');
    loader.AddVariable("w_theta_layer2", 'F');
    loader.AddVariable("w_theta_layer3", 'F');
    loader.AddVariable("w_theta_layer4", 'F');
    loader.AddVariable("w_theta_layer5", 'F');
    loader.AddVariable("w_theta_layer6", 'F');
    loader.AddVariable("w_theta_layer7", 'F');
    loader.AddVariable("w_theta_layer8", 'F');
    loader.AddVariable("w_theta_layer9", 'F');
    loader.AddVariable("w_theta_layer10", 'F');
    loader.AddVariable("w_theta_layer11", 'F');
    loader.AddVariable("w_module_layer0", 'F');
    loader.AddVariable("w_module_layer1", 'F');
    loader.AddVariable("w_module_layer2", 'F');
    loader.AddVariable("w_module_layer3", 'F');
    loader.AddVariable("w_module_layer4", 'F');
    loader.AddVariable("w_module_layer5", 'F');
    loader.AddVariable("w_module_layer6", 'F');
    loader.AddVariable("w_module_layer7", 'F');
    loader.AddVariable("w_module_layer8", 'F');
    loader.AddVariable("w_module_layer9", 'F');
    loader.AddVariable("w_module_layer10", 'F');
    loader.AddVariable("w_module_layer11", 'F');

    // Open your TTree files and get the TTree objects
    TFile *signalFile     = TFile::Open("../80_100_GeV/gamma_variable.root");
    TFile *backgroundFile = TFile::Open("../80_100_GeV/pi0_variable.root");
    TTree *signalTree     = dynamic_cast<TTree*>(signalFile->Get("variable"));
    TTree *backgroundTree = dynamic_cast<TTree*>(backgroundFile->Get("variable"));

    // Prepare the training and testing trees
    loader.AddSignalTree(signalTree, 1.0); // Signal weight is set to 1.0
    loader.AddBackgroundTree(backgroundTree, 1.0); // Background weight is set to 1.0
    loader.PrepareTrainingAndTestTree("", "SplitMode=Random");

    // Book a BDT method
    factory.BookMethod(&loader, TMVA::Types::kBDT, "BDT");

    // Train the BDT
    factory.TrainAllMethods();
    // Evaluate the BDT
    factory.TestAllMethods();
    factory.EvaluateAllMethods();
    // Save the output
    //factory.Close();
    // Close your TFile objects
    signalFile->Close();
    backgroundFile->Close();
}

void runBDT() {
    bdt();
}
