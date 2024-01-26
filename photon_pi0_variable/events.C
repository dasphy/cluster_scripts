#define events_cxx
#include <TFile.h>
#include <TTree.h>
#include "events.h"
#include <TH2.h>
#include <TStyle.h>

void events::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

void events::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();
}

Bool_t events::Process(Long64_t entry)
{
   fReader.SetLocalEntry(entry);

   if (entry % 1000 == 0) {
     std::cout << "------ Processing event : " << entry << " ------" << std::endl;
   }

   // loop over clusters
   for (unsigned int i = 0; i < CorrectedCaloTopoClusters_position_x.GetSize(); i++) {
     E_cluster = 1000*CorrectedCaloTopoClusters_energy[i];  // switch E to MeV
     if (E_cluster < 1000*ClusterEnergyThreshold)  continue;

     //if (E_cluster < 60000 || E_cluster > 80000)  continue;
     //if (E_cluster < 80000)    continue;

     // initialize them for a new cluster
     std::vector<float> E_layer(12, 0.0);
     // for calculation of width in theta/module
     std::vector<float> theta2_E_layer(12, 0.0);
     std::vector<float> theta_E_layer(12, 0.0);
     std::vector<float> module2_E_layer(12, 0.0);
     std::vector<float> module_E_layer(12, 0.0);
     // E of cell in each layer
     std::vector<std::vector<float>> vec_E_cell_layer(12, std::vector<float>(1, -1.0));
     std::vector<std::vector<float>> vec_theta_cell_layer(12, std::vector<float>(1, -1.0));
     std::vector<std::vector<float>> vec_module_cell_layer(12, std::vector<float>(1, -1.0));

     // loop over all positioned cells
     for (unsigned int j = 0; j < PositionedCaloTopoClusterCells_position_x.GetSize(); j ++) {
       int layer  = Layer(PositionedCaloTopoClusterCells_cellID[j]);
       int theta  = ThetaBin(PositionedCaloTopoClusterCells_cellID[j]);
       int module = Module(PositionedCaloTopoClusterCells_cellID[j]);

       vec_E_cell_layer[layer].push_back(1000*PositionedCaloTopoClusterCells_energy[j]);
       vec_theta_cell_layer[layer].push_back(theta);
       vec_module_cell_layer[layer].push_back(module);

     }  // end of loop over all positioned cells

     vec_E_cluster.push_back(E_cluster);
     std::vector<float> E_cell_Max(12, 0.);
     std::vector<float> E_cell_secMax(12, 0.);
     std::vector<float> E_cell_Min(12, std::numeric_limits<float>::max());
     // module IDs are from 0 to 1535
     std::vector<float> module_cell_Max(12, -1);
     std::vector<float> module_cell_Min(12, 1536);

     // loop over 12 layers
     for (unsigned int layer = 0; layer < numLayers; layer ++) {
       // loop over cells
       for (unsigned int i = 0; i < vec_E_cell_layer[layer].size(); i ++) {
         // find the max and secmax E in this layer
         if (vec_E_cell_layer[layer][i] > E_cell_Max[layer]) {
           E_cell_secMax[layer] = E_cell_Max[layer];
           E_cell_Max[layer] = vec_E_cell_layer[layer][i];
         } else if (vec_E_cell_layer[layer][i] > E_cell_secMax[layer] && vec_E_cell_layer[layer][i] < E_cell_Max[layer]) {
           E_cell_secMax[layer] = vec_E_cell_layer[layer][i];
         }
         // find the min E in this layer
         if (vec_E_cell_layer[layer][i] != -1.0 && vec_E_cell_layer[layer][i] < E_cell_Min[layer]) {
           E_cell_Min[layer] = vec_E_cell_layer[layer][i];
         }
         // find the max and min module ID in this layer (module IDs are from 0 to 1535)
         if (vec_module_cell_layer[layer][i] > module_cell_Max[layer]) {
           module_cell_Max[layer] = vec_module_cell_layer[layer][i];
         }
         if (vec_module_cell_layer[layer][i] < module_cell_Min[layer] && vec_module_cell_layer[layer][i] != -1.0) {
           module_cell_Min[layer] = vec_module_cell_layer[layer][i];
         }
       }  // end of loop over cells in this layer
       // check min E_cell
       if (E_cell_Min[layer] > 1e15)    E_cell_Min[layer] = 0.;

       // reset module ID due to module periodicity
       if ((module_cell_Max[layer] - module_cell_Min[layer]) > 1500) {
         for (unsigned int i = 0; i < vec_module_cell_layer[layer].size(); i ++) {
           if (vec_module_cell_layer[layer][i] > 1536/2.)    vec_module_cell_layer[layer][i] = vec_module_cell_layer[layer][i] - 1535;
         }
         std::cout << "------ Reset module ID in layer : " << layer << std::endl;
       }

       // sum these vars over cells in each layer
       for (unsigned int i = 0; i < vec_E_cell_layer[layer].size(); i ++) {
         if (vec_E_cell_layer[layer][i] == -1) {
           vec_E_cell_layer[layer][i] = 1e-9;
           vec_theta_cell_layer[layer][i] = 0.;
           vec_theta_cell_layer[layer][i] = 0.;
         }
         E_layer[layer]         = E_layer[layer]         + vec_E_cell_layer[layer][i];
         theta2_E_layer[layer]  = theta2_E_layer[layer]  + vec_theta_cell_layer[layer][i]  * vec_theta_cell_layer[layer][i]  * vec_E_cell_layer[layer][i];
         theta_E_layer[layer]   = theta_E_layer[layer]   + vec_theta_cell_layer[layer][i]  * vec_E_cell_layer[layer][i];
         module2_E_layer[layer] = module2_E_layer[layer] + vec_module_cell_layer[layer][i] * vec_module_cell_layer[layer][i] * vec_E_cell_layer[layer][i];
         module_E_layer[layer]  = module_E_layer[layer]  + vec_module_cell_layer[layer][i] * vec_E_cell_layer[layer][i];
       }

       vec_E_layer[layer].push_back(E_layer[layer]);
       vec_E_frac_layer[layer].push_back(E_layer[layer]/E_cluster);

      //  std::cout << "E_layer in layer " << layer << " is : " << E_layer[layer] << std::endl;
      //  std::cout << "E_cluster is : " << E_cluster << std::endl;
      //  std::cout << "vec_E_frac_layer in layer " << layer << " is : " << E_layer[layer]/E_cluster << std::endl;
      //  std::cout << ":::::::::::::::::: " << std::endl;

       vec_w_theta_layer[layer].push_back(sqrt(fabs(theta2_E_layer[layer]/E_layer[layer]-std::pow(theta_E_layer[layer]/E_layer[layer],2)))>50 ? 50 : sqrt(fabs(theta2_E_layer[layer]/E_layer[layer]-std::pow(theta_E_layer[layer]/E_layer[layer],2))));
       vec_w_module_layer[layer].push_back(sqrt(fabs(module2_E_layer[layer]/E_layer[layer]-std::pow(module_E_layer[layer]/E_layer[layer],2)))>50 ? 50 : sqrt(fabs(module2_E_layer[layer]/E_layer[layer]-std::pow(module_E_layer[layer]/E_layer[layer],2))));

       // find the index of max and secmax E
       auto it_Max      = std::find(vec_E_cell_layer[layer].begin(), vec_E_cell_layer[layer].end(), E_cell_Max[layer]);
       int index_Max    = std::distance(vec_E_cell_layer[layer].begin(), it_Max);
       auto it_secMax   = std::find(vec_E_cell_layer[layer].begin(), vec_E_cell_layer[layer].end(), E_cell_secMax[layer]);
       int index_secMax = std::distance(vec_E_cell_layer[layer].begin(), it_secMax);

       vec_theta_diff_layer[layer].push_back(fabs(vec_theta_cell_layer[layer][index_Max]-vec_theta_cell_layer[layer][index_secMax])>10 ? 10 : fabs(vec_theta_cell_layer[layer][index_Max]-vec_theta_cell_layer[layer][index_secMax]));
       vec_E_ratio_layer[layer].push_back((E_cell_Max[layer]+E_cell_secMax[layer]==0) ? 1. : (E_cell_Max[layer]-E_cell_secMax[layer])/(E_cell_Max[layer]+E_cell_secMax[layer]));
       vec_delta_E_layer[layer].push_back((E_cell_Max[layer]-E_cell_Min[layer]));

       vec_E_cell_layer[layer].clear();
       vec_theta_cell_layer[layer].clear();
       vec_module_cell_layer[layer].clear();

     }  // end of loop over 12 layers
     //std::cout << "===========================" << std::endl;

   }  // end of loop over clusters

   return kTRUE;
}

void events::SlaveTerminate()
{
}

void events::Terminate()
{
   TFile* fOut = new TFile("gamma_pi0_variable.root", "RECREATE");
   TTree *myTree = new TTree("variable", "variables for MVA");

   myTree->Branch("E_cluster",        &_E_cluster,        "_E_cluster/F");
   myTree->Branch("E_frac_layer0",    &_E_frac_layer0,    "_E_frac_layer0/F");
   myTree->Branch("E_frac_layer1",    &_E_frac_layer1,    "_E_frac_layer1/F");
   myTree->Branch("E_frac_layer2",    &_E_frac_layer2,    "_E_frac_layer2/F");
   myTree->Branch("E_frac_layer3",    &_E_frac_layer3,    "_E_frac_layer3/F");
   myTree->Branch("E_frac_layer4",    &_E_frac_layer4,    "_E_frac_layer4/F");
   myTree->Branch("E_frac_layer5",    &_E_frac_layer5,    "_E_frac_layer5/F");
   myTree->Branch("E_frac_layer6",    &_E_frac_layer6,    "_E_frac_layer6/F");
   myTree->Branch("E_frac_layer7",    &_E_frac_layer7,    "_E_frac_layer7/F");
   myTree->Branch("E_frac_layer8",    &_E_frac_layer8,    "_E_frac_layer8/F");
   myTree->Branch("E_frac_layer9",    &_E_frac_layer9,    "_E_frac_layer9/F");
   myTree->Branch("E_frac_layer10",   &_E_frac_layer10,   "_E_frac_layer10/F");
   myTree->Branch("E_frac_layer11",   &_E_frac_layer11,   "_E_frac_layer11/F");
   myTree->Branch("w_theta_layer0",   &_w_theta_layer0,   "_w_theta_layer0/F");
   myTree->Branch("w_theta_layer1",   &_w_theta_layer1,   "_w_theta_layer1/F");
   myTree->Branch("w_theta_layer2",   &_w_theta_layer2,   "_w_theta_layer2/F");
   myTree->Branch("w_theta_layer3",   &_w_theta_layer3,   "_w_theta_layer3/F");
   myTree->Branch("w_theta_layer4",   &_w_theta_layer4,   "_w_theta_layer4/F");
   myTree->Branch("w_theta_layer5",   &_w_theta_layer5,   "_w_theta_layer5/F");
   myTree->Branch("w_theta_layer6",   &_w_theta_layer6,   "_w_theta_layer6/F");
   myTree->Branch("w_theta_layer7",   &_w_theta_layer7,   "_w_theta_layer7/F");
   myTree->Branch("w_theta_layer8",   &_w_theta_layer8,   "_w_theta_layer8/F");
   myTree->Branch("w_theta_layer9",   &_w_theta_layer9,   "_w_theta_layer9/F");
   myTree->Branch("w_theta_layer10",  &_w_theta_layer10,  "_w_theta_layer10/F");
   myTree->Branch("w_theta_layer11",  &_w_theta_layer11,  "_w_theta_layer11/F");
   myTree->Branch("w_module_layer0",  &_w_module_layer0,  "_w_module_layer0/F");
   myTree->Branch("w_module_layer1",  &_w_module_layer1,  "_w_module_layer1/F");
   myTree->Branch("w_module_layer2",  &_w_module_layer2,  "_w_module_layer2/F");
   myTree->Branch("w_module_layer3",  &_w_module_layer3,  "_w_module_layer3/F");
   myTree->Branch("w_module_layer4",  &_w_module_layer4,  "_w_module_layer4/F");
   myTree->Branch("w_module_layer5",  &_w_module_layer5,  "_w_module_layer5/F");
   myTree->Branch("w_module_layer6",  &_w_module_layer6,  "_w_module_layer6/F");
   myTree->Branch("w_module_layer7",  &_w_module_layer7,  "_w_module_layer7/F");
   myTree->Branch("w_module_layer8",  &_w_module_layer8,  "_w_module_layer8/F");
   myTree->Branch("w_module_layer9",  &_w_module_layer9,  "_w_module_layer9/F");
   myTree->Branch("w_module_layer10", &_w_module_layer10, "_w_module_layer10/F");
   myTree->Branch("w_module_layer11", &_w_module_layer11, "_w_module_layer11/F");
   myTree->Branch("theta_diff_layer0",   &_theta_diff_layer0,   "_theta_diff_layer0/F");
   myTree->Branch("theta_diff_layer1",   &_theta_diff_layer1,   "_theta_diff_layer1/F");
   myTree->Branch("theta_diff_layer2",   &_theta_diff_layer2,   "_theta_diff_layer2/F");
   myTree->Branch("theta_diff_layer3",   &_theta_diff_layer3,   "_theta_diff_layer3/F");
   myTree->Branch("theta_diff_layer4",   &_theta_diff_layer4,   "_theta_diff_layer4/F");
   myTree->Branch("theta_diff_layer5",   &_theta_diff_layer5,   "_theta_diff_layer5/F");
   myTree->Branch("theta_diff_layer6",   &_theta_diff_layer6,   "_theta_diff_layer6/F");
   myTree->Branch("theta_diff_layer7",   &_theta_diff_layer7,   "_theta_diff_layer7/F");
   myTree->Branch("theta_diff_layer8",   &_theta_diff_layer8,   "_theta_diff_layer8/F");
   myTree->Branch("theta_diff_layer9",   &_theta_diff_layer9,   "_theta_diff_layer9/F");
   myTree->Branch("theta_diff_layer10",  &_theta_diff_layer10,  "_theta_diff_layer10/F");
   myTree->Branch("theta_diff_layer11",  &_theta_diff_layer11,  "_theta_diff_layer11/F");
   myTree->Branch("E_ratio_layer0",      &_E_ratio_layer0,      "_E_ratio_layer0/F");
   myTree->Branch("E_ratio_layer1",      &_E_ratio_layer1,      "_E_ratio_layer1/F");
   myTree->Branch("E_ratio_layer2",      &_E_ratio_layer2,      "_E_ratio_layer2/F");
   myTree->Branch("E_ratio_layer3",      &_E_ratio_layer3,      "_E_ratio_layer3/F");
   myTree->Branch("E_ratio_layer4",      &_E_ratio_layer4,      "_E_ratio_layer4/F");
   myTree->Branch("E_ratio_layer5",      &_E_ratio_layer5,      "_E_ratio_layer5/F");
   myTree->Branch("E_ratio_layer6",      &_E_ratio_layer6,      "_E_ratio_layer6/F");
   myTree->Branch("E_ratio_layer7",      &_E_ratio_layer7,      "_E_ratio_layer7/F");
   myTree->Branch("E_ratio_layer8",      &_E_ratio_layer8,      "_E_ratio_layer8/F");
   myTree->Branch("E_ratio_layer9",      &_E_ratio_layer9,      "_E_ratio_layer9/F");
   myTree->Branch("E_ratio_layer10",     &_E_ratio_layer10,     "_E_ratio_layer10/F");
   myTree->Branch("E_ratio_layer11",     &_E_ratio_layer11,     "_E_ratio_layer11/F");
   myTree->Branch("delta_E_layer0",      &_delta_E_layer0,      "_delta_E_layer0/F");
   myTree->Branch("delta_E_layer1",      &_delta_E_layer1,      "_delta_E_layer1/F");
   myTree->Branch("delta_E_layer2",      &_delta_E_layer2,      "_delta_E_layer2/F");
   myTree->Branch("delta_E_layer3",      &_delta_E_layer3,      "_delta_E_layer3/F");
   myTree->Branch("delta_E_layer4",      &_delta_E_layer4,      "_delta_E_layer4/F");
   myTree->Branch("delta_E_layer5",      &_delta_E_layer5,      "_delta_E_layer5/F");
   myTree->Branch("delta_E_layer6",      &_delta_E_layer6,      "_delta_E_layer6/F");
   myTree->Branch("delta_E_layer7",      &_delta_E_layer7,      "_delta_E_layer7/F");
   myTree->Branch("delta_E_layer8",      &_delta_E_layer8,      "_delta_E_layer8/F");
   myTree->Branch("delta_E_layer9",      &_delta_E_layer9,      "_delta_E_layer9/F");
   myTree->Branch("delta_E_layer10",     &_delta_E_layer10,     "_delta_E_layer10/F");
   myTree->Branch("delta_E_layer11",     &_delta_E_layer11,     "_delta_E_layer11/F");

   for (unsigned int i = 0; i < vec_E_layer[0].size(); i++) {
     _E_cluster      = vec_E_cluster[i];
     _E_frac_layer0  = vec_E_frac_layer[0][i];
     _E_frac_layer1  = vec_E_frac_layer[1][i];
     _E_frac_layer2  = vec_E_frac_layer[2][i];
     _E_frac_layer3  = vec_E_frac_layer[3][i];
     _E_frac_layer4  = vec_E_frac_layer[4][i];
     _E_frac_layer5  = vec_E_frac_layer[5][i];
     _E_frac_layer6  = vec_E_frac_layer[6][i];
     _E_frac_layer7  = vec_E_frac_layer[7][i];
     _E_frac_layer8  = vec_E_frac_layer[8][i];
     _E_frac_layer9  = vec_E_frac_layer[9][i];
     _E_frac_layer10 = vec_E_frac_layer[10][i];
     _E_frac_layer11 = vec_E_frac_layer[11][i];
     _w_theta_layer0   = vec_w_theta_layer[0][i];
     _w_theta_layer1   = vec_w_theta_layer[1][i];
     _w_theta_layer2   = vec_w_theta_layer[2][i];
     _w_theta_layer3   = vec_w_theta_layer[3][i];
     _w_theta_layer4   = vec_w_theta_layer[4][i];
     _w_theta_layer5   = vec_w_theta_layer[5][i];
     _w_theta_layer6   = vec_w_theta_layer[6][i];
     _w_theta_layer7   = vec_w_theta_layer[7][i];
     _w_theta_layer8   = vec_w_theta_layer[8][i];
     _w_theta_layer9   = vec_w_theta_layer[9][i];
     _w_theta_layer10  = vec_w_theta_layer[10][i];
     _w_theta_layer11  = vec_w_theta_layer[11][i];
     _w_module_layer0  = vec_w_module_layer[0][i];
     _w_module_layer1  = vec_w_module_layer[1][i];
     _w_module_layer2  = vec_w_module_layer[2][i];
     _w_module_layer3  = vec_w_module_layer[3][i];
     _w_module_layer4  = vec_w_module_layer[4][i];
     _w_module_layer5  = vec_w_module_layer[5][i];
     _w_module_layer6  = vec_w_module_layer[6][i];
     _w_module_layer7  = vec_w_module_layer[7][i];
     _w_module_layer8  = vec_w_module_layer[8][i];
     _w_module_layer9  = vec_w_module_layer[9][i];
     _w_module_layer10 = vec_w_module_layer[10][i];
     _w_module_layer11 = vec_w_module_layer[11][i];
     _theta_diff_layer0   = vec_theta_diff_layer[0][i];
     _theta_diff_layer1   = vec_theta_diff_layer[1][i];
     _theta_diff_layer2   = vec_theta_diff_layer[2][i];
     _theta_diff_layer3   = vec_theta_diff_layer[3][i];
     _theta_diff_layer4   = vec_theta_diff_layer[4][i];
     _theta_diff_layer5   = vec_theta_diff_layer[5][i];
     _theta_diff_layer6   = vec_theta_diff_layer[6][i];
     _theta_diff_layer7   = vec_theta_diff_layer[7][i];
     _theta_diff_layer8   = vec_theta_diff_layer[8][i];
     _theta_diff_layer9   = vec_theta_diff_layer[9][i];
     _theta_diff_layer10  = vec_theta_diff_layer[10][i];
     _theta_diff_layer11  = vec_theta_diff_layer[11][i];
     _E_ratio_layer0   = vec_E_ratio_layer[0][i];
     _E_ratio_layer1   = vec_E_ratio_layer[1][i];
     _E_ratio_layer2   = vec_E_ratio_layer[2][i];
     _E_ratio_layer3   = vec_E_ratio_layer[3][i];
     _E_ratio_layer4   = vec_E_ratio_layer[4][i];
     _E_ratio_layer5   = vec_E_ratio_layer[5][i];
     _E_ratio_layer6   = vec_E_ratio_layer[6][i];
     _E_ratio_layer7   = vec_E_ratio_layer[7][i];
     _E_ratio_layer8   = vec_E_ratio_layer[8][i];
     _E_ratio_layer9   = vec_E_ratio_layer[9][i];
     _E_ratio_layer10  = vec_E_ratio_layer[10][i];
     _E_ratio_layer11  = vec_E_ratio_layer[11][i];
     _delta_E_layer0   = vec_delta_E_layer[0][i];
     _delta_E_layer1   = vec_delta_E_layer[1][i];
     _delta_E_layer2   = vec_delta_E_layer[2][i];
     _delta_E_layer3   = vec_delta_E_layer[3][i];
     _delta_E_layer4   = vec_delta_E_layer[4][i];
     _delta_E_layer5   = vec_delta_E_layer[5][i];
     _delta_E_layer6   = vec_delta_E_layer[6][i];
     _delta_E_layer7   = vec_delta_E_layer[7][i];
     _delta_E_layer8   = vec_delta_E_layer[8][i];
     _delta_E_layer9   = vec_delta_E_layer[9][i];
     _delta_E_layer10  = vec_delta_E_layer[10][i];
     _delta_E_layer11  = vec_delta_E_layer[11][i];

     myTree->Fill();
   }

   myTree->Write();
   fOut->Close();

   // loop over layers
   for (unsigned int layer = 0; layer < numLayers; layer ++) {
     vec_layer[layer] = layer * 1.0;
     // get mean and std_dev
     mean_E_layer[layer] = std::accumulate(vec_E_layer[layer].begin(), vec_E_layer[layer].end(), 0) / vec_E_layer[layer].size();
     float vec_E_layer_sumSquaredDifferences = 0.0;
     for (float value : vec_E_layer[layer]) {
       float difference = value - mean_E_layer[layer];
       vec_E_layer_sumSquaredDifferences += difference * difference;
     }
     // std_dev as error
     mean_E_layer_error[layer] = std::sqrt(vec_E_layer_sumSquaredDifferences / vec_E_layer[layer].size());

     float sum_E_frac = 0.;
     for (float _E_frac : vec_E_frac_layer[layer]) {
       sum_E_frac = sum_E_frac + _E_frac;
       //std::cout << "vec_E_frac_layer in layer " << layer << " is : " << rt << std::endl;
     }
     mean_E_frac_layer[layer] = sum_E_frac / vec_E_frac_layer[layer].size();
     // ?? mean_E_frac_layer[layer] = std::accumulate(vec_E_frac_layer[layer].begin(), vec_E_frac_layer[layer].end(), 0) / vec_E_frac_layer[layer].size();
     float vec_E_frac_layer_sumSquaredDifferences = 0.0;
     for (float value : vec_E_frac_layer[layer]) {
       float difference = value - mean_E_frac_layer[layer];
       vec_E_frac_layer_sumSquaredDifferences += difference * difference;
     }
     mean_E_frac_layer_error[layer] = std::sqrt(vec_E_frac_layer_sumSquaredDifferences / vec_E_frac_layer[layer].size());

     mean_w_theta_layer[layer] = std::accumulate(vec_w_theta_layer[layer].begin(), vec_w_theta_layer[layer].end(), 0) / vec_w_theta_layer[layer].size();
     float vec_w_theta_layer_sumSquaredDifferences = 0.0;
     for (float value : vec_w_theta_layer[layer]) {
       float difference = value - mean_w_theta_layer[layer];
       vec_w_theta_layer_sumSquaredDifferences += difference * difference;
     }
     mean_w_theta_layer_error[layer] = std::sqrt(vec_w_theta_layer_sumSquaredDifferences / vec_w_theta_layer[layer].size());

     mean_w_module_layer[layer] = std::accumulate(vec_w_module_layer[layer].begin(), vec_w_module_layer[layer].end(), 0) / vec_w_module_layer[layer].size();
     float vec_w_module_layer_sumSquaredDifferences = 0.0;
     for (float value : vec_w_module_layer[layer]) {
       float difference = value - mean_w_module_layer[layer];
       vec_w_module_layer_sumSquaredDifferences += difference * difference;
     }
     mean_w_module_layer_error[layer] = std::sqrt(vec_w_module_layer_sumSquaredDifferences / vec_w_module_layer[layer].size());

     std::cout << "mean_E_layer in layer_" << layer << " is : " << mean_E_layer[layer] << std::endl;
     std::cout << "mean_E_layer_error in layer_" << layer << " is : " << mean_E_layer_error[layer] << std::endl;
     std::cout << "mean_E_frac_layer in layer_" << layer << " is : " << mean_E_frac_layer[layer] << std::endl;
     std::cout << "mean_E_frac_layer_error in layer_" << layer << " is : " << mean_E_frac_layer_error[layer] << std::endl;
     std::cout << "mean_w_theta_layer in layer_" << layer << " is : " << mean_w_theta_layer[layer] << std::endl;
     std::cout << "mean_w_theta_layer_error in layer_" << layer << " is : " << mean_w_theta_layer_error[layer] << std::endl;
     std::cout << "mean_w_module_layer in layer_" << layer << " is : " << mean_w_module_layer[layer] << std::endl;
     std::cout << "mean_w_module_layer_error in layer_" << layer << " is : " << mean_w_module_layer_error[layer] << std::endl;
     std::cout << "===========================" << std::endl;
   }

   float max_mean_E_layer        = std::numeric_limits<float>::min();
   float max_mean_E_frac_layer   = std::numeric_limits<float>::min();
   float max_mean_w_theta_layer  = std::numeric_limits<float>::min();
   float max_mean_w_module_layer = std::numeric_limits<float>::min();
   for (float E : mean_E_layer) {
     if (E > max_mean_E_layer) {  max_mean_E_layer = E;  }
   }
   for (float E_f : mean_E_frac_layer) {
     if (E_f > max_mean_E_frac_layer) {  max_mean_E_frac_layer = E_f;  }
   }
   for (float w_t : mean_w_theta_layer) {
     if (w_t > max_mean_w_theta_layer) {  max_mean_w_theta_layer = w_t;  }
   }
   for (float w_m : mean_w_module_layer) {
     if (w_m > max_mean_w_module_layer) {  max_mean_w_module_layer = w_m;  }
   }

   TCanvas* c1 = new TCanvas("canvas1", "canvas1", 1000, 600);
   c1->SetLeftMargin(0.13);
   c1->SetRightMargin(0.06);
   TGraphErrors gr1(vec_layer.size(), &vec_layer[0], &mean_E_layer[0], 0, &mean_E_layer_error[0]);
   gr1.SetMarkerStyle(8);
   gr1.SetMarkerColor(kRed);
   TH2D* axor1 = new TH2D("axor1", "", 100, -0.5, 11.5, 100, -max_mean_E_layer*.1, max_mean_E_layer*1.5);
   axor1->SetDirectory(0);
   axor1->GetXaxis()->SetNdivisions(120);
   axor1->GetXaxis()->SetTitle("Layer");
   axor1->GetYaxis()->SetTitle("Energy in Layer [MeV]");
   axor1->SetStats(0);
   axor1->Draw();
   gr1.DrawClone("SAME P");
   c1->Print("./c1.png");

   TCanvas* c4 = new TCanvas("canvas4", "canvas4", 1000, 600);
   c4->SetLeftMargin(0.13);
   c4->SetRightMargin(0.06);
   TGraphErrors gr4(vec_layer.size(), &vec_layer[0], &mean_E_frac_layer[0], 0, &mean_E_frac_layer_error[0]);
   gr4.SetMarkerStyle(8);
   gr4.SetMarkerColor(kRed);
   TH2D* axor4 = new TH2D("axor4", "", 100, -0.5, 11.5, 100, -max_mean_E_frac_layer*.1, max_mean_E_frac_layer*1.5);
   axor4->SetDirectory(0);
   axor4->GetXaxis()->SetNdivisions(120);
   axor4->GetXaxis()->SetTitle("Layer");
   axor4->GetYaxis()->SetTitle("Energy Fraction in Layer");
   axor4->SetStats(0);
   axor4->Draw();
   gr4.DrawClone("SAME P");
   c4->Print("./c4.png");

   TCanvas* c2 = new TCanvas("canvas2", "canvas2", 1000, 600);
   c2->SetLeftMargin(0.13);
   c2->SetRightMargin(0.06);
   TGraphErrors gr2(vec_layer.size(), &vec_layer[0], &mean_w_theta_layer[0], 0, &mean_w_theta_layer_error[0]);
   gr2.SetMarkerStyle(8);
   gr2.SetMarkerColor(kRed);
   TH2D* axor2 = new TH2D("axor2", "", 100, -0.5, 11.5, 100, -max_mean_w_theta_layer*.5, max_mean_w_theta_layer*2.5);
   axor2->SetDirectory(0);
   axor2->GetXaxis()->SetNdivisions(120);
   axor2->GetXaxis()->SetTitle("Layer");
   axor2->GetYaxis()->SetTitle("Width in Theta");
   axor2->SetStats(0);
   axor2->Draw();
   gr2.DrawClone("SAME P");
   c2->Print("./c2.png");

   TCanvas* c3 = new TCanvas("canvas3", "canvas3", 1000, 600);
   c3->SetLeftMargin(0.13);
   c3->SetRightMargin(0.06);
   TGraphErrors gr3(vec_layer.size(), &vec_layer[0], &mean_w_module_layer[0], 0, &mean_w_module_layer_error[0]);
   gr3.SetMarkerStyle(8);
   gr3.SetMarkerColor(kRed);
   TH2D* axor3 = new TH2D("axor3", "", 100, -0.5, 11.5, 100, -max_mean_w_module_layer*.5, max_mean_w_module_layer*2.5);
   axor3->SetDirectory(0);
   axor3->GetXaxis()->SetNdivisions(120);
   axor3->GetXaxis()->SetTitle("Layer");
   axor3->GetYaxis()->SetTitle("Width in Module");
   axor3->SetStats(0);
   axor3->Draw();
   gr3.DrawClone("SAME P");
   c3->Print("./c3.png");

}
