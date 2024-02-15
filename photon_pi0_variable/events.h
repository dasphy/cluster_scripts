#ifndef events_h
#define events_h
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
// Headers needed by this particular selector
#include <vector>

class events : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<ULong_t> CaloClusterCells_cellID = {fReader, "CaloClusterCells.cellID"};
   TTreeReaderArray<Float_t> CaloClusterCells_energy = {fReader, "CaloClusterCells.energy"};
   TTreeReaderArray<Float_t> CaloClusterCells_energyError = {fReader, "CaloClusterCells.energyError"};
   TTreeReaderArray<Float_t> CaloClusterCells_time = {fReader, "CaloClusterCells.time"};
   TTreeReaderArray<Float_t> CaloClusterCells_position_x = {fReader, "CaloClusterCells.position.x"};
   TTreeReaderArray<Float_t> CaloClusterCells_position_y = {fReader, "CaloClusterCells.position.y"};
   TTreeReaderArray<Float_t> CaloClusterCells_position_z = {fReader, "CaloClusterCells.position.z"};
   TTreeReaderArray<Int_t> CaloClusters_type = {fReader, "CaloClusters.type"};
   TTreeReaderArray<Float_t> CaloClusters_energy = {fReader, "CaloClusters.energy"};
   TTreeReaderArray<Float_t> CaloClusters_energyError = {fReader, "CaloClusters.energyError"};
   TTreeReaderArray<Float_t> CaloClusters_position_x = {fReader, "CaloClusters.position.x"};
   TTreeReaderArray<Float_t> CaloClusters_position_y = {fReader, "CaloClusters.position.y"};
   TTreeReaderArray<Float_t> CaloClusters_position_z = {fReader, "CaloClusters.position.z"};
   TTreeReaderArray<Float_t> CaloClusters_directionError_x = {fReader, "CaloClusters.directionError.x"};
   TTreeReaderArray<Float_t> CaloClusters_directionError_y = {fReader, "CaloClusters.directionError.y"};
   TTreeReaderArray<Float_t> CaloClusters_directionError_z = {fReader, "CaloClusters.directionError.z"};
   TTreeReaderArray<Int_t> _CaloClusters_clusters_index = {fReader, "_CaloClusters_clusters.index"};
   TTreeReaderArray<unsigned int> _CaloClusters_clusters_collectionID = {fReader, "_CaloClusters_clusters.collectionID"};
   TTreeReaderArray<Int_t> _CaloClusters_hits_index = {fReader, "_CaloClusters_hits.index"};
   TTreeReaderArray<unsigned int> _CaloClusters_hits_collectionID = {fReader, "_CaloClusters_hits.collectionID"};
   TTreeReaderArray<Int_t> _CaloClusters_particleIDs_index = {fReader, "_CaloClusters_particleIDs.index"};
   TTreeReaderArray<unsigned int> _CaloClusters_particleIDs_collectionID = {fReader, "_CaloClusters_particleIDs.collectionID"};
   TTreeReaderArray<float> _CaloClusters_shapeParameters = {fReader, "_CaloClusters_shapeParameters"};
   TTreeReaderArray<float> _CaloClusters_subdetectorEnergies = {fReader, "_CaloClusters_subdetectorEnergies"};
   TTreeReaderArray<ULong_t> CaloTopoClusterCells_cellID = {fReader, "CaloTopoClusterCells.cellID"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_energy = {fReader, "CaloTopoClusterCells.energy"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_energyError = {fReader, "CaloTopoClusterCells.energyError"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_time = {fReader, "CaloTopoClusterCells.time"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_position_x = {fReader, "CaloTopoClusterCells.position.x"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_position_y = {fReader, "CaloTopoClusterCells.position.y"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_position_z = {fReader, "CaloTopoClusterCells.position.z"};
   TTreeReaderArray<Int_t> CaloTopoClusters_type = {fReader, "CaloTopoClusters.type"};
   TTreeReaderArray<Float_t> CaloTopoClusters_energy = {fReader, "CaloTopoClusters.energy"};
   TTreeReaderArray<Float_t> CaloTopoClusters_energyError = {fReader, "CaloTopoClusters.energyError"};
   TTreeReaderArray<Float_t> CaloTopoClusters_position_x = {fReader, "CaloTopoClusters.position.x"};
   TTreeReaderArray<Float_t> CaloTopoClusters_position_y = {fReader, "CaloTopoClusters.position.y"};
   TTreeReaderArray<Float_t> CaloTopoClusters_position_z = {fReader, "CaloTopoClusters.position.z"};
   TTreeReaderArray<Float_t> CaloTopoClusters_directionError_x = {fReader, "CaloTopoClusters.directionError.x"};
   TTreeReaderArray<Float_t> CaloTopoClusters_directionError_y = {fReader, "CaloTopoClusters.directionError.y"};
   TTreeReaderArray<Float_t> CaloTopoClusters_directionError_z = {fReader, "CaloTopoClusters.directionError.z"};
   TTreeReaderArray<Int_t> _CaloTopoClusters_clusters_index = {fReader, "_CaloTopoClusters_clusters.index"};
   TTreeReaderArray<unsigned int> _CaloTopoClusters_clusters_collectionID = {fReader, "_CaloTopoClusters_clusters.collectionID"};
   TTreeReaderArray<Int_t> _CaloTopoClusters_hits_index = {fReader, "_CaloTopoClusters_hits.index"};
   TTreeReaderArray<unsigned int> _CaloTopoClusters_hits_collectionID = {fReader, "_CaloTopoClusters_hits.collectionID"};
   TTreeReaderArray<Int_t> _CaloTopoClusters_particleIDs_index = {fReader, "_CaloTopoClusters_particleIDs.index"};
   TTreeReaderArray<unsigned int> _CaloTopoClusters_particleIDs_collectionID = {fReader, "_CaloTopoClusters_particleIDs.collectionID"};
   TTreeReaderArray<float> _CaloTopoClusters_shapeParameters = {fReader, "_CaloTopoClusters_shapeParameters"};
   TTreeReaderArray<float> _CaloTopoClusters_subdetectorEnergies = {fReader, "_CaloTopoClusters_subdetectorEnergies"};
   TTreeReaderArray<Int_t> CorrectedCaloClusters_type = {fReader, "CorrectedCaloClusters.type"};
   TTreeReaderArray<Float_t> CorrectedCaloClusters_energy = {fReader, "CorrectedCaloClusters.energy"};
   TTreeReaderArray<Float_t> CorrectedCaloClusters_energyError = {fReader, "CorrectedCaloClusters.energyError"};
   TTreeReaderArray<Float_t> CorrectedCaloClusters_position_x = {fReader, "CorrectedCaloClusters.position.x"};
   TTreeReaderArray<Float_t> CorrectedCaloClusters_position_y = {fReader, "CorrectedCaloClusters.position.y"};
   TTreeReaderArray<Float_t> CorrectedCaloClusters_position_z = {fReader, "CorrectedCaloClusters.position.z"};
   TTreeReaderArray<Float_t> CorrectedCaloClusters_directionError_x = {fReader, "CorrectedCaloClusters.directionError.x"};
   TTreeReaderArray<Float_t> CorrectedCaloClusters_directionError_y = {fReader, "CorrectedCaloClusters.directionError.y"};
   TTreeReaderArray<Float_t> CorrectedCaloClusters_directionError_z = {fReader, "CorrectedCaloClusters.directionError.z"};
   TTreeReaderArray<Int_t> _CorrectedCaloClusters_clusters_index = {fReader, "_CorrectedCaloClusters_clusters.index"};
   TTreeReaderArray<unsigned int> _CorrectedCaloClusters_clusters_collectionID = {fReader, "_CorrectedCaloClusters_clusters.collectionID"};
   TTreeReaderArray<Int_t> _CorrectedCaloClusters_hits_index = {fReader, "_CorrectedCaloClusters_hits.index"};
   TTreeReaderArray<unsigned int> _CorrectedCaloClusters_hits_collectionID = {fReader, "_CorrectedCaloClusters_hits.collectionID"};
   TTreeReaderArray<Int_t> _CorrectedCaloClusters_particleIDs_index = {fReader, "_CorrectedCaloClusters_particleIDs.index"};
   TTreeReaderArray<unsigned int> _CorrectedCaloClusters_particleIDs_collectionID = {fReader, "_CorrectedCaloClusters_particleIDs.collectionID"};
   TTreeReaderArray<float> _CorrectedCaloClusters_shapeParameters = {fReader, "_CorrectedCaloClusters_shapeParameters"};
   TTreeReaderArray<float> _CorrectedCaloClusters_subdetectorEnergies = {fReader, "_CorrectedCaloClusters_subdetectorEnergies"};
   TTreeReaderArray<Int_t> CorrectedCaloTopoClusters_type = {fReader, "CorrectedCaloTopoClusters.type"};
   TTreeReaderArray<Float_t> CorrectedCaloTopoClusters_energy = {fReader, "CorrectedCaloTopoClusters.energy"};
   TTreeReaderArray<Float_t> CorrectedCaloTopoClusters_energyError = {fReader, "CorrectedCaloTopoClusters.energyError"};
   TTreeReaderArray<Float_t> CorrectedCaloTopoClusters_position_x = {fReader, "CorrectedCaloTopoClusters.position.x"};
   TTreeReaderArray<Float_t> CorrectedCaloTopoClusters_position_y = {fReader, "CorrectedCaloTopoClusters.position.y"};
   TTreeReaderArray<Float_t> CorrectedCaloTopoClusters_position_z = {fReader, "CorrectedCaloTopoClusters.position.z"};
   TTreeReaderArray<Float_t> CorrectedCaloTopoClusters_directionError_x = {fReader, "CorrectedCaloTopoClusters.directionError.x"};
   TTreeReaderArray<Float_t> CorrectedCaloTopoClusters_directionError_y = {fReader, "CorrectedCaloTopoClusters.directionError.y"};
   TTreeReaderArray<Float_t> CorrectedCaloTopoClusters_directionError_z = {fReader, "CorrectedCaloTopoClusters.directionError.z"};
   TTreeReaderArray<Int_t> _CorrectedCaloTopoClusters_clusters_index = {fReader, "_CorrectedCaloTopoClusters_clusters.index"};
   TTreeReaderArray<unsigned int> _CorrectedCaloTopoClusters_clusters_collectionID = {fReader, "_CorrectedCaloTopoClusters_clusters.collectionID"};
   TTreeReaderArray<Int_t> _CorrectedCaloTopoClusters_hits_index = {fReader, "_CorrectedCaloTopoClusters_hits.index"};
   TTreeReaderArray<unsigned int> _CorrectedCaloTopoClusters_hits_collectionID = {fReader, "_CorrectedCaloTopoClusters_hits.collectionID"};
   TTreeReaderArray<Int_t> _CorrectedCaloTopoClusters_particleIDs_index = {fReader, "_CorrectedCaloTopoClusters_particleIDs.index"};
   TTreeReaderArray<unsigned int> _CorrectedCaloTopoClusters_particleIDs_collectionID = {fReader, "_CorrectedCaloTopoClusters_particleIDs.collectionID"};
   TTreeReaderArray<float> _CorrectedCaloTopoClusters_shapeParameters = {fReader, "_CorrectedCaloTopoClusters_shapeParameters"};
   TTreeReaderArray<float> _CorrectedCaloTopoClusters_subdetectorEnergies = {fReader, "_CorrectedCaloTopoClusters_subdetectorEnergies"};
   TTreeReaderArray<ULong_t> ECalBarrelCells_cellID = {fReader, "ECalBarrelCells.cellID"};
   TTreeReaderArray<Float_t> ECalBarrelCells_energy = {fReader, "ECalBarrelCells.energy"};
   TTreeReaderArray<Float_t> ECalBarrelCells_energyError = {fReader, "ECalBarrelCells.energyError"};
   TTreeReaderArray<Float_t> ECalBarrelCells_time = {fReader, "ECalBarrelCells.time"};
   TTreeReaderArray<Float_t> ECalBarrelCells_position_x = {fReader, "ECalBarrelCells.position.x"};
   TTreeReaderArray<Float_t> ECalBarrelCells_position_y = {fReader, "ECalBarrelCells.position.y"};
   TTreeReaderArray<Float_t> ECalBarrelCells_position_z = {fReader, "ECalBarrelCells.position.z"};
   TTreeReaderArray<ULong_t> ECalBarrelCells2_cellID = {fReader, "ECalBarrelCells2.cellID"};
   TTreeReaderArray<Float_t> ECalBarrelCells2_energy = {fReader, "ECalBarrelCells2.energy"};
   TTreeReaderArray<Float_t> ECalBarrelCells2_energyError = {fReader, "ECalBarrelCells2.energyError"};
   TTreeReaderArray<Float_t> ECalBarrelCells2_time = {fReader, "ECalBarrelCells2.time"};
   TTreeReaderArray<Float_t> ECalBarrelCells2_position_x = {fReader, "ECalBarrelCells2.position.x"};
   TTreeReaderArray<Float_t> ECalBarrelCells2_position_y = {fReader, "ECalBarrelCells2.position.y"};
   TTreeReaderArray<Float_t> ECalBarrelCells2_position_z = {fReader, "ECalBarrelCells2.position.z"};
   TTreeReaderArray<ULong_t> ECalBarrelCellsMerged_cellID = {fReader, "ECalBarrelCellsMerged.cellID"};
   TTreeReaderArray<Float_t> ECalBarrelCellsMerged_energy = {fReader, "ECalBarrelCellsMerged.energy"};
   TTreeReaderArray<Float_t> ECalBarrelCellsMerged_position_x = {fReader, "ECalBarrelCellsMerged.position.x"};
   TTreeReaderArray<Float_t> ECalBarrelCellsMerged_position_y = {fReader, "ECalBarrelCellsMerged.position.y"};
   TTreeReaderArray<Float_t> ECalBarrelCellsMerged_position_z = {fReader, "ECalBarrelCellsMerged.position.z"};
   TTreeReaderArray<Int_t> _ECalBarrelCellsMerged_contributions_index = {fReader, "_ECalBarrelCellsMerged_contributions.index"};
   TTreeReaderArray<unsigned int> _ECalBarrelCellsMerged_contributions_collectionID = {fReader, "_ECalBarrelCellsMerged_contributions.collectionID"};
   TTreeReaderArray<ULong_t> ECalBarrelPositionedCells_cellID = {fReader, "ECalBarrelPositionedCells.cellID"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_energy = {fReader, "ECalBarrelPositionedCells.energy"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_energyError = {fReader, "ECalBarrelPositionedCells.energyError"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_time = {fReader, "ECalBarrelPositionedCells.time"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_position_x = {fReader, "ECalBarrelPositionedCells.position.x"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_position_y = {fReader, "ECalBarrelPositionedCells.position.y"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_position_z = {fReader, "ECalBarrelPositionedCells.position.z"};
   TTreeReaderArray<ULong_t> ECalBarrelPositionedCells2_cellID = {fReader, "ECalBarrelPositionedCells2.cellID"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells2_energy = {fReader, "ECalBarrelPositionedCells2.energy"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells2_energyError = {fReader, "ECalBarrelPositionedCells2.energyError"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells2_time = {fReader, "ECalBarrelPositionedCells2.time"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells2_position_x = {fReader, "ECalBarrelPositionedCells2.position.x"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells2_position_y = {fReader, "ECalBarrelPositionedCells2.position.y"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells2_position_z = {fReader, "ECalBarrelPositionedCells2.position.z"};
   TTreeReaderArray<ULong_t> ECalEndcapCells_cellID = {fReader, "ECalEndcapCells.cellID"};
   TTreeReaderArray<Float_t> ECalEndcapCells_energy = {fReader, "ECalEndcapCells.energy"};
   TTreeReaderArray<Float_t> ECalEndcapCells_energyError = {fReader, "ECalEndcapCells.energyError"};
   TTreeReaderArray<Float_t> ECalEndcapCells_time = {fReader, "ECalEndcapCells.time"};
   TTreeReaderArray<Float_t> ECalEndcapCells_position_x = {fReader, "ECalEndcapCells.position.x"};
   TTreeReaderArray<Float_t> ECalEndcapCells_position_y = {fReader, "ECalEndcapCells.position.y"};
   TTreeReaderArray<Float_t> ECalEndcapCells_position_z = {fReader, "ECalEndcapCells.position.z"};
   TTreeReaderArray<ULong_t> ECalEndcapHits_cellID = {fReader, "ECalEndcapHits.cellID"};
   TTreeReaderArray<Float_t> ECalEndcapHits_energy = {fReader, "ECalEndcapHits.energy"};
   TTreeReaderArray<Float_t> ECalEndcapHits_position_x = {fReader, "ECalEndcapHits.position.x"};
   TTreeReaderArray<Float_t> ECalEndcapHits_position_y = {fReader, "ECalEndcapHits.position.y"};
   TTreeReaderArray<Float_t> ECalEndcapHits_position_z = {fReader, "ECalEndcapHits.position.z"};
   TTreeReaderArray<Int_t> _ECalEndcapHits_contributions_index = {fReader, "_ECalEndcapHits_contributions.index"};
   TTreeReaderArray<unsigned int> _ECalEndcapHits_contributions_collectionID = {fReader, "_ECalEndcapHits_contributions.collectionID"};
   TTreeReaderArray<Int_t> genParticles_PDG = {fReader, "genParticles.PDG"};
   TTreeReaderArray<Int_t> genParticles_generatorStatus = {fReader, "genParticles.generatorStatus"};
   TTreeReaderArray<Int_t> genParticles_simulatorStatus = {fReader, "genParticles.simulatorStatus"};
   TTreeReaderArray<Float_t> genParticles_charge = {fReader, "genParticles.charge"};
   TTreeReaderArray<Float_t> genParticles_time = {fReader, "genParticles.time"};
   TTreeReaderArray<Double_t> genParticles_mass = {fReader, "genParticles.mass"};
   TTreeReaderArray<Double_t> genParticles_vertex_x = {fReader, "genParticles.vertex.x"};
   TTreeReaderArray<Double_t> genParticles_vertex_y = {fReader, "genParticles.vertex.y"};
   TTreeReaderArray<Double_t> genParticles_vertex_z = {fReader, "genParticles.vertex.z"};
   TTreeReaderArray<Double_t> genParticles_endpoint_x = {fReader, "genParticles.endpoint.x"};
   TTreeReaderArray<Double_t> genParticles_endpoint_y = {fReader, "genParticles.endpoint.y"};
   TTreeReaderArray<Double_t> genParticles_endpoint_z = {fReader, "genParticles.endpoint.z"};
   TTreeReaderArray<Float_t> genParticles_momentum_x = {fReader, "genParticles.momentum.x"};
   TTreeReaderArray<Float_t> genParticles_momentum_y = {fReader, "genParticles.momentum.y"};
   TTreeReaderArray<Float_t> genParticles_momentum_z = {fReader, "genParticles.momentum.z"};
   TTreeReaderArray<Float_t> genParticles_momentumAtEndpoint_x = {fReader, "genParticles.momentumAtEndpoint.x"};
   TTreeReaderArray<Float_t> genParticles_momentumAtEndpoint_y = {fReader, "genParticles.momentumAtEndpoint.y"};
   TTreeReaderArray<Float_t> genParticles_momentumAtEndpoint_z = {fReader, "genParticles.momentumAtEndpoint.z"};
   TTreeReaderArray<Float_t> genParticles_spin_x = {fReader, "genParticles.spin.x"};
   TTreeReaderArray<Float_t> genParticles_spin_y = {fReader, "genParticles.spin.y"};
   TTreeReaderArray<Float_t> genParticles_spin_z = {fReader, "genParticles.spin.z"};
   TTreeReaderArray<Int_t> genParticles_colorFlow_a = {fReader, "genParticles.colorFlow.a"};
   TTreeReaderArray<Int_t> genParticles_colorFlow_b = {fReader, "genParticles.colorFlow.b"};
   TTreeReaderArray<Int_t> _genParticles_parents_index = {fReader, "_genParticles_parents.index"};
   TTreeReaderArray<unsigned int> _genParticles_parents_collectionID = {fReader, "_genParticles_parents.collectionID"};
   TTreeReaderArray<Int_t> _genParticles_daughters_index = {fReader, "_genParticles_daughters.index"};
   TTreeReaderArray<unsigned int> _genParticles_daughters_collectionID = {fReader, "_genParticles_daughters.collectionID"};
   TTreeReaderArray<ULong_t> PositionedCaloClusterCells_cellID = {fReader, "PositionedCaloClusterCells.cellID"};
   TTreeReaderArray<Float_t> PositionedCaloClusterCells_energy = {fReader, "PositionedCaloClusterCells.energy"};
   TTreeReaderArray<Float_t> PositionedCaloClusterCells_energyError = {fReader, "PositionedCaloClusterCells.energyError"};
   TTreeReaderArray<Float_t> PositionedCaloClusterCells_time = {fReader, "PositionedCaloClusterCells.time"};
   TTreeReaderArray<Float_t> PositionedCaloClusterCells_position_x = {fReader, "PositionedCaloClusterCells.position.x"};
   TTreeReaderArray<Float_t> PositionedCaloClusterCells_position_y = {fReader, "PositionedCaloClusterCells.position.y"};
   TTreeReaderArray<Float_t> PositionedCaloClusterCells_position_z = {fReader, "PositionedCaloClusterCells.position.z"};
   TTreeReaderArray<ULong_t> PositionedCaloTopoClusterCells_cellID = {fReader, "PositionedCaloTopoClusterCells.cellID"};
   TTreeReaderArray<Float_t> PositionedCaloTopoClusterCells_energy = {fReader, "PositionedCaloTopoClusterCells.energy"};
   TTreeReaderArray<Float_t> PositionedCaloTopoClusterCells_energyError = {fReader, "PositionedCaloTopoClusterCells.energyError"};
   TTreeReaderArray<Float_t> PositionedCaloTopoClusterCells_time = {fReader, "PositionedCaloTopoClusterCells.time"};
   TTreeReaderArray<Float_t> PositionedCaloTopoClusterCells_position_x = {fReader, "PositionedCaloTopoClusterCells.position.x"};
   TTreeReaderArray<Float_t> PositionedCaloTopoClusterCells_position_y = {fReader, "PositionedCaloTopoClusterCells.position.y"};
   TTreeReaderArray<Float_t> PositionedCaloTopoClusterCells_position_z = {fReader, "PositionedCaloTopoClusterCells.position.z"};

   float E_cluster;  // in MeV
   float ClusterEnergyThreshold = 70.0;  // min cluster energy (in GeV)
   //float ClusterEnergyThreshold = 2.0;  // min cluster energy (in GeV)
   const int numLayers = 12;

   // E_cluster, length is the number of clusters
   std::vector<float> vec_E_cluster;

   // each has 12 sub-vectors
   std::vector<std::vector<float>> vec_E_layer        = { {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {} };
   std::vector<std::vector<float>> vec_E_frac_layer   = { {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {} };
   std::vector<std::vector<float>> vec_w_theta_layer  = { {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {} };
   std::vector<std::vector<float>> vec_w_module_layer = { {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {} };
   std::vector<std::vector<int>> vec_theta_diff_layer = { {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {} };
   std::vector<std::vector<float>> vec_E_ratio_layer  = { {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {} };
   std::vector<std::vector<float>> vec_delta_E_layer  = { {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {} };

   // longitudinal energy and widths profile, length is 12
   // to plot E, E_frac, w_theta, w_module profile in TGraph
   std::vector<float> mean_E_layer              = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   std::vector<float> mean_E_layer_error        = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   std::vector<float> mean_E_frac_layer         = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   std::vector<float> mean_E_frac_layer_error   = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   std::vector<float> mean_w_theta_layer        = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   std::vector<float> mean_w_theta_layer_error  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   std::vector<float> mean_w_module_layer       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   std::vector<float> mean_w_module_layer_error = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

   // fill TTree
   float _E_cluster;
   float _E_frac_layer0;
   float _E_frac_layer1;
   float _E_frac_layer2;
   float _E_frac_layer3;
   float _E_frac_layer4;
   float _E_frac_layer5;
   float _E_frac_layer6;
   float _E_frac_layer7;
   float _E_frac_layer8;
   float _E_frac_layer9;
   float _E_frac_layer10;
   float _E_frac_layer11;
   float _w_theta_layer0;
   float _w_theta_layer1;
   float _w_theta_layer2;
   float _w_theta_layer3;
   float _w_theta_layer4;
   float _w_theta_layer5;
   float _w_theta_layer6;
   float _w_theta_layer7;
   float _w_theta_layer8;
   float _w_theta_layer9;
   float _w_theta_layer10;
   float _w_theta_layer11;
   float _w_module_layer0;
   float _w_module_layer1;
   float _w_module_layer2;
   float _w_module_layer3;
   float _w_module_layer4;
   float _w_module_layer5;
   float _w_module_layer6;
   float _w_module_layer7;
   float _w_module_layer8;
   float _w_module_layer9;
   float _w_module_layer10;
   float _w_module_layer11;
   int _theta_diff_layer0;
   int _theta_diff_layer1;
   int _theta_diff_layer2;
   int _theta_diff_layer3;
   int _theta_diff_layer4;
   int _theta_diff_layer5;
   int _theta_diff_layer6;
   int _theta_diff_layer7;
   int _theta_diff_layer8;
   int _theta_diff_layer9;
   int _theta_diff_layer10;
   int _theta_diff_layer11;
   float _E_ratio_layer0;
   float _E_ratio_layer1;
   float _E_ratio_layer2;
   float _E_ratio_layer3;
   float _E_ratio_layer4;
   float _E_ratio_layer5;
   float _E_ratio_layer6;
   float _E_ratio_layer7;
   float _E_ratio_layer8;
   float _E_ratio_layer9;
   float _E_ratio_layer10;
   float _E_ratio_layer11;
   float _delta_E_layer0;
   float _delta_E_layer1;
   float _delta_E_layer2;
   float _delta_E_layer3;
   float _delta_E_layer4;
   float _delta_E_layer5;
   float _delta_E_layer6;
   float _delta_E_layer7;
   float _delta_E_layer8;
   float _delta_E_layer9;
   float _delta_E_layer10;
   float _delta_E_layer11;

   // vec A is theta_cell, vec B is E_cell
   // for a given theta, sum up the E_cell over modules
   // then sort the theta vec, update the E vec simultaneously
   std::pair<std::vector<int>, std::vector<float>> MergeSumAndSort(std::vector<int>& A, std::vector<float>& B) {
     std::unordered_map<int, float> elementSum;
     // traverse vec A, merge the same elements (theta ID) in vec A and sum up the corresponding elements (E_cell) in vec B
     for (size_t i = 0; i < A.size(); i ++) {
       elementSum[A[i]] += B[i];
     }
     std::vector<int> A_new;
     std::vector<float> B_new;
     // re-build vec A and vec B from elementSum
     for (const auto& entry : elementSum) {
       A_new.push_back(entry.first);
       B_new.push_back(entry.second);
     }
     std::vector<int> vec_1(A_new.size());
     std::vector<float> vec_2(B_new.size());
     std::vector<size_t> indices(A_new.size());
     for (size_t i = 0; i < indices.size(); ++ i) {
       indices[i] = i;
     }
     // sort the theta vec, update the E vec simultaneously
     std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
       return A_new[a] < A_new[b];
     });
     for (size_t i = 0; i < indices.size(); ++ i) {
       vec_1[i] = A_new[indices[i]];
       vec_2[i] = B_new[indices[i]];
     }
     return std::make_pair(vec_1, vec_2);
   }

//   // too slow !
//   // get mean of a vector
//   float get_Mean(std::vector<float> vec_data) {
//     return std::accumulate(vec_data.begin(), vec_data.end(), 0) / vec_data.size();
//   }
//   // get RMS of a vector
//   float get_RMS(std::vector<float> vec_data) {
//     return std::sqrt(std::inner_product(vec_data.begin(), vec_data.end(), vec_data.begin(), 0.0) / vec_data.size());
//   }
//   // get std_dev of a vector
//   float get_StdDev(std::vector<float> vec_data) {
//     float sumSquaredDifferences = 0.0;
//     for (const auto& value : vec_data) {
//       float difference = value - std::accumulate(vec_data.begin(), vec_data.end(), 0) / vec_data.size();
//       sumSquaredDifferences += difference * difference;
//     }
//     return std::sqrt(sumSquaredDifferences / vec_data.size());
//   }

   // extract layer number from cellID
   ULong_t Layer(ULong_t cellID) {
     const ULong_t mask = (1<<8) -1;
     return (cellID >> 11) & mask;
   }
   // extract module number from cellID
   ULong_t Module(ULong_t cellID) {
     const ULong_t mask = (1<<11) -1;
     return (cellID >> 19) & mask;
   }
   // extract theta bin from cellID
   ULong_t ThetaBin(ULong_t cellID) {
     const ULong_t mask = (1<<10) -1;
     return (cellID >> 30) & mask;
   }

   events(TTree * /*tree*/ =0) { }
   virtual ~events() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   ClassDef(events,0);
};

#endif

// draw E-theta TGraph
int count_cluster = 0;
std::vector<std::vector<int>> theta_bin_clu(100, {0});
std::vector<std::vector<float>> E_theta_bin_clu(100, {0});

#ifdef events_cxx
void events::Init(TTree *tree)
{
   fReader.SetTree(tree);
}

Bool_t events::Notify()
{
   return kTRUE;
}

#endif // #ifdef events_cxx
