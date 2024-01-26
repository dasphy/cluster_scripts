#ifndef theta_phi_events_h
#define theta_phi_events_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "edm4hep/CalorimeterHitData.h"

#include <vector>

#include "edm4hep/ClusterData.h"

#include "podio/ObjectID.h"

#include "edm4hep/SimCalorimeterHitData.h"

#include "edm4hep/MCParticleData.h"

#include "podio/GenericParameters.h"



class theta_phi_events : public TSelector {
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
   TTreeReaderArray<ULong_t> ECalBarrelPositionedHits_cellID = {fReader, "ECalBarrelPositionedHits.cellID"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedHits_energy = {fReader, "ECalBarrelPositionedHits.energy"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedHits_position_x = {fReader, "ECalBarrelPositionedHits.position.x"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedHits_position_y = {fReader, "ECalBarrelPositionedHits.position.y"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedHits_position_z = {fReader, "ECalBarrelPositionedHits.position.z"};
   TTreeReaderArray<Int_t> _ECalBarrelPositionedHits_contributions_index = {fReader, "_ECalBarrelPositionedHits_contributions.index"};
   TTreeReaderArray<unsigned int> _ECalBarrelPositionedHits_contributions_collectionID = {fReader, "_ECalBarrelPositionedHits_contributions.collectionID"};
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

   theta_phi_events(TTree * /*tree*/ =0) { }
   virtual ~theta_phi_events() { }
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

   ClassDef(theta_phi_events,0);

   // Need to be modified before process the Tree in ROOT file
   //bool fixed_phi = true;
   bool fixed_phi = false;

   TString label;
   int nbins = 100;
   float hist_min = -0.005;
   float hist_max = 0.005;
   std::vector<double> theta_range;
   std::vector<double> phi_range;
   std::vector<double> resol_theta;
   std::vector<double> resol_phi;
   std::vector<double> x_axis_theta;
   std::vector<double> x_axis_phi;

   TH1 *P_theta = new TH1D("Particle Theta", "", nbins, 0., 3.2);
   TH1 *P_phi   = new TH1D("Particle Phi", "", nbins, -3.4, 3.4);
   TH1 *C_theta_1 = new TH1D("Theta resolution 1", "", nbins, hist_min, hist_max);
   TH1 *C_theta_2 = new TH1D("Theta resolution 2", "", nbins, hist_min, hist_max);
   TH1 *C_theta_3 = new TH1D("Theta resolution 3", "", nbins, hist_min, hist_max);
   TH1 *C_theta_4 = new TH1D("Theta resolution 4", "", nbins, hist_min, hist_max);
   TH1 *C_theta_5 = new TH1D("Theta resolution 5", "", nbins, hist_min, hist_max);
   TH1 *C_phi_1   = new TH1D("Phi resolution 1", "", nbins, hist_min, hist_max);
   TH1 *C_phi_2   = new TH1D("Phi resolution 2", "", nbins, hist_min, hist_max);
   TH1 *C_phi_3   = new TH1D("Phi resolution 3", "", nbins, hist_min, hist_max);
   TH1 *C_phi_4   = new TH1D("Phi resolution 4", "", nbins, hist_min, hist_max);
   TH1 *C_phi_5   = new TH1D("Phi resolution 5", "", nbins, hist_min, hist_max);

};

#endif

#ifdef theta_phi_events_cxx
void theta_phi_events::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t theta_phi_events::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef theta_phi_events_cxx
