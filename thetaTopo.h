#ifndef thetaTopo_h
#define thetaTopo_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>

class thetaTopo : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
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
   ///TTreeReaderArray<Int_t> genParticles#0_index = {fReader, "genParticles#0.index"};
   ///TTreeReaderArray<Int_t> genParticles#0_collectionID = {fReader, "genParticles#0.collectionID"};
   ///TTreeReaderArray<Int_t> genParticles#1_index = {fReader, "genParticles#1.index"};
   ///TTreeReaderArray<Int_t> genParticles#1_collectionID = {fReader, "genParticles#1.collectionID"};
   TTreeReaderArray<ULong_t> ECalBarrelPositionedHits_cellID = {fReader, "ECalBarrelPositionedHits.cellID"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedHits_energy = {fReader, "ECalBarrelPositionedHits.energy"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedHits_position_x = {fReader, "ECalBarrelPositionedHits.position.x"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedHits_position_y = {fReader, "ECalBarrelPositionedHits.position.y"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedHits_position_z = {fReader, "ECalBarrelPositionedHits.position.z"};
   ///TTreeReaderArray<Int_t> ECalBarrelPositionedHits#0_index = {fReader, "ECalBarrelPositionedHits#0.index"};
   ///TTreeReaderArray<Int_t> ECalBarrelPositionedHits#0_collectionID = {fReader, "ECalBarrelPositionedHits#0.collectionID"};
   TTreeReaderArray<ULong_t> ECalEndcapHits_cellID = {fReader, "ECalEndcapHits.cellID"};
   TTreeReaderArray<Float_t> ECalEndcapHits_energy = {fReader, "ECalEndcapHits.energy"};
   TTreeReaderArray<Float_t> ECalEndcapHits_position_x = {fReader, "ECalEndcapHits.position.x"};
   TTreeReaderArray<Float_t> ECalEndcapHits_position_y = {fReader, "ECalEndcapHits.position.y"};
   TTreeReaderArray<Float_t> ECalEndcapHits_position_z = {fReader, "ECalEndcapHits.position.z"};
   ///TTreeReaderArray<Int_t> ECalEndcapHits#0_index = {fReader, "ECalEndcapHits#0.index"};
   ///TTreeReaderArray<Int_t> ECalEndcapHits#0_collectionID = {fReader, "ECalEndcapHits#0.collectionID"};
   TTreeReaderArray<ULong_t> ECalBarrelCellsStep1_cellID = {fReader, "ECalBarrelCellsStep1.cellID"};
   TTreeReaderArray<Float_t> ECalBarrelCellsStep1_energy = {fReader, "ECalBarrelCellsStep1.energy"};
   TTreeReaderArray<Float_t> ECalBarrelCellsStep1_energyError = {fReader, "ECalBarrelCellsStep1.energyError"};
   TTreeReaderArray<Float_t> ECalBarrelCellsStep1_time = {fReader, "ECalBarrelCellsStep1.time"};
   TTreeReaderArray<Float_t> ECalBarrelCellsStep1_position_x = {fReader, "ECalBarrelCellsStep1.position.x"};
   TTreeReaderArray<Float_t> ECalBarrelCellsStep1_position_y = {fReader, "ECalBarrelCellsStep1.position.y"};
   TTreeReaderArray<Float_t> ECalBarrelCellsStep1_position_z = {fReader, "ECalBarrelCellsStep1.position.z"};
   TTreeReaderArray<ULong_t> ECalBarrelCellsStep2_cellID = {fReader, "ECalBarrelCellsStep2.cellID"};
   TTreeReaderArray<Float_t> ECalBarrelCellsStep2_energy = {fReader, "ECalBarrelCellsStep2.energy"};
   TTreeReaderArray<Float_t> ECalBarrelCellsStep2_position_x = {fReader, "ECalBarrelCellsStep2.position.x"};
   TTreeReaderArray<Float_t> ECalBarrelCellsStep2_position_y = {fReader, "ECalBarrelCellsStep2.position.y"};
   TTreeReaderArray<Float_t> ECalBarrelCellsStep2_position_z = {fReader, "ECalBarrelCellsStep2.position.z"};
   ///TTreeReaderArray<Int_t> ECalBarrelCellsStep2#0_index = {fReader, "ECalBarrelCellsStep2#0.index"};
   ///TTreeReaderArray<Int_t> ECalBarrelCellsStep2#0_collectionID = {fReader, "ECalBarrelCellsStep2#0.collectionID"};
   TTreeReaderArray<ULong_t> ECalBarrelCells_cellID = {fReader, "ECalBarrelCells.cellID"};
   TTreeReaderArray<Float_t> ECalBarrelCells_energy = {fReader, "ECalBarrelCells.energy"};
   TTreeReaderArray<Float_t> ECalBarrelCells_energyError = {fReader, "ECalBarrelCells.energyError"};
   TTreeReaderArray<Float_t> ECalBarrelCells_time = {fReader, "ECalBarrelCells.time"};
   TTreeReaderArray<Float_t> ECalBarrelCells_position_x = {fReader, "ECalBarrelCells.position.x"};
   TTreeReaderArray<Float_t> ECalBarrelCells_position_y = {fReader, "ECalBarrelCells.position.y"};
   TTreeReaderArray<Float_t> ECalBarrelCells_position_z = {fReader, "ECalBarrelCells.position.z"};
   TTreeReaderArray<ULong_t> ECalBarrelPositionedCells_cellID = {fReader, "ECalBarrelPositionedCells.cellID"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_energy = {fReader, "ECalBarrelPositionedCells.energy"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_energyError = {fReader, "ECalBarrelPositionedCells.energyError"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_time = {fReader, "ECalBarrelPositionedCells.time"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_position_x = {fReader, "ECalBarrelPositionedCells.position.x"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_position_y = {fReader, "ECalBarrelPositionedCells.position.y"};
   TTreeReaderArray<Float_t> ECalBarrelPositionedCells_position_z = {fReader, "ECalBarrelPositionedCells.position.z"};
   TTreeReaderArray<ULong_t> ECalEndcapCells_cellID = {fReader, "ECalEndcapCells.cellID"};
   TTreeReaderArray<Float_t> ECalEndcapCells_energy = {fReader, "ECalEndcapCells.energy"};
   TTreeReaderArray<Float_t> ECalEndcapCells_energyError = {fReader, "ECalEndcapCells.energyError"};
   TTreeReaderArray<Float_t> ECalEndcapCells_time = {fReader, "ECalEndcapCells.time"};
   TTreeReaderArray<Float_t> ECalEndcapCells_position_x = {fReader, "ECalEndcapCells.position.x"};
   TTreeReaderArray<Float_t> ECalEndcapCells_position_y = {fReader, "ECalEndcapCells.position.y"};
   TTreeReaderArray<Float_t> ECalEndcapCells_position_z = {fReader, "ECalEndcapCells.position.z"};
   TTreeReaderArray<ULong_t> emptyCaloCells_cellID = {fReader, "emptyCaloCells.cellID"};
   TTreeReaderArray<Float_t> emptyCaloCells_energy = {fReader, "emptyCaloCells.energy"};
   TTreeReaderArray<Float_t> emptyCaloCells_energyError = {fReader, "emptyCaloCells.energyError"};
   TTreeReaderArray<Float_t> emptyCaloCells_time = {fReader, "emptyCaloCells.time"};
   TTreeReaderArray<Float_t> emptyCaloCells_position_x = {fReader, "emptyCaloCells.position.x"};
   TTreeReaderArray<Float_t> emptyCaloCells_position_y = {fReader, "emptyCaloCells.position.y"};
   TTreeReaderArray<Float_t> emptyCaloCells_position_z = {fReader, "emptyCaloCells.position.z"};
   TTreeReaderArray<Int_t> CaloTopoClusters_type = {fReader, "CaloTopoClusters.type"};
   TTreeReaderArray<Float_t> CaloTopoClusters_energy = {fReader, "CaloTopoClusters.energy"};
   TTreeReaderArray<Float_t> CaloTopoClusters_energyError = {fReader, "CaloTopoClusters.energyError"};
   TTreeReaderArray<Float_t> CaloTopoClusters_position_x = {fReader, "CaloTopoClusters.position.x"};
   TTreeReaderArray<Float_t> CaloTopoClusters_position_y = {fReader, "CaloTopoClusters.position.y"};
   TTreeReaderArray<Float_t> CaloTopoClusters_position_z = {fReader, "CaloTopoClusters.position.z"};
   TTreeReaderArray<Float_t> CaloTopoClusters_directionError_x = {fReader, "CaloTopoClusters.directionError.x"};
   TTreeReaderArray<Float_t> CaloTopoClusters_directionError_y = {fReader, "CaloTopoClusters.directionError.y"};
   TTreeReaderArray<Float_t> CaloTopoClusters_directionError_z = {fReader, "CaloTopoClusters.directionError.z"};
   ///TTreeReaderArray<Int_t> CaloTopoClusters#0_index = {fReader, "CaloTopoClusters#0.index"};
   ///TTreeReaderArray<Int_t> CaloTopoClusters#0_collectionID = {fReader, "CaloTopoClusters#0.collectionID"};
   ///TTreeReaderArray<Int_t> CaloTopoClusters#1_index = {fReader, "CaloTopoClusters#1.index"};
   ///TTreeReaderArray<Int_t> CaloTopoClusters#1_collectionID = {fReader, "CaloTopoClusters#1.collectionID"};
   ///TTreeReaderArray<Int_t> CaloTopoClusters#2_index = {fReader, "CaloTopoClusters#2.index"};
   ///TTreeReaderArray<Int_t> CaloTopoClusters#2_collectionID = {fReader, "CaloTopoClusters#2.collectionID"};
   TTreeReaderArray<float> CaloTopoClusters_0 = {fReader, "CaloTopoClusters_0"};
   TTreeReaderArray<float> CaloTopoClusters_1 = {fReader, "CaloTopoClusters_1"};
   TTreeReaderArray<ULong_t> CaloTopoClusterCells_cellID = {fReader, "CaloTopoClusterCells.cellID"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_energy = {fReader, "CaloTopoClusterCells.energy"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_energyError = {fReader, "CaloTopoClusterCells.energyError"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_time = {fReader, "CaloTopoClusterCells.time"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_position_x = {fReader, "CaloTopoClusterCells.position.x"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_position_y = {fReader, "CaloTopoClusterCells.position.y"};
   TTreeReaderArray<Float_t> CaloTopoClusterCells_position_z = {fReader, "CaloTopoClusterCells.position.z"};

   thetaTopo(TTree * /*tree*/ =0) { }
   virtual ~thetaTopo() { }
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

   ClassDef(thetaTopo,0);

////
   TString label;

   std::vector<double> theta_range;
   std::vector<double> phi_range;

   std::vector<double> resol_theta;
   std::vector<double> resol_phi;

   std::vector<double> x_axis_theta;
   std::vector<double> x_axis_phi;

   TH1 *P_theta = new TH1D("particle theta", "", 50, 0., 3.2);
   TH1 *P_phi   = new TH1D("particle phi", "", 50, -3.4, 3.4);

   TH1 *C_theta_1 = new TH1D("topo_theta_reso_1", "", 50, -0.3, 0.3);
   TH1 *C_theta_2 = new TH1D("topo_theta_reso_2", "", 50, -0.3, 0.3);
   TH1 *C_theta_3 = new TH1D("topo_theta_reso_3", "", 50, -0.3, 0.3);
   TH1 *C_theta_4 = new TH1D("topo_theta_reso_4", "", 50, -0.3, 0.3);
   TH1 *C_theta_5 = new TH1D("topo_theta_reso_5", "", 50, -0.3, 0.3);

   TH1 *C_phi_1   = new TH1D("topo_phi_reso_1", "", 50, -0.3, 0.3);
   TH1 *C_phi_2   = new TH1D("topo_phi_reso_2", "", 50, -0.3, 0.3);
   TH1 *C_phi_3   = new TH1D("topo_phi_reso_3", "", 50, -0.3, 0.3);
   TH1 *C_phi_4   = new TH1D("topo_phi_reso_4", "", 50, -0.3, 0.3);
   TH1 *C_phi_5   = new TH1D("topo_phi_reso_5", "", 50, -0.3, 0.3);

};

#endif

#ifdef thetaTopo_cxx
void thetaTopo::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t thetaTopo::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef thetaTopo_cxx
