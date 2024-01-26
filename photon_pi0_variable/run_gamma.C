{
  TChain ch("events");
  ch.Add("10_GeV/output_evts_20000_pdg_22_10_GeV_ThetaMinMax_90_90_PhiMinMax_0_0.root");
  ch.Add("10_GeV/output_evts_30000_pdg_22_10_GeV_ThetaMinMax_90_90_PhiMinMax_0_0.root");
  ch.Add("10_GeV/output_evts_40000_pdg_22_10_GeV_ThetaMinMax_90_90_PhiMinMax_0_0.root");
  ch.Process("events.C");
}
