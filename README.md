# thetaphi_script
to produce (topo)cluster responses and resolutions plots over theta and phi

"fixed_phi" need to be assigned in theta_phi_events.h before processing:
```
bool fixed_phi = true/false;
```

## usage
```
root -l output.root
events->Process("theta_phi_events.C")
```

the .h file need to be regenerated once the contents of TTree is changed:
```
events->MakeSelector()
```

"events" is the TTree name.

## TChain
to enlarge the statistics, there might be more than one root files.

no need to "hadd" the root files if they too large. one can use "TChain":
```
TChain ch("events");
ch.Add("./output_evts_5000_pdg_11_50_GeV_ThetaMinMax_90_90_PhiMinMax_0_6.28318_1.root");
ch.Add("./output_evts_5000_pdg_11_50_GeV_ThetaMinMax_90_90_PhiMinMax_0_6.28318_2.root");
ch.Process("theta_phi_events.C");
```
