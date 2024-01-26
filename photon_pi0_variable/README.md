# photon_pi0_variable
To calculate variables regarding shower shape of photon and pi0 clustering

## usage
```
root output_evts_20000_pdg_22_5_100_GeV_ThetaMinMax_80_100_PhiMinMax_0_0.5.root
```
Then process with
```
events->Process("events.C")
```
```events.h``` need to be regenerated if the contents of TTree is changed

```
events->MakeSelector()
```

## Train BDT
```
root runBDT.C
```

## Plot the ROC curve of BDT
```
root get_ROC_curve.C
```
Then do
```
root plot_roc.C
```
