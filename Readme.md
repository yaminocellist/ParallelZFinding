The z-hit-extension trick is only applied to x/y vertex finding, but not onto the ad (Angular Distance) cut finding;

# Usage
For maximum display resolution check:
```bash
system_profiler SPDisplaysDataType | grep Resolution
```

For cpp files:
```bash
cd ./build
cmake .. && make && ./Para_dPhi <METHOD> <TARGET> <CEN_LOW> <CEN_HIGH> <Z_LOW> <Z_HIGH>
```

<METHOD> is divided into 3 parts by delimit '_', 
part one: 
    nomix
    perCen
    perZ
    mix
part two:
    w   
    wo  

Arguments definition for "dPhi_Fitter.C":
```bash
root 'dPhi_Fitter.C("OPT1_OPT2_OPT3")'
```
    OPT1: option 1, could be:
        single: process only one .root file;
        all   : process all found .root files;
        A     : analyze the multiplicity densities (v.s. different centralities)

    OPT2: option 2, could be:
        root:   using ROOT default Fit;
        Chi2:   using Chi2 fit method;
        LS  :   using Least Square fit method;
        TLS :   using Total Least Square fit method.
    OPT3: option 3, could be:
        integer number, to indicate which .root file to be processed.