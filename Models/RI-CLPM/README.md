
# Experiments: CLPM vs RI-CLPM in lavaan vs brms

In this folder I'm experimenting with CLPM vs RI-CLPM in lavaan vs brms.

Wrote a function to simulate a 3-wave RI-CLPM: `simulate-RICLPM-data.R`  

- [x] Fit CLPM and RICLPM with lavaan
- [x] Simulate DGP for CLPM: `CLPM-lavaan-simulation.Rmd` 
- [x] Simulate DGP for RICLPM: `RICLPM-lavaan-simulation.Rmd`
- [] Recover CLPM using brms/multilevel model instead: `CLPM-brms.Rmd` 
- [] Recover RICLPM using brms/multilevel model instead: `RICLPM-brms-simulation.Rmd`

Rendered documents are in `html/` 

Current conclusions:
- RI-CLPM is a decent baseline model
  - or any model that can deal with time-invariant confounding
- Easy to fit with lavaan
- Difficult to fit with brms... ?
- Assume no measurement error and stationarity?
- Could probaly fit multilevel version of it in Stan?
- Could add measurement error model?
  - Would be easy to simulate ME

/ Kristoffer