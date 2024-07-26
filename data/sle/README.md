# COSMOPharm Data Sources

This document details the sources for the free-volume-related parameters (`v298`, `v_hc`) of APIs that are not covered in [our *Mol. Pharm.* paper](https://doi.org/10.1021/acs.molpharmaceut.4c00342). These specific parameters are included in the `table_params.xlsx` file. For the sources of other parameters, please refer to the above *Mol. Pharm.* paper or the publication by [Klajmon (2022)](https://doi.org/10.1021/acs.molpharmaceut.2c00573).

## Sources for `v298`

- **CBZ (Carbamazepine)**
  - Liquid densities from [molecular dynamics simulations](https://doi.org/10.1021/acs.molpharmaceut.2c00573) extrapolated to 298 K.

- **RBV (Ribavirin) and VST (Valsartan)**
  - As in case of GSF and SIM, estimated using a quantitative structure–property relationship (QSPR) approach implemented in the [Amsterdam Modeling Suite](https://www.scm.com/), version 2022.101. Should be taken with caution.
 
## Sources for `v_hc`

- **All APIs**
  - `v_hc` based on Bondi atomic radii and calculated using the fast method proposed by [Zhao *et al.*](https://doi.org/10.1021/jo034808o)

Please make sure to refer to the `table_params.xlsx` for the detailed parameters and their respective values.
