# GENESIS 1.1.6 modified to include NMMD.

## Installation: 
*See https://www.r-ccs.riken.jp/labs/cbrt/installation/ for installation requirements*

Clone the NMMD repository to `/path/to/NMMD`, then :
```
cd /path/to/NMMD
./configure
make install
cd ElNemo
make
```

## Usage:
```
bin/atdyn pathToNMMD boolNMMD numberModes dtNMFactor
```

- `pathToNMMD` : path to NMMD ending with a `/` character (example : `/path/to/NMMD/`)
- `boolNMMD` : 1 to perform NMMD, 0 for standard MD
- `numberModes` : The number of normal modes to use after skipping the 6 first modes (example : 4 will use modes from 7-10)
- `dtNMFactor` : Factor of time step of normal mode integration compared to the standard time step `NM_dt = dtNMFactor * dt` (example : 1.0 `NM_dt = dt`)

## Limitations:
- NMMD is available only for the VVER integrator
- SHAKE algorithm have to be turn off
