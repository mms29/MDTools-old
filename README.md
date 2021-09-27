# NMMD in GENESIS 1.4

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
/path/to/NMMD/bin/atdyn INP pathToNMMD boolNMMD numberModes dtNMFactor
```
- `INP` : path to INP file (see how to define the INP file in https://www.r-ccs.riken.jp/labs/cbrt/usage/)
- `pathToNMMD` : path to NMMD ending with a `/` character (example : `/path/to/NMMD/`)
- `numberModes` : The number of normal modes to use after skipping the 6 first modes (example : 4 will use modes from 7-10)
- `massModes` : The normal mode mass (example : 1.0)

## Limitations:
- NMMD is available only for the VVER integrator
- NMMD is available only for LANGEVIN temperature control
- SHAKE algorithm have to be turn off
