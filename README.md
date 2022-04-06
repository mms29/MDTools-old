# NMMD in GENESIS 1.7.1

## Installation: 
*See https://www.r-ccs.riken.jp/labs/cbrt/installation/ for installation requirements*

Clone the NMMD repository to `/path/to/NMMD`. The repository contains a modified version of GENESIS 1.7.1 to perform NMMD.
To build GENESIS :
```
cd /path/to/NMMD
./configure
make install
```

## Usage:
```
/path/to/NMMD/bin/atdyn INP
```

The INP file is used as in GENESIS (see how to define the INP file in https://www.r-ccs.riken.jp/labs/cbrt/usage/). To use NMMD, you must chose the `integrator=NMMD` in the `DYNAMICS` section and add a section `NMMD` to the INP file as follow :

```
...
[DYNAMICS]
integrator = NMMD
...
[NMMD]
nm_number = 5
nm_mass = 10.0
nm_file = /NMfile/for/NMMD
nm_init = 0.0 0.0 0.0 0.0 0.0
nm_dt = 0.001
```
- `nm_number` : The number of normal modes to use after skipping the 6 first modes (example : 10 will use modes from 7-16), default 10
- `nm_mass` : Mass value of Normal modes for NMMD, default : 10.0
- `nm_file` : Normal mode file, text file with normal mode vectors. The vectors contained in the file will be used (if you want to skip some modes, e.g. first 6 modes, you must remove them from the file manually). The syntax of the file for a system with 2 normal modes and N atoms is as follow :
```
 VECTOR    1
 -------------------------------------
 1.460600e-02   -1.566100e-02   1.267500e-02
 1.129600e-02   -1.015400e-02   1.141000e-02
 5.195700e-03   -6.018800e-03   1.088200e-02
...
-1.047400e-02    6.799600e-03   1.267000e-02
-4.738500e-03    2.462300e-03   1.760800e-02
-4.374800e-03   -1.732600e-05   2.199900e-02
 VECTOR    2
 -----------------------------------
 -2.198000e-02   3.097700e-02   -1.289900e-02
 -1.971500e-02   2.571300e-02   -9.875500e-03
 -1.413200e-02   2.378900e-02   -7.625900e-03
...
 -3.139400e-02   4.586600e-02   -4.489600e-03
 -3.235900e-02   5.275000e-02   -5.051600e-03
 -3.891800e-02   6.174600e-02   -6.280400e-03
```
- `nm_dt` : Normal mode integration time step in picoseconds. I most cases, it should be the same as the MD time step, default : 0.001
- `nm_init` : (Optional) Initial values for normal mode amplitudes, must match the number of modes, default : None

## Limitations:
- NMMD is available only for LANGEVIN temperature control in the NVT enemble
- SHAKE/RATTLE algorithm have to be turn off
