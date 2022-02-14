# EMFIT on 2D Images in GENESIS 1.4

## Installation: 
*See https://www.r-ccs.riken.jp/labs/cbrt/installation/ for installation requirements*

Clone the NMMD repository to `/path/to/NMMD`. The repository contains a modified version of GENESIS 1.4.0 to perform NMMD and a NMA program called ElNemo.
To build GENESIS :
```
cd /path/to/NMMD
./configure
make install
```
To build ElNemo :
```
cd /path/to/NMMD/ElNemo
make
```

## Usage NMMD:

In the INP file : 
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
nm_number = 10
nm_mass = 10.0
nm_limit = 1000.0
elnemo_cutoff = 8.0 
elnemo_rtb_block = 10
elnemo_path = /path/to/NMMD/ElNemo
nm_prefix = /prefix/for/NMMD
```
- `nm_number` : The number of normal modes to use after skipping the 6 first modes (example : 10 will use modes from 7-16), default 10
- `nm_mass` : Mass value of Normal modes for NMMD, default : 10.0
- `nm_limit` : Threshold of normal mode amplitude above/below which the normal modes are updated, default : 1000.0
- `elnemo_cutoff` : Cutoff distance for elastic network model (Angstrom) default : 8.0
- `elnemo_rtb_block` : Number of residue per RTB block in the NMA computation : 10
- `elnemo_path` : Path to ElNemo installation, default : None
- `nm_prefix` : Path to a temporary directory to store NMMD files during the simulation, default : None

 

## Limitations:
- NMMD is available only for LANGEVIN temperature control in the NVT enemble
- SHAKE/RATTLE algorithm have to be turn off

## Usage Emfit images:

```

[EXPERIMENTS]
emfit = YES                            # YES/NO
emfit_type = IMAGE                     # VOLUME/IMAGE 
emfit_target = path/to/image/file.spi  # 2D SPIDER file
emfit_sigma = 2.0                      # sigma of 2D gaussian
emfit_tolerance = 0.01                 # Gaussian truncation threshold
emfit_period = 1                       # not used
emfit_roll_angle = 0.0                 # Euler roll angle in degrees
emfit_tilt_angle = 0.0                 # Euler tilt angle in degrees
emfit_yaw_angle = 0.0                  # Euler yaw angle in degrees
emfit_shift_x = 0.0                    # Shift in x direction in pixels
emfit_shift_y = 0.0                    # Shift in y direction in pixels
emfit_pixel_size = 1.0                 # Size of a pixel in Angstrom
```
