# EMFIT on 2D Images in GENESIS 1.4

## Installation: 
*See https://www.r-ccs.riken.jp/labs/cbrt/installation/ for installation requirements*


## Usage:

In the INP file : 
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
