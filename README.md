#GENESIS 1.1.6 modified to include NMMD.

##Installation: 
```
cd /path/to/NMMD
./configure
make install
cd ElNemo
make
```

##Usage:
`bin/atdyn pathToNMMD boolNMMD numberModes dtNMFactor`

-pathToNMMD : path to NMMD ending with a "/" character (example : /path/to/NMMD/)
-boolNMMD : 1 to perform NMMD, 0 for standard MD
-numberModes : The number of normal modes to use after skipping the 6 first modes (example : 4 will use modes from 7-10)
-dtNMFactor : Factor of time step of normal mode integration compared to the standard time step NM_dt = dtNMFactor * dt (example : 1.0 NM_dt = dt)
  
*NMMD is available only for the VVER integrator *

