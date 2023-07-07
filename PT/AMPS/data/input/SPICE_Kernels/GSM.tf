Definition of the GSM frame

The magnetic latitude of the Earth magnetic north pole is
79.8 degrees North, the longitude is 288.2 degrees East.

This frame is defined based on SPK data: different planetary
ephemerides for Earth, Sun and the Sun Barycenter
will lead to a different frame orientation at a given time.

It is strongly recommended to indicate what data have been used
in the evaluation of this frame when referring to it, e.g.
BC_GSE using de405 ephemerides.

The kernel is found at http://spiftp.esac.esa.int/data/SPICE/BEPICOLOMBO/kernels/fk/bc_sci_v03.tf

The definition of the BepiColombo Geocentric Solar Magnetospheric
frame is as follows:

\begindata
   FRAME_GSM                       =  -236966
   FRAME_-236966_NAME           = 'GSM'
   FRAME_-236966_CLASS          =  5
   FRAME_-236966_CLASS_ID       =  -236966
   FRAME_-236966_CENTER         =  399
   FRAME_-236966_RELATIVE       = 'J2000'
   FRAME_-236966_DEF_STYLE      = 'PARAMETERIZED'
   FRAME_-236966_FAMILY         = 'TWO-VECTOR'
   FRAME_-236966_PRI_AXIS       = 'X'
   FRAME_-236966_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
   FRAME_-236966_PRI_OBSERVER   = 'EARTH'
   FRAME_-236966_PRI_TARGET     = 'SUN'
   FRAME_-236966_PRI_ABCORR     = 'NONE'
   FRAME_-236966_SEC_AXIS       = 'Z'
   FRAME_-236966_SEC_VECTOR_DEF = 'CONSTANT'
   FRAME_-236966_SEC_SPEC       = 'LATITUDINAL'
   FRAME_-236966_SEC_UNITS      = 'DEGREES'
   FRAME_-236966_SEC_LONGITUDE  =  288.2
   FRAME_-236966_SEC_LATITUDE   =   79.8
   FRAME_-236966_SEC_FRAME      = 'EARTH_FIXED'
\begintext
