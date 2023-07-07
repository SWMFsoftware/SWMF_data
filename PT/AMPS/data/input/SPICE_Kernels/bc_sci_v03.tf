KPL/FK

Frame (FK) SPICE kernel file for BepiColombo-specific Science frames
=============================================================================

   These frames are sorted in two groups: those that are BEPICOLOMBO mission 
   specific and those that are Mercury generic. The first group contains the 
   frames defined by and for the BEPICOLOMBO mission, while the second 
   provides the frames that are commonly accepted by the scientific
   community for Mercury.
   
   The IAU body-fixed rotational frames for Mercury and Venus is an exception 
   to this grouping, as they are provided in a separate PCK kernel file. 


Version and Date
-------------------------------------------------------------------------------

  Version 0.3 -- February 12, 2018 -- Marc Costa Sitja, ESAC/ESA

    General update following the indications of the Hermean Environment
    Working Group (HEWG).

  Version 0.2 -- June 20, 2017 -- Marc Costa Sitja, ESAC/ESA

    Added Mercury Solar Magnetospheric frame and Mercury Orbital generic 
    frame. Updated some definitions and file outline.

  Version 0.1 -- July 04, 2016 -- Marc Costa Sitja, ESAC/ESA

     Updated BEPICOLOMBO MPO IDs from -69 to -121. 

  Version 0.0 -- February 26, 2015 -- Boris Semenov, NAIF

     Initial version.


References
-------------------------------------------------------------------------------

   1. ``Frames Required Reading''
 
   2. ``Kernel Pool Required Reading''
 
   3.  Russell, C.T., "Geophysical Coordinate Transformations,"
       Cosmic Electrodynamics, 2, 184-196, 1971
 
   4. pck00009.tpc, based on IAU 2006 constants
 
   5. pck00010.tpc, based on IAU 2009 constants
 
   6. MERCURY_MME/999999 frame definition provided by Dr. Jonathan
      McAuliffe, ESA, 20 Feb 2015

   7. ``MESSENGER Dynamic Frame Definitions Kernel'', msgr_dyn_v500.tf,
      extracted from the MESS-E/V/H-SPICE-6-V1.0 Archive Dataset,
      Scott Turner, July 2008.

   8. Russell, C.T., "Geophysical Coordinate Transformations,"
      Cosmic Electrodynamics, 2, 184-196, 1971

   9. Planetary Fact Sheets - Mercury Fact Sheet,
      David R. Williams, NASA Goddard Space Flight Center
      https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
      Accessed on 12th February 2018


Contact Information
-----------------------------------------------------------------------------

   If you have any questions regarding this file contact the ESA SPICE
   Service at ESAC:

           Marc Costa Sitja
           (+34) 91-8131-457
           mcosta@sciops.esa.int, esa_spice@sciops.esa.int


Implementation Notes
-----------------------------------------------------------------------------

   This file is used by the SPICE system as follows: programs that make use
   of this frame kernel must "load" the kernel normally during program
   initialization. Loading the kernel associates the data items with
   their names in a data structure called the "kernel pool".  The SPICELIB
   routine FURNSH loads a kernel into the pool as shown below:
 
     FORTRAN: (SPICELIB)
 
       CALL FURNSH ( frame_kernel_name )
 
     C: (CSPICE)
 
       furnsh_c ( frame_kernel_name );
 
     IDL: (ICY)
 
       cspice_furnsh, frame_kernel_name
 
     MATLAB: (MICE)
       
          cspice_furnsh ( 'frame_kernel_name' )
 
     PYTHON: (SPICEYPY)*
 
          furnsh( frame_kernel_name )
 
   In order for a program or routine to extract data from the pool, the
   SPICELIB routines GDPOOL, GIPOOL, and GCPOOL are used.  See [2] for
   more details.
 
   This file was created and may be updated with a text editor or word
   processor.
 
   * SPICEPY is a non-official, community developed Python wrapper for the
     NAIF SPICE toolkit. Its development is managed on Github.
     It is available at: https://github.com/AndrewAnnex/SpiceyPy


BepiColombo Science Frames names and NAIF ID Codes
-------------------------------------------------------------------------------

  The following BepiColombo science frames are defined in this kernel file:

      SPICE Frame Name             Long-name
      --------------------------   ------------------------------------------

  Mercury based frames for the in-orbit phase(**):

      BC_MSO                       Mercury-Centric Solar Orbital
      BC_MSO_AB                    Mercury-Centric Solar Orbital Aberrated
      BC_MSM                       Mercury-Centric Solar Magnetospheric
      BC_MBF                       Mercury Body Fixed
      BC_MME_IAU2006_OF_DATE       Mean Mercury Equator and IAU vector of date
                                   using IAU 2006 Mercury rotation constants 
      BC_MME_IAU2006_J2000         BC_MME_IAU2006_OF_DATE frame frozen at 
                                   J2000 TDB
      BC_MME_IAU2009_OF_DATE       Mean Mercury Equator and IAU vector of date
                                   using IAU 2009 Mercury rotation constants 
      BC_MME_IAU2009_J2000         BC_MME_IAU2009_OF_DATE frame frozen at 
                                   J2000 TDB

  Sun based frames for the interplanetary cruise phase:

      BC_MPO_RTN                   MPO Radial-Tangential-Normal Heliocentric
      BC_MMO_RTN                   MMO Radial-Tangential-Normal Heliocentric

  Earth based frames for the interplanetary cruise flybys:

      BC_GSE                       Geocentric Solar Ecliptic
      BC_GSM                       Geocentric Solar Magnetospheric

  Venus based frames for the interplanetary cruise flybys:

      BC_VSO                       Venus Solar Orbital


   (**) These frames are commonly used by other missions for data analysis
        and scientific research. In the future NAIF may include some of them
        in their official generic frames kernel for the Jupiter system.
        When this happens the frames will be removed from this kernel.


   These frames have the following centers, frame class and NAIF
   IDs:
   
      SPICE Frame Name              Center      Class    NAIF ID
      --------------------------    ----------  -------  ----------
      BC_MSO                        MERCURY     DYNAMIC     -121901
      BC_MSO_AB                     MERCURY     FIXED       -121902
      BC_MSM                        MERCURY     DYNAMIC     -121903
      BC_MBF                        MERCURY     FIXED       -121904
      BC_MME_IAU2006_OF_DATE        MERCURY     DYNAMIC     -121911
      BC_MME_IAU2006_J2000          MERCURY     FIXED       -121912
      BC_MME_IAU2009_OF_DATE        MERCURY     DYNAMIC     -121913
      BC_MME_IAU2009_J2000          MERCURY     FIXED       -121914
      BC_MPO_RTN                    SUN         DYNAMIC     -121921
      BC_MMO_RTN                    SUN         DYNAMIC     -121922
      BC_GSE                        EARTH       DYNAMIC     -121931
      BC_GSM                        EARTH       DYNAMIC     -121932
      BC_VSO                        VENUS       DYNAMIC     -121941

   The keywords implementing these frame definitions are located in the 
   section "Frame Definitions."


General Notes About This File
-------------------------------------------------------------------------------

   About Required Data:
   --------------------
   Most of the dynamic frames defined in this file require at least one
   of the following kernels to be loaded prior to their evaluation, 
   normally during program initialization:

     - Planetary ephemeris data (SPK), i.e. de403, de405, etc;
     - Planetary constants data (PCK);
     - Earth generic frames definitions (FK).

   Note that loading different kernels will lead to different
   orientations of the same frame at a given epoch, providing different
   results from each other, in terms of state vectors referred to these 
   frames.
 
 
   About Implementation:
   ---------------------
   The SPICE frames defined within this file and their corresponding 
   references in literature might not be equivalent, both due to 
   variations in the SPICE kernels on which the SPICE frame depends, 
   and due to possible differences in both the frame's definition and 
   implementation (e.g. GSE can be defined using the instantaneous 
   orbital plane or mean ecliptic; the mean ecliptic is a function of 
   the ecliptic model). Please refer to each applicable frame 
   description section for particular details on the current SPICE
   kernel implementation.   


Frame Definitions
-------------------------------------------------------------------------------

   This section contains the definitions of the Mercury Science frames.

Mercury Based Frames
-------------------------------------------------------------------------------

   These dynamic frames are used for analyzing data in a reference
   frame tied to the dynamics of Mercury.


BepiColombo Mercury-centric Solar Orbital frame (BC_MSO)
------------------------------------------------------------------------

   Definition:
   -----------

   The Mercury-centric solar orbital frame is defined as follows:

      -  +X axis is the position of the Sun relative to
         Mercury; it's the primary vector and points
         from Mercury to Sun;

      -  +Y axis is the component of the inertially referenced
         velocity of Sun relative to Mercury orthogonal
         to the +X axis;

      -  +Z axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         Mercury.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------

   This frame is defined as a two-vector frame using two different
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position'
   vector. Therefore, the ephemeris data required to compute the
   Mercury-Sun vector in J2000 reference frame have
   to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target velocity'
   vector. Therefore, the ephemeris data required to compute the
   Mercury-Sun velocity vector in the J2000 reference frame
   have to be loaded before using this frame.


   Remarks:
   --------

   This frame is defined based on SPK data: different planetary
   ephemerides for Mercury, Sun and the Sun Barycenter
   will lead to a different frame orientation at a given time.

   It is strongly recommended to indicate what data have been used
   in the evaluation of this frame when referring to it, e.g.
   MERCURY_SUN_ORB using de405 ephemerides.


  \begindata

      FRAME_BC_MSO                  =  -121901
      FRAME_-121901_NAME            = 'BC_MSO'
      FRAME_-121901_CLASS           =  5
      FRAME_-121901_CLASS_ID        =  -121901
      FRAME_-121901_CENTER          =  199
      FRAME_-121901_RELATIVE        = 'J2000'
      FRAME_-121901_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_-121901_FAMILY          = 'TWO-VECTOR'
      FRAME_-121901_PRI_AXIS        = 'X'
      FRAME_-121901_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_-121901_PRI_OBSERVER    = 'MERCURY'
      FRAME_-121901_PRI_TARGET      = 'SUN'
      FRAME_-121901_PRI_ABCORR      = 'NONE'
      FRAME_-121901_SEC_AXIS        = 'Y'
      FRAME_-121901_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
      FRAME_-121901_SEC_OBSERVER    = 'MERCURY'
      FRAME_-121901_SEC_TARGET      = 'SUN'
      FRAME_-121901_SEC_ABCORR      = 'NONE'
      FRAME_-121901_SEC_FRAME       = 'J2000'

  \begintext


BepiColombo Mercury-centric Solar Orbital Aberrated frame (BC_MSO_AB)
---------------------------------------------------------------------

   Definition:
   -----------

   The Mercury-centric solar orbital aberrated frame is defined is defined,
   based on the definition in [8], section 3.3.7) for the Geocentric Solar Wind
   frame, as follows:

      -  +Z is perpendicular to the orbital plane of Mercury and is
         positive towards North;

      -  +X axis is the component of the solar wind direction vector
         that is orthogonal to the +Z-axis. This axis lies in the
         orbital plane of Mercury and is positive in the direction
         opposite to the solar wind;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the center of mass of Mercury.


   Uses and applications:
   ----------------------

   The Jupiter Solar Wind frame is useful in analyzing the impact
   of solar wind on hemispheric events (from [8], section 3.3.7).


   Remarks:
   --------

   A critical issue to consider in the definition of this
   frame is aberration - the bow shock of Mercury is rotated in
   the plane of the orbit of Mercury by angle Vplanet/Vsolar_wind,
   where Vplanet is the velocity of Mercury in its orbit and
   Vsolar_wind is the solar wind speed. The sense of rotation is
   such the bow shock lags its un-rotated location on the anti-sunward
   side of Mercury. This is the +Y direction for Mercury-centric solar
   orbital frame (BC_MSO). Strictly speaking Vplanet should
   be the component perpendicular to the solar wind.

   Th solar wind streams off of the Sun in all directions at speeds
   of about 400 km/s. Nevertheless, the solar wind is not uniform.
   Although it is always directed away from the Sun, it changes speed
   and carries with it magnetic clouds, interacting regions where high
   speed wind catches up with slow speed wind, and composition
   variations. The solar wind speed can range from high (800 km/s) over
   coronal holes to low velocities (300 km/s) over streamers
   (see [13]).

   Mercury's orbital velocity ranges from 58.98 km/s to 38.86 km/s with an
   average of 47.36 km/s (see [9]).

   For Mercury under normal solar wind conditions (v ~400 km/s) the
   angle of aberration ranges from 8.4483 degrees for the maximum
   orbital velocity to 5.5663 degrees for the minimum orbital
   velocity of Mercury.

   By convention, the average velocity for Mercury in its orbit and for
   the solar wind at Mercury will be used in the definition of this
   frame, which is highly significant for determining magnetopause and
   bow shock locations. This value is determined to be 6.7838 degrees.

   Since the minimum aberration is of approximately 2.7831 degrees
   (maximum solar wind and minimum Mercury orbital velocities) and
   the maximum aberration is of approximately 11.2643 degrees (minimum
   solar wind and maximum Mercury orbital velocities), the error in
   the definition of the +X-axis direction is between -4.0007 to
   4.4805 degrees.

   This frame is defined relative to BC_MSO, which is dynamic. This aspect
   shall be taken into account if/when using this frame to define other
   dynamic frames. For further details, please refer to [1].

  \begindata

      FRAME_BC_MSO_AB            = -121902
      FRAME_-121902_NAME         = 'BC_MSO_AB'
      FRAME_-121902_CLASS        =  4
      FRAME_-121902_CLASS_ID     = -121902
      FRAME_-121902_CENTER       =  599
      TKFRAME_-121902_SPEC       = 'ANGLES'
      TKFRAME_-121902_RELATIVE   = 'BC_MSO'
      TKFRAME_-121902_ANGLES     = (   -6.7838,  0.0,  0.0 )
      TKFRAME_-121902_AXES       = (         3,    2,    3 )
      TKFRAME_-121902_UNITS      = 'DEGREES'

  \begintext


BepiColombo Mercury-centric Solar Magnetospheric (BC_MSM)
------------------------------------------------------------------------

   Definition:
   -----------

   The Mercury-centric solar magnetospheric frame is defined as follows:

      - +Z axis is the projection of the Mercury magnetic dipole axis
        (positive north) on to the plane perpendicular to the +X axis.

      -  +X axis is the position of the Sun relative to Mercury;
         it's the primary vector and points from Mercury to Sun;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         Mercury.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------

   This frame is defined as a two-vector frame using two different
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position'
   vector. Therefore, the ephemeris data required to compute the
   Mercury-Sun vector in J2000 reference frame have
   to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target velocity'
   vector. Therefore, the ephemeris data required to compute the
   Mercury-Sun velocity vector in the J2000 reference frame
   have to be loaded before using this frame.


   Remarks:
   --------

   This frame is defined based on SPK data: different planetary
   ephemerides for Mercury, Sun and the Sun Barycenter
   will lead to a different frame orientation at a given time.

   It is strongly recommended to indicate what data have been used
   in the evaluation of this frame when referring to it, e.g.
   MERCURY_SUN_ORB using de405 ephemerides.


  \begindata

      FRAME_BC_MSM                  =  -121903
      FRAME_-121903_NAME            = 'BC_MSM'
      FRAME_-121903_CLASS           =  5
      FRAME_-121903_CLASS_ID        =  -121903
      FRAME_-121903_CENTER          =  199
      FRAME_-121903_RELATIVE        = 'J2000'
      FRAME_-121903_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_-121903_FAMILY          = 'TWO-VECTOR'
      FRAME_-121903_PRI_AXIS        = 'X'
      FRAME_-121903_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_-121903_PRI_OBSERVER    = 'MERCURY'
      FRAME_-121903_PRI_TARGET      = 'SUN'
      FRAME_-121903_PRI_ABCORR      = 'NONE'
      FRAME_-121903_SEC_AXIS        = 'Z'
      FRAME_-121903_SEC_VECTOR_DEF  = 'CONSTANT'
      FRAME_-121903_SEC_SPEC        = 'LATITUDINAL'
      FRAME_-121903_SEC_UNITS       = 'DEGREES'
      FRAME_-121903_SEC_LONGITUDE   = 0.0
      FRAME_-121903_SEC_LATITUDE    = 90.0
      FRAME_-121903_SEC_FRAME       = 'IAU_MERCURY'

  \begintext


BepiColombo Mercury Body Fixed (BC_MBF)
------------------------------------------------------------------------

   Definition:
   -----------

   The Mercury Body Fixed frame is equivalent to the predefined Mercury Body
   Fixed frame (IAU_MERCURY) which is defined as a PCK (body-fixed) frame.
   PCK frames are reference frames whose orientation with respect to inertial
   frames is supplied through either binary or text PCK files.

   \begindata

       FRAME_BC_MBF             = -121904
       FRAME_-121904_NAME       = 'BC_MBF'
       FRAME_-121904_CLASS      = 4
       FRAME_-121904_CLASS_ID   = -121904
       FRAME_-121904_CENTER     = 199

       TKFRAME_-121904_SPEC     = 'MATRIX'
       TKFRAME_-121904_RELATIVE = 'IAU_MERCURY'
       TKFRAME_-121904_MATRIX   = ( 1, 0, 0,
                                   0, 1, 0,
                                   0, 0, 1 )

   \begintext


BepiColombo Mercury Mean Equator of Date Frame based on IAU2006 Constants
-------------------------------------------------------------------------

  The BC_MME_IAU2006_OF_DATE frame is based on Mean Mercury Equator
  and IAU vector of date computed using IAU 2006 Mercury rotation
  constants

  The BC_MME_IAU2006_OF_DATE frame is implemented as an Euler frame
  mathematically identical to the PCK frame IAU_MERCURY based on IAU
  2006 Mercury rotation constants but without prime meridian rotation
  terms.

  The IAU 2006 PCK data from [4] defining the IAU_MERCURY frame are:

     BODY199_POLE_RA          = (  281.01     -0.033      0. )
     BODY199_POLE_DEC         = (   61.45     -0.005      0. )
     BODY199_PM               = (  329.548     6.1385025  0. )

  Here pole RA/Dec terms in the PCK are in degrees and degreeS/Century;
  the rates here have been converted to degrees/sec ( = (rate deg/cen)/
  (86400.0 * 36525.0)). Prime meridia terms from the PCK are disregarded.

  The 3x3 transformation matrix M defined by the angles is

     M = [    0.0]   [angle_2]   [angle_3]
                  3           1           3

  Vectors are mapped from the J2000 base frame to the
  BC_MME_IAU2006_OF_DATE frame via left multiplication by M.

  The relationship of these Euler angles to RA/Dec for the
  J2000-to-IAU Mercury Mean Equator and IAU vector of date
  transformation is as follows:

     angle_1 is      0.0
     angle_2 is 90 - Dec
     angle_3 is 90 + RA, mapped into the range 0 < angle_3 < 2*pi

  Since when we define the BC_MME_IAU2006_OF_DATE frame we're defining
  the *inverse* of the above transformation, the angles for our Euler
  frame definition are reversed and the signs negated:

     angle_1 is -90 - RA, mapped into the range 0 < angle_3 < 2*pi
     angle_2 is -90 + Dec
     angle_3 is       0.0

  Then our frame definition is:

  \begindata

     FRAME_BC_MME_IAU2006_OF_DATE =  -121911
     FRAME_-121911_NAME            = 'BC_MME_IAU2006_OF_DATE'
     FRAME_-121911_CLASS           =  5
     FRAME_-121911_CLASS_ID        =  -121921
     FRAME_-121911_CENTER          =  199
     FRAME_-121911_RELATIVE        = 'J2000'
     FRAME_-121911_DEF_STYLE       = 'PARAMETERIZED'
     FRAME_-121911_FAMILY          = 'EULER'
     FRAME_-121911_EPOCH           =  @2000-JAN-1/12:00:00
     FRAME_-121911_AXES            =  ( 3  1  3 )
     FRAME_-121911_UNITS           = 'DEGREES'
     FRAME_-121911_ANGLE_1_COEFFS  = (  -11.01
                                         1.0457068978629554e-11 )
     FRAME_-121911_ANGLE_2_COEFFS  = (  -28.55
                                       -1.5844043907014476e-12 )
     FRAME_-121911_ANGLE_3_COEFFS  = (    0.0                     )
     FRAME_-121911_ROTATION_STATE  = 'INERTIAL'

  \begintext


BepiColombo Mercury Mean Equator at J2000 Frame based on IAU2006 Constants
--------------------------------------------------------------------------

  The BC_MME_IAU2006_J2000 frame is the BC_MME_IAU2006_OF_DATE frame
  frozen at J2000 TDB. For efficiency it is defined as a offset frame
  using pre-computed rotation matrix transforming vectors from
  BC_MME_IAU2006_OF_DATE to J2000 at J2000 TDB.

  For consistency with earlier BepiColombo usage the matrix below is
  from [6].

  \begindata

     FRAME_BC_MME_IAU2006_J2000    = -121912
     FRAME_-121912_NAME            = 'BC_MME_IAU2006_J2000'
     FRAME_-121912_CLASS           = 4
     FRAME_-121912_CLASS_ID        = -121912
     FRAME_-121912_CENTER          = 199
     TKFRAME_-121912_RELATIVE      = 'J2000'

     TKFRAME_-121912_SPEC          = 'MATRIX'

     TKFRAME_-121912_MATRIX        = (

        0.981593866044678,      0.190980318733265,      1.45064023353692E-15,
       -0.167757184264224,      0.862232423481673,      0.477925491080635,
        0.0912743626173337,    -0.469128730471140,      0.878400378515027

                                    )

  \begintext


BepiColombo Mercury Mean Equator of Date Frame based on IAU2009 Constants
-------------------------------------------------------------------------

  The BC_MME_IAU2009_OF_DATE frame is based on Mean Mercury Equator
  and IAU vector of date computed using IAU 2009 Mercury rotation
  constants

  The BC_MME_IAU2009_OF_DATE frame is implemented as as Euler frame
  mathematically identical to the PCK frame IAU_MERCURY based on IAU
  2009 Mercury rotation constants but without prime meridian rotation
  terms.

  The IAU 2009 PCK data from [5] defining the IAU_MERCURY frame are:

     BODY199_POLE_RA          = (  281.0097   -0.0328     0. )
     BODY199_POLE_DEC         = (   61.4143   -0.0049     0. )
     BODY199_PM               = (  329.5469    6.1385025  0. )

  Here pole RA/Dec terms in the PCK are in degrees and degreeS/Century;
  the rates here have been converted to degrees/sec ( = (rate deg/cen)/
  (86400.0 * 36525.0)). Prime meridia terms from the PCK are disregarded.

  The 3x3 transformation matrix M defined by the angles is

     M = [    0.0]   [angle_2]   [angle_3]
                  3           1           3

  Vectors are mapped from the J2000 base frame to the
  BC_MME_IAU2009_OF_DATE frame via left multiplication by M.

  The relationship of these Euler angles to RA/Dec for the
  J2000-to-IAU Mercury Mean Equator and IAU vector of date
  transformation is as follows:

     angle_1 is      0.0
     angle_2 is 90 - Dec
     angle_3 is 90 + RA, mapped into the range 0 < angle_3 < 2*pi

  Since when we define the BC_MME_IAU2009_OF_DATE frame we're defining
  the *inverse* of the above transformation, the angles for our Euler
  frame definition are reversed and the signs negated:

     angle_1 is -90 - RA, mapped into the range 0 < angle_3 < 2*pi
     angle_2 is -90 + Dec
     angle_3 is       0.0

  Then our frame definition is:

  \begindata

     FRAME_BC_MME_IAU2009_OF_DATE =  -121913
     FRAME_-121913_NAME            = 'BC_MME_IAU2009_OF_DATE'
     FRAME_-121913_CLASS           =  5
     FRAME_-121913_CLASS_ID        =  -121913
     FRAME_-121913_CENTER          =  199
     FRAME_-121913_RELATIVE        = 'J2000'
     FRAME_-121913_DEF_STYLE       = 'PARAMETERIZED'
     FRAME_-121913_FAMILY          = 'EULER'
     FRAME_-121913_EPOCH           =  @2000-JAN-1/12:00:00
     FRAME_-121913_AXES            =  ( 3  1  3 )
     FRAME_-121913_UNITS           = 'DEGREES'
     FRAME_-121913_ANGLE_1_COEFFS  = (  -11.0097
                                         1.0393692803001496e-11 )
     FRAME_-121913_ANGLE_2_COEFFS  = (  -28.5857
                                        -1.5527163028874185e-12 )
     FRAME_-121913_ANGLE_3_COEFFS  = (    0.0                     )
     FRAME_-121911_ROTATION_STATE  = 'INERTIAL'

  \begintext


Sun Based Frames
-------------------------------------------------------------------------------


BepiColombo MPO Heliocentric Radial-Tangential-Normal (BC_MPO_RTN)
------------------------------------------------------------------------

   Definition:
   -----------

   The BepiColombo MPO Heliocentric Radial-Tangential-Normal frame is defined
   as follows (from [4]):

      -  the position of MPO relative to the Sun is the
         primary vector: +X axis points from the Sun to
         MPO;

      -  the projection of the solar rotational axis perpendicular to
         the +X axis defines the frame's +Z axis;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the center of mass of the Sun.

   All vectors are geometric: no aberration corrections are used.


   Uses and applications (from [6]):
   ----------------------------------

   This frame is used to define the velocity and field direction of the
   plasma environment that the spacecraft finds itself in.


   Required Data:
   --------------

   This frame is defined as a two-vector frame using two different
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position'
   vector, therefore, the ephemeris data required to compute the
   Sun-MPO position vector in the J2000 reference frame must be
   loaded before using this frame.

   The secondary vector is defined as a constant vector in the IAU_SUN
   frame, which is a PCK-based frame, therefore a PCK file containing
   the orientation constants for the Sun must be loaded before using
   this frame.


   Remarks:
   --------

   This frame is defined based on SPK data: different planetary
   ephemerides for Mercury, the Sun, the Solar System Barycenter and
   MPO spacecraft will lead to different frame orientation at a given time.

   This frame is also defined based on the IAU_SUN frame, whose
   evaluation is based on the data included in the loaded PCK file:
   different orientation constants for the Sun's spin axis will lead to
   different frame orientation at a given time.

   It is strongly recommended to indicate what data have been used
   in the evaluation of this frame when referring to it, e.g.
   BC_RTN using the IAU 2009 constants, the DE432 ephemeris and
   the MPO ephemeris version N.


   \begindata

      FRAME_BC_MPO_RTN               = -121921
      FRAME_-121921_NAME             = 'BC_MPO_RTN'
      FRAME_-121921_CLASS            =  5
      FRAME_-121921_CLASS_ID         = -121921
      FRAME_-121921_CENTER           =  10
      FRAME_-121921_RELATIVE         = 'J2000'
      FRAME_-121921_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_-121921_FAMILY           = 'TWO-VECTOR'
      FRAME_-121921_PRI_AXIS         = 'X'
      FRAME_-121921_PRI_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_-121921_PRI_OBSERVER     = 'SUN'
      FRAME_-121921_PRI_TARGET       = 'MPO'
      FRAME_-121921_PRI_ABCORR       = 'NONE'
      FRAME_-121921_SEC_AXIS         = 'Z'
      FRAME_-121921_SEC_VECTOR_DEF   = 'CONSTANT'
      FRAME_-121921_SEC_FRAME        = 'IAU_SUN'
      FRAME_-121921_SEC_SPEC         = 'RECTANGULAR'
      FRAME_-121921_SEC_VECTOR       = ( 0, 0, 1 )

   \begintext


BepiColombo MMO Heliocentric Radial-Tangential-Normal (BC_MMO_RTN)
------------------------------------------------------------------------

   Definition:
   -----------

   The BepiColombo MMO Heliocentric Radial-Tangential-Normal frame is defined
   as follows (from [4]):

      -  the position of MMO relative to the Sun is the
         primary vector: +X axis points from the Sun to
         MPO;

      -  the projection of the solar rotational axis perpendicular to
         the +X axis defines the frame's +Z axis;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the center of mass of the Sun.

   All vectors are geometric: no aberration corrections are used.


   Uses and applications (from [6]):
   ----------------------------------

   This frame is used to define the velocity and field direction of the
   plasma environment that the spacecraft finds itself in.


   Required Data:
   --------------

   This frame is defined as a two-vector frame using two different
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position'
   vector, therefore, the ephemeris data required to compute the
   Sun-MMO position vector in the J2000 reference frame must be
   loaded before using this frame.

   The secondary vector is defined as a constant vector in the IAU_SUN
   frame, which is a PCK-based frame, therefore a PCK file containing
   the orientation constants for the Sun must be loaded before using
   this frame.


   Remarks:
   --------

   This frame is defined based on SPK data: different planetary
   ephemerides for Mercury, the Sun, the Solar System Barycenter and
   MPO spacecraft will lead to different frame orientation at a given time.

   This frame is also defined based on the IAU_SUN frame, whose
   evaluation is based on the data included in the loaded PCK file:
   different orientation constants for the Sun's spin axis will lead to
   different frame orientation at a given time.

   It is strongly recommended to indicate what data have been used
   in the evaluation of this frame when referring to it, e.g.
   BC_RTN using the IAU 2009 constants, the DE432 ephemeris and
   the MMO ephemeris version N.


   \begindata

      FRAME_BC_MMO_RTN               = -121922
      FRAME_-121922_NAME             = 'BC_MMO_RTN'
      FRAME_-121922_CLASS            =  5
      FRAME_-121922_CLASS_ID         = -121922
      FRAME_-121922_CENTER           =  10
      FRAME_-121922_RELATIVE         = 'J2000'
      FRAME_-121922_DEF_STYLE        = 'PARAMETERIZED'
      FRAME_-121922_FAMILY           = 'TWO-VECTOR'
      FRAME_-121922_PRI_AXIS         = 'X'
      FRAME_-121922_PRI_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
      FRAME_-121922_PRI_OBSERVER     = 'SUN'
      FRAME_-121922_PRI_TARGET       = 'MMO'
      FRAME_-121922_PRI_ABCORR       = 'NONE'
      FRAME_-121922_SEC_AXIS         = 'Z'
      FRAME_-121922_SEC_VECTOR_DEF   = 'CONSTANT'
      FRAME_-121922_SEC_FRAME        = 'IAU_SUN'
      FRAME_-121922_SEC_SPEC         = 'RECTANGULAR'
      FRAME_-121922_SEC_VECTOR       = ( 0, 0, 1 )

   \begintext


Earth based Frames
-------------------------------------------------------------------------------


BepiColombo Geocentric Solar Ecliptic (BC_GSE)
------------------------------------------------------------------------

   Definition:
   -----------

   The Bepicolombo Geoentric solar ecliptic frame is defined as follows:

      -  +Z axis is perpendicular to the plane of the Earth's orbit around
         the Sun (positive North);

      -  +X axis is the position of the Sun relative to
         the Earth; it's the primary vector and points
         from Earth to Sun;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         the Earth.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------

   This frame is defined as a two-vector frame using two different
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position'
   vector. Therefore, the ephemeris data required to compute the
   Earth-Sun vector in J2000 reference frame have
   to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target velocity'
   vector. Therefore, the ephemeris data required to compute the
   Earth-Sun velocity vector in the J2000 reference frame
   have to be loaded before using this frame.


   Remarks:
   --------

   This frame is defined based on SPK data: different planetary
   ephemerides for Earth, Sun and the Sun Barycenter
   will lead to a different frame orientation at a given time.

   It is strongly recommended to indicate what data have been used
   in the evaluation of this frame when referring to it, e.g.
   BC_GSE using de405 ephemerides.

   The definition of the BepiColombo Geocentric Solar Ecliptic frame is as
   follows:


   \begindata

      FRAME_BC_GSE                  = -121931
      FRAME_-121931_NAME            = 'BC_GSE'
      FRAME_-121931_CLASS           = 5
      FRAME_-121931_CLASS_ID        = -121931
      FRAME_-121931_CENTER          = 399
      FRAME_-121931_RELATIVE        = 'J2000'
      FRAME_-121931_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_-121931_FAMILY          = 'TWO-VECTOR'
      FRAME_-121931_PRI_AXIS        = 'X'
      FRAME_-121931_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_-121931_PRI_OBSERVER    = 'EARTH'
      FRAME_-121931_PRI_TARGET      = 'SUN'
      FRAME_-121931_PRI_ABCORR      = 'NONE'
      FRAME_-121931_SEC_AXIS        = '-Y'
      FRAME_-121931_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
      FRAME_-121931_SEC_OBSERVER    = 'SUN'
      FRAME_-121931_SEC_TARGET      = 'EARTH'
      FRAME_-121931_SEC_FRAME       = 'J2000'
      FRAME_-121931_SEC_ABCORR      = 'NONE

   \begintext


BepiColombo Geocentric Solar Magnetospheric (BC_GSM)
------------------------------------------------------------------------

   Definition:
   -----------

   The Bepicolombo Geoventric solar magnetospheric frame is defined as follows:

      -  +Z axis is the projection of the Earth's magnetic dipole axis
         (positive North) on to the plane perpendicular to the X
         axis;

      -  +X axis is the position of the Sun relative to
         the Earth; it's the primary vector and points
         from Earth to Sun;

      -  +Y axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         the Earth.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------

   This frame is defined as a two-vector frame using two different
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position'
   vector. Therefore, the ephemeris data required to compute the
   Earth-Sun vector in J2000 reference frame have
   to be loaded before using this frame.


   Remarks:
   --------

   The magnetic latitude of the Earth magnetic north pole is
   79.8 degrees North, the longitude is 288.2 degrees East.

   This frame is defined based on SPK data: different planetary
   ephemerides for Earth, Sun and the Sun Barycenter
   will lead to a different frame orientation at a given time.

   It is strongly recommended to indicate what data have been used
   in the evaluation of this frame when referring to it, e.g.
   BC_GSE using de405 ephemerides.

   The definition of the BepiColombo Geocentric Solar Magnetospheric
   frame is as follows:


   \begindata

      FRAME_BC_GSM                  = -121932
      FRAME_-121932_NAME            = 'BC_GSM'
      FRAME_-121932_CLASS           = 5
      FRAME_-121932_CLASS_ID        = -121932
      FRAME_-121932_CENTER          = 399
      FRAME_-121932_RELATIVE        = 'J2000'
      FRAME_-121932_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_-121932_FAMILY          = 'TWO-VECTOR'
      FRAME_-121932_PRI_AXIS        = 'X'
      FRAME_-121932_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_-121932_PRI_OBSERVER    = 'EARTH'
      FRAME_-121932_PRI_TARGET      = 'SUN'
      FRAME_-121932_PRI_ABCORR      = 'NONE'
      FRAME_-121932_SEC_AXIS        = 'Z'
      FRAME_-121932_SEC_VECTOR_DEF  = 'CONSTANT'
      FRAME_-121932_SEC_SPEC        = 'LATITUDINAL'
      FRAME_-121932_SEC_UNITS       = 'DEGREES'
      FRAME_-121932_SEC_LONGITUDE   = 288.2
      FRAME_-121932_SEC_LATITUDE    = 79.8
      FRAME_-121932_SEC_FRAME       = 'EARTH_FIXED'

   \begintext


Venus based Frames
-------------------------------------------------------------------------------


BepiColombo Venus-centric Solar Orbital frame (BC_VSO)
------------------------------------------------------------------------

   Definition:
   -----------

   The Venus-centric solar orbital frame is defined as follows:

      -  +X axis is the position of the Sun relative to
         Venus; it's the primary vector and points
         from Mercury to Sun;

      -  +Y axis is the component of the inertially referenced
         velocity of Sun relative to Venus orthogonal
         to the +X axis;

      -  +Z axis completes the right-handed system;

      -  the origin of this frame is the center of mass of
         Mercury.

   All vectors are geometric: no corrections are used.


   Required Data:
   --------------

   This frame is defined as a two-vector frame using two different
   types of specifications for the primary and secondary vectors.

   The primary vector is defined as an 'observer-target position'
   vector. Therefore, the ephemeris data required to compute the
   Venus-Sun vector in J2000 reference frame have
   to be loaded before using this frame.

   The secondary vector is defined as an 'observer-target velocity'
   vector. Therefore, the ephemeris data required to compute the
   Venus-Sun velocity vector in the J2000 reference frame
   have to be loaded before using this frame.


   Remarks:
   --------

   This frame is defined based on SPK data: different planetary
   ephemerides for Venus, Sun and the Sun Barycenter
   will lead to a different frame orientation at a given time.

   It is strongly recommended to indicate what data have been used
   in the evaluation of this frame when referring to it, e.g.
   BC_VSO using de405 ephemerides.


  \begindata

      FRAME_BC_VSO                  =  -121941
      FRAME_-121941_NAME            = 'BC_VSO'
      FRAME_-121941_CLASS           =  5
      FRAME_-121941_CLASS_ID        =  -121941
      FRAME_-121941_CENTER          =  299
      FRAME_-121941_RELATIVE        = 'J2000'
      FRAME_-121941_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_-121941_FAMILY          = 'TWO-VECTOR'
      FRAME_-121941_PRI_AXIS        = 'X'
      FRAME_-121941_PRI_VECTOR_DEF  = 'OBSERVER_TARGET_POSITION'
      FRAME_-121941_PRI_OBSERVER    = 'VENUS'
      FRAME_-121941_PRI_TARGET      = 'SUN'
      FRAME_-121941_PRI_ABCORR      = 'NONE'
      FRAME_-121941_SEC_AXIS        = 'Y'
      FRAME_-121941_SEC_VECTOR_DEF  = 'OBSERVER_TARGET_VELOCITY'
      FRAME_-121941_SEC_OBSERVER    = 'VENUS'
      FRAME_-121941_SEC_TARGET      = 'SUN'
      FRAME_-121941_SEC_ABCORR      = 'NONE'
      FRAME_-121941_SEC_FRAME       = 'J2000'

  \begintext


End of FK file.