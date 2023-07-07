#CODEVERSION
 9.20                        CodeVersion

#USERMODULE
Europa Single-Fluid MHD Pe, Xianzhe Jia, June 2013
 0.20                  VersionUserModule

#COMPONENT
GM                              NameComp

#PRECISION
8                              nByteReal

#EQUATION
MHD with electron pressure                    NameEquation
       9                            nVar

#CHECKGRIDSIZE
       8                              nI
       8                              nJ
       8                              nK
      64                     MinBlockALL

#NEWRESTART
F                         DoRestartBFace

#RESTARTINFILE
one                              TypeRestartInFile

#NSTEP
    3000                           nStep

#STARTTIME
    2000                           iYear
       1                          iMonth
       1                            iDay
       0                           iHour
       0                         iMinute
       0                         iSecond

#TIMESIMULATION
 0.000000000000000E+00       tSimulation

#GRIDGEOMETRY
spherical_lnr               TypeGeometry

#GRID
       4                     nRootBlockX
       4                     nRootBlockY
       4                     nRootBlockZ
-6.400000000000000E+01              xMin
 6.400000000000000E+01              xMax
-6.400000000000000E+01              yMin
 6.400000000000000E+01              yMax
-6.400000000000000E+01              zMin
 6.400000000000000E+01              zMax

#LIMITRADIUS
 1.000000000000000E+00         RadiusMin
 1.120000000000000E+02         RadiusMax

#COORDSYSTEM
GSE                      TypeCoordSystem

#SOLARWIND
 1.500000000000000E+01            SwNDim
 1.044000000000000E+06            SwTDim
-1.000000000000000E+02           SwUxDim
 0.000000000000000E+00           SwUyDim
 0.000000000000000E+00           SwUzDim
 0.000000000000000E+00          SwBxDdim
 0.000000000000000E+00           SwByDim
-4.250000000000000E+02           SwBzDim

#IOUNITS
PLANETARY                     TypeIoUnit

#NORMALIZATION
READ
 1.569000000000000E+06        No2SiUnitX
 1.569000000000000E+06        No2SiUnitU
 1.672600000000000E-21      No2SiUnitRho

#PLOTFILENAME
hour                     NameMaxTimeUnit

#BODY
T                                UseBody
 1.000000000000000E+00             rBody
 0.000000000000000E+00         rCurrents
 1.000000000000000E+04          BodyNDim
 6.000000000000000E+02          BodyTDim

#END

Additional info

       F time_accurate

Io2Si_V=
  1.56900E+06  1.00000E+03  1.67260E-21  1.00000E+00  1.00000E+06  1.00000E-09  1.00000E-09  1.00000E+00  1.00000E+00  1.00000E+00  1.00000E-06  1.00000E-03  1.00000E+00  6.37349E-16  5.72958E+01
No2Io_V=
  1.00000E+00  1.56900E+03  1.00000E+00  1.00000E+00  1.00000E+00  4.11754E+00  7.19323E+01  2.62431E-15  4.11754E-09  6.46042E-03  3.64830E-02  1.12862E+02  2.98221E+08  7.19323E+01  1.74533E-02
