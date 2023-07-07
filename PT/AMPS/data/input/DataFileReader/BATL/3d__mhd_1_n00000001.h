#HEADFILE
GM/IO2/3d__mhd_1_n00000001.h
10           nProc
T            save_binary
8            nByteReal

#NDIM
3            nDim

#GRIDBLOCKSIZE
       4        BlockSize1
       4        BlockSize2
       4        BlockSize3

#ROOTBLOCK
      16       nRootBlock1
       4       nRootBlock2
       4       nRootBlock3

#GRIDGEOMETRYLIMIT
cartesian                 TypeGeometry

#PERIODIC
       F       IsPeriodic1
       F       IsPeriodic2
       F       IsPeriodic3

#NSTEP
       1             nStep

#TIMESIMULATION
  5.0000000000E-02    TimeSimulation

#NCELL
650752        nCellPlot

#CELLSIZE
  5.0000000000E-01      CellSizeMin1
  5.0000000000E-01      CellSizeMin2
  5.0000000000E-01      CellSizeMin3

#PLOTRANGE
 -2.4000000000E+02         Coord1Min
  1.6000000000E+01         Coord1Max
 -3.2000000000E+01         Coord2Min
  3.2000000000E+01         Coord2Max
 -3.2000000000E+01         Coord3Min
  3.2000000000E+01         Coord3Max

#PLOTRESOLUTION
 -1.0000000000E+00       DxSavePlot1
 -1.0000000000E+00       DxSavePlot2
 -1.0000000000E+00       DxSavePlot3


#SCALARPARAM
      10            nParam
  1.66667E+00            Param1
  2.99790E+08            Param2
 -0.00000E+00            Param3
  1.60000E+01            Param4
  4.00000E+00            Param5
  4.00000E+00            Param6
  4.00000E+00            Param7
  4.00000E+00            Param8
  4.00000E+00            Param9
  1.00000E+00            Param10
  2.99790E+08            cLight
 -0.00000E+00         ThetaTild
  3.00000E+00             rBody

#PLOTVARIABLE
11                     nPlotVar
rho ux uy uz bx by bz p jx jy jz g c th p1 p2 p3 NX NY NZ R
R R R Mp/cc km/s km/s km/s nT nT nT nPa uA/m2 uA/m2 uA/m2

#OUTPUTFORMAT
ascii
