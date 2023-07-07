#!/bin/sh

# ----------------------------------------------------
# 1D tests
# ----------------------------------------------------

# Plain

cp run_1d_mars/UA/data/log00000002.dat srcData/log00000002.1d.dat.mars
#cp run_1d_mars/UA/data/run_information.txt srcData/run_information.1d.txt
#cp run_1d_mars/UA/data/v15.1d.ps srcData
#cp run_1d_mars/UA/data/v25.1d.ps srcData

# ----------------------------------------------------
# 3D tests
# ----------------------------------------------------

cp run_3d_mars/UA/data/log00000002.dat srcData/log00000002.3d.dat.mars
#cp run_3d_mars/UA/data/run_information.txt srcData/run_information.3d.txt
#cp run_3d_mars/UA/data/3d.ps srcData

#rm srcData/*~
#cvs commit srcData/*.1d.* srcData/*.eclipse.* srcData/*3d.*
