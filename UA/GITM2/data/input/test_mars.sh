#!/bin/sh

MPI=/usr/local/bin/mpirun

# ----------------------------------------------------
# 1D tests
# ----------------------------------------------------

./Config.pl -mars -g=1,1,120,4
make

# Plain 1D

rm -rf run_1d_mars
make rundir
mv run run_1d_mars
cp srcData/UAM.mars.in.1d run_1d_mars/UAM.in
cd run_1d_mars
$MPI -np 1 ./GITM.exe
../share/Scripts/DiffNum.pl -b -r=1e-5 UA/data/log0000000?.dat UA/DataIn/log00000002.1d.dat.mars >& test_1d_mars.diff
#diff UA/data/run_information.txt UA/DataIn/run_information.txt >> ../test_1d_mars.diff
#cd UA ; $MPI -np 1 ./pGITM ; cd ..
#cd UA/data ; idl < ../DataIn/idl_input.1d ; cd ../..
#cd ..
ls -l test_1d_mars.diff
cd ..

# ----------------------------------------------------
# 3D tests
# ----------------------------------------------------

./Config.pl -mars -g=8,4,120,4
make

rm -rf run_3d_mars
make rundir
mv run run_3d_mars
cp srcData/UAM.mars.in.3d run_3d_mars/UAM.in
cd run_3d_mars
$MPI -np 4 ./GITM.exe
../share/Scripts/DiffNum.pl -b -r=1e-5 UA/data/log0000000?.dat UA/DataIn/log00000002.3d.dat.mars >& test_3d_mars.diff
#diff UA/data/run_information.txt UA/DataIn/run_information.3d.txt >> ../test_3d_mars.diff
#$cd UA ; $MPI pGITM ; cd ..
#cd UA/data ; idl < ../DataIn/idl_input.3d ; cd ../..
#cd ..
ls -l test_3d_mars.diff
cd ..

