#!/bin/bash

# take snapshot of time for the uniqueness in filenames further.
timestamp=$(date +"%Y_%m_%d_%H_%M_%S")

# Begin with this .sh file location
home_path=$(dirname "$0")
cd $home_path

# home path assigned for future reference
home_path=$(pwd)

# input file should also be secured for the tee_the_data.f code it needs after completion of the executable
cd $home_path
cp ../run.machine.dat ../trash/broken/input_${timestamp}_Dipole.dat
cp ../output_files.dat ../trash/broken/input2_${timestamp}_Dipole.dat

# Compile the file to update summary data
cd ../summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}_Dipole.o"

#Everything is ready now start the Executable
cd $home_path
cd ../dipole_standalone
make
./runDipole | tee "../trash/broken/output_${timestamp}.Dipole"

# Results are now ready we can combine the data of tee to the output files.
cd $home_path
cd ../summary/compare
./tee_the_data_${timestamp}_Dipole.o "output_${timestamp}.Dipole" 'dipole' "input_${timestamp}_Dipole.dat" "input2_${timestamp}_Dipole.dat"

# This executable and input data is no more required
rm -f "tee_the_data_${timestamp}_Dipole.o" 
cd $home_path
cd ../trash/broken
rm -f input_${timestamp}_Dipole.dat
rm -f input2_${timestamp}_Dipole.dat
