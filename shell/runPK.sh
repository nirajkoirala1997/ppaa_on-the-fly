#!/bin/bash

# take snapshot of time for the uniqueness in filenames further.
timestamp=$(date +"%Y_%m_%d_%H_%M_%S")

# Begin with this .sh file location
home_path=$(dirname "$0")
cd $home_path

# home path assighed for future reference
home_path=$(pwd)

# input file should also be secured for the tee_the_data.f code it needs after completion of the executable
cd $home_path
cp ../run.machine.dat ../trash/broken/input_${timestamp}_PK.dat
cp ../output_files.dat ../trash/broken/input2_${timestamp}_PK.dat


# Compile the file to update summary data
cd ../summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}_PK.o"

#Everything is ready now start the Executable
cd $home_path
cd ../PK_integral_standalone
make
./runPK | tee "../trash/broken/output_${timestamp}.PK"

# Results are now ready we can combine the data of tee to the output files.
cd $home_path
cd ../summary/compare
./tee_the_data_${timestamp}_PK.o "output_${timestamp}.PK" 'PK' "input_${timestamp}_PK.dat" "input2_${timestamp}_PK.dat"

# This executable and input data is no more required
rm -f "tee_the_data_${timestamp}_PK.o" 
cd $home_path
cd ../trash/broken
rm -f input_${timestamp}_PK.dat
rm -f input2_${timestamp}_PK.dat
