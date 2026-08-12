#!/bin/bash

home_path=$(dirname "$0")
cd $home_path
home_path=$(pwd)
cd ../../

# Define the commands to run in parallel, and add a unique identifier to each command
commands=(
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
    "(./shell/runDipole.sh)"
)

# Define the new values for the "max # of distribution increment step_size from xq"
distribution_steps=(100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500)

# Function to modify the input files
modify_input_files() {
    local step_size=$1
    local index=$2

#    echo "Modifying files for index: $index with step size: $step_size"

    # Modify run.machine.dat
    sed -i "5s/.*/${step_size}/" run.machine.dat
#    echo "Updated run.machine.dat with step size: $step_size"

    # Modify output_files.dat
    sed -i "s/real[0-9]*.dat/real${index}.dat/" output_files.dat
#    echo "Updated output_files.dat for index: $index"
}

# Create a temporary directory to hold copies of the modified input files
tmp_dir=$(mktemp -d)
tmp_dir2=$(mktemp -d)

cp run.machine.dat output_files.dat $tmp_dir
cp run.machine.dat output_files.dat $tmp_dir2
#
## Export the modify_input_files function for parallel
export -f modify_input_files
#
## Run the commands in parallel using GNU Parallel
parallel --line-buffer ::: \
    "(sleep  0; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[ 0]}  0; ${commands[ 0]})" \
    "(sleep  1; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[ 1]}  1; ${commands[ 1]})" \
    "(sleep  2; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[ 2]}  2; ${commands[ 2]})" \
    "(sleep  3; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[ 3]}  3; ${commands[ 3]})" \
    "(sleep  4; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[ 4]}  4; ${commands[ 4]})" \
    "(sleep  5; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[ 5]}  5; ${commands[ 5]})" \
    "(sleep  6; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[ 6]}  6; ${commands[ 6]})" \
    "(sleep  7; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[ 7]}  7; ${commands[ 7]})" \
    "(sleep  8; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[ 8]}  8; ${commands[ 8]})" \
    "(sleep  9; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[ 9]}  9; ${commands[ 9]})" \
    "(sleep 10; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[10]} 10; ${commands[10]})" \
    "(sleep 11; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[11]} 11; ${commands[11]})" \
    "(sleep 12; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[12]} 12; ${commands[12]})" \
    "(sleep 13; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[13]} 13; ${commands[13]})" \
    "(sleep 14; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[14]} 14; ${commands[14]})" \
    "(sleep 15; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[15]} 15; ${commands[15]})" \
    "(sleep 16; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[16]} 16; ${commands[16]})" \
    "(sleep 17; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[17]} 17; ${commands[17]})" \
    "(sleep 18; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[18]} 18; ${commands[18]})" \
    "(sleep 19; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[19]} 19; ${commands[19]})" \
    "(sleep 20; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[20]} 20; ${commands[20]})" \
    "(sleep 21; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[21]} 21; ${commands[21]})" \
    "(sleep 22; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[22]} 22; ${commands[22]})" \
    "(sleep 23; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[23]} 23; ${commands[23]})" \
    "(sleep 24; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[24]} 24; ${commands[24]})" 


# Clean up the temporary directory
cp $tmp_dir2/run.machine.dat $tmp_dir2/output_files.dat .
rm -rf $tmp_dir
rm -rf $tmp_dir2


#~~~~~~~~[Combine Data generated from Parallel Computation]
cd $home_path
cd ../../

output_dir=$(sed -n '7p' run.machine.dat | awk '{print $1}')

output_summary_file="summary/${output_dir}/real_all.dat"
temp_file="temp_summary.dat"

# Ensure the summary file is empty before we start
> $output_summary_file
> $temp_file

# Define the Q values
#distribution_steps=(100 400 700 1000 1300 1600 1900 2200 2500)
distribution_steps=(100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500)

# Iterate over the distribution steps and process the corresponding output files
for i in ${!distribution_steps[@]}; do
    q_value=${distribution_steps[$i]}
    file="summary/${output_dir}/real${i}.dat"
    
    # Check if the file exists
    if [[ -f $file ]]; then
        # Read the first line and extract the required values
        first_line=$(head -n 1 $file)
        # Add the first line's data to the temporary summary file
        echo $first_line | awk -v q_value=$q_value '{printf "%10d %25.15e %25.15e\n", q_value, $2, $3}' >> $temp_file
        # Append the rest of the file's content to the output summary file
        tail -n +2 $file >> $output_summary_file
    else
        echo "Warning: File $file not found."
    fi
done

# Prepend the processed first lines to the output summary file
cat $temp_file $output_summary_file > "temp_all.dat" && mv "temp_all.dat" $output_summary_file

# Remove the temporary file
rm -f $temp_file

echo "  "
echo "  "
echo "Processing complete. Summary saved to $output_summary_file."
echo "  "
