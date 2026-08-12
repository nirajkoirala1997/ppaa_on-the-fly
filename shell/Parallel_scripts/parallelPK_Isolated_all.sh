#!/bin/bash

#home_path=$(dirname "$0")
home_path=$(dirname "$(readlink -f "$0")")
cd $home_path
home_path=$(pwd)
cd ../../
echo $home_path

# Define the commands to run in parallel, and add a unique identifier to each command
commands=(
    "(./shell/parallel_PK_Plus_Del_Reg.sh)")

# Define the new values for the "max # of distribution increment step_size from xq"
distribution_steps=(170 180 190 200 210)

# Function to modify the input files
modify_input_files() {
    local step_size=$1
    local index=$2

#    echo "Modifying files for index: $index with step size: $step_size"

    # Modify run.machine.dat
    sed -i "5s/.*/${step_size}/" run.machine.dat
#    echo "Updated run.machine.dat with step size: $step_size"

    # Modify output_files.dat
    sed -i "s/PK_Plus[0-9]*.dat/PK_Plus${index}.dat/" output_files.dat
    sed -i "s/PK_Regular[0-9]*.dat/PK_Regular${index}.dat/" output_files.dat
    sed -i "s/PK_Delta[0-9]*.dat/PK_Delta${index}.dat/" output_files.dat
#    echo "Updated output_files.dat for index: $index"
}

# Create a temporary directory to hold copies of the modified input files
tmp_dir=$(mktemp -d)
tmp_dir2=$(mktemp -d)
cp run.machine.dat output_files.dat $tmp_dir
cp run.machine.dat output_files.dat $tmp_dir2

# Export the modify_input_files function for parallel
export -f modify_input_files

# Run the commands in parallel using GNU Parallel
parallel --line-buffer ::: \
    "(sleep 0; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[0]} 0; ${commands[0]})" \
    "(sleep 1; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[1]} 1; ${commands[0]})" \
    "(sleep 2; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[2]} 2; ${commands[0]})" \
    "(sleep 3; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[3]} 3; ${commands[0]})" \
    "(sleep 4; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[4]} 4; ${commands[0]})" 
#    "(sleep 5; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[5]} 5; ${commands[0]})" \
#    "(sleep 6; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[6]} 6; ${commands[0]})" 
#    "(sleep 7; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[7]} 7; ${commands[0]})" \
#    "(sleep 8; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[8]} 8; ${commands[0]})" 
#    "(sleep 9; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[9]} 9; ${commands[0]})" \
#    "(sleep 10; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[10]} 10; ${commands[0]})" 
#    "(sleep 11; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[11]} 11; ${commands[0]})" \
#    "(sleep 12; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[12]} 12; ${commands[0]})" \
#    "(sleep 13; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files ${distribution_steps[13]} 13; ${commands[0]})" 

# Clean up the temporary directory
cp $tmp_dir2/run.machine.dat $tmp_dir2/output_files.dat .
rm -rf $tmp_dir
rm -rf $tmp_dir2

#~~~~~~~~[Combine Data generated from Parallel Computation]
cd $home_path
cd ../../

output_dir=$(sed -n '7p' run.machine.dat | awk '{print $1}')

output_summary_file1="PK_Isolated/summary/${output_dir}/PK_Plus_all.dat"
output_summary_file2="PK_Isolated/summary/${output_dir}/PK_Regular_all.dat"
output_summary_file3="PK_Isolated/summary/${output_dir}/PK_Delta_all.dat"
#output_summary_bk_file="PK_Isolated/summary/${output_dir}/PK_Plus_all_backup.dat"
temp_file1="temp_summary1.dat"
temp_file2="temp_summary2.dat"
temp_file3="temp_summary3.dat"

# Ensure the summary file is empty before we start
> $output_summary_file1
> $output_summary_file2
> $output_summary_file3
> $temp_file1
> $temp_file2
> $temp_file3

# Iterate over the distribution steps and process the corresponding output files
for i in ${!distribution_steps[@]}; do
    q_value=${distribution_steps[$i]}

    file1="PK_Isolated/summary/${output_dir}/PK_Plus${i}.dat"
    file2="PK_Isolated/summary/${output_dir}/PK_Regular${i}.dat"
    file3="PK_Isolated/summary/${output_dir}/PK_Delta${i}.dat"
    
    # Check if the file exists
    if [[ -f $file1 ]]; then
        # Read the first line and extract the required values
        first_line=$(head -n 1 $file1)
        # Add the first line's data to the temporary summary file
        echo $first_line | awk -v q_value=$q_value '{printf "%10d %25.15e %25.15e\n", q_value, $2, $3}' >> $temp_file1
        # Append the rest of the file's content to the output summary file
        tail -n +2 $file1 >> $output_summary_file1
	rm -f $file1
    else
        echo "Warning: File $file1 not found."
    fi
    # Check if the file exists
    if [[ -f $file2 ]]; then
        # Read the first line and extract the required values
        first_line=$(head -n 1 $file2)
        # Add the first line's data to the temporary summary file
        echo $first_line | awk -v q_value=$q_value '{printf "%10d %25.15e %25.15e\n", q_value, $2, $3}' >> $temp_file2
        # Append the rest of the file's content to the output summary file
        tail -n +2 $file2 >> $output_summary_file2
	rm -f $file2
    else
        echo "Warning: File $file2 not found."
    fi
    # Check if the file exists
    if [[ -f $file3 ]]; then
        # Read the first line and extract the required values
        first_line=$(head -n 1 $file3)
        # Add the first line's data to the temporary summary file
        echo $first_line | awk -v q_value=$q_value '{printf "%10d %25.15e %25.15e\n", q_value, $2, $3}' >> $temp_file3
        # Append the rest of the file's content to the output summary file
        tail -n +2 $file3 >> $output_summary_file3
	rm -f $file3
    else
        echo "Warning: File $file3 not found."
    fi

done

# Prepend the processed first lines to the output summary file
    cat $temp_file1 $output_summary_file1 > "temp_all1.dat" && mv "temp_all1.dat" "$output_summary_file1"
    cat $temp_file2 $output_summary_file2 > "temp_all2.dat" && mv "temp_all2.dat" "$output_summary_file2"
    cat $temp_file3 $output_summary_file3 > "temp_all3.dat" && mv "temp_all3.dat" "$output_summary_file3"
    echo "  "
    echo "  "
    echo "Processing complete. Summary saved to $output_summary_file1."
    echo "Processing complete. Summary saved to $output_summary_file2."
    echo "Processing complete. Summary saved to $output_summary_file3."
    echo "  "
    echo "  "


#cat $temp_file $output_summary_file > "temp_all.dat" && mv "temp_all.dat" $output_summary_file

# Remove the temporary file
rm -f $temp_file1
rm -f $temp_file2
rm -f $temp_file3

