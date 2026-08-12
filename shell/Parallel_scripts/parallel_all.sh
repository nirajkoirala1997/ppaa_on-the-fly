#!/bin/bash

home_path=$(dirname "$0")
cd $home_path
home_path=$(pwd)

cd ../../
# Create a temporary directory to hold copies of the modified input files
tmp_dir=$(mktemp -d)
cp run.machine.dat output_files.dat $tmp_dir

cd $home_path
# Define the commands to run in parallel, and add a unique identifier to each command
commands=(
    "(./parallelDipole.sh)"
    "(./parallelPK.sh)"
    "(./parallelVirtual.sh)"
)

# Run the commands in parallel using GNU Parallel
parallel --line-buffer ::: \
    "(sleep 0; ${commands[0]})" \
    "(sleep 1; ${commands[1]})" \
    "(sleep 2; ${commands[2]})" 

cd ../../
# Restore input files to the original state.
cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .
