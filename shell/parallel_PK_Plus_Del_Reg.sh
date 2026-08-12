#!/bin/bash

home_path=$(dirname "$0")
cd $home_path
home_path=$(pwd)
cd ../

# Define the commands to run in parallel, and add a unique identifier to each command
command1="(./shell/runPK_Plus.sh)"
command2="(./shell/runPK_Regular.sh)"
command3="(./shell/runPK_Delta.sh)"

# Run the commands in parallel using GNU Parallel and capture the output
parallel --line-buffer ::: "$command1" "$command2" "$command3" 
#parallel --line-buffer ::: "$command1" "$command3" 
