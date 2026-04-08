#!/bin/sh

pathnamedata="/home/clas12/MaxenceTest/juin2023/CodeDamienNew/"
pathnamerays="/mnt/nas_clas12/DATA/CosmicBench/2022/root_files/"
pathworkdir="workdir/"

declare -a acos=(
"CosTb_TPOT_Z3_M430_D100_RD3_M470_D600_datrun_211117_17H34")

# get length of the ped array, assumed that the cos array is the same length
length=${#acos[@]}

# use for loop to read all values and indexes
# Feu file number in _000_0[56].fdf
for (( i=1; i<${length}+1; i++ ));
do
    echo "copying stuff"
    echo "Data file: $pathnamedata"/output_"${acos[$i-1]}"_000.root
    cp -v "$pathnamedata"/output_"${acos[$i-1]}"_000.root "$pathworkdir"/output.root
    cp -v "$pathnamedata"/RMSPed_"${acos[$i-1]}".dat "$pathworkdir"/RMSPed.dat
    echo "Ray file : $pathnamerays"/rays_"${acos[$i-1]}".root
    cp -v "$pathnamerays"/rays_"${acos[$i-1]}".root "$pathworkdir"/run_rays.root

    root -l -q TBanalysisRD3.C++
    mv -v output.txt output_"${acos[$i-1]}".txt
done