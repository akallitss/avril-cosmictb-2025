#!/bin/sh

pathnameped="/mnt/nas_clas12/DATA/CosmicBench/2022/W02/"
pathname="/mnt/nas_clas12/DATA/CosmicBench/2022/W02/"

declare -a aped=("CosTb_TPOT_ZZZ_PS2_z1_430_z2_460_long_100fC_pedthr_220215_17H33")

declare -a acos=(
"CosTb_TPOT_ZZZ_PS2_z1_430_z2_460_long_100fC_datrun_220215_17H33")

#make pedestals/RMS

for num in $(seq -f "%03g" 0 5 44)
do
    ((first=num))
    ((last=num+4))
    echo "Cosmics data files: $pathname${acos[$i-1]}"_"$num"_01.fdf
    echo    "\"data_file_first\" : "$first","   
    echo    "\"data_file_last\" : "$last","   
done
