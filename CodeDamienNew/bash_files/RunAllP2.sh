#!/bin/sh
pathnameped="/data/cosmic_data/2025/W43/test_Alex/"
pathname="/data/cosmic_data/2025/W43/test_Alex/"

declare -a aped=("CosTb_det_480_800_pedthr_260114_09H58") 
declare -a acos=("CosTb_det_480_800_datrun_260114_09H58")

nfiles=1
#make pedestals/RMS

# get length of the ped array, assumed that the cos array is the same length
length=${#aped[@]}

# use for loop to read all values and indexes
# Feu file number in _000_0[56].fdf
for (( i=1; i<${length}+1; i++ ));
do
    echo "Pedestal files: $pathname${aped[$i-1]}"_000_03.fdf
    ./DreamDataReader "$pathname${aped[$i-1]}"_000_03.fdf #"$pathname${aped[$i-1]}"_000_04.fdf "$pathname${aped[$i-1]}"_000_07.fdf
    root -l -q 'Pedestal.C++("output.root")'
    root -l -q 'CommonNoisePedSubstr.C++("output.root")'
    root -l -q 'RMSPedestaux.C++("output.root")'
    cp -v RMSPed.dat RMSPed_"${acos[$i-1]}".dat

    for num in $(seq -f "%03g" 0 $nfiles)
    do
        echo "Set $num en cours d'analyse"
        echo "-> $pathname${acos[$i-1]}"_"$num"_03.fdf
        if [ -f "$pathname${acos[$i-1]}"_"$num"_03.fdf ];
        then
            echo "file here"
            ./DreamDataReader "$pathname${acos[$i-1]}"_"$num"_03.fdf #"$pathname${acos[$i-1]}"_"$num"_04.fdf "$pathname${acos[$i-1]}"_"$num"_07.fdf
            root -l -q 'CommonNoisePedSubstr.C++("output.root")'
            mv -v output.root/data/cosmic_data/2026/W2/decoded_data/workdir_tmp/output_"${acos[$i-1]}"_"$num".root  
        fi
    
        echo "file not here"
    done
    hadd/data/cosmic_data/2026/W2/decoded_data/workdir_tmp/output_"${acos[$i-1]}"_all.root/data/cosmic_data/2026/W2/decoded_data/workdir_tmp/output_"${acos[$i-1]}"_*.root
done
