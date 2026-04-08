#!/bin/sh

pathnameped="/mnt/cosmic_data/P2/W49/"
pathname="/mnt/cosmic_data/P2/W49/"

nfiles=3

declare -a aped=("CosTb_cosmiclong_P2_M400_D900_pedthr_241204_13H46")
declare -a acos=("CosTb_cosmiclong_P2_M400_D900_datrun_241204_13H46")

#make pedestals/RMS

# get length of the ped array, assumed that the cos array is the same length
length=${#aped[@]}

# use for loop to read all values and indexes
# Feu file number in _000_0[56].fdf
for (( i=1; i<${length}+1; i++ ));
do
    echo "Pedestal files: $pathname${aped[$i-1]}"_000_03.fdf
    ./DreamDataReader "$pathname${aped[$i-1]}"_000_03.fdf "$pathname${aped[$i-1]}"_000_04.fdf
    root -l -q 'Pedestal.C++("output.root")'
    root -l -q 'CommonNoisePedSubstr.C++("output.root")'
    root -l -q 'RMSPedestaux.C++("output.root")'
    cp -v RMSPed.dat RMSPed_"${acos[$i-1]}".dat

    for num in $(seq -f "%03g" 0 $nfiles)
    do
        echo "Set $num en cours d'analyse"
        echo "-> $pathname${acos[$i-1]}"_"$num"_"$feu".fdf
        if [ -f "$pathname${acos[$i-1]}"_"$num"_"$feu".fdf ];
        then
            echo "file here"
            ./DreamDataReader "$pathname${acos[$i-1]}"_"$num"_03.fdf "$pathname${acos[$i-1]}"_"$num"_04.fdf
            root -l -q 'CommonNoisePedSubstr.C++("output.root")'
            mv -v output.root output_"${acos[$i-1]}"_"$num".root
        fi
        echo "file not here"
    done
    hadd  output_"${acos[$i-1]}"_all.root  output_"${acos[$i-1]}"_*.root
done