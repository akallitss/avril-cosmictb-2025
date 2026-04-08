#!/bin/sh

pathnameped="/mnt/nas_clas12/DATA/CosmicBench/2023/W08/"
pathname="/mnt/nas_clas12/DATA/CosmicBench/2023/W08/"

nfiles=5

feu="05"

declare -a aped=("selfTPOTFe_URWELL_D680_M320_R0_pedthr_230419_13H58")
declare -a acos=("selfTPOTFe_URWELL_D680_M320_R0_LOWACT5_datrun_230419_17H12")



#make pedestals/RMS

# get length of the ped array, assumed that the cos array is the same length
length=${#aped[@]}

# use for loop to read all values and indexes
# Feu file number in _000_0[56].fdf
for (( i=1; i<${length}+1; i++ ));
do
    echo "Pedestal files: $pathname${aped[$i-1]}"_000_"$feu".fdf
    ./DreamDataReader "$pathname${aped[$i-1]}"_000_"$feu".fdf #"$pathname${aped[$i-1]}"_000_02.fdf "$pathname${aped[$i-1]}"_000_03.fdf
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
            ./DreamDataReader "$pathname${acos[$i-1]}"_"$num"_"$feu".fdf #"$pathname${acos[$i-1]}"_"$num"_02.fdf "$pathname${acos[$i-1]}"_"$num"_03.fdf
            root -l -q 'CommonNoisePedSubstr.C++("output.root")'
            mv -v output.root output_"${acos[$i-1]}"_"$num".root
        fi
        echo "file not here"
    done
done