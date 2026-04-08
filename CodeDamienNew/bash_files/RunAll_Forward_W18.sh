#!/bin/sh

make DreamDataReader
filenameped="/media/mcubedisk1/data/2015W18/RawPed_CR6C1_Test_150430_16H17"
filename="/media/mcubedisk1/data/2015W18/RawCos_CR6C1_Test_150430_16H18"

#make pedestals/RMS
echo "Pedestal files: $filenameped"_000_0[789].fdf "$filenameped"_000_10.fdf
./DreamDataReader "$filenameped"_000_0[789].fdf "$filenameped"_000_10.fdf
root -l -q 'Pedestal.C++("output.root")'
root -l -q 'CommonNoisePedSubstr.C++("output.root")'
root -l -q 'RMSPedestaux.C++("output.root")'

#[78910]case 
#nb of data file 
for num in $(seq -f "%03g" 0 13)
do
    echo "Set $num en cours d'analyse"
    echo "-> $filename"_"$num"_0[789].fdf "$filename"_"$num"_10.fdf
    ./DreamDataReader "$filename"_"$num"_0[789].fdf "$filename"_"$num"_10.fdf
    root -l -q 'CommonNoisePedSubstrZS.C++("output.root",4.5)'
    mv -v signal.root signal_"$num".root
    rm output.root
done


