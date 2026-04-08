#!/bin/sh

make DreamDataReader
filenameped="/media/data/Clas12/CosmicBench/2015/w24/Fwd4Ped_150612_18H26"

#[7810]case 
#make pedestals/RMS
echo "Pedestal files: $filenameped"_000_02.fdf  "$filenameped"_000_10.fdf 
./DreamDataReader "$filenameped"_000_02.fdf  "$filenameped"_000_10.fdf 
root -l -q 'Pedestal.C++("output.root")'
root -l -q 'RMSPedestauxNZS.C++("output.root")'
root -l -q 'CommonNoisePedSubstr.C++("output.root")'
root -l -q 'RMSPedestaux.C++("output.root")'


#filename="/media/mcubedisk1/data/2015W18/RawCos_CR6C1_440_150430_18H16"

#./DreamDataReader "$filename"_000_0[78].fdf "$filename"_000_10.fdf
#root -l -q 'CommonNoisePedSubstr.C++("output.root")'



#sudo ~/bin/FeuUdpControl -a 0x7FE -w 1024 -c ../Tcm4Fwd.cfg -f Fwd4Ped
