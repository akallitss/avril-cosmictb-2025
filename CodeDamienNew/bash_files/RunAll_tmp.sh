#!/bin/sh

make DreamDataReader
filename="/media/compassextdisk/banc_cosmiques/data/2015W15/RawCos_Fwd460_MGoutbench_150410_17H42"
#[2456]
for num in $(seq -f "%03g" 0 54)
do
echo "Set $num en cours d'analyse"
echo "-> $filename"_"$num"_0[24].fdf
./DreamDataReader "$filename"_"$num"_0[24].fdf
root -l -q 'CommonNoisePedSubstrZS.C++("output.root",4.5)'
mv -v signal.root signal_"$num".root
rm output.root
done


