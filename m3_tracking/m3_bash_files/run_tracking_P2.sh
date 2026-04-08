#!/bin/sh
pathnameped="/data/cosmic_data/2025/W43/test_Alex/"
pathname="/data/cosmic_data/2025/W43/test_Alex/"

declare -a aped=("CosTb_det_480_800_pedthr_260114_09H58") 
declare -a acos=("CosTb_det_480_800_datrun_260114_09H58")

first_file_num=0
last_file_num=1

#make pedestals/RMS
echo "Pedestal files: $pathnameped${aped}"_000_01.fdf  
#make ending of config file
echo    "\"data_file_first\" : 0," > end.txt
echo    "\"data_file_last\" : 0," >> end.txt
echo	"\"data_file_basename\" : \"$pathnameped${aped}_\"}"  >> end.txt
cat config_ref_2022.json end.txt > config.json
  
#clean previous ped
rm root_files/test_signal.root root_files/test_RMSPed.dat root_files/test_Ped.dat
./DataReader config.json read
./DataReader config.json ped

#make cosmic config file
#clean previous cos
rm config_cos.json end.txt
num=$(printf "%03d" "$first_file_num")
echo "Cosmics data files: $pathname${acos}"_"$num"_01.fdf
echo    "\"data_file_first\" : "$first_file_num"," > end.txt  
echo    "\"data_file_last\" : "$last_file_num"," >> end.txt  
echo	"\"data_file_basename\" : \"$pathname${acos}_\"}"  >> end.txt
cat config_ref_2022.json end.txt > config_cos.json
./DataReader config_cos.json analyse
./tracking config_cos.json rays output_"$num".root
