
#!/bin/sh

pathnameped="/mnt/nas_clas12/DATA/CosmicBench/2024/W05/"
pathname="/mnt/nas_clas12/DATA/CosmicBench/2024/W05/"
nfiles=1

declare -a aped=("CosTb_380V_stats_pedthr_240212_11H42")
declare -a acos=("CosTb_380V_stats_datrun_240212_11H42")

#make pedestals/RMS

# get length of the ped array, assumed that the cos array is the same length
length=${#aped[@]}

# use for loop to read all values and indexes
for (( i=1; i<${length}+1; i++ ));
do

    echo "Pedestal files: $pathnameped${aped[$i-1]}"_000_01.fdf  
    #make ending of config file
    echo    "\"data_file_first\" : 0," > end.txt
    echo    "\"data_file_last\" : 0," >> end.txt
    echo	"\"data_file_basename\" : \"$pathnameped${aped[$i-1]}_\"}"  >> end.txt
    cat config_ref_2022.json end.txt > config.json
  
    #clean previous ped
    rm -v test_ped_signal.root test_RMSPed.dat test_Ped.dat
    ./DataReader config.json read
    ./DataReader config.json ped

    #make cosmic config file
    #clean previous cos
    rm -v config_cos.json end.txt
    for num in $(seq -f "%03g" 0 5 $nfiles)
    do
	((first=num))
	((last=num+4))
	echo "Cosmics data files: $pathname${acos[$i-1]}"_000_01.fdf
	echo    "\"data_file_first\" : "$first"," > end.txt  
	echo    "\"data_file_last\" : "$last"," >> end.txt  
	echo	"\"data_file_basename\" : \"$pathname${acos[$i-1]}_\"}"  >> end.txt
	cat config_ref_2022.json end.txt > config_cos.json
	#./DataReader config_cos.json analyse
	#./tracking config_cos.json rays root_files/rays_"${acos[$i-1]}"_"$num".root
	./tracking config_cos.json rays output_"$num".root
    done
    # hadd root_files/rays_"${acos[$i-1]}"_all.root root_files/rays_"${acos[$i-1]}"_0*.root
    # mv -v root_files/rays_"${acos[$i-1]}"_all.root /mnt/nas_clas12/DATA/CosmicBench/2022/root_files/
    # rm -fv root_files/rays_"${acos[$i-1]}"_0*.root
    
done
