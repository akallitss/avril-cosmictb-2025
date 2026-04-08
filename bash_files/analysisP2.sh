#!/bin/sh
pathnameped="/data/cosmic_data/2025/W43/test_P2/cosmic_data_raw/"
pathname="/data/cosmic_data/2025/W43/test_P2/cosmic_data_raw/"

declare -a aped=("CosTb_p2_m400v_d600v_cosmic_longrun_pedthr_260205_09H11") 
declare -a acos=("CosTb_p2_m400v_d600v_cosmic_longrun_datrun_260205_09H11")


nfiles=16 #last file num
#make pedestals/RMS

# get length of the ped array, assumed that the cos array is the same length
length=${#aped[@]}

# use for loop to read all values and indexes
# Feu file number in _000_0[56].fdf
for (( i=1; i<${length}+1; i++ ));
do

    #-------- ped for tracking ---------#
    cd /local/home/usernsw/avril/oct2025_code_reference_ctb/m3_tracking
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

    #-------- ped for P2 ---------#
    cd /local/home/usernsw/avril/oct2025_code_reference_ctb/CodeDamienNew
    echo "Pedestal files: $pathname${aped[$i-1]}"_000_03.fdf
    ./DreamDataReader "$pathname${aped[$i-1]}"_000_03.fdf "$pathname${aped[$i-1]}"_000_04.fdf "$pathname${aped[$i-1]}"_000_07.fdf
    root -l -q 'Pedestal.C++("output.root")'
    root -l -q 'CommonNoisePedSubstr.C++("output.root")'
    root -l -q 'RMSPedestaux.C++("output.root")'
    cp -v RMSPed.dat /data/cosmic_data/2026/W2/decoded_data/workdir_tmp/RMSPed.dat

    #------ loop on file -----#
    for num in $(seq -f "%03g" 0 $nfiles)
    do
        #-------- decode tracking --------#
        cd /local/home/usernsw/avril/oct2025_code_reference_ctb/m3_tracking
        #make cosmic config file
        #clean previous cos
        rm config_cos.json end.txt
        file_num=$((10#$num))
        echo "Cosmics data files: $pathname${acos}"_"$file_num"_01.fdf
        echo    "\"data_file_first\" : "$file_num"," > end.txt  
        echo    "\"data_file_last\" : "$file_num"," >> end.txt  
        echo	"\"data_file_basename\" : \"$pathname${acos}_\"}"  >> end.txt
        cat config_ref_2022.json end.txt > config_cos.json
        ./DataReader config_cos.json analyse
        ./tracking config_cos.json rays /data/cosmic_data/2026/W2/decoded_data/workdir_tmp/run_rays.root
        #-------- decode P2 ------#
        cd /local/home/usernsw/avril/oct2025_code_reference_ctb/CodeDamienNew
        echo "Set $num en cours d'analyse"
        echo "-> $pathname${acos[$i-1]}"_"$num"_03.fdf
        if [ -f "$pathname${acos[$i-1]}"_"$num"_03.fdf ];
        then
            echo "file here"
            ./DreamDataReader "$pathname${acos[$i-1]}"_"$num"_03.fdf "$pathname${acos[$i-1]}"_"$num"_04.fdf "$pathname${acos[$i-1]}"_"$num"_07.fdf
            root -l -q 'CommonNoisePedSubstr.C++("output.root")'
            mv -v output.root /data/cosmic_data/2026/W2/decoded_data/workdir_tmp/output.root  
        fi
        echo "file not here"

        cd /local/home/usernsw/avril/oct2025_code_reference_ctb
        root -l -q 'TBanalysisP2_2025_mod.C++("/data/cosmic_data/2026/W2/decoded_data/workdir_tmp/")'
        rm -f outtree_temp_"$num".root
        mv -v outtree.root outtree_temp_"$num".root
        rm -f /data/cosmic_data/2026/W2/decoded_data/workdir_tmp/output.root
        rm -f /data/cosmic_data/2026/W2/decoded_data/workdir_tmp/run_rays.root
    done
    rm -f outtree_final.root
    hadd outtree_final.root outtree_temp_*.root
    rm -f /data/cosmic_data/2026/W2/decoded_data/workdir_tmp/RMSPed.dat
    #rm -f outtree_temp_*.root
done
