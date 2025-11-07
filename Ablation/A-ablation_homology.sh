#!/bin/bash


WorkDir="/home/yagui/Busqueda_antigenos/ALL/"
RScripts="/home/yagui/AvP-Analysis"
AgDir="/home/yagui/Busqueda_antigenos/AgProtect"

homology_analysis() {
    local strain_dir=$1
    local subdir_name=$(basename "$strain_dir")

    Rscript "$RScripts"/A.1-Ablation_homology_AvP.R "$WorkDir" "$subdir_name" $evalue $bits $ide || {
        log_message  "ERROR: R script failed for $subdir_name"
        exit 1
    } 
}

declare -a evalues=("1e-3" "1e-5" "1e-10")
declare -a bitss=("30" "50" "70")
declare -a ides=("15" "25" "35")

for strain_dir in "$WorkDir"/*; do
	for evalue in "${evalues[@]}"; do
		for bits in "${bitss[@]}"; do
			for ide in "${ides[@]}"; do
        			homology_analysis "$strain_dir" >> Ablacion.txt
        		done
        	done
        done 
done



for evalue in "${evalues[@]}"; do
	for bits in "${bitss[@]}"; do
		for ide in "${ides[@]}"; do
       			Rscript "$RScripts"/A.2-Ablation_homology_AgProtect.R $AgDir $evalue $bits $ide  >> Ablacion_ag.txt
       		done
       	done
done 


