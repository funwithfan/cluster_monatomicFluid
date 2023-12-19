#!/bin/bash

#programDir="/home/users/jcfan/MyPrograms/cluster_energy_SingleAtomMolecule_1D_neigh/build/cluster_sherlock"
programDir="../build/cluster_sherlock"

P_list=(65) #(240 300 350 460 540)
numFrame=1
save_every=50
resultDir="/public/home/ustc_fanjc/cluster_energy_SingleAtomMolecule_new"
simulationDir="/public/home/ustc_fanjc/cluster_energy_SingleAtomMolecule_new"
forceField="ReaxFF"

# 使用 seq 命令生成 T_list
T_list=(155)


num_P=${#P_list[@]}
num_T=${#T_list[@]}

for ((j=0; j<$num_P; j++)); do
    P=${P_list[j]}
    for ((i=0; i<$num_T; i++)); do
        T=${T_list[i]}
        job_name="${P}_${T}"  
        output_file="${P}_${T}_${forceField}_%j.out" 
        error_file="${P}_${T}_${forceField}_%j.err" 
        time_limit="1:00:00"
        partition="hfacnormal02"
        cpus="1"
        #memory="5GB"

        # 编写 SLURM 脚本文件
        echo "#!/bin/bash" > slurm_script_${j}_${i}.sh
        echo "#SBATCH --job-name=${job_name}" >> slurm_script_${j}_${i}.sh
        echo "#SBATCH --output=${output_file}" >> slurm_script_${j}_${i}.sh
        echo "#SBATCH --error=${error_file}" >> slurm_script_${j}_${i}.sh
        echo "#SBATCH --time=${time_limit}" >> slurm_script_${j}_${i}.sh
        echo "#SBATCH -p ${partition}" >> slurm_script_${j}_${i}.sh
        echo "#SBATCH -c ${cpus}" >> slurm_script_${j}_${i}.sh
        #echo "#SBATCH --mem=${memory}" >> slurm_script_${j}_${i}.sh

        echo "${programDir} ${P} ${T} ${numFrame} ${save_every} ${resultDir} ${simulationDir} ${forceField}" >> slurm_script_${j}_${i}.sh
        
        # 提交作业
        sbatch slurm_script_${j}_${i}.sh

        # 删除生成的 SLURM 脚本文件
        rm slurm_script_${j}_${i}.sh
    done
done

