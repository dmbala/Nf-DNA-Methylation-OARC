#!/bin/bash
#SBATCH --job-name=m-seq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=03:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
module use /projects/community/modulefiles
module load nextflow 
export NXF_OPTS='-Xms1g -Xmx4g'
export NF_Work_Dir="/scratch/${USER}/NFWorkDir/${PWD}/work"
mkdir -p $NF_Work_Dir
srun nextflow run bb4.nf --SingleEnd="false" --reads="/projects/oarc/NF-Seq/SampleData/chrX_data/samples/*_{1,2}.fastq.gz" -w $NF_Work_Dir -with-trace -with-report DNA-Meth-nf.html  -with-timeline DNA-Meth-nf-timeline.html  -with-dag dna-methyl-dag.png -resume 
#srun nextflow run bb4.nf --SingleEnd="false" --reads="/projects/oarc/NF-Seq/SampleData/chrX_data/samples/*_{1,2}.fastq.gz" --task_cpus=$SLURM_CPUS_PER_TASK -w $NF_Work_Dir -with-trace -resume 

#pid=$!
#wait $pid
#if  [ $pid  -eq  0 ];
#then
#    echo " Pipeline completed. Removing WorkDir files"
#    rm -rf $NF_Work_Dir
#else
#    echo "Pipeline not completed." 
#fi


