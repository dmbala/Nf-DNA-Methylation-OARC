
# DNA methylation analysis pipeline using nextflow. 
<img src="https://github.com/dmbala/Nf-DNA-Methylation-OARC/blob/master/Fig/dna-methyl-dag.png" width="450px" height="350px" />

## Files in the repo

 * DNA-Meth-nf-timeline.html: Pipeline Execution Summary
 * DNA-Meth-nf.html:  Pipeline MultiQC results
 * dna-methylation-analysis.nf: Nextflow pipeline (need to include more on analysis)
 * job-submit.s: slurm job description file
 * nextflow.config: Nextflow config file


## Major steps 

 * STEP 1 - Quality control with FastQC
 * STEP 2A - Adaptor trimming and filtering (Trimgalore) 
 * STEP 2B - Adaptor trimming and filtering (bbduk) 
 * STEP 2C - Adaptor trimming and filtering (fastp) 
 * STEP 3 - Quality control with FastQC after Trim
 * STEP 4A - Build Bismark index if not provided
 * STEP 4 - align with Bismark
 * STEP 5 - Bismark deduplicate
 * STEP 6 - Bismark methylation extraction
 * STEP 7 - Bismark Sample Report
 * STEP 8 - Bismark Summary Report
 * STEP 9 - dmrcov
 * STEP 10 - dmrcsv
 * STEP 11 - MultiQC

## Storage
Produces large temporary fastq files that contain methylation info (e.g C_to_T.fastq.gz or G_to_A.fastq.gz)

## Single End vs Pairs
For single end reads
```
params.SingleEnd = "true"
```
For pair reads
```
params.SingleEnd = "false" 
```
