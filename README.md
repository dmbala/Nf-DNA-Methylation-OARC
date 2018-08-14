
# DNA methylation analysis pipeline using nextflow. 
<img src="https://github.com/dmbala/Nf-DNA-Methylation-OARC/blob/master/Fig/dna-methyl-dag.png" width="450px" height="350px" />

## Files in the repo

 * DNA-Meth-nf-timeline.html: Pipeline Execution Summary
 * DNA-Meth-nf.html:  Pipeline MultiQC results
 * dna-methylation-analysis.nf: Nextflow pipeline (need to include more on analysis)
 * job-submit.s: slurm job description file
 * nextflow.config: Nextflow config file


## Major steps 
 * STEP 1 Quality control check with FastQC
 * STEP 2 Remove adaptors with trimgalore
 * STEP 3 Align with HISAT2
 * STEP 4 Generate BED from gtf file
 * STEP 5 Perform RSeQC analysis after alighment
 * STEP 6 Mark duplicates
 * STEP 7 Feature counts
 * STEP 8 Merge featurecounts
 * STEP 9 edgeR MDS and heatmap
 * STEP 10 Deseq2 analysis 
 * STEP 11 MultiQC

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
