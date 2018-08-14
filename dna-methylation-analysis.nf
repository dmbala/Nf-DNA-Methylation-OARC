
/*
===============================================================
 DNA Methylation analysis pipeline
===============================================================
 ##  DNA-Meth Analysis Pipeline. Started July 2018.
 ##   https://github.com/dmbala
 ##   Bala Desinghu <dmbala@gmail.com>
 ##   Major steps for DNA-meth are outlined here at https://github.com/SciLifeLab and https://github.com/renyiwu/bioseq/blob/master/Methyl-seq_general_workflow.txt. More detail about nextflow can be found at https://www.nextflow.io
 ##   ---------------------------------------------------------------
 ##   This pipeline is prepared for 
 ##   the Office of Advanced Research Computing (OARC), Rutgers.  
 ##   The pipeline is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. All required packages are available in OARC cluster. 
 ##   ---------------------------------------------------------------
 #
*/

Program_Dir = "/projects/oarc/NF-Seq/"

params.reads = "/projects/oarc/NF-Seq/SampleData/chrX_data/samples/*_{1,2}.fastq.gz"
params.SingleEnd=false
params.genome = "/projects/community/genomics_references/Mus_musculus/NCBI/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
params.gtf = "/projects/community/genomics_references/Mus_musculus/NCBI/GRCm38/Annotation/Genes/genes.gtf"
params.bismark_index = "/projects/community/genomics_references/BismarkIndex/GRCm38/"
params.outdir = "${PWD}/results"
params.task_cpus = 1
params.do_trimgalore = false 
params.do_bbduk = false 
params.do_fastp = true 

Channel
    .fromFilePairs( params.reads, size: params.SingleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB" }
    .into {read_files_fastqc; read_files_trim_galore; read_files_trim_bbduk; read_files_trim_fastp; ch_print }
ch_print.subscribe {println "read_pairs: $it"}

/*
 * STEP 1 - Quality control with FastQC
 */

process fastqc {
    tag "FastQC"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    $Program_Dir/FastQC/fastqc -q $reads
    """
}

/*
 * STEP 2A - Adaptor trimming and filtering (Trimgalore) 
 */


if(params.do_trimgalore){
    process trim_galore {
        tag "Trimgalore"
        module "python/2.7.12"
        afterScript = 'rm -rf *.clumped*.fastq.gz'
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else null
            }

        input:
        set val(name), file(reads) from read_files_trim_galore

        output:
        set val(name), file ("*trimmed.fq.gz") into trimmed_reads, check_trim_reads, trimmed_reads_fastqc
        file "*trimming_report.txt" into trim_galore_results
        file "*_fastqc.{zip,html}" into trim_galore_fastqc_reports

        script:
        if (params.SingleEnd) {
            """
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in=$reads out=${name}.clumped.fastq.gz markduplicates subs=0
            $Program_Dir/TrimGalore/trim_galore \\
                --fastqc -q 20 --gzip ${name}.clumped.fastq.gz
            rename _val_1.fq.gz _R1_trimmed.fq.gz *
            """
        } else {
            """
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in1=${reads[0]} in2=${reads[1]}   \\
                out1=${name}.clumped_R1.fastq.gz out2=${name}.clumped_R2.fastq.gz \\
                markduplicates subs=0
            $Program_Dir/TrimGalore/trim_galore \\
                --paired -q 20 --fastqc \\
                --gzip ${name}.clumped_R1.fastq.gz ${name}.clumped_R2.fastq.gz  
            rename _val_1.fq.gz _R1_trimmed.fq.gz *
            rename _val_2.fq.gz _R2_trimmed.fq.gz *
            """
        }
      
    }
}

/*
 * STEP 2B - Adaptor trimming and filtering (bbduk) 
 */

if(params.do_bbduk){
    process trim_bbduk {
        tag "BBDUK"
        module "java/1.8.0_161"
        afterScript = 'rm -rf *.clumped*.fastq.gz'
        publishDir "${params.outdir}/trim_bbduk", mode: 'copy', 
            saveAs: {filename ->
                if (filename.indexOf(".out") > 0) "$filename"
                else null
            }

        input:
        set val(name), file(reads) from read_files_trim_bbduk

        output:
        set val(name), file ("*trimmed.fq.gz") into trimmed_reads, check_trim_reads, trimmed_reads_fastqc
        set val(name), file ("*.out") into trim_bbduk_results

        script:
        if (params.SingleEnd) {
            """
            echo "SingleEnd"
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in=${reads} out=${name}.clumped.fastq.gz \\
                markduplicates subs=0
            $Program_Dir/bbmap_38.20/bbduk.sh -Xmx8g qin=auto \\
                in=${name}.clumped.fastq.gz \\
                out=${name}_trimmed.fq.gz \\
                ref=$Program_Dir/bbmap_38.20/resources/adapters.fa \\
                stats=${name}.stats.out aqhist=${name}_aqhist.out \\
                ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=20 ziplevel=9  
            """
        } else {
            """
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in1=${reads[0]} in2=${reads[1]}   \\
                out1=${name}.clumped_R1.fastq.gz out2=${name}.clumped_R2.fastq.gz \\
                markduplicates subs=0
            $Program_Dir/bbmap_38.20/bbduk.sh -Xmx8g qin=auto \\
                in1=${name}.clumped_R1.fastq.gz in2=${name}.clumped_R2.fastq.gz \\
                out1=${name}_R1_trimmed.fq.gz out2=${name}_R2_trimmed.fq.gz \\
                ref=$Program_Dir/bbmap_38.20/resources/adapters.fa \\
                stats=${name}_stats.out aqhist=${name}_aqhist.out \\
                ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=20 ziplevel=9  
            """
        }
    }
check_trim_reads.subscribe {println "trimed read: $it"}
}

/*
 * STEP 2C - Adaptor trimming and filtering (fastp) 
 */

if(params.do_fastp){
    process trim_fastp {
        tag "fastp-trimm"
        afterScript = 'rm -rf *.clumped*.fastq.gz'
        module "gcc/5.4"
        publishDir "${params.outdir}/trim_fastp", mode: 'copy', 
            saveAs: {filename ->
                if (filename.indexOf("_out.html") > 0) "$filename"
                else if (filename.indexOf("_out.json") > 0) "$filename"
                else null
            }

        input:
        set val(name), file(reads) from read_files_trim_fastp


        output:
        set val(name), file ("*trimmed.fq.gz") into trimmed_reads, check_trim_reads, trimmed_reads_fastqc
        set val(name), file ("*_out.*") into trim_fastp_results

        script:
        afterScript = 'rm -rf ${name}.clumped*.fastq.gz'
        if (params.SingleEnd) {
            """
            echo "SingleEnd"
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in=${reads} out=${name}.clumped.fastq.gz \\
                markduplicates subs=0
            $Program_Dir/fastp/fastp -i ${name}.clumped.fastq.gz \\
                -o ${name}_trimmed.fq.gz \\
                -h ${name}_fastp_out.html \\
                -j ${name}_fastp_out.json -z 9
            """
        } else {
            """
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in1=${reads[0]} in2=${reads[1]}   \\
                out1=${name}.clumped_R1.fastq.gz out2=${name}.clumped_R2.fastq.gz \\
                markduplicates subs=0
            $Program_Dir/fastp/fastp -i ${name}.clumped_R1.fastq.gz \\
                -I ${name}.clumped_R2.fastq.gz -c \\
                -o ${name}_R1_trimmed.fq.gz -O ${name}_R2_trimmed.fq.gz \\
                -h ${name}_fastp_out.html \\
                -j ${name}_fastp_out.json -z 9
            """
        }
    }
check_trim_reads.subscribe {println "trimed read: $it"}
}

/*
 * STEP 3 - Quality control with FastQC after Trim
 */

process trimmed_fastqc {
    tag "trimmed_FastQC"
    publishDir "${params.outdir}/trimmed_fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into trimmed_fastqc_results

    script:
    """
    $Program_Dir/FastQC/fastqc -q $reads
    """
}

if (!params.do_trimgalore && !params.do_bbduk && !params.do_fastp) {
    read_files_trim_galore.into{trimmed_reads; check_trim_reads}
    check_trim_reads.subscribe {println "trim avoided read: $it"}
}

/*
 * STEP 4A - Build Bismark index if not provided
 */

Channel
    .fromPath(params.genome)
    .ifEmpty { exit 1, "Cannot find any genome file in: ${params.genome}\nNB" }
    .into {read_file_genome; ch_print }
ch_print.subscribe {println "read_file_genome: $it"}

if(params.bismark_index){
    Channel
        .fromPath(params.bismark_index)
        .ifEmpty { exit 1, "Bismark index not found: ${params.bismark_index}" }
        .into { bismark_index; bismark_index1; bismark_index2 }
}
else {
    process generateBismarkIndex {
        tag "BismarkIndex"
        publishDir "${params.outdir}/reference_genome", mode: 'copy' 
        module "bowtie2/2.2.9"

        input:
        file genome_file from read_file_genome

        output:
        file "BismarkIndex" into bismark_index, bismark_index1, bismark_index2

        script:
        """
        mkdir BismarkIndex
        cp $genome_file BismarkIndex/.
        ${Program_Dir}/Bismark/bismark_genome_preparation BismarkIndex
        """
    }
}

/*
 * STEP 4 - align with Bismark
 */

process bismark_alignment {
    tag "Bismark Align"
    module "samtools/1.3.1:bowtie2/2.2.9"
    publishDir "${params.outdir}/bismark_alignments", mode: 'copy',
          saveAs: {filename ->
                if (filename.indexOf("report.txt") > 0) "$filename"
                else null
           }	

    input:
    set val(name), file(reads) from trimmed_reads
    file index from bismark_index.collect()

    output:
    file "*.bam" into bismark_bam, bismark_bam_1, bismark_bam_2, bismark_bam_3
    file "*report.txt" into bismark_align_log_1, bismark_align_log_2, bismark_align_log_3
    script:
    if (params.SingleEnd) {
        """
        ${Program_Dir}/Bismark/bismark \\
            --multicore 1 \\
            --genome $index \\
            --path_to_bowtie ${Program_Dir}/bowtie2-2.3.4.2-linux-x86_64/ \\
            //--bowtie2 -D 5 -R 1 -N 0 -L 22 \\
            --bowtie2 --gzip\\
            $reads
        """
    } else {
        """
        ${Program_Dir}/Bismark/bismark \\
            --multicore 1 \\
            --genome $index \\
            --path_to_bowtie ${Program_Dir}/bowtie2-2.3.4.2-linux-x86_64/ \\
            //--bowtie2  -D 5 -R 1 -N 0 -L 22 \\
            --bowtie2 --gzip\\
            -1 ${reads[0]} \\
            -2 ${reads[1]}
        """
    }
}

/*
 * STEP 5 - Bismark deduplicate
 */

process bismark_deduplicate {
    tag "Bismark Deduplicate"
    module "samtools/1.3.1:bowtie2/2.2.9"
    publishDir "${params.outdir}/bismark_deduplicated", 
        saveAs: {filename -> filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename"}

    input:
    file bamfile from bismark_bam

    output:
    file "${bamfile.baseName}.deduplicated.bam" into bismark_bam_dedup, bismark_dedupbam_dmr
    file "${bamfile.baseName}.deduplication_report.txt" into bismark_dedup_log_1, bismark_dedup_log_2, bismark_dedup_log_3

    script:
    if (params.SingleEnd) {
        """
        ${Program_Dir}/Bismark/deduplicate_bismark -s --bam $bamfile
        """
    } else {
        """
        ${Program_Dir}/Bismark/deduplicate_bismark -p --bam $bamfile
        """
    }
}




/*
 * STEP 6 - Bismark methylation extraction
 */

process bismark_methXtract {
    tag "Bismark MethXtract"
    module "samtools/1.3.1:bowtie2/2.2.9"
    publishDir "${params.outdir}/bismark_methylation_calls", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("splitting_report.txt") > 0) "logs/$filename"
            else if (filename.indexOf("M-bias") > 0) "m-bias/$filename"
            else if (filename.indexOf(".cov") > 0) "methylation_coverage/$filename"
            else if (filename.indexOf("bedGraph") > 0) "bedGraph/$filename"
            else "methylation_calls/$filename"
        }

    input:
    file bam from bismark_bam_dedup
    file index from bismark_index1.collect()

    output:
    file "${bam.baseName}_splitting_report.txt" into bismark_splitting_report_1, bismark_splitting_report_2, bismark_splitting_report_3
    file "${bam.baseName}.M-bias.txt" into bismark_mbias_1, bismark_mbias_2, bismark_mbias_3
    file '*.{png,gz}' into bismark_methXtract_results

    script:
    if (params.SingleEnd) {
        """
        ${Program_Dir}/Bismark/bismark_methylation_extractor  \\
            --genome_folder $index \\
            --comprehensive --merge_non_CpG \\
            --bedGraph \\
            --counts \\
            --gzip \\
            -s \\
            --report \\
            $bam
        """
    } else {
        """
        ${Program_Dir}/Bismark/bismark_methylation_extractor \\
            --genome_folder $index \\
            --comprehensive --merge_non_CpG \\
            --bedGraph \\
            --counts \\
            --gzip \\
            -p \\
            --no_overlap \\
            --report \\
            $bam
        """
    }
}


/*
 * STEP 7 - Bismark Sample Report
 */

process bismark_report {
    tag "Bismark Sample Report"
    publishDir "${params.outdir}/bismark_reports", mode: 'copy'

    input:
    file bismark_align_log_1
    file bismark_dedup_log_1
    file bismark_splitting_report_1
    file bismark_mbias_1

    output:
    file '*{html,txt}' into bismark_reports_results

    script:
    name = bismark_align_log_1.toString() - ~/(_R1)?(_trimmed|_val_1).+$/
    """
    ${Program_Dir}/Bismark/bismark2report \\
        --alignment_report $bismark_align_log_1 \\
        --dedup_report $bismark_dedup_log_1 \\
        --splitting_report $bismark_splitting_report_1 \\
        --mbias_report $bismark_mbias_1
    """
}

/*
 * STEP 8 - Bismark Summary Report
 */
process bismark_summary {
    publishDir "${params.outdir}/bismark_summary", mode: 'copy'

    input:
    file ('*') from bismark_bam_2.collect()
    file ('*') from bismark_align_log_2.collect()
    file ('*') from bismark_dedup_log_2.collect()
    file ('*') from bismark_splitting_report_2.collect()
    file ('*') from bismark_mbias_2.collect()

    output:
    file '*{html,txt}' into bismark_summary_results

    script:
    """
    ${Program_Dir}/Bismark/bismark2summary
    """
}

/*
 * STEP 9 - dmrcov
 */
process dmr_covfiles {
    tag "dmr"
    module "samtools/1.3.1:bowtie2/2.2.9:intel/17.0.2:python/2.7.12"
    publishDir "${params.outdir}/DMR", mode: 'copy'

    input:
    file deduped_bam from bismark_dedupbam_dmr

    output:
    file "${deduped_bam.baseName}*.cov" into dmr_cov_files, dmr_cov_check 

    script:
    """
    samtools view -h $deduped_bam > ${deduped_bam}.samviewout
    python ${Program_Dir}/DMRfinder/extract_CpG_data.py \
           -i ${deduped_bam}.samviewout -o ${deduped_bam.baseName}.cov 
    rm ${deduped_bam}.samviewout
    """
}
dmr_cov_check.subscribe {println "dmr_cov_check: $it"}

/*
 * STEP 10 - dmrcsv
*/
process dmrcsv {
    tag "dmrcsv"
    module "python/2.7.12"
    publishDir "${params.outdir}/DMR", mode: 'copy'

    input:
    file cov_files from dmr_cov_files.collect()

    output:
    file "combined_covs.csv" into dmr_csv_file, dmr_csv_check
    script:
    """
    python ${Program_Dir}/DMRfinder/combine_CpG_sites.py -o combined_covs.csv $cov_files
    """
}


/*
 * STEP 11 - MultiQC
 */
process multiqc {
    tag "MultiQC report"
    module "intel/17.0.2:python/2.7.12"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    if (params.do_trimgalore) "file ('trim_galore/*') from trim_galore_results.collect()"
    if (params.do_bbduk) "file ('trim_bbduk/*') from trim_bbduk_results.collect()"
    if (params.do_fastp) "file ('trim_fastp/*') from trim_fastp_results.collect()"
    file ('trimmed_fastqc/*') from trimmed_fastqc_results.collect()
    file ('bismark/*') from bismark_align_log_3.collect()
    file ('bismark/*') from bismark_dedup_log_3.collect()
    file ('bismark/*') from bismark_splitting_report_3.collect()
    file ('bismark/*') from bismark_mbias_3.collect()
    file ('bismark/*') from bismark_reports_results.collect()
    file ('bismark/*') from bismark_summary_results.collect()

    output:
    file "*_report.html" into multiqc_report
    file "*_data"
    file '.command.err' into multiqc_stderr

    script:
    rtitle = "DNA-Methylation analysis on OARC HPC Cluster"
    """
    set +u
    source ${Program_Dir}/PyProgVirt/mqc/bin/activate
    ${Program_Dir}/PyProgVirt/mqc/bin/multiqc -h
    ${Program_Dir}/PyProgVirt/mqc/bin/multiqc -d ${params.outdir} \\
         -i "DNA-Meth" \\
         -b "DNA Methylation MultiQC report"
    deactivate
    set -u
    """

}


workflow.onComplete {
    println "WORKFLOW SUMMARY"
    println "Pipeline completed at: $workflow.complete"
    println "Duration: $workflow.duration"
    println "WorkDir: $workflow.workDir"
    println "Exit Status: $workflow.exitStatus"
    println ( workflow.complete ? "Okay" : "Oops ..NOT OKAY" )
}

