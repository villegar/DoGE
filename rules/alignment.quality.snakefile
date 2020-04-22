rule alignment_quality:
    input:
        rules.alignment.output
    output:
        ALIGNMENT_QC + "{raw_reads}_stats.txt"
    log:
        ALIGNMENT_QC + "{raw_reads}_stats.log"
    message:
        "Assessing alignment quality"
    shell:
        "SAMstats --sorted_sam_file {input} --outf {output} 2> {log}"
