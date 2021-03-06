rule fastqc_trimmed:
    input:
        rules.trim_reads.output
    output:
        html = TRIMMED_READS_FASTQC + "{raw_reads}_fastqc.html",
        zip  = TRIMMED_READS_FASTQC + "{raw_reads}_fastqc.zip"
    message:
        "FastQC on trimmed data"
    log:
        TRIMMED_READS_FASTQC + "{raw_reads}.log"
    threads:
        CPUS_FASTQC
    shell:
        "fastqc -o " + TRIMMED_READS_FASTQC + " -t {threads} {input} 2> {log}"
