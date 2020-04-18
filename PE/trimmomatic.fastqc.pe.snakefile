rule fastqc_trimmed:
    input:
        rules.trim_reads.output
    output:
        html = "3.QC.TRIMMED/{raw_reads}{raw_ends}_fastqc.html",
        zip  = "3.QC.TRIMMED/{raw_reads}{raw_ends}_fastqc.zip"
    message:
        "FastQC on trimmed data"
    log:
        "3.QC.TRIMMED/{raw_reads}{raw_ends}.log"
    threads:
        CPUS_FASTQC
    shell:
        "fastqc -o 3.QC.TRIMMED -t {threads} {input} 2> {log}"
