rule fastqc_raw:
    input:
        reads = READS + "/{raw_reads}." + EXTENSION
    output:
        html = RAW_FASTQC + "{raw_reads}_fastqc.html",
        zip  = RAW_FASTQC + "{raw_reads}_fastqc.zip"
    message:
        "FastQC on raw data"
    log:
        RAW_FASTQC + "{raw_reads}.log"
    threads:
        CPUS_FASTQC
    run:
        "fastqc -o 1.QC.RAW -t {threads} {input.reads} 2> {log}"
