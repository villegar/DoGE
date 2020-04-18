rule fastqc_raw:
    input:
        reads = READS + "/{raw_reads}." + EXTENSION
    output:
        html = "1.QC.RAW/{raw_reads}_fastqc.html",
        zip  = "1.QC.RAW/{raw_reads}_fastqc.zip"
    message:
        "FastQC on raw data"
    log:
        "1.QC.RAW/{raw_reads}.log"
    threads:
        CPUS_FASTQC
    run:
        "fastqc -o 1.QC.RAW -t {threads} {input.reads} 2> {log}"
