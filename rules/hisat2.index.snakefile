rule hisat2_index:
    input:
        expand(REF_GENOME + "{file}", file = GENOME_FILENAMES["FA"])
    output:
        directory(REF_GENOME + "HISAT2_INDEX")
    log:
        "hisat2_index.log"
    message:
        "Creating HISAT2 index"
    threads:
        CPUS_HISAT2_INDEX
    shell:
        "mkdir -p " + REF_GENOME + "HISAT2_INDEX && hisat2-build -p {threads} {input} {output}/idx 2> {log}"
