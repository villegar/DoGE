rule hisat2_index:
    input:
        expand("GENOME/{file}", hisat2 = ["HISAT2_INDEX"], file = GENOME_FILENAMES[1])
    output:
        directory("GENOME/HISAT2_INDEX")
    log:
        "hisat2_index.log"
    message:
        "Creating HISAT2 index"
    threads:
        CPUS_HISAT2_INDEX
    shell:
        "mkdir -p GENOME/HISAT2_INDEX && hisat2-build -p {threads} {input} {output}/idx 2> {log}"
