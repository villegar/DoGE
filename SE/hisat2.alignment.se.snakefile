rule alignment:
    input:
        index = rules.hisat2_index.output,
        reads = rules.trim_reads.output
    output:
        "4.ALIGNMENT/{raw_reads}_sorted.sam"
    log:
        "4.ALIGNMENT/{raw_reads}_sam.log"
    message:
        "Genome alignment"
    threads:
        CPUS_ALIGNMENT
    shell:
        "hisat2 --phred33 -p {threads} --qc-filter -x {input.index}/idx \
        -U {input.reads} -S {output} 2> {log}"
