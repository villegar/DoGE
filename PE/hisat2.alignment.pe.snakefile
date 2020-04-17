rule alignment:
    input:
        index = rules.hisat2_index.output,
        r1    = rules.trim_reads.output.r1,
        r2    = rules.trim_reads.output.r2
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
         -1 {input.r1} -2 {input.r2} -S {output} 2> {log}"
