rule feature_counts:
    input:
        gtf = expand(REF_GENOME + "{file}", file = GENOME_FILENAMES["GTF"]),
        bam = el([ALIGNMENT],el(LIBS,["_sorted.bam"]))
    output:
        COUNTS + "counts.txt"
    threads:
        CPUS_READCOUNTS
    shell:
        "mkdir -p " + COUNTS + " && \
        featureCounts -T {threads} -t exon -g gene_id -Q 30 -F GTF \
        -a {input.gtf} \
        -o {output} \
        {input.bam}"
