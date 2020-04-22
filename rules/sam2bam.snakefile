rule sam2bam:
    input:
        rules.alignment.output
    output:
        ALIGNMENT + "{raw_reads}_sorted.bam"
    log:
        ALIGNMENT + "{raw_reads}_bam.log"
    message:
        "Converting SAM to BAM"
    threads:
        CPUS_ALIGNMENT
    shell:
        "samtools view -@ {threads} -bS {input} | samtools sort -@ {threads} -o {output} 2> {log}"
