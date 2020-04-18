rule trim_reads:
    input:
        adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters"),
        reads   = READS + "/{raw_reads}." + EXTENSION
    output:
        TRIMMED_READS + "{raw_reads}." + EXTENSION
    params:
        options = TRIMMOMATIC_OPTIONS
    log:
        TRIMMED_READS + "{raw_reads}.log"
    message:
        "Using Single End Trimming"
    threads:
        CPUS_TRIMMING
    shell:
        "trimmomatic SE -threads {threads} {input.reads} {output} \
        {params.options} 2> {log}"
