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
        "trimmomatic SE -threads {threads} {params.options} " +
        "{input.reads} {output} 2> {log}"
