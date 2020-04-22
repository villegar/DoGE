rule trim_reads:
    input:
        adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters"),
        r1      = READS + "/{raw_reads}_" + RAW_ENDS[0] + "." + EXTENSION,
        r2      = READS + "/{raw_reads}_" + RAW_ENDS[1] + "." + EXTENSION
    output:
        r1    = TRIMMED_READS + "{raw_reads}_" + RAW_ENDS[0] + "." + EXTENSION,
        r1_un = TRIMMED_READS + "{raw_reads}_" + RAW_ENDS[0] + "_un." + EXTENSION,
        r2    = TRIMMED_READS + "{raw_reads}_" + RAW_ENDS[1] + "." + EXTENSION,
        r2_un = TRIMMED_READS + "{raw_reads}_" + RAW_ENDS[1] + "_un." + EXTENSION
    params:
        options = TRIMMOMATIC_OPTIONS
    log:
        TRIMMED_READS + "{raw_reads}.log"
    message:
        "Using Paired End Trimming"
    threads:
        CPUS_TRIMMING
    shell:
        "trimmomatic PE -threads {threads} {params.options} {input.r1} "+
        "{input.r2} {output}  2> {log}"
