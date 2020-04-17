rule trim_reads:
    input:
        adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters"),
        r1      = READS + "/{raw_reads}_" + RAW_ENDS[0] + "." + EXTENSION,
        r2      = READS + "/{raw_reads}_" + RAW_ENDS[1] + "." + EXTENSION
    output:
        r1    = "2.TRIMMED/{raw_reads}_" + RAW_ENDS[0] + "." + EXTENSION,
        r1_un = "2.TRIMMED/{raw_reads}_" + RAW_ENDS[0] + "_un." + EXTENSION,
        r2    = "2.TRIMMED/{raw_reads}_" + RAW_ENDS[1] + "." + EXTENSION,
        r2_un = "2.TRIMMED/{raw_reads}_" + RAW_ENDS[1] + "_un." + EXTENSION
    params:
        options = TRIMMOMATIC_OPTIONS
    log:
        "2.TRIMMED/{raw_reads}.log"
    message:
        "Using Paired End Trimming"
    threads:
        CPUS_TRIMMING
    shell:
        "trimmomatic PE -threads {threads} {input.r1} {input.r2} \
        {output} {params.options} 2> {log}"
