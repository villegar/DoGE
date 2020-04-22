rule all:
    input:
        expand(RAW_FASTQC + "{raw_reads}_fastqc.{format}",
            raw_reads = LIBS, format = ["html","zip"]),
        expand(TRIMMED_READS_FASTQC + "{raw_reads}_fastqc.{format}",
            raw_reads = LIBS, format = ["html","zip"]),
        expand(ALIGNMENT_QC + "{raw_reads}_stats.txt", raw_reads = LIBS),
        RMD + "doge_report_simple.html"
    output:
        expand(RMD + "Report_{step}.html",
        step = ["FastQC_Raw","Trimming","FastQC_Trimmed","Alignment"])
    run:
        shell("multiqc -f -o {params.reports} -n Report_FastQC_Raw.html -d " + RAW_FASTQC)
        shell("multiqc -f -o {params.reports} -n Report_Trimming.html -d " + TRIMMED_READS)
        shell("multiqc -f -o {params.reports} -n Report_FastQC_Trimmed.html -d " + TRIMMED_READS_FASTQC)
        shell("multiqc -f -o {params.reports} -n Report_Alignment.html -d " + ALIGNMENT)
        shell("mkdir -p {params.logs} && mv *.log {params.logs}")
        shell("tar -zvcf 7.RMD.tar.gz 7.RMD")
