rule all:
    input:
        expand("1.QC.RAW/{raw_reads}{raw_ends}_fastqc.{format}",
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"]),
        expand("3.QC.TRIMMED/{raw_reads}{raw_ends}_fastqc.{format}",
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"]),
        expand("5.QC.ALIGNMENT/{raw_reads}_stats.txt", raw_reads = LIBS),
        "7.RMD/doge_report.html"
    output:
        logs 	= directory("0.LOGS"),
        reports	= directory("10.MULTIQC")
    run:
        shell("multiqc -o {output.reports} -n 1.Report_FastQC_Raw.html -d 1.QC.RAW")
        shell("multiqc -o {output.reports} -n 2.Report_Trimming.html -d 2.TRIMMED")
        shell("multiqc -o {output.reports} -n 3.Report_FastQC_Trimmed.html -d 3.QC.TRIMMED")
        shell("multiqc -o {output.reports} -n 4.Report_Alignment.html -d 4.ALIGNMENT")
        #shell("multiqc -o {output.reports} -n 5.Report_AlignmentQC.html -d 5.QC.ALIGNMENT")
        shell("mkdir -p {output.logs} && mv *.log {output.logs}")
