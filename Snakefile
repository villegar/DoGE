####### Libraries #######
from utils import extractFilenames, findLibraries, which
from utils import expand_list as el

####### Global variables #######
EXTENSION = config["reads"]["extension"]
PREFIX = config["reads"]["prefix"]
READS = config["reads"]["path"]
FORWARD_READ_ID = config["reads"]["forward_read_id"]
REVERSE_READ_ID = config["reads"]["reverse_read_id"]
PAIRED_END = [True if config["reads"]["end_type"] == "pe" else False][0]
TRIMMOMATIC_OPTIONS = config["trimmomatic"]["options"]
if PAIRED_END:
    ENDS = el(["_"],[FORWARD_READ_ID,REVERSE_READ_ID])
    ENDS = [FORWARD_READ_ID,REVERSE_READ_ID]
    FORWARD_READ_ID = [FORWARD_READ_ID]
    REVERSE_READ_ID = [REVERSE_READ_ID]
    SUFFIX = "_" + FORWARD_READ_ID[0] + "." + EXTENSION
else:
    ENDS = []
    FORWARD_READ_ID = []
    REVERSE_READ_ID = []
    SUFFIX = "." + EXTENSION

LIBS = findLibraries(READS,PREFIX,SUFFIX)

###### Multithread configuration #####
CPUS_FASTQC = 4
CPUS_TRIMMING = 5
CPUS_HISAT2_INDEX = 40
CPUS_ALIGNMENT = 10
CPUS_READCOUNTS = 20

ADAPTER = which("trimmomatic")

####### Output directories #######
REF_GENOME = "GENOME/"
RAW_FASTQC = "1.QC.RAW/"
TRIMMED_READS = "2.TRIMMED/"
TRIMMED_READS_FASTQC = "3.QC.TRIMMED/"
ALIGNMENT = "4.ALIGNMENT/"
ALIGNMENT_QC = "5.QC.ALIGNMENT/"
COUNTS = "6.COUNTS/"
RMD = "7.RMD/"

####### Reference datasets #######
GENOME = config["genome"]
GENOME_FILENAMES = extractFilenames(GENOME.keys(),".gz")

RAW_ENDS = [""]
if PAIRED_END:
    RAW_ENDS = el(["_"],ENDS)

####### Rules #######
rule all:
    input:
        expand(RAW_FASTQC + "{raw_reads}{raw_ends}_fastqc.{format}",
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"]),
        expand(TRIMMED_READS_FASTQC + "{raw_reads}{raw_ends}_fastqc.{format}",
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"]),
        expand(ALIGNMENT_QC + "{raw_reads}_stats.txt", raw_reads = LIBS),
        #expand(COUNTS + "counts.{format}", format = ["txt","matrix"]),
        RMD + "doge_report_simple.html"
    output:
        expand(RMD + "Report_{step}.html",
        step = ["FastQC_Raw","Trimming","FastQC_Trimmed","Alignment"])
    params:
        logs 	= directory("0.LOGS"),
        reports	= directory(RMD)
    run:
        shell("multiqc -f -o {params.reports} -n Report_FastQC_Raw.html -d " + RAW_FASTQC)
        shell("multiqc -f -o {params.reports} -n Report_Trimming.html -d " + TRIMMED_READS)
        shell("multiqc -f -o {params.reports} -n Report_FastQC_Trimmed.html -d " + TRIMMED_READS_FASTQC)
        shell("multiqc -f -o {params.reports} -n Report_Alignment.html -d " + ALIGNMENT)
        shell("mkdir -p {params.logs} && mv *.log {params.logs}")
        shell("tar -zvcf 7.RMD.tar.gz 7.RMD")

rule fastqc_raw:
    input:
        reads = READS + "/{raw_reads}{raw_ends}." + EXTENSION
    output:
        html = RAW_FASTQC + "{raw_reads}{raw_ends}_fastqc.html",
        zip  = RAW_FASTQC + "{raw_reads}{raw_ends}_fastqc.zip"
    message:
        "FastQC on raw data"
    log:
        RAW_FASTQC + "{raw_reads}{raw_ends}.log"
    threads:
        CPUS_FASTQC
    shell:
        "fastqc -o " + RAW_FASTQC + " -t {threads} {input.reads} 2> {log}"

if PAIRED_END:
    rule trim_reads:
        input:
            adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters"),
            r1      = READS + "/{raw_reads}_" + ENDS[0] + "." + EXTENSION,
            r2      = READS + "/{raw_reads}_" + ENDS[1] + "." + EXTENSION
        output:
            r1    = TRIMMED_READS + "{raw_reads}_" + ENDS[0] + "." + EXTENSION,
            r1_un = TRIMMED_READS + "{raw_reads}_" + ENDS[0] + "_un." + EXTENSION,
            r2    = TRIMMED_READS + "{raw_reads}_" + ENDS[1] + "." + EXTENSION,
            r2_un = TRIMMED_READS + "{raw_reads}_" + ENDS[1] + "_un." + EXTENSION
        params:
            options = TRIMMOMATIC_OPTIONS
        log:
            TRIMMED_READS + "{raw_reads}.log"
        message:
            "Using Paired End Trimming"
        threads:
            CPUS_TRIMMING
        shell:
            "trimmomatic PE -threads {threads} {input.r1} {input.r2} \
            {output} {params.options} 2> {log}"

else:
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
            "trimmomatic SE -threads {threads} {input.reads} {output} {params.options} 2> {log}"

if PAIRED_END:
    rule fastqc_trimmed:
        input:
            rules.trim_reads.output
        output:
            html = TRIMMED_READS_FASTQC + "{raw_reads}{raw_ends}_fastqc.html",
            zip  = TRIMMED_READS_FASTQC + "{raw_reads}{raw_ends}_fastqc.zip"
        message:
            "FastQC on trimmed data"
        log:
            TRIMMED_READS_FASTQC + "{raw_reads}{raw_ends}.log"
        threads:
            CPUS_FASTQC
        shell:
            "fastqc -o " + TRIMMED_READS_FASTQC + " -t {threads} {input} 2> {log}"
else:
    rule fastqc_trimmed:
        input:
            rules.trim_reads.output
        output:
            html = TRIMMED_READS_FASTQC + "{raw_reads}_fastqc.html",
            zip  = TRIMMED_READS_FASTQC + "{raw_reads}_fastqc.zip"
        message:
            "FastQC on trimmed data"
        log:
            TRIMMED_READS_FASTQC + "{raw_reads}.log"
        threads:
            CPUS_FASTQC
        shell:
            "fastqc -o " + TRIMMED_READS_FASTQC + " -t {threads} {input} 2> {log}"

rule hisat2_index:
    input:
        expand(REF_GENOME + "{file}", file = GENOME_FILENAMES[1])
    output:
        directory(REF_GENOME + "HISAT2_INDEX")
    log:
        "hisat2_index.log"
    message:
        "Creating HISAT2 index"
    threads:
        CPUS_HISAT2_INDEX
    shell:
        "mkdir -p " + REF_GENOME + "HISAT2_INDEX && hisat2-build -p {threads} {input} {output}/idx 2> {log}"

if PAIRED_END:
    rule alignment:
        input:
            index = rules.hisat2_index.output,
            r1    = rules.trim_reads.output.r1,
            r2    = rules.trim_reads.output.r2,
        output:
            ALIGNMENT + "{raw_reads}_sorted.sam"
        log:
            ALIGNMENT + "{raw_reads}_sam.log"
        message:
            "Genome alignment"
        threads:
            CPUS_ALIGNMENT
        shell:
            "hisat2 --phred33 -p {threads} --qc-filter -x {input.index}/idx -1 {input.r1} -2 {input.r2} -S {output} 2> {log}"

else:
    rule alignment:
        input:
            index = rules.hisat2_index.output,
            reads = rules.trim_reads.output
        output:
            ALIGNMENT + "{raw_reads}_sorted.sam"
        log:
            ALIGNMENT + "{raw_reads}_sam.log"
        message:
            "Genome alignment"
        threads:
            CPUS_ALIGNMENT
        shell:
            "hisat2 --phred33 -p {threads} --qc-filter -x {input.index}/idx -U {input.reads} -S {output} 2> {log}"

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

rule alignment_quality:
    input:
        rules.alignment.output
    output:
        ALIGNMENT_QC + "{raw_reads}_stats.txt"
    log:
        ALIGNMENT_QC + "{raw_reads}_stats.log"
    message:
        "Assessing alignment quality"
    shell:
        "SAMstats --sorted_sam_file {input} --outf {output} 2> {log}"

rule feature_counts:
    input:
        gtf = expand(REF_GENOME + "{file}", file = GENOME_FILENAMES[0]),
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

rule quantification_table:
    input:
        counts = rules.feature_counts.output
    params:
        libs = "".join(el(["\\t"],LIBS)),
        cols = "1," + ",".join([str(i) for i in range(7, 7 + len(LIBS))])
    output:
        COUNTS + "counts.matrix"
    shell:
        "cat {input.counts} | grep -v '^#' | cut -f {params.cols} | \
        sed '1d' | sed '1i\Geneid{params.libs}' > {output}"

rule annotation_table:
    input:
        gtf = expand(REF_GENOME + "{file}", file = GENOME_FILENAMES[0])
    output:
        RMD + "gene_annotation.txt"
    shell:
        "sed '/^[[:blank:]]*#/d;s/#.*//' {input.gtf} | awk '($3 == \"gene\")' | \
        awk -F';' '$1=$1' OFS='\\t' > {output}"

rule rmd_report:
    input:
        annotation = rules.annotation_table.output,
        counts = rules.quantification_table.output,
        experiment = "exp_design.csv"
    output:
        RMD + "doge_report_simple.html"
    shell:
        "cp html/* rmd/* " + RMD + " && " +
        "cp {input.counts} " + RMD + " && " +
        "cp {input.experiment} " + RMD + " && cd " + RMD + " && " +
        "Rscript -e \"rmarkdown::render(\'doge_report_simple.Rmd\'," +
        "output_dir=\'.\', clean = TRUE, quiet = TRUE)\" && " +
        "cd ../"
