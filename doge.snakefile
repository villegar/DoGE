####### Libraries #######
from utils import extractFilenames, findLibraries, loadGenome, verifyGenome, which
from utils import expand_list as el

####### Global variables #######
EXTENSION = config["reads"]["extension"]
PREFIX = config["reads"]["prefix"]
READS = config["reads"]["path"]
PAIRED_END = [True if config["reads"]["end_type"] == "pe" else False][0]
try:
    TRIMMOMATIC_OPTIONS = config["trimmomatic"]["options"]
except:
    raise ValueError("trimmomatic > options not found in the configuration file")
if PAIRED_END:
    FORWARD_READ_ID = [config["reads"]["forward_read_id"]]
    REVERSE_READ_ID = [config["reads"]["reverse_read_id"]]
    #ENDS = el(["_"],[FORWARD_READ_ID,REVERSE_READ_ID])
    ENDS = [FORWARD_READ_ID,REVERSE_READ_ID]
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
FA,GTF = loadGenome(config["genome"])
GENOME_FILENAMES = {"FA":FA,"GTF":GTF}
verifyGenome(config["genome"],REF_GENOME + FA, REF_GENOME + GTF)

RAW_ENDS = [""]
if PAIRED_END:
    RAW_ENDS = el(["_"],ENDS)

####### Rules #######
if PAIRED_END:
    include:"./PE/all.pe.snakefile"
else:
    include:"./SE/all.se.snakefile"

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


rule hisat2_index:
    input:
        expand(REF_GENOME + "{file}", file = GENOME_FILENAMES["FA"])
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
        gtf = expand(REF_GENOME + "{file}", file = GENOME_FILENAMES["GTF"])
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
