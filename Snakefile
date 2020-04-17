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
TRIMMOMATIC_OPTIONS = config["trimmomatic"]
if PAIRED_END:
    ENDS = el(["_"],[FORWARD_READ_ID,REVERSE_READ_ID])
    ENDS = [FORWARD_READ_ID,REVERSE_READ_ID]
    FORWARD_READ_ID = [FORWARD_READ_ID]
    MODE = ["paired","unpaired"]
    MODE = ["","_un"]
    #MODE = ["1","2"]
    REVERSE_READ_ID = [REVERSE_READ_ID]
    SUFFIX = "_" + FORWARD_READ_ID[0] + "." + EXTENSION
    MODE = [""]

else:
    ENDS = []
    FORWARD_READ_ID = []
    MODE = []
    REVERSE_READ_ID = []
    SUFFIX = "." + EXTENSION

LIBS = findLibraries(READS,PREFIX,SUFFIX)

###### Multithread configuration #####
CPUS_FASTQC = 3
CPUS_PHIX = 15
CPUS_TRIMMING = 5
CPUS_HISAT2_INDEX = 20
CPUS_ALIGNMENT = 10
CPUS_STAR = 20
CPUS_ARIA = 16
CPUS_KRAKEN = 20
CPUS_READCOUNTS = 5
CPUS_RNA = 20

ADAPTER = which("trimmomatic")

####### Reference datasets #######
GENOME = config["genome"]
GENOME_FILENAMES = extractFilenames(GENOME.keys(),".gz")
#GENOME4STAR = config["genome4star"]
#GENOME4STAR_FILENAMES = extractFilenames(GENOME4STAR.keys(),".gz")
#GENOME4PHIX = config["genome4phiX"]
#KRAKEN_DB = config["krakenDB"]
#KRAKEN_DB_FILENAMES = extractFilenames(KRAKEN_DB.keys(),".tgz")
#rRNA = config["rRNAref"]
#rRNA_FILES = list(rRNA.keys())

RAW_ENDS = [""]
if PAIRED_END:
    RAW_ENDS = el(["_"],ENDS)
    RAW_ENDS = ENDS

####### Rules #######
rule all:
    input:
        expand("1.QC.RAW/{raw_reads}{raw_ends}_fastqc.{format}",
            raw_reads = LIBS, raw_ends = el(["_"],ENDS), format = ["html","zip"]),
        # expand("2.TRIMMED/{raw_reads2}{raw_ends2}.{format}",
        #     raw_reads2 = LIBS, raw_ends2 = RAW_ENDS, format = EXTENSION),
        expand("3.QC.TRIMMED/{raw_reads}{raw_ends}_fastqc.{format}",
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"]),
        expand("4.ALIGNMENT/{raw_reads}_sorted.bam", raw_reads = LIBS),
        expand("5.QC.ALIGNMENT/{raw_reads}_stats.txt", raw_reads = LIBS),
        expand("6.COUNTS/counts.{format}", format = ["txt","matrix"])
        #READS + "/{raw_reads}_sorted.bam"
    output:
        logs 	= directory("0.LOGS"),
        reports	= directory("10.MULTIQC")
    run:
        shell("multiqc -o {output.reports} -n 1.Report_FastQC_Raw.html -d 1.QC.RAW")
        shell("multiqc -o {output.reports} -n 2.Report_Trimming.html -d 2.TRIMMED")
        shell("multiqc -o {output.reports} -n 3.Report_FastQC_Trimmed.html -d 3.QC.TRIMMED")
        shell("multiqc -o {output.reports} -n 4.Report_Alignment.html -d 4.ALIGNMENT")
        shell("multiqc -o {output.reports} -n 5.Report_AlignmentQC.html -d 5.QC.ALIGNMENT")
        shell("mkdir -p {output.logs} && mv *.log {output.logs}")

rule fastqc_raw:
    input:
        reads = READS + "/{raw_reads}{raw_ends}." + EXTENSION
    output:
        html = "1.QC.RAW/{raw_reads}{raw_ends}_fastqc.html",
        zip  = "1.QC.RAW/{raw_reads}{raw_ends}_fastqc.zip"
    message:
        "FastQC on raw data"
    log:
        "1.QC.RAW/{raw_reads}{raw_ends}.log"
    threads:
        CPUS_FASTQC
    run:
        shell("fastqc -o 1.QC.RAW -t {threads} {input.reads} 2> {log}")

if PAIRED_END:
    rule trim_reads:
        input:
            adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters")
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

else:
    rule trim_reads:
        input:
            adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters"),
            reads   = READS + "/{raw_reads}." + EXTENSION
        output:
            "2.TRIMMED/{raw_reads}." + EXTENSION
        params:
            options = TRIMMOMATIC_OPTIONS
        log:
            "2.TRIMMED/{raw_reads}.log"
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
            html = "3.QC.TRIMMED/{raw_reads}{raw_ends}_fastqc.html",
            zip  = "3.QC.TRIMMED/{raw_reads}{raw_ends}_fastqc.zip"
        message:
            "FastQC on trimmed data"
        log:
            "3.QC.TRIMMED/{raw_reads}{raw_ends}.log"
        threads:
            CPUS_FASTQC
        shell:
            "fastqc -o 3.QC.TRIMMED -t {threads} {input} 2> {log}"
else:
    rule fastqc_trimmed:
        input:
            rules.trim_reads.output
        output:
            html = "3.QC.TRIMMED/{raw_reads}_fastqc.html",
            zip  = "3.QC.TRIMMED/{raw_reads}_fastqc.zip"
        message:
            "FastQC on trimmed data"
        log:
            "3.QC.TRIMMED/{raw_reads}.log"
        threads:
            CPUS_FASTQC
        shell:
            "fastqc -o 3.QC.TRIMMED -t {threads} {input} 2> {log}"

rule hisat2_index:
    input:
        expand("GENOME/{file}", hisat2 = ["HISAT2_INDEX"], file = GENOME_FILENAMES[1])
        # primary_assembly = GENOME_FILENAMES[0]
    output:
        directory("GENOME/HISAT2_INDEX")
    #    expand("GENOME/{hisat2}/{file}", hisat2 = ["HISAT2_INDEX"], file = GENOME_FILENAMES)
    log:
        "hisat2_index.log"
    message:
        "Creating HISAT2 index"
    threads:
        CPUS_HISAT2_INDEX
    shell:
        "mkdir -p GENOME/HISAT2_INDEX && hisat2-build -p {threads} {input} {output}/idx 2> {log}"

if PAIRED_END:
    rule alignment:
        input:
            index = rules.hisat2_index.output,
            r1    = rules.trim_reads.output.r1,
            r2    = rules.trim_reads.output.r2,
            # r1    = READS + "/{raw_reads}" + ENDS[0] + "." + EXTENSION,
            # r2    = READS + "/{raw_reads}" + ENDS[1] + "." + EXTENSION
        output:
            "4.ALIGNMENT/{raw_reads}_sorted.sam"
            # "4.ALIGNMENT/{raw_reads}{raw_ends}_sorted.bam"
        log:
            "4.ALIGNMENT/{raw_reads}_sam.log"
            # "4.ALIGNMENT/{raw_reads}{raw_ends}.log"
        message:
            "Genome alignment"
        threads:
            CPUS_ALIGNMENT
        shell:
            "hisat2 --phred33 -p {threads} --qc-filter -x {input.index}/idx -1 {input.r1} -2 {input.r2} -S {output} 2> {log}"
            #"hisat2 --phred33 -p {threads} --qc-filter -x {input.index}/idx -1 {input.r1} -2 {input.r2} | samtools view -@ {threads} -bS - | samtools sort -@ {threads} -o {output}"

else:
    rule alignment:
        input:
            index = rules.hisat2_index.output,
            reads = rules.trim_reads.output
            # reads = READS + "/{raw_reads}." + EXTENSION
        output:
            "4.ALIGNMENT/{raw_reads}_sorted.sam"
        log:
            "4.ALIGNMENT/{raw_reads}_sam.log"
        message:
            "Genome alignment"
        threads:
            CPUS_ALIGNMENT
        shell:
            "hisat2 --phred33 -p {threads} --qc-filter -x {input.index}/idx -U {input.reads} -S {output} 2> {log}"
            #"hisat2 --phred33 -p {threads} --qc-filter -x {input.index}/idx -U {input.reads} | samtools view -@ {threads} -bS - | samtools sort -@ {threads} -o {output}"

rule sam2bam:
    input:
        rules.alignment.output
    output:
        "4.ALIGNMENT/{raw_reads}_sorted.bam"
    log:
        "4.ALIGNMENT/{raw_reads}_bam.log"
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
        "5.QC.ALIGNMENT/{raw_reads}_stats.txt"
    log:
        "5.QC.ALIGNMENT/{raw_reads}_stats.log"
    message:
        "Assessing alignment quality"
    threads:
        CPUS_ALIGNMENT
    shell:
        "SAMstatsParallel --sorted_sam_file {input} --outf {output} --threads {threads} 2> {log}"

rule feature_counts:
    input:
        gtf = expand("GENOME/{file}", file = GENOME_FILENAMES[0]),
        bam = el(["4.ALIGNMENT/"],el(LIBS,["_sorted.bam"]))
    output:
        "6.COUNTS/counts.txt"
    threads:
        20
    shell:
        "mkdir -p 6.COUNTS && \
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
        "6.COUNTS/counts.matrix"
    shell:
        "cat {input.counts} | grep -v '^#' | cut -f {params.cols} | \
        sed '1d' | sed '1i\Geneid{params.libs}' > {output}"
