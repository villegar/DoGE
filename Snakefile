####### Libraries #######
from utils import extractFilenames, fastqc, findLibraries, which
from utils import expand_list as el

####### Global variables #######
EXTENSION = config["reads"]["extension"]
PREFIX = config["reads"]["prefix"]
READS = config["reads"]["path"]
FORWARD_READ_ID = config["reads"]["forward_read_id"]
REVERSE_READ_ID = config["reads"]["reverse_read_id"]
PAIRED_END = [True if config["reads"]["end_type"] == "pe" else False][0]
if PAIRED_END:
    DIRECTION = ["_forward_","_reverse_"]
    # #DIRECTION = ["_1_","_2_"]
    # DIRECTION = el(["_"],el([FORWARD_READ_ID,REVERSE_READ_ID],["_"]))
    # DIRECTION = el(["_"],[FORWARD_READ_ID,REVERSE_READ_ID])
    ENDS = [FORWARD_READ_ID,REVERSE_READ_ID]
    FORWARD_READ_ID = [FORWARD_READ_ID]
    MODE = ["paired","unpaired"]
    MODE = ["","_un"]
    #MODE = ["1","2"]
    REVERSE_READ_ID = [REVERSE_READ_ID]
    SUFFIX = "_" + FORWARD_READ_ID[0] + "." + EXTENSION
    DIRECTION = [""]
    MODE = [""]

else:
    DIRECTION = []
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
CPUS_STAR = 20
CPUS_ARIA = 16
CPUS_KRAKEN = 20
CPUS_READCOUNTS = 5
CPUS_RNA = 20

ADAPTER = which("trimmomatic")

####### Reference datasets #######
GENOME = config["genome"]
#GENOME4STAR = config["genome4star"]
#GENOME4STAR_FILENAMES = extractFilenames(GENOME4STAR.keys(),".gz")
GENOME4PHIX = config["genome4phiX"]
#KRAKEN_DB = config["krakenDB"]
#KRAKEN_DB_FILENAMES = extractFilenames(KRAKEN_DB.keys(),".tgz")
#rRNA = config["rRNAref"]
#rRNA_FILES = list(rRNA.keys())

print("Raw Names")
RAW_ENDS = [""]
RAW_LIBS = LIBS
TRM_LIBS = LIBS
TRM_LIBS_OUT = [""]
RAW_LIBS_R1 = LIBS
RAW_LIBS_R2 = LIBS
if (len(ENDS) > 0):
    RAW_ENDS = el(["_"],ENDS)
    RAW_LIBS = el(LIBS,el(["_"],ENDS))
    TRM_LIBS = el(LIBS,el(DIRECTION,MODE))
    RAW_LIBS_R1 = el(LIBS,el(["_"],FORWARD_READ_ID))
    RAW_LIBS_R2 = el(LIBS,el(["_"],REVERSE_READ_ID))
    TRM_LIBS_OUT = el(DIRECTION,MODE)

print(el(["3.QC.TRIMMED/"],el(LIBS,el([""],TRM_LIBS_OUT))))
print("ENDS")
print(ENDS)

####### Rules #######
rule all:
    input:
        expand("1.QC.RAW/{raw_reads}{raw_ends}_fastqc.{format}",
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"]),
        expand("3.QC.TRIMMED/{raw_reads}{raw_ends}_fastqc.{format}",
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"])
    output:
        logs 	= directory("0.LOGS"),
        reports	= directory("10.MULTIQC")
    run:
        shell("multiqc -o {output.reports} -n 1.Report_FastQC_Raw.html -d 1.QC.RAW")
        shell("multiqc -o {output.reports} -n 2.Report_Trimming.html -d 2.TRIMMED")
        shell("multiqc -o {output.reports} -n 3.Report_FastQC_Trimmed.html -d 3.QC.TRIMMED")
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
    print("Using Paired End Trimming")
    rule trim_reads:
        input:
            # reads   = READS + "/{raw_reads}{raw_ends}." + EXTENSION,
            r1      = READS + "/{raw_reads}" + ENDS[0] + "." + EXTENSION,
            r2      = READS + "/{raw_reads}" + ENDS[1] + "." + EXTENSION,
            adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters")
        output:
            # "2.TRIMMED/{raw_reads}{raw_ends}." + EXTENSION,
            # "2.TRIMMED/{raw_reads}{raw_ends}_un." + EXTENSION,
            "2.TRIMMED/{raw_reads}" + ENDS[0] + "." + EXTENSION,
            "2.TRIMMED/{raw_reads}" + ENDS[0] + "_un." + EXTENSION,
            "2.TRIMMED/{raw_reads}" + ENDS[1] + "." + EXTENSION,
            "2.TRIMMED/{raw_reads}" + ENDS[1] + "_un." + EXTENSION
        log:
            #"2.TRIMMED/{raw_reads}{raw_ends}.log"
            "2.TRIMMED/{raw_reads}.log"
        threads:
            CPUS_TRIMMING
        shell:
            #"trimmomatic SE -threads {threads} {input.reads} {output} ILLUMINACLIP:{input.adapter}/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:20 TRAILING:3 MINLEN:36 2> {log}"
            "trimmomatic PE -threads {threads} {input.r1} {input.r2} {output} ILLUMINACLIP:{input.adapter}/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:20 TRAILING:3 MINLEN:36 2> {log}"
            #"trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} ILLUMINACLIP:{input.adapter}/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:20 TRAILING:3 MINLEN:36 2> {log}"

else:
    print("Using Single End Trimming")
    rule trim_reads:
        input:
            adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters"),
            reads   = READS + "/{raw_reads}." + EXTENSION
            #reads = rules.reads.input.reads
        output:
            "2.TRIMMED/{raw_reads}." + EXTENSION
        log:
            "2.TRIMMED/{raw_reads}.log"
        threads:
            CPUS_TRIMMING
        shell:
            "trimmomatic SE -threads {threads} {input.reads} {output} ILLUMINACLIP:{input.adapter}/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:20 TRAILING:3 MINLEN:36 2> {log}"

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
        run:
            shell("fastqc -o 3.QC.TRIMMED -t {threads} {input} 2> {log}")
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
        run:
            shell("fastqc -o 3.QC.TRIMMED -t {threads} {input} 2> {log}")
