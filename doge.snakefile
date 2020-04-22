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
    include: "PE/all.pe.snakefile"
    include: "PE/fastqc.pe.snakefile"
    include: "PE/trimmomatic.pe.snakefile"
    include: "PE/trimmomatic.fastqc.pe.snakefile"
    include: "PE/hisat2.alignment.pe.snakefile"
else:
    include: "SE/all.se.snakefile"
    include: "SE/fastqc.se.snakefile"
    include: "SE/trimmomatic.se.snakefile"
    include: "SE/trimmomatic.fastqc.se.snakefile"
    include: "SE/hisat2.alignment.se.snakefile"

include: "rules/hisat2.index.snakefile"
include: "rules/sam2bam.snakefile"
include: "rules/alignment.quality.snakefile"
include: "rules/quantification.table.snakefile"
include: "rules/annotation.table.snakefile"
include: "rules/rmd.report.snakefile"
