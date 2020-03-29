####### Libraries #######
import glob
import os

####### Util functions #######
def extractFilenames(fullnames,suffix):
        names = []
        for file in fullnames:
            names.append(os.path.basename(file).split(suffix)[0])
        return sorted(names)

def findLibraries(path,prefix,suffix):
	filenames_path = glob.glob(os.path.join(path,prefix) + "*" + suffix)
	names = []
	for file in filenames_path:
	    library = os.path.basename(file).split(suffix)[0]
	    if(library not in names):
		    names.append(library)
	return sorted(names)

def which(file):
        for path in os.environ["PATH"].split(os.pathsep):
                if os.path.exists(os.path.join(path, file)):
                        return path
        return None

####### Global variables #######
EXTENSION = config["reads"]["extension"]
PREFIX = config["reads"]["prefix"]
READS = config["reads"]["path"]
FORWARD_READ_ID = config["reads"]["forward_read_id"]
SUFFIX = "_" + FORWARD_READ_ID + "." + EXTENSION
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
GENOME4STAR = config["genome4star"]
#GENOME4STAR_FILENAMES = extractFilenames(GENOME4STAR.keys(),".gz")
GENOME4PHIX = config["genome4phiX"]
#KRAKEN_DB = config["krakenDB"]
#KRAKEN_DB_FILENAMES = extractFilenames(KRAKEN_DB.keys(),".tgz")
#rRNA = config["rRNAref"]
#rRNA_FILES = list(rRNA.keys())

####### Rules #######
rule all:
	input:
		expand("1.QC.RAW/{library}_{end}_fastqc.{format}", library=LIBS, end=[1, 2], format=["html","zip"]),
        expand("3.QC.TRIMMED/{library}_{direction}_{mode}_fastqc.{format}",
			library=LIBS, direction=["forward","reverse"], mode=["paired","unpaired"], format=["html","zip"])
	output:
		logs 	= directory("0.LOGS"),
		reports	= directory("10.MULTIQC")
	run:
		shell("multiqc -o {output.reports} -n 1.Report_FastQC_Raw.html -d 1.QC.RAW")
		shell("mkdir -p {output.logs} && mv *.log {output.logs}")

rule reads:
	input:
		reads = READS + "/{library}_{end}." + EXTENSION,
		r1    = READS + "/{library}_1." + EXTENSION,
		r2    = READS + "/{library}_2." + EXTENSION
	message:
		"Gathering reads"

rule fastqc_raw:
	input:
		reads = rules.reads.input.reads
	output:
		html = "1.QC.RAW/{library}_{end}_fastqc.html",
		zip  = "1.QC.RAW/{library}_{end}_fastqc.zip"
	message:
		"FastQC on raw data"
	log:
		"1.QC.RAW/{library}_{end}.log"
	threads:
		CPUS_FASTQC
	shell:
		"fastqc -o 1.QC.RAW -t {threads} {input} 2> {log}"

rule trim_reads:
	input:
		adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters"),
		r1 = rules.reads.input.r1,
		r2 = rules.reads.input.r2
	output:
		forward_paired   = "2.TRIMMED/{library}_forward_paired.fastq.gz",
		forward_unpaired = "2.TRIMMED/{library}_forward_unpaired.fastq.gz",
		reverse_paired   = "2.TRIMMED/{library}_reverse_paired.fastq.gz",
		reverse_unpaired = "2.TRIMMED/{library}_reverse_unpaired.fastq.gz"
	message:
		"Trimming reads"
	log:
		"2.TRIMMED/{library}.log"
	threads:
		CPUS_TRIMMING
	shell:
		"trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} ILLUMINACLIP:{input.adapter}/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:20 TRAILING:3 MINLEN:36 2> {log}"

rule fastqc_trimmed:
	input:
		rules.trim_reads.output.forward_paired,
		rules.trim_reads.output.forward_unpaired,
		rules.trim_reads.output.reverse_paired,
		rules.trim_reads.output.reverse_unpaired
	output:
		html = "3.QC.TRIMMED/{library}_{direction}_{mode}_fastqc.html",
		zip  = "3.QC.TRIMMED/{library}_{direction}_{mode}_fastqc.zip"
	message:
		"FastQC on trimmed data"
	log:
		"3.QC.TRIMMED/{library}_{direction}_{mode}.log"
	threads:
		CPUS_FASTQC
	shell:
                "fastqc -o 3.QC.TRIMMED -t {threads} {input} 2> {log}"
