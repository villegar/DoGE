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


####### Reference datasets #######
GENOME4STAR = config["genome4star"]
GENOME4STAR_FILENAMES = extractFilenames(GENOME4STAR.keys(),".gz")
GENOME4PHIX = config["genome4phiX"]
KRAKEN_DB = config["krakenDB"]
KRAKEN_DB_FILENAMES = extractFilenames(KRAKEN_DB.keys(),".tgz")
rRNA = config["rRNAref"]
rRNA_FILES = list(rRNA.keys())

####### Rules #######
rule all:
	input:
		expand("1.QC.RAW/{library}_{end}_fastqc.{format}", library=LIBS, end=[1, 2], format=["html","zip"])

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
