# D<sub>o</sub>GE: Differential Gene Expression Analysis Pipeline
![doge](images/logo.jpg)
 
## BIO792: Next Generation Sequencing Data Analysis
![Rule Graph](images/rule-graph.png?raw=true "Rule Graph")

Table of Contents
=================

* [Requirements](#requirements)
	* [R libraries:](#r-libraries)
	* [Setup](#setup)
	* [Genome file](#genome-file)
* [Execution](#execution)
	* [Single node](#single-node)
	* [Multi-node](#multi-node)
	* [Alternatively](#alternatively)
	* [Cluster configuration (cluster.json)](#cluster-configuration-clusterjson)
	* [Pipeline configuration (config.json)](#pipeline-configuration-configjson)
* [Study Case](#study-case)
	* [Data set](#data-set)
	
## Requirements
-	FastQC 0.11.9+
-	FeatureCounts 2.0.0+
-	MultiQC 1.7+
-	Python 3.6+ (using [Ana](https://anaconda.org)([mini](https://docs.conda.io/en/latest/miniconda.html))conda)
-	Samtools 1.3+
-	Snakemake 5.7+
-	SRA Toolkit 2.9.6+
-	Trimmomatic 0.39+

### R libraries:
-	calibrate [http://cran.r-project.org/package=calibrate]
-	DESeq2 [https://bioconductor.org/packages/3.10/bioc/html/DESeq2.html]
-	dplyr [http://cran.r-project.org/package=dplyr]
-	GGally [http://cran.r-project.org/package=GGally]
-	ggplot2 [http://cran.r-project.org/package=ggplot2]
-	ggrepel [http://cran.r-project.org/package=ggrepel]
-	gplots [http://cran.r-project.org/package=gplots]
-	gridExtra [http://cran.r-project.org/package=gridExtra]
-	kableExtra [http://cran.r-project.org/package=kableExtra]
- 	knitr [http://cran.r-project.org/package=knitr]
-	latex2exp [http://cran.r-project.org/package=latex2exp]
-	pander [http://cran.r-project.org/package=pander]
-	RColorBrewer [http://cran.r-project.org/package=RColorBrewer]
-	rmarkdown [http://cran.r-project.org/package=rmarkdown]

## Setup
```bash
git clone https://github.com/villegar/doge
cd doge
conda env create -f environment.yml
conda activate DoGE or source activate DoGE
python download.genome.py genomes/X-genome.json
```

### Genome file
A good place to get some reference genomes and gene annotations is http://uswest.ensembl.org/info/data/ftp/index.html. The reference must be stored in JSON format (see below template),`X-genome.json`
```bash
{
	"X.fa.gz":
            "ftp://ftp.ensembl.org/pub/some/path/to/X.fa.gz",
        "X.gtf.gz":
            "ftp://ftp.ensembl.org/pub/some/path/to/X.gtf.gz"
}
```

## Execution
### Single node
```bash
snakemake -j CPUS \ # maximum number of CPUs available to Snakemake
	  --configfile config.json # configuration file
```

### Multi-node
```bash
snakemake -j JOBS  \ # maximum number of simultaneous jobs to spawn
	  --configfile config.json # configuration file
          --latency-wait 1000 \ # files latency in seconds
          --cluster-config cluster.json \ # cluster configuration file
          --cluster "sbatch --job-name={cluster.name} 
                            --nodes={cluster.nodes} 
                            --ntasks-per-node={cluster.ntasks} 
                            --output={cluster.log} 
                            --partition={cluster.partition} 
                            --time={cluster.time}"
```
#### Alternatively
```bash
bash run_cluster config.json &> log &
```

#### Cluster configuration (cluster.json)
```bash
{
    "__default__" :
    {
        "time" : "1-00:00:00",
        "nodes" : 1,
        "partition" : "compute",
	"ntasks": "{threads}",
	"name": "DoGE-{rule}",
	"log": "DoGE-{rule}-%J.log"
    }
}
```

#### Pipeline configuration (config.json)
```bash
{
    "genome":
    {
        "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz":
            "ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        "Homo_sapiens.GRCh38.99.gtf.gz":
            "ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz"
    },

    "_comment": "Below READS section shows the configuration for a directory containing reads in the format: /{PATH}/SRR{LIBRARY}_1.fastq",
    "reads": {
        "extension": "fastq",
        "end_type": "se",
        "forward_read_id": "1",
        "reverse_read_id": "2",
        "path": "/gpfs/scratch/Classes/bio792/reads",
        "prefix": "SRR"
    },

    "trimmomatic":{
      "options": "ILLUMINACLIP:{input.adapter}/TruSeq3-SE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:20 TRAILING:3 MINLEN:24"
    }
}
```


# Study Case
## Data set
For this study case the following article title [`LncRNA DEANR1 facilitates human endoderm differentiation by activating FOXA2 expression`](https://www.ncbi.nlm.nih.gov/gds/?term=(SRP019241)%20AND%20gds_sra[filter]
) was consulted
