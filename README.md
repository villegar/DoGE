# D<sub>o</sub>GE: Differential Gene Expression Analysis Pipeline <img src="images/logo.jpg" alt="doge-log" align="right" height=150px/>
 
## BIO792: Next Generation Sequencing Data Analysis
<p align="center">
	<a href="https://raw.githubusercontent.com/villegar/DoGE/master/images/rule-graph.png">
		<img src=images/rule-graph.png width=500px />
	</a>
</p>


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
	* [Execution](#execution)
	
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
conda env create -f environment.yml -n DoGE
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
-	The `genome` section __MUST__ point to the path for the `X-genome.json` file.
-	The `reads` section points the pipeline to the location (`path`), format (`extension`), type (`end_type`), and prefix (`prefix`) of the raw reads. Optionally, if `end_type = pe` (paired-end), both the forward (`forward_read_id`) and reverse (`reverse_read_id`) reads identifier (e.g. 1, R1, 2, R2, etc.) should be specified.
-	The `trimmomatic` section should contain a sub-key called `options` with the parameters for trimming, excluding the input and ouput names, which will be set up by the pipeline.

```bash
{
    "genome": "/path/to/X-genome.json",
    "reads": {
        "extension": "fastq",
        "end_type": "se",
        "forward_read_id": "1",
        "reverse_read_id": "2",
        "path": "/path/to/raw/reads",
        "prefix": "SRR"
    },
    "trimmomatic":{
      "options": "ILLUMINACLIP:{input.adapter}/TruSeq3-SE-2.fa:2:30:10:2:keepBothReads"
    }
}
```


# Study Case
## Data set
For this study case the following article title [`LncRNA DEANR1 facilitates human endoderm differentiation by activating FOXA2 expression`](https://www.ncbi.nlm.nih.gov/gds/?term=(SRP019241)%20AND%20gds_sra[filter]
) was consulted.
https://doi.org/10.1016/j.celrep.2015.03.008

### Accession numbers
```
SRR1958165
SRR1958166
SRR1958167
SRR1958168
SRR1958169
SRR1958170
```

### Configuration file
```
{
    "genome": "genomes/human-genome.json",
    "reads": {
        "extension": "fastq",
        "end_type": "se",
        "path": "/path/to/reads",
        "prefix": "SRR"
    },
    "trimmomatic":{
      "options": "ILLUMINACLIP:{input.adapter}/TruSeq3-SE-2.fa:2:30:10:2:keepBothReads TRAILING:3 MINLEN:24"
    }
}
```

## Execution
It is a good practice to perform a `dry-run` of the workflow before submitting for execution. This can be done by appending the `-n` option to the `snakemake` command:

```
snakemake --configfile config.json -n
```

The output will display a summary of each job that will be processed and a final summary that should look like:
```
Job counts:
        count   jobs
        6       alignment
        6       alignment_quality
        1       all
        1       annotation_table
        6       fastqc_raw
        6       fastqc_trimmed
        1       feature_counts
        1       hisat2_index
        1       quantification_table
        1       rmd_report
        6       sam2bam
        6       trim_reads
        42
```

For a graphical summary of above jobs, check the directed acyciclic graph: https://raw.githubusercontent.com/villegar/DoGE/master/images/dag.png

### Single node execution
```
snakemake -j CPUS \ # maximum number of CPUs available to Snakemake
	  --configfile config.json # configuration file
```
