# D<sub>o</sub>GE: Differential Gene Expression Analysis Pipeline
![doge](images/logo.jpg)
 
## BIO792: Next Generation Sequencing Data Analysis
![Rule Graph](images/rule-graph.png?raw=true "Rule Graph")

## Requirements
-	Aria2c 1.34.0+
-	FastQC 0.11.8+
-	FeatureCounts 1.6.4+
-	MultiQC 1.7+
-	Python 3.6+ (using [Ana](https://anaconda.org)([mini](https://docs.conda.io/en/latest/miniconda.html))conda)
-	Samtools 1.9+
-	Snakemake 5.7+
-	SRA Toolkit 2.9.6+

## Setup
```bash
git clone https://github.com/villegar/doge
cd doge
conda env create -f environment.yml
conda activate DoGE or source activate DoGE
python download.genome.py human-genome.json
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
bash run_cluster &> log &
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
    "genome4phiX":
    {
        "PhiX": 
            "ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz"
    },

    "genome4star":
    {
        "GRCm38.primary_assembly.genome.fa.gz": 
            "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/GRCm38.primary_assembly.genome.fa.gz",
        "gencode.vM22.annotation.gtf.gz": 
            "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.annotation.gtf.gz"
    },

    "_comment" :
        "Below READS section shows the configuration for a directory containing reads in the format:
	/{PATH}/SRR{LIBRARY}_1.fastq",
    "reads":
    {
        "extension": "fastq.gz",
        "path": "/gpfs/scratch/Classes/bio792/reads",
        "prefix": "SRR"
    }
}
```
