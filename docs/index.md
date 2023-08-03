## Metagenomics analysis procedure for APDC project

### Required tools

## **Required tools**

For this analysis, we use the following tools. Most of the tools can be installed from `conda` via `bioconda` channel. As such, we recommend creating a `conda` environment, installing all the tools in that particular environment and adding the environment `bin` to the system `PATH`. Below is a list of bioinformatics tools required for this analysis.

* `fastqc` for checking the quality of raw read data and `multiqc` is used to aggregate `fastqc` results for multiple samples into a single file for easier inspection.
* `trim_galore` for trimming off adapters and low quality reads.
* `kraken2` for taxonomic classification of raw reads.
* `bowtie2` and `tanoti` for mapping reads onto reference genomes.
* `tablet` for visualising  alignment files.
* `samtools` for manipulating and inspecting sequence alignment files.
* `metaspades` for de novo assembly of the reads.
* `diamond` for comparing contigs reconstructed from denovo assemby to the NCBI non-redundant protein database. 
* `mafft` for multiple sequence alignment.
* `IQ-tree` and `iTOL` for constructing and visualising phylogenetic trees respectively.

## **Required Databases**

Pre-built databases: [https://benlangmead.github.io/aws-indexes/bowtie](https://benlangmead.github.io/aws-indexes/)
* Bowtie2 index of the Human genome: [https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip](https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip)
* Kraken2 viral DB: [https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230605.tar.gz](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230605.tar.gz)

For `DIAMOND` and `BLASTn`, we build on-prem databases using the RVDB database [https://rvdb.dbi.udel.edu/](https://rvdb.dbi.udel.edu/).
* Unclustered nucleotide sequences: [https://rvdb.dbi.udel.edu/download/U-RVDBv26.0.fasta.gz](https://rvdb.dbi.udel.edu/download/U-RVDBv26.0.fasta.gz)
* Proteic version of RVDB: [https://rvdb-prot.pasteur.fr/files/U-RVDBv26.0-prot.fasta.xz](https://rvdb-prot.pasteur.fr/files/U-RVDBv26.0-prot.fasta.xz)

## **Configuration**

The analysis procedure presented in this document assumes a computing infrastructure that uses `SLURM` for Job scheduling and management, in particular the Uganda Medical Informatics Centre. For clarity, a `slurm` script is created for each step of the analysis. For a particular job (or set of jobs), a user has upto a maximum of 4 nodes and 32 CPUs, this is reflected in the scripts below, so feel free to change according to resources that maybe available to you.

## **Quality control and trimming**

Here we assess the quality of the raw reads using `Fastqc` program. This gives the basic statistics for each of the forward and backward reads.

```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
#!/bin/bash
#SBATCH --partition=all_2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=FastQC
#SBATCH --mail-type=ALL
#SBATCH --mail-user=assekagiri@gmail.com
#SBATCH --output=slurm_%j.out

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

##  Run FastQC
mkdir fastqc_results
fastqc /mnt/lustre01/projects/viral_discovery/users/alfred/data/Metagenomics/Fastq/*.fastq -o fastqc_results
```

After inspecting the quality of the reads, we use `trim_galore` to trim reads by length and quality. Note that our input is the fastq reads as specified earlier.

```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
#!/bin/bash
#SBATCH --partition=all_2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=FastQC
#SBATCH --mail-type=ALL
#SBATCH --mail-user=assekagiri@gmail.com
#SBATCH --output=slurm_%j.out

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

mkdir trimmed
trim_galore -q 30 --dont_gzip --paired data/sample1_R1.fq data/sample1_R2.fq -o trimmed
```

## **Taxonomic classification of short reads**

```
#!/bin/bash
#SBATCH --partition=all_2
#SBATCH --nodes=4
#SBATCH --cpus-per-task=32
#SBATCH --job-name=kraken2-viral
#SBATCH --mail-type=ALL
#SBATCH --mail-user=assekagiri@gmail.com
#SBATCH --output=slurm_%j.out

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES

##  Run kraken2 on just the viral DB
for f in $(ls /mnt/lustre01/projects/viral_discovery/users/alfred/data/Metagenomics/Fastq/*R1*.fastq); do echo $f;    
            sample=$(basename $f '_L001_R1_001.fastq'); echo $sample; \
            kraken2 --paired --report ${sample}.report \
            --output ${sample}.txt ${sample}_L001_R1_001.fastq \
            ${sample}_L001_R2_001.fastq \
            --db /mnt/lustre01/projects/viral_discovery/users/alfred/databases/krakenDB/standard; 
done
```

## **Reference based mapping**

At this point, we have clean reads and we are ready to map the reads onto reference genomes of interest using `bowtie2`.

```
#!/bin/bash
#SBATCH --partition=all_2
#SBATCH --nodes=4
#SBATCH --cpus-per-task=32
#SBATCH --job-name=Mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=assekagiri@gmail.com
#SBATCH --output=slurm_%j.out

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES

##  Align host-free reads to a reference sequence

An example for SARS-CoV-2 as the reference genome

data_path="/mnt/lustre01/projects/viral_discovery/users/alfred/data/Fastq"
ref_path="/mnt/lustre01/projects/viral_discovery/users/alfred/analysis/mapping/bowtie2-ref/sars"

for f in $(ls $data_path/*_L001_R1_001.fastq); do
  fname=$(basename $f "_L001_R1_001.fastq");
  bowtie2 -1 $data_path/${fname}_L001_R1_001.fastq\
          -2 $data_path/${fname}_L001_R2_001.fastq \
          -x $ref_path -S ${fname}.sam
  samtools view -b ${fname}.sam > ${fname}.bam
  samtools sort ${fname}.bam > ${fname}_sorted.bam
  samtools index ${fname}_sorted.bam
done
```

## **De novo assembly**
This section comprises of three main stages; removing host reads from the read data, assembling short reads into longer contiguous sequences using two different assemblers,
merging contigs constructed by two assemblers and lastly blasting the merged contigs using
`diamond`.

### **Remove host reads**
Prior to performing de novo assembly, we remove the host reads from the 
data. Here we use `bowtie2` to map reads to the human genome and retain  the reads that do not map to the human genome for de novo assembly.

```
#!/bin/bash
#SBATCH --partition=all_2
#SBATCH --nodes=4
#SBATCH --cpus-per-task=32
#SBATCH --job-name=Mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=assekagiri@gmail.com
#SBATCH --output=slurm_%j.out

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES

##  Align reads to the Human genome

data_path="/mnt/lustre01/projects/viral_discovery/users/alfred/data/Fastq"
ref_path="/mnt/lustre01/projects/viral_discovery/users/databases/humangenome/humangenome"

for f in $(ls $data_path/*_L001_R1_001.fastq); do
  fname=$(basename $f "_L001_R1_001.fastq");
  bowtie2 -1 $data_path/${fname}_L001_R1_001.fastq\
          -2 $data_path/${fname}_L001_R2_001.fastq \
          -x $ref_path -S ${fname}.sam \
          --un-conc-gz ${fname}_clean > ${fname}.sam;
  mv ${fname}_clean.1 ${fname}_clean.1.gz
  mv ${fname}_clean.2 ${fname}_clean.2.gz
  gunzip ${fname}_clean.1.gz ${fname}_clean.2.gz
done
```

### **De novo assembly using spades**

Assemble the short reads into longer contigs using `spades`.

```
#!/bin/bash
#SBATCH --partition=all_2
#SBATCH --nodes=4
#SBATCH --cpus-per-task=32
#SBATCH --job-name=MetasPAdes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=assekagiri@gmail.com
#SBATCH --output=slurm_%j.out

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES

##  Run  metaspades
mkdir spades_output
data_path="/mnt/lustre01/projects/viral_discovery/users/alfred/analysis/hostfree"
for f in $(ls $data_path/*1.fastq); \
 do echo $f; sample=$(basename $f '_clean.1.fastq'); echo $sample; \
  metaspades.py -1 $data_path/${sample}_clean.1.fastq \
                -2 $data_path/${sample}_clean.2.fastq \
                -o ${sample}_spades_output; 
done
```

### **Blast de-novo constructed contigs using diamond**

```
mkdir diamond_out
#!/bin/bash
#SBATCH --partition=all_2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=DIAMOND
#SBATCH --mail-type=ALL
#SBATCH --mail-user=assekagiri@gmail.com
#SBATCH --output=slurm_%j.out

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES

##  Run  DIAMOND blastx
data_path="/mnt/lustre01/projects/viral_discovery/users/alfred/analysis/hostfree/scaffolds"
db_path="/mnt/lustre01/projects/viral_discovery/users/alfred/databases/dmnd/viral"
for f in $(ls $data_path/*.fasta); \
 do echo $f; sample=$(basename $f 'scaffold.fasta'); echo $sample; \
  diamond blastx --db $db_path --query $f --out ${sample}_dmnd.tsv --outfmt 6 
done
```

### **Classify assembled contigs using Kraken2**

```
#!/bin/bash
#SBATCH --partition=all_2
#SBATCH --nodes=4
#SBATCH --cpus-per-task=32
#SBATCH --job-name=kraken2-viral
#SBATCH --mail-type=ALL
#SBATCH --mail-user=assekagiri@gmail.com
#SBATCH --output=slurm_%j.out

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES

##  Run kraken2 on just the viral DB

for f in $(ls /mnt/lustre01/projects/viral_discovery/users/alfred/data/Metagenomics/assembly/*.fasta); do echo $f;    
            sample=$(basename $f '.fasta'; echo $sample; \
            kraken2 --report ${sample}.report \
            --output ${sample}.txt  \
            ${sample}.fasta \
            --db /mnt/lustre01/projects/viral_discovery/users/alfred/databases/krakenDB; 
done
```

## **Phylogenetic analyses**

Multiple sequence alignment using `mafft`.
```
mafft all_contigs.fa > all_contigs_aln.fa
```

Phylogenetic tree construction using `iqtree`
```
iqtree -s all_contigs_aln.fa -m MF
```

Uploaded the tree to [iTol](https://itol.embl.de) for visualisation and further inspection. 

## **Contact**
Principal Investigator, APDC project: Dr. Nicholas Bbosa at [nicholas.bbosa@mrcuganda.org](nicholas.bbosa@mrcuganda.org)
