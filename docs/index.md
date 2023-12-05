## Metagenomics analysis procedure for APDC project

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

Pre-built databases: Mantained by [BenLangmead at John Hopkins University](https://langmead-lab.org/).
* `Bowtie2` index of the Human genome: [https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip](https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip)
* `Kraken2` viral DB: [https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230605.tar.gz](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230605.tar.gz)

For `DIAMOND` and `BLASTn`, we build custom databases using the  [RVDB database](https://rvdb.dbi.udel.edu/). This resource is mantained by the Arifa Khan's group at CBER, FDA for enhancing virus detection using NGS.

* Unclustered nucleotide sequences: [https://rvdb.dbi.udel.edu/download/U-RVDBv26.0.fasta.gz](https://rvdb.dbi.udel.edu/download/U-RVDBv26.0.fasta.gz)
* Proteic version of RVDB: [https://rvdb-prot.pasteur.fr/files/U-RVDBv26.0-prot.fasta.xz](https://rvdb-prot.pasteur.fr/files/U-RVDBv26.0-prot.fasta.xz)

## **Configuration**

The analysis procedure presented in this document assumes a computing infrastructure that uses `SLURM` for Job scheduling and management, in particular the Uganda Medical Informatics Centre. For clarity, a `slurm` script is created for each step of the analysis. For a particular job (or set of jobs), a user has upto a maximum of 4 nodes and 32 CPUs, this is reflected in the scripts below, so feel free to change according to resources that maybe available to you.

## **Quality control and trimming**

Here we assess the quality of the raw reads using `Fastqc` and `MultiQC`. This gives the basic statistics for each of the forward and backward reads. 

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
data_path='/mnt/lustre01/projects/viral_discovery/users/alfred/data/Metagenomics'

mkdir fastqc_results
fastqc $data_path/Fastq/*.fastq -o fastqc_results
multiqc fastqc_results/* -o fastqc_results 
```

At this point, we use `scp` to download the MultiQC report and have a look at it, according to our assessement of the report, we can choose whether or not to do some adaptor/quality trimming. After inspecting the quality of the reads, we use `trim_galore` for adaptor and quality trimming. Choice of cut-offs on quality scores, length, e.t.c is guided by the assessment made on the QC plots generated above.

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

data_path="/mnt/lustre01/projects/viral_discovery/users/alfred/data/Metagenomics"
phred_score=30

mkdir trimmed
for f in $(ls $data_path/Fastq/*R1*.fastq); do echo $f;
sample=$(basename $f '_L001_R1_001.fastq'); echo $sample;  
trim_galore -q $phred_score --dont_gzip --paired $data_path/${sample}_L001_R1_001.fastq $data_path/${sample}_L001_R2_001.fastq -o trimmed
done
```

## **Read-based taxonomic identification using Kraken2**

To get an idea of the pathogens pathogens that could be in this sample, we do screen the read data using kraken2. Basically, the taxonomic identification is done in comparison/reference to a kraken-compatible pre-built database. We use two databases, the standard one (for viral, bacterial and Human reads classification) and the viral database for only virus classification. Change the `db_path` in the script below to chose a desired database.

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

##  Run kraken2 on the standard Kraken2 DB
data_path="/mnt/lustre01/projects/viral_discovery/users/alfred/data/Metagenomics"
db_path="/mnt/lustre01/projects/viral_discovery/users/alfred/databases/krakenDB/standard"
#db_path="/mnt/lustre01/projects/viral_discovery/users/alfred/databases/krakenDB/viral"

for f in $(ls $data_path/Fastq/*R1*.fastq); do echo $f;    
            sample=$(basename $f '_L001_R1_001.fastq'); echo $sample; \
            kraken2 --paired --report ${sample}.report \
            --output ${sample}.txt ${sample}_L001_R1_001.fastq \
            ${sample}_L001_R2_001.fastq \
            --db $db_path; 
done
```

Process Kraken output using [KrakenTools](https://github.com/jenniferlu717/KrakenTools) and visualise results using KronaTools.

```
# convert to mpa format
for f in `ls *.report`; do b=$(basename $f '.report'); echo $b;
python /mnt/lustre01/projects/viral_discovery/users/alfred/analysis/scripts/KrakenTools/kreport2mpa.py -r $f -o ${b}.report.mpa  --display-header;
done
python /mnt/lustre01/projects/viral_discovery/users/alfred/analysis/scripts/KrakenTools/combine_mpa.py -i *.mpa -o kraken-combined.txt

# generate krona plots
for f in `ls *.report`; do b=$(basename $f '.report'); echo $b;
python /mnt/lustre01/projects/viral_discovery/users/alfred/analysis/scripts/KrakenTools/kreport2krona.py -r $f -o ${b}.krona ;
done
for f in `ls *.krona`; do b=$(basename $f '.krona'); echo $b; ktImportText $f -o ${b}.krona.html ; done
```

**Extracting reads classified as a certain taxon**

- We need the taxon ID, in the example below we use taxon ID 1392.
- The FASTQ files from whoch to extract classified reads
- The `Kraken` TXT file   

```
data_path='path/to/FASTQ/files'
kraken_res='path/to/Kraken/results'

for f in `ls $kraken_res/*.txt`; do b=$(basename $f '.txt'); echo $b;
python extract_kraken_reads.py -k ${kraken_res}/${b}.txt -1 ${data_path}/${b}_L001_R1_001.fastq -2 ${data_path}/${b}_L001_R2_001.fastq -t 1392 -o ${b}_1.fasta -o2 ${b}_2.fasta ;
cat ${b}*.fasta > ${b}.fa
done
```

## **Mapping mNGS data onto reference genomes**

After inspecting the Kraken2 output, we hand pick taxa of interest e.g viruses, bacteria e.t.c. We then download corresponding reference genomes (we usually download reference sequences of the genus to which a particular pathogen of interest belongs). We use two tools for this purpose; `bowtie2` and `tanoti`. 

#### Mapping mNGS data onto reference genomes using Bowtie2 

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


#### Mapping mNGS data onto reference genomes using Tanoti 

Map the short reads onto the reference genome, obtain mapping statistics and generate consensus sequences from alignment maps

```
data_path="/mnt/lustre01/projects/viral_discovery/users/alfred/data/Fastq"
ref_path="/mnt/lustre01/projects/viral_discovery/users/alfred/analysis/mapping/bowtie2-ref/refs"

for f in $(ls $data_path/*_L001_R1_001.fastq); do
            fname=$(basename $f "_L001_R1_001.fastq");
            mkdir $fname;
            tanoti -r $ref_path/ref.fasta -i $data_path/${fname}_L001_R1_001.fastq $data_path/${fname}_L001_R2_001.fastq -p 1;
            mv FinalAssembly.sam $fname/ ;
            SAM_STATS  $fname/FinalAssembly.sam
            #generate consensus
            SAM2CONSENSUS -i $fname/FinalAssembly.sam -o $fname/${fname}_consensus.fa
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

### **Blast de-novo constructed contigs using BLASTn**

We compare contigs to a database of reference nucleotide sequences from Refseq. This can be adjusted using the `db_path` variable in the script below. The reference sequences can be the entire NCBI `nt`, sequences of particular `genus` or `family` of viruses e.t.c, as long as it is in `.fasta` format.

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
db_path="/mnt/lustre01/projects/viral_discovery/users/alfred/databases/blastn/virus.fa"
for f in $(ls $data_path/*.fasta); \
 do echo $f; sample=$(basename $f 'scaffold.fasta'); echo $sample; \
  blastn -query $f -subject $db_path -outfmt 6 > ${sample}_blastn.tsv
done

# Add sample names to blastn output
for f in $(ls *_blastn.tsv); do bn=$(basename $f 'blastn.tsv'); sed "s/$/\t$bn/g" $f > ${bn}_blastn_labelled.tsv; done
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
  diamond blastx --db $db_path --query $f --out ${sample}_dmnd.tsv --outfmt 6 > ${sample}_dmnd.tsv
done

# Add sample names to dmnd output
for f in $(ls *_dmnd.tsv); do bn=$(basename $f '_scaffolds.fasta_dmnd.tsv'); sed "s/$/\t$bn/g" $f > ${bn}_blastx.tsv; done
```

## **Phylogenetic analyses**

After analysing the Blastn and Blastx results, identify representative reference genomes from the Family/Genus of the virus of interest - refer to NCBI RefSeq and/or [ICTV](https://ictv.global/). Construct a base multiple sequence alignment for that particular family/genus. Add sequences (contigs) and construct a phylogeny using IQtree/RAxML.

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
