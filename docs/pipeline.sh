# quality control
mkdir fastqc_results
fastqc /mnt/lustre01/projects/viral_discovery/users/alfred/data/Metagenomics/Fastq/*.fastq -o fastqc_results

# adaptor and quality trimming
mkdir trimmed
trim_galore -q 30 --dont_gzip --paired data/sample1_R1.fq data/sample1_R2.fq -o trimmed

# Taxonomic classification using kraken
for f in $(ls /mnt/lustre01/projects/viral_discovery/users/alfred/data/Metagenomics/Fastq/*R1*.fastq); do echo $f;    
            sample=$(basename $f '_L001_R1_001.fastq'); echo $sample; \
            kraken2 --paired --report ${sample}.report \
            --output ${sample}.txt ${sample}_L001_R1_001.fastq \
            ${sample}_L001_R2_001.fastq \
            --db /mnt/lustre01/projects/viral_discovery/users/alfred/databases/krakenDB/standard; 
done

# Mapping reads onto reference genomes
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

# Remove host reads by aligning reads onto the Human genome and retain the unmapped reads
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


# De novo assembly of host depleted reads
mkdir spades_output
data_path="/mnt/lustre01/projects/viral_discovery/users/alfred/analysis/hostfree"
for f in $(ls $data_path/*1.fastq); \
 do echo $f; sample=$(basename $f '_clean.1.fastq'); echo $sample; \
  metaspades.py -1 $data_path/${sample}_clean.1.fastq \
                -2 $data_path/${sample}_clean.2.fastq \
                -o ${sample}_spades_output; 
done

# Blast de-novo constructed contigs using diamond
mkdir diamond
data_path="/mnt/lustre01/projects/viral_discovery/users/alfred/analysis/hostfree/scaffolds"
db_path="/mnt/lustre01/projects/viral_discovery/users/alfred/databases/dmnd/viral"
for f in $(ls $data_path/*.fasta); \
 do echo $f; sample=$(basename $f 'scaffold.fasta'); echo $sample; \
  diamond blastx --db $db_path --query $f --out ${sample}_dmnd.tsv --outfmt 6 
done

# Phylogenetic analysis
mafft all_contigs.fa > all_contigs_aln.fa
iqtree -s all_contigs_aln.fa -m MF
