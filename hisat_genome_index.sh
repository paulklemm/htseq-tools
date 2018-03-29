#!/bin/bash

# Create Splicesites
# Docker image Bioinf
# hisat2_extract_splice_sites.py Homo_sapiens.GRCh38.91.gtf > splicesites.tsv
# hisat2_extract_exons.py Homo_sapiens.GRCh38.91.gtf > exons.tsv
# hisat2-build -p 32 --ss splicesites.tsv --exon exons.tsv Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly

# Make Hisat Index
ENSEMBL_BUILD=91
SPECIES=Homo_sapiens
GENOME=GRCh38
OUTPUT=$(readlink -f ../genome)
logs=$(readlink -f ../slurm_logs)/

SHIFTER="/beegfs/bin/shifter/latest/bin/shifter --image=mpgagebioinformatics/bioinformatics_software:v1.1.3 bash"

mkdir $OUTPUT
mkdir $logs

cd $OUTPUT

wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_BUILD}/fasta/${SPECIES,,}/dna/${SPECIES}.${GENOME}.dna.primary_assembly.fa.gz
gunzip ${SPECIES}.${GENOME}.dna.primary_assembly.fa.gz
# Download transcriptome
wget ftp://ftp.ensembl.org/pub/release-91/gtf/${SPECIES,,}/${SPECIES}.${GENOME}.${ENSEMBL_BUILD}.gtf.gz
gunzip ${SPECIES}.${GENOME}.${ENSEMBL_BUILD}.gtf.gz

sbatch << EOF
#!/bin/bash
#SBATCH --output ${logs}hisat_generate_genome_${ENSEMBL_BUILD}_${SPECIES}.%j.out
#SBATCH --partition=blade-b
#SBATCH -c 32
#SBATCH --job-name='hstgen'

${SHIFTER} << SHI
#!/bin/bash
source /beegfs/scratch/bruening_scratch/pklemm/shifter/home/.bashrc
module load hisat
cd ${OUTPUT}
hisat2_extract_splice_sites.py ${SPECIES}.${GENOME}.${ENSEMBL_BUILD}.gtf > splicesites.tsv
hisat2_extract_exons.py ${SPECIES}.${GENOME}.${ENSEMBL_BUILD}.gtf > exons.tsv
hisat2-build -p 32 --ss splicesites.tsv --exon exons.tsv ${SPECIES}.${GENOME}.dna.primary_assembly.fa ${SPECIES}.${GENOME}.dna.primary_assembly
SHI
EOF
