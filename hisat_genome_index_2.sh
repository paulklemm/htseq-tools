#!/bin/bash

# Example usage: 
# ./hisat_genome_index.sh --ensemblbuild 91 --species "Homo_sapiens" --genome GRCh38 --output "../genome" --logs "../slurm_logs"
# ./hisat_genome_index.sh --genome ftp://ftp.ensembl.org/pub/release-91/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz
# ./hisat_genome_index_2.sh --genome "ftp://ftp.ensembl.org/pub/release-91/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz" --output "../genome" --logs "../slurm_logs"

# Parameter read-in from https://stackoverflow.com/a/14203146/2274058
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -g|--genome)
    GENOME_URL="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--logs)
    LOGS="$2"
    shift # past argument
    shift # past value
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo "*** Run Make HISAT Genome Index ***"

OUTPUT=$(readlink -f ${OUTPUT})
LOGS=$(readlink -f ${LOGS})/

mkdir $OUTPUT
mkdir $LOGS

ENSEMBL_BUILD=$(echo "${GENOME_URL}" | sed 's/.*release-\([^/]*\).*/\1/')
SPECIES=$(echo "${GENOME_URL}" | sed 's/.*dna\/\([^\.]*\).*/\1/')
GENOME=$(echo "${GENOME_URL}" | sed 's/.*dna\/[^\.]*\.\([^\.]*\).*/\1/')
GENOME_COMPLETE_NAME=$(echo "${GENOME_URL}" | sed 's/.*\///')
GENOME_COMPLETE_NAME_WITHOUT_EXTENSION=$(echo "${GENOME_COMPLETE_NAME}" | sed 's/\.fa.gz//')

echo ENSEMBL BUILD                          = "${ENSEMBL_BUILD}"
echo SPECIES                                = "${SPECIES}"
echo GENOME                                 = "${GENOME}"
echo OUTPUT                                 = "${OUTPUT}"
echo LOGS                                   = "${LOGS}"
echo GENOME_COMPLETE_NAME                   = "${GENOME_COMPLETE_NAME}"
echo GENOME_COMPLETE_NAME_WITHOUT_EXTENSION = "${GENOME_COMPLETE_NAME_WITHOUT_EXTENSION}"

# Create Splicesites
# Docker image Bioinf
# hisat2_extract_splice_sites.py Homo_sapiens.GRCh38.91.gtf > splicesites.tsv
# hisat2_extract_exons.py Homo_sapiens.GRCh38.91.gtf > exons.tsv
# hisat2-build -p 32 --ss splicesites.tsv --exon exons.tsv Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly

# Make Hisat Index
# ENSEMBL_BUILD=91
# SPECIES=Homo_sapiens
# GENOME=GRCh38

SHIFTER="/beegfs/bin/shifter/latest/bin/shifter --image=mpgagebioinformatics/bioinformatics_software:v1.1.3 bash"

sbatch << EOF
#!/bin/bash
#SBATCH --output ${LOGS}hisat_generate_genome_${ENSEMBL_BUILD}_${SPECIES}.%j.out
#SBATCH --partition=blade-b
#SBATCH -c 32
#SBATCH --job-name='hstgen'

${SHIFTER} << SHI
#!/bin/bash
source /beegfs/scratch/bruening_scratch/pklemm/shifter/home/.bashrc
module load hisat
cd ${OUTPUT}
# Download FA and GTFs
# wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_BUILD}/fasta/${SPECIES,,}/dna/${SPECIES}.${GENOME}.dna.primary_assembly.fa.gz
# gunzip ${SPECIES}.${GENOME}.dna.primary_assembly.fa.gz
wget ${GENOME_URL}
gunzip ${GENOME_COMPLETE_NAME}
# Download transcriptome
wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_BUILD}/gtf/${SPECIES,,}/${SPECIES}.${GENOME}.${ENSEMBL_BUILD}.gtf.gz
gunzip ${SPECIES}.${GENOME}.${ENSEMBL_BUILD}.gtf.gz

# Make index
hisat2_extract_splice_sites.py ${SPECIES}.${GENOME}.${ENSEMBL_BUILD}.gtf > splicesites.tsv
hisat2_extract_exons.py ${SPECIES}.${GENOME}.${ENSEMBL_BUILD}.gtf > exons.tsv
# hisat2-build -p 32 --ss splicesites.tsv --exon exons.tsv ${SPECIES}.${GENOME}.dna.primary_assembly.fa ${SPECIES}.${GENOME}.dna.primary_assembly
hisat2-build -p 32 --ss splicesites.tsv --exon exons.tsv $(echo "${GENOME_COMPLETE_NAME_WITHOUT_EXTENSION}.fa") ${GENOME_COMPLETE_NAME_WITHOUT_EXTENSION}
SHI
EOF
