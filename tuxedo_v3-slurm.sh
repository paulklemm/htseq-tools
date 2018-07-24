#!/bin/bash

"This script needs to run form inside the folder scripts in a working project with the following structure:
project/scripts
project/raw_data


Raw data needs to be labeled in the following fashion:

Sample_serial-Folder-Line-Time_point/day_of_life(day_post_Treament)-treament-REPlicate-READ

and with the exact number of characters as in the example bellow:

S_001-F_MAFATS-L____NHa_0-____-REP_1-READ_1.fastq.gz

S_XXX-F_MAFAXX-L_XXXXXXXX-XXXX-REP_X-READ_x.fastq.gz

Please notice that for paired samples, the S_XXX is the same.


Make sure you have edited the last section of this script - cuffdiff - before you execute this script." > /dev/null 2>&1


#############################################################################

# Define series as SE or PE and stranded or unstranded
# TODO: Set labels
SE_unstr=()
SE_str=()
PE_str=("MAFA")
PE_uns=()
mix=()

unstr=()
str=("MAFA")
#mix=("Yid3")


# Which series do you which to work on:

series="MAFA"

# Reference genome

# TODO: Set annotation files
#ann=/beegfs/common/genomes/caenorhabditis_elegans/89/
ori_GTF=/beegfs/scratch/bruening_scratch/pklemm/2018-04-Mafalda-RiboSeq/genome/Mus_musculus.GRCm38.92.gtf
hisat_index=/beegfs/scratch/bruening_scratch/pklemm/2018-04-Mafalda-RiboSeq/genome/Mus_musculus.GRCm38.92.dna.toplevel
#adapters_file=/beegfs/group_bit/home/JBoucas/documents/TruSeqAdapters.txt
genome=${hisat_index}.fa

## Custom commands for the dataset
##

SHIFTER="/beegfs/bin/shifter/latest/bin/shifter --image=mpgagebioinformatics/bioinformatics_software:v1.1.3 bash"

#############################################################################


echo "Creating required folders"
mkdir -p ../slurm_logs
mkdir -p ../fastqc_output
mkdir -p ../tmp
mkdir -p ../flexbar_output
mkdir -p ../hisat_output
mkdir -p ../stringtie_output
mkdir -p ../cuffmerge_output
mkdir -p ../cuffdiff_output
mkdir -p ../cuffquant_output
mkdir -p ../cufflinks_output

top=$(readlink -f ../)/
tmp=$(readlink -f ../tmp)/
raw=$(readlink -f ../raw_data)/
rawt=$(readlink -f ../flexbar_output)/
merg=$(readlink -f ../cuffmerge_output)/ 
qua=$(readlink -f ../cuffquant_output)/ 
logs=$(readlink -f ../slurm_logs)/

# Required function
function contains() {
    local n=$#
    local value=${!n}
    for ((i=1;i < $#;i++)) {
        if [ "${!i}" == "${value}" ]; then
            echo "y"
            return 0
        fi
    }
    echo "n"
    return 1
}

#############################################################################

echo "Starting FASTQC"

cd ${raw}
for serie in $series; do
    cd ${raw}
    
    for file in $(ls *${serie}*.fastq.gz); do 
    echo "Starting FastQC for file $file"
sbatch << EOF
#!/bin/bash
#SBATCH --output ${logs}${file%..fastq.gz}_fastqc_.%j.out
#SBATCH --error ${logs}${file%..fastq.gz}_fastqc_.%j.err
#SBATCH --partition=blade-b
#SBATCH -c 4
#SBATCH --job-name='fastqc'

${SHIFTER} << SHI
#!/bin/bash
source /beegfs/scratch/bruening_scratch/pklemm/shifter/home/.bashrc
module load fastqc
cd ${raw}
# FASTQC call
fastqc -t 4 -o ../fastqc_output ${file}
SHI
EOF

    done
done

############################################################################

ids=

cd ${rawt}
for serie in $series; do
    cd ${raw}
    for file in $(ls *${serie}*1.fastq.gz); do
        
        # Libraries and paired end vs. single end settings for HISAT    
 
        if [[ $(contains "${SE_unstr[@]}" "$serie") == "y" ]]; then
            lib=
            files="-U ${file}"
        elif [[ $(contains "${PE_uns[@]}" "$serie") == "y" ]]; then
            lib=
            files="-1 ${file} -2 ${file::(-10)}2.fastq.gz"
        elif [[ $(contains "${SE_str[@]}" "$serie") == "y" ]]; then
            lib="--rna-strandness R"
            files="-U ${file}"
        elif [[ $(contains "${PE_str[@]}" "$serie") == "y" ]]; then
            lib="--rna-strandness RF"
            files="-1 ${file} -2 ${file::(-10)}2.fastq.gz"
        elif [[ $(contains "${mix[@]}" "$serie") == "y" ]]; then
            files=-U ${file}
            REP=${file:30:5}
            if [[ ${REP} == REP_3 ]]; then
                lib="--rna-strandness R"
            else
                lib=
            fi
        fi

rm -rf ${logs}hisat_${file::(-16)}.*.out 

ids=${ids}:$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH --output ${logs}hisat_${file::(-16)}.%j.out
#SBATCH --error ${logs}hisat_${file::(-16)}.%j.err
#SBATCH --partition=blade-b
#SBATCH --cpus-per-task=14 
#SBATCH --job-name='HS_ST'

${SHIFTER} << SHI
#!/bin/bash
source /beegfs/scratch/bruening_scratch/pklemm/shifter/home/.bashrc
cd ${raw}
module load bowtie
module load hisat

# HISAT call 

echo "`date '+%Y-%m-%d_%H-%M-%S'` Run Hisat"
echo "hisat2 -p 14 ${lib} --dta-cufflinks --met-file ${top}hisat_output/${file::(-16)}.stats -x ${hisat_index} -S ${top}hisat_output/${file::(-16)}.sam ${files}"

hisat2 -p 14 ${lib} --dta-cufflinks --met-file ${top}hisat_output/${file::(-16)}.stats \
-x ${hisat_index} -S ${top}hisat_output/${file::(-16)}.sam \
${files}

echo "`date '+%Y-%m-%d_%H-%M-%S'` Hisat done"

cd ${top}hisat_output
module load samtools

# Use samtools to select mapped reads and sort them

echo "`date '+%Y-%m-%d_%H-%M-%S'` Start samtools"
echo "samtools view -@ 14 -bhS -F 4 ${file::(-16)}.sam | samtools sort -@ 14 -o ${file::(-16)}.bam -"

samtools view -@ 14 -bhS -F 4 ${file::(-16)}.sam | samtools sort -@ 14 -o ${file::(-16)}.bam -
# rm -rf ${file::(-16)}.sam
mkdir -p ${top}stringtie_output/${file::(-16)}

echo "`date '+%Y-%m-%d_%H-%M-%S'` samtools done"

module load stringtie

# StringTie call

echo "`date '+%Y-%m-%d_%H-%M-%S'` Start stringtie"
echo "stringtie ${file::(-16)}.bam -o ${top}stringtie_output/${file::(-16)}.gtf -p 14 -G ${ori_GTF} -f 0.99 -C ${top}stringtie_output/${file::(-16)}_full_cov.gtf -b ${top}stringtie_output/${file::(-16)}"

stringtie ${file::(-16)}.bam -o ${top}stringtie_output/${file::(-16)}.gtf \
-p 14 -G ${ori_GTF} -f 0.99 \
-C ${top}stringtie_output/${file::(-16)}_full_cov.gtf \
-b ${top}stringtie_output/${file::(-16)} 

echo "`date '+%Y-%m-%d_%H-%M-%S'` stringtie done"

SHI
EOF
)

    done
done

echo "Waiting for HISAT and StringTie jobs${ids} to complete"
srun -p blade-b -d afterok${ids} echo "HiSat and StringTie done. Starting cuffmerge"
 
#############################################################################

# # Custom cufflinks call

# ids=

# for serie in $series; do

#     # Library settings for cuffquant

#     if [[ $(contains "${unstr[@]}" "$serie") == "y" ]]; then
#         lib="fr-unstranded"
#     elif [[ $(contains "${str[@]}" "$serie") == "y" ]]; then
#         lib="fr-firststrand"
#     elif [[ $(contains "${mix[@]}" "$serie") == "y" ]]; then
#         lib="fr-unstranded"
#     fi

#     cd ${top}hisat_output
#     for file in $(ls *${serie}*.bam); do 

# ids=${ids}:$(sbatch --parsable << EOF
# #!/bin/bash
# #SBATCH -p blade-b
# #SBATCH -o ${logs}cufflinks.${serie}.${file::(-4)}.%j.out
# #SBATCH -o ${logs}cufflinks.${serie}.${file::(-4)}.%j.err
# #SBATCH --cpus-per-task=18 
# #SBATCH --job-name='cfflnks'

# ${SHIFTER} << SHI
# #!/bin/bash
# source /beegfs/scratch/bruening_scratch/pklemm/shifter/home/.bashrc
# cd ${top}

# # Cufflinks call

# module load cufflinks

# cd ${top}/hisat_output
# echo "`date '+%Y-%m-%d_%H-%M-%S'` Start cufflinks"
# echo $(pwd)
# echo "cufflinks -p 14 --library-type fr-firststrand -o cufflinks_output/${file::(-4)} -g ${ori_GTF} ${file}"
# cufflinks -p 14 --library-type fr-firststrand -o cufflinks_output/${file::(-4)} -g ${ori_GTF} ${file}

# echo "`date '+%Y-%m-%d_%H-%M-%S'` cufflinks done"

# SHI
# EOF
# )
#     done
# done

# srun -p blade-b -d afterok:${ids} echo "Done with cufflinks"

# #############################################################################


for serie in $series; do
    rm -rf ${tmp}assemblies_${serie}.txt

    cd ${top}stringtie_output
    mkdir -p full_coverage
    mv *_full_cov.gtf full_coverage
    
    # Select only transcripts which have full coverage

    cd full_coverage
    for gtf in $(ls *${serie}*.gtf); do
        readlink -f ${gtf} >> ${tmp}assemblies_${serie}.txt
    done

    cd ${top}
    mkdir -p cuffmerge_output/${serie}
    cmout=$(readlink -f cuffmerge_output/${serie})/
    echo ${serie}
id=$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH -p blade-b
#SBATCH -o ${logs}cuffmerge.${serie}.%j.out
#SBATCH --job-name='cffmrg'


${SHIFTER} << SHI
#!/bin/bash
source /beegfs/scratch/bruening_scratch/pklemm/shifter/home/.bashrc
cd ${top}

module load cufflinks
    
# Cuffmerge call

cuffmerge -p 2 \
-o ${cmout} --min-isoform-fraction 1.0 \
-g ${ori_GTF} -s ${genome} ${tmp}assemblies_${serie}.txt

SHI
EOF
)
done

srun -p blade-b -d afterok:${id} echo "Done with cuffmerge"

cd ${tmp}

#############################################################################

echo "Starting cuffquant"

ids=

for serie in $series; do

    # Library settings for cuffquant

    if [[ $(contains "${unstr[@]}" "$serie") == "y" ]]; then
        lib="fr-unstranded"
    elif [[ $(contains "${str[@]}" "$serie") == "y" ]]; then
        lib="fr-firststrand"
    elif [[ $(contains "${mix[@]}" "$serie") == "y" ]]; then
        lib="fr-unstranded"
    fi

    cd ${top}hisat_output
    for file in $(ls *${serie}*.bam); do 
rm -rf ${logs}quant_${file::(-4)}.*.out
ids=${ids}:$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH -p blade-b
#SBATCH --cpus-per-task=18 
#SBATCH -o ${logs}quant_${file::(-4)}.%j.out
#SBATCH --job-name='cffqnt'

${SHIFTER} << SHI       
#!/bin/bash
source /beegfs/scratch/bruening_scratch/pklemm/shifter/home/.bashrc
cd ${top}cuffquant_output
mkdir ${serie}
cd ${serie}
module load cufflinks

# Cuffquant call

cuffquant -p 18 --library-type ${lib} \
-o ${file::(-4)} \
${top}cuffmerge_output/${serie}/merged.gtf \
${top}hisat_output/${file}

SHI
EOF
)
    done
done


echo "Waiting for cuffquant jobs${ids} to complete"
srun -p blade-b -d afterok${ids} echo "Cuffquant done. Starting cuffdiff."


#############################################################################

#### cuff diff >>>> one section per serie ######

# TODO: Set serie
serie=MAFA
mkdir -p ${top}cuffdiff_output/${serie}
dout=$(readlink -f ${top}cuffdiff_output/${serie})
lib="fr-firststrand"

rm -rf ${logs}cuffdiff.${serie}.*.out
sbatch --parsable << EOF
#!/bin/bash
#SBATCH -p blade-b
#SBATCH --cpus-per-task=18 
#SBATCH -o ${logs}cuffdiff.${serie}.%j.out
#SBATCH --job-name='cffdff'

${SHIFTER} << SHI
#!/bin/bash
source /beegfs/scratch/bruening_scratch/pklemm/shifter/home/.bashrc

#!/bin/bash
cd ${qua}${serie}

module load cufflinks

# Cuffdiff call

# TODO: Set cxb files and -L labels
cuffdiff -p 18 --library-type ${lib} \
-L IN_Control,IN_Experiment,IP_Control,IP_Experiment \
-o ${dout} --dispersion-method per-condition \
${top}cuffmerge_output/${serie}/merged.gtf \
S_001-F_MAFA-L____IN-1C-____-REP_1/abundances.cxb,S_005-F_MAFA-L____IN-2C-____-REP_1/abundances.cxb,S_007-F_MAFA-L____IN-3C-____-REP_1/abundances.cxb,S_013-F_MAFA-L____IN-4C-____-REP_1/abundances.cxb \
S_003-F_MAFA-L____IN-1E-____-REP_1/abundances.cxb,S_009-F_MAFA-L____IN-2E-____-REP_1/abundances.cxb,S_011-F_MAFA-L____IN-3E-____-REP_1/abundances.cxb,S_015-F_MAFA-L____IN-4E-____-REP_1/abundances.cxb \
S_002-F_MAFA-L____IP-1C-____-REP_1/abundances.cxb,S_006-F_MAFA-L____IP-2C-____-REP_1/abundances.cxb,S_008-F_MAFA-L____IP-3C-____-REP_1/abundances.cxb,S_014-F_MAFA-L____IP-4C-____-REP_1/abundances.cxb \
S_004-F_MAFA-L____IP-1E-____-REP_1/abundances.cxb,S_010-F_MAFA-L____IP-2E-____-REP_1/abundances.cxb,S_012-F_MAFA-L____IP-3E-____-REP_1/abundances.cxb,S_016-F_MAFA-L____IP-4E-____-REP_1/abundances.cxb

SHI
EOF

#### END section

exit
