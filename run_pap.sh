#!/bin/bash
: '
██████  █████  ██████
██  ██ ██   ██ ██  ██
██████ ███████ ██████
██     ██   ██ ██
██     ██   ██ ██
                        

POLARtron Analysis Pipeline (PAP)
By default, PAP will...
> Use BWA MEM to align data against a custom SARS-CoV-2 reference sequence.
> Use SAMtools...
  > fixmate to fill in mate coordinates and insert size fields.
  > sort to sort alignments based on position.
  > merg to merge alignments into a singular alignments file.
  > markdup to mark reads believed to be PCR duplicates.
> Use an AWK script to seperate control from viral alignments based on the presence of specific SNPs.
> Use SAMtools addreplacerg to add tags to each read specifying if they are control or viral.
> Use SAMtools depth to calculate the depth per base of the viral alignments.
> Use SAMtools to calculate statistics about viral and total alignments.
> Use Samtools merge to merge the tagged control and viral alignments back into a singular alignments file.
'

# Define variables
DEFAULT_THREADS=16
DEFAULT_CLEAN="True"
PATHOGEN_NAME="SARS-CoV-2"
DEPTH_MAPQ=4
QCUTOFF=0
GOODBASECHANGE=1
PE=0
SEQSPLIT=1

# Define paths
DEFAULT_TOP_DIR=$(pwd)
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Define scripts
REMRECOMBO="${PIPELINE_DIR}/scripts/accugenomics/remRecombo.sh"
COMPILE_RESULT="${PIPELINE_DIR}/scripts/compile_pap_results.py"

# Define reference files
NT_TO_IS="${PIPELINE_DIR}/scripts/accugenomics/NT_IS_LOOKUP_TABLE_v0.4.2_seperate.txt"
REFERENCE="${PIPELINE_DIR}/sars_cov_2_accukit_ISv0.4.1/sars_cov_2_accukit_ISv0.4.1.fasta"

# Define help function
printHelpAndExit() {
    cat <<PRINTHELPANDEXIT
Program: POLARtron Analysis Pipeline (PAP)
Version: 1.0
Contact: Per A. Adastra <adastra.aspera.per@gmail.com>
Usage:
        $0 [options]
Options:

   -d  Top level directory which must contain a subdirectory (fastq/) with fastq files (Default: $DEFAULT_TOP_DIR)
   -t  Number of threads for BWA alignment (Default: $DEFAULT_THREADS)
   -c  Clean up after pipeline completes (Default: $DEFAULT_CLEAN)
   -h  Print this help and exit

PRINTHELPANDEXIT
exit
}

# Parse user input
while getopts "d:t:c:h" opt;
do
    case $opt in
        d) TOP_DIR=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        c) CLEAN=$OPTARG ;;
        h) printHelpAndExit ;;
        *) printHelpAndExit ;;
    esac
done

# Define threads, top directory and clean options based on user input
if [ -z "${THREADS}" ]; then THREADS=$DEFAULT_THREADS; fi
if [ -z "${CLEAN}" ]; then CLEAN=$DEFAULT_CLEAN; fi
if [ -z "${TOP_DIR}" ]; then TOP_DIR=$DEFAULT_TOP_DIR; fi

# Define library name as the directory name
LIB_NAME=$(echo $TOP_DIR | awk -F "/" '{print $NF}')

# Define work directory
export WORK_DIR="${TOP_DIR}/pap"

# Check basic dependencies
echo "ʕ·ᴥ·ʔ : Checking dependencies..."

command -v bwa >/dev/null 2>&1 || { echo >&2 "ʕ·ᴥ·ʔ : BWA required but it's not installed!"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "ʕ·ᴥ·ʔ : Samtools required but it's not installed!"; exit 1; }
command -v python >/dev/null 2>&1 || { echo >&2 "ʕ·ᴥ·ʔ : Python required but it's not installed!"; exit 1; }

# Check for data (FASTQ) files
READ1_STR="R1"
READ2_STR="R2"
FASTQ_DIR="${TOP_DIR}/fastq"

if ls $FASTQ_DIR 1> /dev/null 2>&1;
then
    echo "ʕ·ᴥ·ʔ : Looking for fastq files...fastq files exist"
    testname=$(ls -l ${FASTQ_DIR})
    if [ "${testname: -3}" == ".gz" ]
    then
      read1=${TOP_DIR}"/fastq/*${READ1_STR}*.fastq.gz"
    else
      read1=${TOP_DIR}"/fastq/*${READ1_STR}*.fastq"
    fi
else
    echo "ʕ·ᴥ·ʔ : Failed to find any files matching ${FASTQ_DIR}"
    exit
fi

# Check to make sure output folders do not already exist
if ! mkdir "${WORK_DIR}" >/dev/null 2>&1; then echo "ʕ·ᴥ·ʔ : Unable to create ${WORK_DIR}! Exiting!"; exit 1; fi
if ! mkdir "${WORK_DIR}/aligned">/dev/null 2>&1; then echo "ʕ·ᴥ·ʔ : Unable to create ${WORK_DIR}/aligned! Exiting!"; exit 1; fi
if ! mkdir "${WORK_DIR}/debug">/dev/null 2>&1; then echo "ʕ·ᴥ·ʔ : Unable to create ${WORK_DIR}/debug! Exiting!"; exit 1; fi
if ! mkdir "${WORK_DIR}/final">/dev/null 2>&1; then echo "ʕ·ᴥ·ʔ : Unable to create ${WORK_DIR}/final! Exiting!"; exit 1; fi

# Create an array comprised of a FASTQ files
declare -a read1files=()
declare -a read2files=()

for i in ${read1}
do
    ext=${i#*$READ1_STR}
    name=${i%$READ1_STR*}
    name1=${name}${READ1_STR}
    name2=${name}${READ2_STR}
    read1files+=($name1$ext)
    read2files+=($name2$ext)
done

#########################################################################################################################
## WORK STARTS BELLOW
#########################################################################################################################

##### First block of work: Alignment of reads to reference
echo "ʕ·ᴥ·ʔ : Aligning files matching $FASTQ_DIR to $PATHOGEN_NAME reference assembly"

for ((i = 0; i < ${#read1files[@]}; ++i)); do
    file1=${read1files[$i]}
    file2=${read2files[$i]}

    FILE=$(basename ${file1%$read1str})
    ALIGNED_FILE=${WORK_DIR}/aligned/${FILE}"_mapped_temp"

    # Align reads to viral reference
    bwa mem -t $THREADS $REFERENCE $file1 $file2 > $ALIGNED_FILE".sam" 2> ${WORK_DIR}/debug/align.out

    # Fill in mate coordinates for deduping
    samtools fixmate -m $ALIGNED_FILE".sam" $ALIGNED_FILE"_matefixd.sam" 2> ${WORK_DIR}/debug/matefix.out

    # Sort reads based on position for deduping
    samtools sort -o $ALIGNED_FILE"_matefixd_sorted.sam" $ALIGNED_FILE"_matefixd.sam" 2> ${WORK_DIR}/debug/sort.out

done

# Merge BAMs if multiple SAMs were generated
samtools merge "${WORK_DIR}/aligned/"*"_matefixd_sorted.sam" -o "${WORK_DIR}/aligned/temp_sorted_merged.sam" 2> ${WORK_DIR}/debug/first_merge.out

# Mark duplicates
samtools markdup "${WORK_DIR}/aligned/temp_sorted_merged.sam" "${WORK_DIR}/aligned/temp_sorted_merged_dupd.sam" 2> ${WORK_DIR}/debug/dedup.out

# Classify alignments
echo "ʕ·ᴥ·ʔ : Removing recombinants reads..."
"${REMRECOMBO}" "${NT_TO_IS}" "${WORK_DIR}/aligned/temp_sorted_merged.sam" "${QCUTOFF}" "${GOODBASECHANGE}" "${PE}" "${SEQSPLIT}" 2> ${WORK_DIR}/debug/recombo.out

# Compress files
for SAM in ${WORK_DIR}/aligned/*sam;
do
    samtools view -hb $SAM > ${SAM%.sam}".bam"
done

# Mark alignments with classification
samtools addreplacerg -r "@RG\tID:GOOD" "${WORK_DIR}/aligned/temp_sorted_merged-good.bam" -o "${WORK_DIR}/aligned/temp_sorted_merged-good_marked.bam" 2> ${WORK_DIR}/debug/good_mark.out
samtools addreplacerg -r "@RG\tID:BAD" "${WORK_DIR}/aligned/temp_sorted_merged-bad.bam" -o "${WORK_DIR}/aligned/temp_sorted_merged-bad_marked.bam" 2> ${WORK_DIR}/debug/bad_mark.out
samtools addreplacerg -r "@RG\tID:CC" "${WORK_DIR}/aligned/temp_sorted_merged-cc.bam" -o "${WORK_DIR}/aligned/temp_sorted_merged-cc_marked.bam" 2> ${WORK_DIR}/debug/cc_mark.out
samtools addreplacerg -r "@RG\tID:UKN" "${WORK_DIR}/aligned/temp_sorted_merged-ukn.bam" -o "${WORK_DIR}/aligned/temp_sorted_merged-ukn_marked.bam" 2> ${WORK_DIR}/debug/ukn_mark.out
samtools addreplacerg -r "@RG\tID:IS" "${WORK_DIR}/aligned/temp_sorted_merged-IS.bam" -o "${WORK_DIR}/aligned/temp_sorted_merged-IS_marked.bam" 2> ${WORK_DIR}/debug/IS_mark.out

# Get coverage of viral reference from read catagories
echo "ʕ·ᴥ·ʔ : Analyzing coverage..."

# Get depth per base
samtools depth -a -Q "${DEPTH_MAPQ}" "${WORK_DIR}/aligned/temp_sorted_merged-good_marked.bam" | awk '$1=="MN908947.3"' > "${WORK_DIR}/aligned/viral_depth_per_base.txt" 2> "${WORK_DIR}/debug/viral_depth.out"

# Gather alignment qc statistics for all reads
echo "ʕ·ᴥ·ʔ :samtools flagstat result" > ${WORK_DIR}/aligned/all_alignment_stats.txt
samtools flagstat "${WORK_DIR}/aligned/temp_sorted_merged.bam"  >> "${WORK_DIR}/aligned/all_alignment_stats.txt" 2> "${WORK_DIR}/debug/all_stats.out"
echo "ʕ·ᴥ·ʔ : samtools stats result " >> ${WORK_DIR}/aligned/all_alignment_stats.txt
samtools stats "${WORK_DIR}/aligned/temp_sorted_merged.bam"  >> "${WORK_DIR}/aligned/all_alignment_stats.txt" 2> "${WORK_DIR}/debug/all_stats.out"

# Gather alignment qc statistics for all viral reads
echo "ʕ·ᴥ·ʔ :samtools flagstat result" > "${WORK_DIR}/aligned/viral_alignment_stats.txt"
samtools flagstat "${WORK_DIR}/aligned/temp_sorted_merged-good_marked.bam"  >> "${WORK_DIR}/aligned/viral_alignment_stats.txt" 2> "${WORK_DIR}/debug/viral_stats.out"
echo "ʕ·ᴥ·ʔ : samtools stats result " >> "${WORK_DIR}/aligned/viral_alignment_stats.txt"
samtools stats "${WORK_DIR}/aligned/temp_sorted_merged-good_marked.bam" >> "${WORK_DIR}/aligned/viral_alignment_stats.txt" 2> "${WORK_DIR}/debug/viral_stats.out"

# Write results to a file
echo "ʕ·ᴥ·ʔ : Compiling results"
python "$COMPILE_RESULT" $LIB_NAME $WORK_DIR


# Combine classified alignments into single file
samtools merge "${WORK_DIR}/aligned/"*"marked.bam" -o "${WORK_DIR}/aligned/classified_alignments.bam" 2> "${WORK_DIR}/debug/final_merge.out"

# Clean up
if [[ "${CLEAN}" == "True" ]];
then
  echo "ʕ·ᴥ·ʔ : Cleaning up files."
  rm "${WORK_DIR}/aligned/"*"temp"*
  rm "${WORK_DIR}/aligned/viral_depth_per_base.txt"
fi

echo "ʕ·ᴥ·ʔ : Pipeline completed, check ${WORK_DIR}/final for result"