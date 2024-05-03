#! /usr/bin/env bash
# Chadarevian_2024_Leukoencephalopathy_RNASeqProcessing.sh version 1.0
# Author: Jonathan Hasselmann (2024)
#
# ======================================= Program Description ====================================== #
# This script the associated pipeline are designed to run on Linux systems and in conjuction with the 
# complete repository located at https://github.com/jhasselmann808/Chadarevian_2024_Leukoencephalopathy
# in addition to the reference files located at https://drive.google.com/drive/folders/1gzf8VM6-btcsj0bY6BxgDQj6_EOcFGLD?usp=sharing. 
# That will provide all of the necessary packages, scripts, metadata, and reference files needed to
# generate the analysis that was used in the Chadarevian et al. (Neuron, 2024) publication. The FASTQ
# files will be automatically downloaded from GEO (GSE253623) and renamed to work with the downstream 
# analysis. If you have previously downloaded the FASTQ files, you will need to ensure that the file 
# names align with the filenames listed in the "./Metadata/Chadarevian_2024_Leukoencephalopathy_Bulk_Metadata.tsv" file
#
# You can set the number of threads (NTHREADS) and the amount of RAM (RAMGIG) to match with your system.
# As noted below, this pipeline takes a substantial amount of time to run on our system, so be aware of
# that before beginning the analysis.
#
# This script will call both Python and R scripts. The required R packages will be installed automatically
# but the script will need to be called from a Python environment with access to the relevant packages.
# See the README for details on which Python packages are required.
#
# This script uses a number of tools made available by other people/groups. If you use parts of this
# script for subsequent analysis, please be sure to properly cite the original creators of the tools.
#
#	FASTQC
#		Andrews, S. (2014). FastQC: A quality control tool for high throughput sequencing data.
#
#	BBDUK
#		Bushnell B. - sourceforge.net/projects/bbmap/
#
#	byPrim.py
#		Song X, Gao HY, Herrup K, Hart RP. Optimized splitting of mixed-species RNA sequencing data. 
#		J Bioinform Comput Biol. 2022 Apr;20(2):2250001. doi: 10.1142/S0219720022500019.
#
#	HISAT2
#		Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and 
#		HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). doi: 10.1038/s41587-019-0201-4
#
#	Kallisto
#		Nicolas L Bray, Harold Pimentel, Pall Melsted and Lior Pachter. Near-optimal
#		probabilistic RNA-seq quantification. Nature Biotechnology. 34, 525-27 (2016)
#		doi:10.1038/nbt.3519
#
#	DESeq2/tximport
#		Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and
#		dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#				    
#		Soneson C, Love MI, Robinson MD (2015). “Differential analyses for RNA-seq: 
#		transcript-level estimates improve gene-level inferences.” F1000Research, 4. 
#		doi: 10.12688/f1000research.7563.1.
# 	
#	Tidyverse
#		Wickham H, Averick M, Bryan J, et al. (2019). “Welcome to the tidyverse.” 
# 		Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686. 
#	
#	ggplot2
#		Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag 
#		New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org. 
#
#
# The steps executed by this script are as follows:
# 	1. Create the required output directories
#	2. Download files from GEO and output md5sum values
# 	3. Perform QC on the original FASTQ files using FASTQC
#	4. Perform FASTQ preprocessing using BBDUK
# 	5. Perform QC on the preprocessed FASTQ files using FASTQC
#	6. Align the preprocessed FASTQ files to the Mouse/Human genome using HISAT2
# 	7. Classify each read as Mouse or Human using the "byPrim.py" script
#	8. Sort the FASTQ files into Mouse and Human specific reads
# 	9. Perform QC on the Mouse and Human-specific FASTQ files using FASTQC
#	10. Pseudoalign the Human FASTQ files to the Human transcriptome using Kallisto
# 	11. Pseudoalign the Mouse FASTQ files to the Mouse transcriptome using Kallisto
#	12. Perform differential genes expression analysis and data plot generation using R/DESeq2
# 
# ================================================================================================== #

# ===============================================================================================================================
# ====================================================== Main Code Runs Below ===================================================
# ===============================================================================================================================

echo "===================== Welcome to the Chadarevian 2024 Data Processor ====================="
echo
echo "WARNING: This process takes about three weeks on our system using 50 threads and 126GB of RAM."
echo "WARNING: The FASTQ sorting steps will use around 30GB of memory per parallel file, and less than" 
echo "4 files analyzed in parallel (~120GB RAM) will substantially increase processing time."
echo
echo "If you are better at coding than I am, please improve upon this before running (lines 408-466 and/or MBJ_FASTQSorter.py scripts)."
echo

# Update the threads and RAM amounts to reflect the maximums available on your system
# Number of threads to use for processing steps
NTHREADS=40

# Gigabytes of RAM to use for BBDuk steps. Other steps will use all available memory
# Also check lines 376-385 for additional parallel processing RAM/file input limitation settings
RAMGIG=120

# ===============================================================================================================================
# ================================================ Mapping All Needed Directories ===============================================
# ===============================================================================================================================

# Collecting the directory for the script location and building input and output directory variables
ScriptDIR=$(dirname $0)

# Subsetting the ScriptDIR variable to point to the base directory
BASEDIR="${ScriptDIR%/*}"
echo "BASEDIR set to: ${BASEDIR}"
echo


# Creating package and script path variables
Fetch_PATH="${BASEDIR}/Scripts_Packages/SRA_Toolkit_v3.1.0/bin/prefetch"
FQ_PATH="${BASEDIR}/Scripts_Packages/SRA_Toolkit_v3.1.0/bin/fasterq-dump"
FASTQC_PATH="${BASEDIR}/Scripts_Packages/FastQC_v0.11.9/fastqc"
Reformat_PATH="${BASEDIR}/Scripts_Packages/BBMap_v38.79/reformat.sh"
BBDuk_PATH="${BASEDIR}/Scripts_Packages/BBMap_v38.79/bbduk.sh"
HISAT2_PATH="${BASEDIR}/Scripts_Packages/HISAT2_v2.2.1/hisat2"
Kallisto_PATH="${BASEDIR}/Scripts_Packages/Kallisto_v0.46.1/kallisto"
byPRIM_PATH="${BASEDIR}/Scripts_Packages/MBJ_byPrim.py"
Sorter_PATH="${BASEDIR}/Scripts_Packages/MBJ_FASTQSorter.py"
RScript_PATH="${BASEDIR}/Scripts_Packages/Chadarevian_2024_Leukoencephalopathy_DESeq2.R"


# Creating variables for input directories
BBDuk_Adapter_Reference="${BASEDIR}/Scripts_Packages/BBMap_v38.79/resources/adapters.fa" 		# BBDuk Adapter Reference
BBDuk_Filt_Reference="${BASEDIR}/References/BBDuk/PhiX_rRNA.fasta"  					# BBDuk Filtering Reference
HISAT2_IDX="${BASEDIR}/References/HISAT2/Index"								# HISAT2 Mixed Species Reference
Kallisto_IDX_Hu="${BASEDIR}/References/Kallisto/GRCh38_99/Homo_sapiens.GRCh38.99.cdna.all.idx" 		# Kallisto Human Reference
Kallisto_IDX_Ms="${BASEDIR}/References/Kallisto/GRCm39_107/Mus_musculus.GRCm39.107.cdna.all.idx" 	# Kallisto Mouse Reference


# Creating variables for output directories
GEO_Output="${BASEDIR}/Original_Fastq/"
FASTQC_Output_1="${BASEDIR}/Original_Fastq/QC"
Preprocessed_Output="${BASEDIR}/Original_Fastq/Preprocessed"
Preprocessed_Stats_Output="${Preprocessed_Output}/Stats"
FASTQC_Output_2="${Preprocessed_Output}/QC"
HISAT2_Output="${BASEDIR}/HISAT2_Output"
byPrim_Output="${BASEDIR}/byPrim_Output"
Sorted_FASTQ_Output="${BASEDIR}/Sorted_Fastq"
FASTQC_Output_3_Hu="${Sorted_FASTQ_Output}/Human/QC"
FASTQC_Output_3_Ms="${Sorted_FASTQ_Output}/Mouse/QC"
Kallisto_Output_Hu="${BASEDIR}/Kallisto_GRCh38_99"
Kallisto_Output_Ms="${BASEDIR}/Kallisto_GRCm39_107"


# Saving output variable names to an array
OUTARR=("GEO_Output" "FASTQC_Output_1" "Preprocessed_Output" "Preprocessed_Stats_Output" "FASTQC_Output_2" "HISAT2_Output" "byPrim_Output" "Sorted_FASTQ_Output" "FASTQC_Output_3_Hu" "FASTQC_Output_3_Ms" "Kallisto_Output_Hu" "Kallisto_Output_Ms")


# Creating folders for each output directory in OUTARR
for outdir in "${OUTARR[@]}"
do
	# Indirect expansion of current array item
	currDir=${!outdir}
	# If the current directory doesn't exist, create it with any necessary parent directories
	if [ ! -d "${currDir}" ]
	then
		echo "${currDir} doesn't exist. Creating directory..."
		echo
		mkdir -p "${currDir}"
	fi

done


# ===============================================================================================================================
# ================================================== Define Semaphore Functions =================================================
# ===============================================================================================================================

# Using semaphores to parallelize portions of the code on a FIFO basis
# initialize a semaphore with a given number of tokens
open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}

# run the given command asynchronously and pop/push tokens
run_with_lock(){
    local x
    # this read waits until there is something to read
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    # push the return code of the command to the semaphore
    printf '%.3d' $? >&3
    )&
}


# ===============================================================================================================================
# ============================================== Download and Rename Files from GEO =============================================
# ===============================================================================================================================

# Read the GEO Accession file into an array
IFS=$'\n' read -d '' -r -a ACCESSION < "${BASEDIR}/GEO_Accession.txt"

# Set the number of semaphores that can be used based on available RAM and THREADS
if (( "${NTHREADS}" < "${#ACCESSION[*]}" ))
then
	N="${NTHREADS}"
else
	N="${#ACCESSION[*]}"
fi

# Adjust the number of threads and RAM each semaphore can use
let "GEO_RAM = ${RAMGIG} / $N"
let "GEO_THREADS = ${NTHREADS} / $N"

# Move into the GEO download directory
cd "${GEO_Output}"

# Initiate semaphores and prefetch the SRA files from GEO
open_sem ${N}
for ((i=0;i<((${#ACCESSION[*]}));++i))
do
	
	# Download the current accession from GEO and split into FASTQ files
	run_with_lock "${Fetch_PATH}" "${ACCESSION[i]}"

done && wait

# Initiate semaphores and split the SRA files into FASTQ files
open_sem ${N}
for ((i=0;i<((${#ACCESSION[*]}));++i))
do
	
	# Split the prefetch SRA files into FASTQ files
	run_with_lock "${FQ_PATH}" "${ACCESSION[i]}" --split-files --mem "${GEO_RAM}" --threads "${GEO_THREADS}" -t "${BASEDIR}/temp/scratch"

done && wait

# Remove the SRA prefetch directories
for ((i=0;i<((${#ACCESSION[*]}));++i))
do
	# Recursively remove each SRA directory
	rm -r "${ACCESSION[i]}"

done

# Compressing the fastq files using gzip
for file in ./*.fastq
do
	pigz --processes "${NTHREADS}" "${file}"
done

# Changing all _1 to _R1 to make file selection easier later
for file in *'_1.'*
do
	[ -f "$file" ] || continue
	
	mv "$file" "${file//_1./_R1.}"
done
	
# Changing all _2 to _R2 to make file selection easier late
for file in *'_2.'*
do
	[ -f "$file" ] || continue
	
	mv "$file" "${file//_2./_R2.}"
done

# Checking md5sum hash values
echo "Checking md5sum 'hash' values..."
echo "Please be patient..."
find . -name "*.gz" | parallel md5sum >> GEO_Download_md5sum.txt
echo


# ===============================================================================================================================
# ================================================ FASTQC on Original FASTQ Files ===============================================
# ===============================================================================================================================

# Running FASTQC on the original FASTQ files
echo "Running FASTQC on the original FASTQ files..."
echo

find . -maxdepth 1 -name "*.gz" | 
xargs "${FASTQC_PATH}" -t ${NTHREADS} --outdir="${FASTQC_Output_1}"


# ===============================================================================================================================
# ========================================== BBDuk Preprocessing on Original FASTQ Files ========================================
# ===============================================================================================================================

# Collecting all Read 1 files in the working directory and saving each to a numbered variable
# This should find all files with a _R1 before the extension then remove the extension and set
# the filename as an item in a numbered list.
# This should also work for both .fastq.gz and .txt.gz
for filename in ./*_R1*.gz
do

	# Removing the extension from each filename
	_name="${filename##*/}"
	_noExt="${_name%%.*}"

	# Saving each forward read filename into an array
	FOR+=("${_name}")
	ForFile+=("${_noExt}")
	
	# Extracting the file extension for this dataset (e.g. fastq.gz, txt.gz. fq.gz, etc.)
	# This only needs to occur once
	if [ -z "{_ext}" ]
	then
		_ext="${filename##*.}"
	fi

done


# Repeating the above for loop but now collecting all Read 2 files
for filename in ./*_R2*.gz
do

	# Removing the extension from each filename
	_name="${filename##*/}"
	_noExt="${_name%%.*}"
	
	# Saving each reverse read filename into an array
	REV+=("${_name}")
	RevFile+=("${_noExt}")

done


# Using the file arrays made in the previous for loops as input for the BBDuk commands
echo "Running BBDUK on the original FASTQ files..."
echo

for ((i=0;i<((${#FOR[*]}));++i))
do

	# Interleaving the files, trimming adapters, filtering low quality and contaminating reads/bases, and splitting the files back into forward and reverse files
	"${Reformat_PATH}" -Xmx${RAMGIG}g in1="${FOR[i]}" in2="${REV[i]}" out=stdout.fastq int=f |
	"${BBDuk_PATH}" -Xmx${RAMGIG}g in=stdin.fastq out=stdout.fastq int=t stats="${Preprocessed_Stats_Output}""/adaptstats_""${RevFile[i]}".txt ref="${BBDuk_Adapter_Reference}" ktrim=r mink=8 k=23 hdist=1 overwrite=true |
	"${BBDuk_PATH}" -Xmx${RAMGIG}g in=stdin.fastq out=stdout.fastq int=t stats="${Preprocessed_Stats_Output}""/filtStats_""${RevFile[i]}".txt ref="${BBDuk_Filt_Reference}" k=31 hdist=1 qtrim=rl trimq=28 minlen=45 overwrite=true |
	"${Reformat_PATH}" -Xmx${RAMGIG}g in=stdin.fastq vint ain out1="${Preprocessed_Output}""/QualTrim_""${FOR[i]}" out2="${Preprocessed_Output}""/QualTrim_""${REV[i]}" int=t overwrite=true

done


# ===============================================================================================================================
# ============================================== FASTQC on Preprocessed FASTQ Files =============================================
# ===============================================================================================================================

# Running FASTQC on the preprocessed FASTQ files
echo "Running FASTQC on the preprocessed FASTQ files..."
echo

cd "${Preprocessed_Output}"
find . -maxdepth 1 -name "*.gz" | 
xargs "${FASTQC_PATH}" -t ${NTHREADS} --outdir="${FASTQC_Output_2}"


# ===============================================================================================================================
# ========================================= HISAT2 Alignment on Preprocessed FASTQ Files ========================================
# ===============================================================================================================================

echo "Running HISAT2 on the preprocessed FASTQ files..."
echo

# Creating the list variables
FOR=()
ForFile=()
REV=()

# Collecting all Forward read files in the working directory
# File names are stored in an array
for filename in ./*_R1*.gz
do
	# Removing the extension from each filename
	_name="${filename##*/}"
	_interim="${_name##*QualTrim_}"
	_noExt="${_interim%%_R1*}"
	FOR+=("${_name}")
	ForFile+=("${_noExt}")
done

# Repeating the above for loop but now collecting all Reverse read files
# File names are stored in an array
for filename in ./*_R2*.gz
do
	# Removing the extension from each filename
	_name="${filename##*/}"
	REV+=("${_name}")
done

# Change directory to the HISAT2 index directory
cd "${HISAT2_IDX}"

# Using length of FOR[] to loop over all forward and reverse files
for ((i=0;i<((${#FOR[*]}));++i))
do
	
	echo "Performing HISAT2 alignment for file ${ForFile[i]}"
	echo
	# Performing HISA2 Alignment to human and mouse mixed species genome
	"${HISAT2_PATH}" -q -t --new-summary --met-file "${HISAT2_Output}/${ForFile[i]}_Metrics.txt" --threads ${NTHREADS} -x "Mixed_Sapiens_Musculus" -1 "${Preprocessed_Output}/${FOR[i]}" -2 "${Preprocessed_Output}/${REV[i]}" -S "${HISAT2_Output}/${ForFile[i]}.sam"	
done	


# ===============================================================================================================================
# ============================================= byPrim Classifying on HISAT2 Output =============================================
# ===============================================================================================================================

# Classifying reads as human or mouse using the MBJ_byPRIM.py script
echo "Running primary classification on the HISAT2 outputs..."
echo

cd "${byPrim_Output}"

# The conditionals below will allocate semaphores based on available memory, number of files, and available threads
# Assuming ~35GB RAM per file
if (( "${RAMGIG}" > 35 ))
then
	N=$(("${RAMGIG}"/35))
else
	N=1
fi

# If memory allows for more threads than there are files, reset number of threads to the number of files
if (( "$N" > "${#ForFile[@]}" ))
then
	N="${#ForFile[@]}"

# If memory or number of files allow for more threads that were assigned on line 92, set to NTHREADS
elif (( "$N" > "${NTHREADS}" ))
then
	N="${NTHREADS}"
fi

# Running the MBJ_byPrim classification
open_sem ${N}
for ((i=0;i<((${#ForFile[*]}));++i))
do

    # Calling the MBJ_byPRIM.py script for each HISAT2 output file
    run_with_lock python "${byPRIM_PATH}" --genome Musculus --altgenome Sapiens --sample ${ForFile[i]} --File "${HISAT2_Output}/${ForFile[i]}.sam"

done && wait


# ===============================================================================================================================
# ===================================== Sorting Preprocessed FASTQ on Species Classification ====================================
# ===============================================================================================================================

# Running FASTQC on the Human FASTQ files
echo
echo "Sorting preprocessed FASTQ files based on species..."

# Running the MBJ_FASTQSorter
# Using semaphores to parallelize the sorting on a FIFO basis
open_sem ${N}
for ((i=0;i<((${#ForFile[*]}));++i))
do

    # Calling the MBJ_FASTQSorter.py script for each preprocessed FASTQ file
    run_with_lock python "${Sorter_PATH}" --sample "${ForFile[i]}" --FastIn "${Preprocessed_Output}" --FastOut "${Sorted_FASTQ_Output}"

done && wait


# ===============================================================================================================================
# ============================================ FASTQC on Human and Mouse FASTQ Files ============================================
# ===============================================================================================================================

# Running FASTQC on the Human FASTQ files
echo "Running FASTQC on the sorted Human FASTQ files..."
echo

cd "${Sorted_FASTQ_Output}/Human"

find . -maxdepth 1 -name "*.gz" | 
xargs "${FASTQC_PATH}" -t ${NTHREADS} --outdir="${FASTQC_Output_3_Hu}"



# Running FASTQC on the Mouse FASTQ files
echo "Running FASTQC on the sorted Mouse FASTQ files..."
echo

cd "${Sorted_FASTQ_Output}/Mouse"
find . -maxdepth 1 -name "*.gz" | 
xargs "${FASTQC_PATH}" -t ${NTHREADS} --outdir="${FASTQC_Output_3_Ms}"


# ===============================================================================================================================
# ======================================== Kallisto Pseudoalignment on Human FASTQ Files ========================================
# ===============================================================================================================================

# Running Kallisto on the Human FASTQ files
echo
echo "Running Kallisto on the sorted Human FASTQ files..."
echo

cd "${Sorted_FASTQ_Output}/Human"

# Resetting the file list variables
FOR=()
ForFile=()
REV=()

# Collecting all Forward read files in the working directory
# File names are stored in an array
for filename in ./*_R1*.gz
do

	# Removing the extension from each filename
	_name="${filename##*/}"
	_interim="${_name##*QualTrim_}"
	_noExt="${_interim%%_R1*}"
	
	# Saving the filenames into a arrays
	FOR+=("${_name}")
	ForFile+=("${_noExt}")

done


# Repeating the above for loop but now collecting all Reverse read files
# File names are stored in an array
for filename in ./*_R2*.gz
do

	# Removing the extension from each filename
	_name="${filename##*/}"
	
	# Saving the filenames into an array
	REV+=("${_name}")

done


# Running Kalliso pseudoalignment without bootstraps
# Using length of FOR[] to loop over all forward and reverse files
for ((i=0;i<((${#FOR[*]}));++i))
do
	echo "Performing pseudoalignment on file ${ForFile[i]}"
	echo
	# Performing Kallisto alignement without bootstraps on all compiled files
	"${Kallisto_PATH}" quant -i "${Kallisto_IDX_Hu}" -o "${Kallisto_Output_Hu}"/"${ForFile[i]}" -t ${NTHREADS} "${FOR[i]}" "${REV[i]}"	

done


# ===============================================================================================================================
# ======================================== Kallisto Pseudoalignment on Mouse FASTQ Files ========================================
# ===============================================================================================================================

# Running Kallisto on the Human FASTQ files
echo "Running Kallisto on the sorted Mouse FASTQ files..."
echo

cd "${Sorted_FASTQ_Output}/Mouse"


# Resetting the file list variables
FOR=()
ForFile=()
REV=()


# Collecting all Forward read files in the working directory
# File names are stored in an array
for filename in ./*_R1*.gz
do

	# Removing the extension from each filename
	_name="${filename##*/}"
	_interim="${_name##*QualTrim_}"
	_noExt="${_interim%%_R1*}"
	
	# Saving the filenames into a arrays
	FOR+=("${_name}")
	ForFile+=("${_noExt}")

done


# Repeating the above for loop but now collecting all Reverse read files
# File names are stored in an array
for filename in ./*_R2*.gz
do

	# Removing the extension from each filename
	_name="${filename##*/}"
	
	# Saving the filenames into an array
	REV+=("${_name}")

done


# Running Kalliso pseudoalignment without bootstraps
# Using length of FOR[] to loop over all forward and reverse files
for ((i=0;i<((${#FOR[*]}));++i))
do
	echo "Performing pseudoalignment on file ${ForFile[i]}"
	echo
	# Performing Kallisto alignement without bootstraps on all compiled files
	"${Kallisto_PATH}" quant -i "${Kallisto_IDX_Ms}" -o "${Kallisto_Output_Ms}"/"${ForFile[i]}" -t ${NTHREADS} "${FOR[i]}" "${REV[i]}"	

done

echo "Finished the alignment pipeline!"
echo "Moving on to the DESeq2 Analysis pipeline..."
echo

# ===============================================================================================================================
# ===================================== Differential Gene Expression Analysis Using DESeq2 ======================================
# ===============================================================================================================================

# Calling the DESeq2 R script and importing the BASEDIR variable
Rscript "${RScript_PATH}" "${BASEDIR}"
